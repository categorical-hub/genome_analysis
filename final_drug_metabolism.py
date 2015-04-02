#!/usr/bin/python
import re
from sys import argv

allele_rating = {}
rsid_dict = {}
allele_dict = {}
drug_dict = {}

class color:
	PURPLE = '\033[95m'
	CYAN = '\033[96m'
	DARKCYAN = '\033[36m'
	BLUE = '\033[94m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
	END = '\033[0m'

# Lists of drugs by their metabolizing CYP enzyme
drugs_1A2 = ["amitriptyline", "caffeine", "clomipramine", "clozapine", "cyclobenzaprine", "estradiol", "fluvoxamine", "haloperidol", "imipramine", "N-DeMe", "mexilletine", "naproxen", "olanzapine", "ondansetron", "phenacetin_", "acetaminophen", "NAPQI", "propranolol", "riluzole", "ropivacaine", "tacrine", "theophylline", "tizanidine", "verapamil", "(R)warfarin", "zileuton", "zolmitriptan"]

drugs_2A6 = []

drugs_2B6 = ["bupropion", "cyclophosphamide", "efavirenz", "ifosfamide", "methadone"]

drugs_2C9 = ["diclofenac", "ibuprofen", "lornoxicam", "meloxicam", "S-naproxen_Nor", "piroxicam", "suprofen", "tolbutamide", "glipizide", "losartan", "irbesartan", "glyburide", "glibenclamide", "glipizide", "glimepiride", "tolbutamide", "amitriptyline", "celecoxib", "fluoxetine", "fluvastatin", "glyburide", "nateglinide", "phenytoin-4-OH2", "rosiglitazone", "tamoxifen", "torsemide", "S-warfarin"]

drugs_2D6 = ["carvedilol", "S-metoprolol", "propafenone", "timolol", "amitriptyline", "clomipramine", "desipramine", "imipramine", "paroxetine", "haloperidol", "perphenazine", "risperidone->9OH", "thioridazine", "zuclopenthixol", "alprenolol", "amphetamine", "aripiprazole", "atomoxetine", "bufuralol", "chlorpheniramine", "chlorpromazine", "codeine(O->desMe)", "debrisoquine", "dexfenfluramine", "dextromethorphan", "duloxetine", "encainide", "flecainide", "fluoxetine", "fluvoxamine", "lidocaine", "metoclopramide", "methoxyamphetamine", "mexilletine", "minaprine", "nebivolol", "nortriptyline", "ondansetron", "oxycodone", "perhexiline", "phenacetin", "phenformin", "promethazine", "propranolol", "sparteine", "tamoxifen", "tramadol", "venlafaxine"]


drugs_2E1 = ["enflurane", "halothane", "isoflurane", "methoxyflurane", "sevoflurane", "acetaminophen->NAPQI", "aniline2", "benzene", "chlorzoxazone", "ethanol", "N,N-dimethylformamide", "theophylline->8-OH"]

drugs_3A4 = ["clarithromycin", "erythromycin(not 3A5)", "telithromycin", "quinidine->3OH(not 3A5)", "alprazolam", "diazepam->3OH", "midazolam", "triazolam", "cyclosporine", "tacrolimus(FK506)", "indinavir", "nelfinavir", "ritonavir", "saquinavir", "cisapride", "astemizole", "chlorpheniramine", "terfenadine", "amlodipine", "diltiazem", "felodipine", "lercanidipine", "nifedipine2", "nisoldipine", "nitrendipine", "verapamil", "atorvastatin", "cerivastatin", "lovastatin", "simvastatin", "estradiol", "hydrocortisone", "progesterone", "testosterone", "alfentanyl", "aprepitant", "aripiprazole", "buspirone", "cafergot", "caffeine_TMU", "cilostazol", "cocaine", "codeineNdemethylation", "dapsone", "dexamethasone", "dextromethorphan", "docetaxel", "domperidone", "eplerenone", "fentanyl", "finasteride", "gleevec", "haloperidol", "irinotecan", "LAAM", "lidocaine", "methadone", "nateglinide", "ondansetron", "pimozide", "propranolol", "quetiapine", "quinine", "risperidone", "salmeterol", "sildenafil", "sirolimus", "tamoxifen", "taxol", "terfenadine", "trazodone", "vincristine", "zaleplon", "ziprasidone", "zolpidem"]

drugs_lookup = {
	'CYP1A2':drugs_1A2,
	'CYP2B6':drugs_2B6,
	'CYP2C9':drugs_2C9,
	'CYP2D6':drugs_2D6,
	'CYP2E1':drugs_2E1,
	'CYP3A4':drugs_3A4,
	'CYP2A6':drugs_2A6
}

# Load a 23andMe report to a global dictionary by rsid
def load_report_to_dictionary():
	if len(argv) < 2:
		print("usage: %s <genome data filename>" % argv[0])
		exit()
	
	f = open(argv[1])
	for line in iter(f):
		if line[0] == '#':
			continue
		splitted = line.split('\t')
		rsid = splitted[0]
		genotype = splitted[3].replace('\r\n', '')
		if rsid not in rsid_dict:
			rsid_dict[rsid] = genotype
		else:
			print('error! rsid encountered twice(??)')
			exit()
	f.close()
	return

# Load snpedia db, and for each known rsid check the user for a match
# if a match is found, add the drug to the "modified metabolizm" drug dictionary
# along with the reason for it (the affected CYP enzyme)
def analyze_report():
	f = open('snpedia.csv')
	for line in iter(f):
		item_arr = line.split(',')
		rsid = item_arr[0]
		desc = item_arr[3]
		if desc[0] == '"':
			desc = line[line.find(desc)+1:line.rfind(',')-1]
		link = item_arr[-1]
		match = re.search('\((.*?)\)$', link)
		if (match):
			genotype = match.group(1).replace(';', '').replace('(','').replace(')','')
			if genotype == '':
				genotype = 'BB'
		else:
			continue

		if 'CYP' in desc:
			# report contains known rsid
			if rsid in rsid_dict:
				#print(rsid_dict[rsid] + '<==>' + genotype)
				if (rsid_dict[rsid] == genotype):
					for search_str, result_dict in drugs_lookup.items():
						if (search_str in desc):
							for drug in result_dict:
								if drug in drug_dict:
									drug_dict[drug] = drug_dict[drug] + "," + desc
								else:
									drug_dict[drug] = desc

def main():
	load_report_to_dictionary()
	analyze_report()

	print(color.BOLD + "Drugs with modified metabolism for: " + argv[1] + color.END)
	for drug, desc in sorted(drug_dict.items()):
		print(color.BOLD + drug + color.END + " modified metabolism due to: " + color.PURPLE + desc + color.END)		

if __name__ == "__main__":
	main()
