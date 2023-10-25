# load MP measurements for low concentration estimation based on detection limits

# following Verbovšek 2011 [A comparison of parameters below the limit of detection in geochemical analyses by substitution methods]

atenolol = readRDS("input_data/RDS_files/Ate.rds") %>% print(., n =134)
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "atenolol")
benzotriazole = readRDS("input_data/RDS_files/Ben.rds") %>% 
  filter(., str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "benzotriazole")
carbamazepine = readRDS("input_data/RDS_files/Car.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "carbamazepine")
clarithromycin = readRDS("input_data/RDS_files/Cla.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "clarithromycin")
diclofenac = readRDS("input_data/RDS_files/Dic.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "diclofenac")
gabapentin = readRDS("input_data/RDS_files/Gab.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "gabapentin")
hydrochlorothiazide = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "hydrochlorothiazide")
Ibuprofen = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "Ibuprofen")
irbesartan = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "irbesartan")
ketoprofen = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "ketoprofen")
lidocaine = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "lidocaine")
metoprolol = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "metoprolol")
propranolol = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "propranolol")
sotalol = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "sotalol")
trimethoprim = readRDS("input_data/RDS_files/Ate.rds") %>% 
  filter(str_detect(measured_value, "<")) %>% 
  mutate(sampling_date = ymd(sampling_date)) %>% 
  select(sampling_date, measured_value) %>% 
  group_by(sampling_date) %>% 
  unique() %>% 
  mutate(DL = as.numeric(str_sub(measured_value, start = 2))) %>% 
  mutate(value = DL/sqrt(2),
         analysis = "trimethoprim")

DLs = 
bind_rows(atenolol,benzotriazole,carbamazepine,clarithromycin,diclofenac,gabapentin,hydrochlorothiazide,Ibuprofen, irbesartan, ketoprofen,lidocaine,metoprolol,propranolol,sotalol,trimethoprim)
