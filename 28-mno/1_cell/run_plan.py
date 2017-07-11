import pickle

if __name__ == '__main__':
  with open('plan.pickle','rb') as fin:
    plan = pickle.load(fin)

  plan.nextstep()

  report = plan.generate_report()

  for rep in report:
    print(rep['id'],rep.keys())
  
  import json
  json.dump(report,open('report.json','w'),indent=1)
