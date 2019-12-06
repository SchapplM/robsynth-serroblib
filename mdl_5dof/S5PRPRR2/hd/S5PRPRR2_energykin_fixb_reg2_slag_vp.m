% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:12
% EndTime: 2019-12-05 15:45:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (135->24), mult. (279->69), div. (0->0), fcn. (172->8), ass. (0->29)
t15 = qJD(2) + qJD(4);
t14 = t15 ^ 2;
t31 = t14 / 0.2e1;
t25 = cos(qJ(2));
t13 = qJD(2) * pkin(2) + t25 * qJD(1);
t18 = sin(pkin(9));
t19 = cos(pkin(9));
t22 = sin(qJ(2));
t29 = qJD(1) * t22;
t11 = t18 * t13 + t19 * t29;
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t10 = t19 * t13 - t18 * t29;
t9 = qJD(2) * pkin(3) + t10;
t6 = t24 * t11 + t21 * t9;
t5 = -t21 * t11 + t24 * t9;
t3 = -t15 * pkin(4) - t5;
t30 = t15 * t3;
t28 = qJD(5) * t15;
t27 = qJD(1) * qJD(2);
t26 = qJD(1) ^ 2;
t23 = cos(qJ(5));
t20 = sin(qJ(5));
t17 = qJD(2) ^ 2 / 0.2e1;
t16 = qJD(3) ^ 2 / 0.2e1;
t4 = t15 * pkin(7) + t6;
t2 = t20 * qJD(3) + t23 * t4;
t1 = t23 * qJD(3) - t20 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t26 / 0.2e1, 0, 0, 0, 0, 0, t17, t25 * t27, -t22 * t27, 0, (t22 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * t26, 0, 0, 0, 0, 0, t17, t10 * qJD(2), -t11 * qJD(2), 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t16, 0, 0, 0, 0, 0, t31, t5 * t15, -t6 * t15, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16, t20 ^ 2 * t31, t20 * t14 * t23, t20 * t28, t23 ^ 2 * t31, t23 * t28, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t23 * t30, -t2 * qJD(5) + t20 * t30, (-t1 * t20 + t2 * t23) * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
