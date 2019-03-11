% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:53
% EndTime: 2019-03-09 02:13:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (418->50), mult. (893->110), div. (0->0), fcn. (555->6), ass. (0->40)
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t46 = sin(qJ(4));
t47 = cos(qJ(4));
t29 = (t43 * t47 + t44 * t46) * qJD(1);
t48 = qJD(1) ^ 2;
t41 = t48 / 0.2e1;
t34 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t50 = -pkin(7) * qJD(1) + t34;
t26 = t50 * t43;
t27 = t50 * t44;
t14 = t47 * t26 + t46 * t27;
t11 = qJD(4) * pkin(8) + t14;
t31 = (-t43 * t46 + t44 * t47) * qJD(1);
t36 = qJD(1) * qJ(2) + qJD(3);
t51 = qJD(1) * t43;
t32 = pkin(3) * t51 + t36;
t12 = t29 * pkin(4) - t31 * pkin(8) + t32;
t45 = sin(qJ(5));
t52 = cos(qJ(5));
t4 = t52 * t11 + t45 * t12;
t3 = -t45 * t11 + t52 * t12;
t13 = -t46 * t26 + t47 * t27;
t10 = -qJD(4) * pkin(4) - t13;
t40 = t44 ^ 2;
t39 = t43 ^ 2;
t37 = -qJD(1) * pkin(1) + qJD(2);
t28 = qJD(5) + t29;
t25 = t28 ^ 2 / 0.2e1;
t21 = t45 * qJD(4) + t52 * t31;
t19 = -t52 * qJD(4) + t45 * t31;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t16 = t21 * t28;
t15 = t19 * t28;
t7 = t21 * t19;
t5 = t19 * pkin(5) + qJD(6) + t10;
t2 = -t19 * qJ(6) + t4;
t1 = t28 * pkin(5) - t21 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t37 * qJD(1), t48 * qJ(2), qJ(2) ^ 2 * t41 + t37 ^ 2 / 0.2e1, t40 * t41, -t44 * t48 * t43, 0, t39 * t41, 0, 0, t36 * t51, t36 * t44 * qJD(1) (-t39 - t40) * t34 * qJD(1), t36 ^ 2 / 0.2e1 + (t39 / 0.2e1 + t40 / 0.2e1) * t34 ^ 2, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * qJD(4), t29 ^ 2 / 0.2e1, -t29 * qJD(4), qJD(4) ^ 2 / 0.2e1, t13 * qJD(4) + t32 * t29, -t14 * qJD(4) + t32 * t31, -t13 * t31 - t14 * t29, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t18, -t7, t16, t17, -t15, t25, t10 * t19 + t3 * t28, t10 * t21 - t4 * t28, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t18, -t7, t16, t17, -t15, t25, t1 * t28 + t5 * t19, -t2 * t28 + t5 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
