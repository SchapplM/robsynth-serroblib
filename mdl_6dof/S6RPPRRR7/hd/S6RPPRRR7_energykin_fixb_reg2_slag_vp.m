% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:46
% EndTime: 2019-03-09 02:33:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (604->52), mult. (1291->130), div. (0->0), fcn. (866->8), ass. (0->41)
t49 = qJD(1) ^ 2;
t41 = t49 / 0.2e1;
t43 = sin(pkin(10));
t33 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t50 = -pkin(7) * qJD(1) + t33;
t26 = t50 * t43;
t44 = cos(pkin(10));
t27 = t50 * t44;
t47 = sin(qJ(4));
t48 = cos(qJ(4));
t15 = -t47 * t26 + t48 * t27;
t30 = (-t43 * t47 + t44 * t48) * qJD(1);
t10 = qJD(4) * pkin(4) - t30 * pkin(8) + t15;
t16 = t48 * t26 + t47 * t27;
t28 = (-t43 * t48 - t44 * t47) * qJD(1);
t11 = t28 * pkin(8) + t16;
t46 = sin(qJ(5));
t53 = cos(qJ(5));
t6 = t46 * t10 + t53 * t11;
t52 = cos(qJ(6));
t51 = qJD(1) * t43;
t36 = qJD(1) * qJ(2) + qJD(3);
t31 = pkin(3) * t51 + t36;
t18 = -t53 * t28 + t46 * t30;
t21 = -t28 * pkin(4) + t31;
t5 = t53 * t10 - t46 * t11;
t45 = sin(qJ(6));
t40 = qJD(4) + qJD(5);
t39 = t44 ^ 2;
t38 = t43 ^ 2;
t37 = -pkin(1) * qJD(1) + qJD(2);
t20 = t46 * t28 + t53 * t30;
t17 = qJD(6) + t18;
t14 = t52 * t20 + t45 * t40;
t12 = t45 * t20 - t52 * t40;
t7 = t18 * pkin(5) - t20 * pkin(9) + t21;
t4 = t40 * pkin(9) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t52 * t4 + t45 * t7;
t1 = -t45 * t4 + t52 * t7;
t8 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t37 * qJD(1), t49 * qJ(2), qJ(2) ^ 2 * t41 + t37 ^ 2 / 0.2e1, t39 * t41, -t44 * t49 * t43, 0, t38 * t41, 0, 0, t36 * t51, t36 * t44 * qJD(1) (-t38 - t39) * t33 * qJD(1), t36 ^ 2 / 0.2e1 + (t38 / 0.2e1 + t39 / 0.2e1) * t33 ^ 2, t30 ^ 2 / 0.2e1, t30 * t28, t30 * qJD(4), t28 ^ 2 / 0.2e1, t28 * qJD(4), qJD(4) ^ 2 / 0.2e1, t15 * qJD(4) - t31 * t28, -t16 * qJD(4) + t31 * t30, -t15 * t30 + t16 * t28, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t40, t18 ^ 2 / 0.2e1, -t18 * t40, t40 ^ 2 / 0.2e1, t21 * t18 + t5 * t40, t21 * t20 - t6 * t40, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
