% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:57
% EndTime: 2019-03-09 01:44:57
% DurationCPUTime: 0.17s
% Computational Cost: add. (342->51), mult. (697->116), div. (0->0), fcn. (387->8), ass. (0->39)
t34 = sin(pkin(10));
t36 = cos(pkin(10));
t39 = sin(qJ(4));
t40 = cos(qJ(4));
t21 = (t34 * t40 + t36 * t39) * qJD(1);
t35 = sin(pkin(9));
t26 = (pkin(1) * t35 + qJ(3)) * qJD(1);
t41 = qJD(1) ^ 2;
t33 = t41 / 0.2e1;
t37 = cos(pkin(9));
t45 = -pkin(1) * t37 - pkin(2);
t24 = qJD(3) + (-pkin(7) + t45) * qJD(1);
t15 = -t39 * qJD(2) + t40 * t24;
t47 = qJ(5) * qJD(1);
t10 = qJD(4) * pkin(4) - t40 * t47 + t15;
t16 = t40 * qJD(2) + t39 * t24;
t11 = -t39 * t47 + t16;
t6 = t34 * t10 + t36 * t11;
t51 = pkin(1) * t41;
t50 = cos(qJ(6));
t49 = t26 * qJD(1);
t48 = t26 ^ 2 / 0.2e1;
t46 = qJD(1) * qJD(4);
t5 = t36 * t10 - t34 * t11;
t20 = t39 * qJD(1) * pkin(4) + qJD(5) + t26;
t38 = sin(qJ(6));
t32 = qJD(2) ^ 2 / 0.2e1;
t31 = qJD(4) ^ 2 / 0.2e1;
t25 = t45 * qJD(1) + qJD(3);
t23 = (-t34 * t39 + t36 * t40) * qJD(1);
t17 = qJD(6) + t21;
t14 = t38 * qJD(4) + t50 * t23;
t12 = -t50 * qJD(4) + t38 * t23;
t7 = t21 * pkin(5) - t23 * pkin(8) + t20;
t4 = qJD(4) * pkin(8) + t6;
t3 = -qJD(4) * pkin(5) - t5;
t2 = t38 * t7 + t50 * t4;
t1 = -t38 * t4 + t50 * t7;
t8 = [0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t37 * t51, -t35 * t51, 0, t32 + (t35 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t41, t33, 0, 0, 0, 0, 0, 0, t25 * qJD(1), t49, t32 + t48 + t25 ^ 2 / 0.2e1, t40 ^ 2 * t33, -t40 * t41 * t39, t40 * t46, t39 ^ 2 * t33, -t39 * t46, t31, t15 * qJD(4) + t39 * t49, -t16 * qJD(4) + t40 * t49 (-t15 * t40 - t16 * t39) * qJD(1), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t48, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * qJD(4), t21 ^ 2 / 0.2e1, -t21 * qJD(4), t31, t5 * qJD(4) + t20 * t21, -t6 * qJD(4) + t20 * t23, -t6 * t21 - t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
