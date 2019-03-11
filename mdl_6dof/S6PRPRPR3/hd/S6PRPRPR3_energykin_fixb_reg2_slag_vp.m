% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:38
% EndTime: 2019-03-08 19:37:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (257->50), mult. (613->113), div. (0->0), fcn. (395->10), ass. (0->48)
t43 = qJD(2) ^ 2;
t32 = t43 / 0.2e1;
t55 = -pkin(4) - pkin(9);
t42 = cos(qJ(2));
t34 = sin(pkin(6));
t54 = qJD(1) * t34;
t21 = qJD(2) * pkin(2) + t42 * t54;
t33 = sin(pkin(11));
t35 = cos(pkin(11));
t39 = sin(qJ(2));
t50 = t39 * t54;
t16 = t33 * t21 + t35 * t50;
t14 = qJD(2) * pkin(8) + t16;
t36 = cos(pkin(6));
t26 = t36 * qJD(1) + qJD(3);
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t9 = t41 * t14 + t38 * t26;
t53 = qJD(2) * t41;
t52 = t38 * qJD(2);
t51 = qJD(2) * qJD(4);
t49 = qJD(2) * t54;
t48 = t38 * t51;
t47 = t41 * t51;
t46 = -qJ(5) * t38 - pkin(3);
t8 = -t38 * t14 + t41 * t26;
t15 = t35 * t21 - t33 * t50;
t45 = qJD(5) - t8;
t6 = -qJD(4) * qJ(5) - t9;
t44 = qJD(1) ^ 2;
t40 = cos(qJ(6));
t37 = sin(qJ(6));
t31 = qJD(4) ^ 2 / 0.2e1;
t29 = t41 ^ 2 * t32;
t28 = t38 ^ 2 * t32;
t27 = qJD(6) + t52;
t25 = t38 * t43 * t41;
t20 = t40 * qJD(4) - t37 * t53;
t18 = t37 * qJD(4) + t40 * t53;
t13 = -qJD(2) * pkin(3) - t15;
t10 = (-pkin(4) * t41 + t46) * qJD(2) - t15;
t7 = (t55 * t41 + t46) * qJD(2) - t15;
t5 = -qJD(4) * pkin(4) + t45;
t4 = pkin(5) * t53 - t6;
t3 = pkin(5) * t52 + t55 * qJD(4) + t45;
t2 = t37 * t3 + t40 * t7;
t1 = t40 * t3 - t37 * t7;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t44 / 0.2e1, 0, 0, 0, 0, 0, t32, t42 * t49, -t39 * t49, 0 (t36 ^ 2 / 0.2e1 + (t39 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * t34 ^ 2) * t44, 0, 0, 0, 0, 0, t32, t15 * qJD(2), -t16 * qJD(2), 0, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t28, t25, t48, t29, t47, t31, t8 * qJD(4) - t13 * t53, -t9 * qJD(4) + t13 * t52 (-t38 * t8 + t41 * t9) * qJD(2), t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t31, -t48, -t47, t28, t25, t29 (t38 * t5 - t41 * t6) * qJD(2), t5 * qJD(4) + t10 * t53, -t6 * qJD(4) - t10 * t52, t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t27, t18 ^ 2 / 0.2e1, -t18 * t27, t27 ^ 2 / 0.2e1, t1 * t27 + t4 * t18, -t2 * t27 + t4 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
