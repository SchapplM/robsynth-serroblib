% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:08
% EndTime: 2019-03-09 09:57:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (583->62), mult. (1357->125), div. (0->0), fcn. (898->6), ass. (0->48)
t48 = qJD(1) ^ 2;
t61 = t48 / 0.2e1;
t46 = sin(qJ(2));
t47 = cos(qJ(2));
t28 = (-pkin(2) * t47 - qJ(3) * t46 - pkin(1)) * qJD(1);
t55 = t47 * qJD(1);
t34 = pkin(7) * t55 + qJD(2) * qJ(3);
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t22 = t44 * t28 - t43 * t34;
t56 = qJD(1) * t46;
t31 = t43 * qJD(2) + t44 * t56;
t10 = -pkin(3) * t55 - t31 * pkin(8) + t22;
t23 = t43 * t28 + t44 * t34;
t29 = -t44 * qJD(2) + t43 * t56;
t13 = -t29 * pkin(8) + t23;
t45 = sin(qJ(4));
t60 = cos(qJ(4));
t7 = t45 * t10 + t60 * t13;
t18 = t60 * t29 + t45 * t31;
t20 = -t45 * t29 + t60 * t31;
t59 = t18 * t20;
t36 = -qJD(4) + t55;
t14 = t20 * t36;
t15 = t36 * t18;
t58 = t47 * t48;
t57 = pkin(4) + qJ(6);
t16 = t18 ^ 2 / 0.2e1;
t17 = t20 ^ 2 / 0.2e1;
t54 = qJD(1) * qJD(2);
t53 = t46 * t54;
t52 = t47 * t54;
t5 = t36 * qJ(5) - t7;
t6 = t60 * t10 - t45 * t13;
t33 = -qJD(2) * pkin(2) + pkin(7) * t56 + qJD(3);
t51 = qJD(5) - t6;
t24 = t29 * pkin(3) + t33;
t50 = -t20 * qJ(5) + t24;
t42 = t47 ^ 2;
t41 = t46 ^ 2;
t38 = t42 * t61;
t35 = t36 ^ 2 / 0.2e1;
t8 = t18 * pkin(4) + t50;
t4 = t36 * pkin(4) + t51;
t3 = t57 * t18 + t50;
t2 = -t18 * pkin(5) + qJD(6) - t5;
t1 = t20 * pkin(5) + t57 * t36 + t51;
t9 = [0, 0, 0, 0, 0, t61, 0, 0, 0, 0, t41 * t61, t46 * t58, t53, t38, t52, qJD(2) ^ 2 / 0.2e1, pkin(1) * t58 - pkin(7) * t53, -t48 * pkin(1) * t46 - pkin(7) * t52 (t41 + t42) * t48 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t48, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t55, t29 ^ 2 / 0.2e1, t29 * t55, t38, -t22 * t55 + t33 * t29, t23 * t55 + t33 * t31, -t22 * t31 - t23 * t29, t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t17, -t59, -t14, t16, t15, t35, t24 * t18 - t6 * t36, t24 * t20 + t7 * t36, -t7 * t18 - t6 * t20, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t35, t14, -t15, t17, -t59, t16, t5 * t18 + t4 * t20, -t8 * t18 - t4 * t36, -t8 * t20 + t5 * t36, t8 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t35, -t15, -t14, t16, t59, t17, t1 * t20 - t2 * t18, -t2 * t36 - t3 * t20, t1 * t36 + t3 * t18, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t9;
