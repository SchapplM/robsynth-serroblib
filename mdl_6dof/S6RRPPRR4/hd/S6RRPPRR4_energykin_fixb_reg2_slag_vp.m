% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:01
% EndTime: 2019-03-09 09:05:01
% DurationCPUTime: 0.24s
% Computational Cost: add. (834->70), mult. (2382->149), div. (0->0), fcn. (1816->10), ass. (0->58)
t48 = sin(pkin(11));
t50 = cos(pkin(11));
t54 = sin(qJ(2));
t55 = cos(qJ(2));
t49 = sin(pkin(6));
t68 = qJD(1) * t49;
t35 = (t48 * t55 + t50 * t54) * t68;
t75 = pkin(3) + pkin(9);
t61 = t55 * t68;
t62 = t54 * t68;
t33 = t48 * t62 - t50 * t61;
t39 = qJD(3) + (-pkin(2) * t55 - pkin(1)) * t68;
t57 = -t35 * qJ(4) + t39;
t12 = t75 * t33 + t57;
t53 = sin(qJ(5));
t74 = cos(qJ(5));
t67 = cos(pkin(6)) * qJD(1);
t46 = qJD(2) + t67;
t63 = pkin(1) * t67;
t45 = t55 * t63;
t26 = t46 * pkin(2) + t45 + (-pkin(8) - qJ(3)) * t62;
t38 = pkin(8) * t61 + t54 * t63;
t29 = qJ(3) * t61 + t38;
t15 = t50 * t26 - t48 * t29;
t59 = qJD(4) - t15;
t9 = t35 * pkin(4) - t75 * t46 + t59;
t6 = t74 * t12 + t53 * t9;
t73 = cos(qJ(6));
t72 = t35 * t33;
t71 = t35 * t46;
t70 = t46 * t33;
t56 = qJD(1) ^ 2;
t69 = t49 ^ 2 * t56;
t16 = t48 * t26 + t50 * t29;
t66 = t33 ^ 2 / 0.2e1;
t65 = t35 ^ 2 / 0.2e1;
t64 = t55 * t69;
t14 = -t46 * qJ(4) - t16;
t60 = t69 / 0.2e1;
t23 = -t74 * t33 + t53 * t46;
t10 = -t33 * pkin(4) - t14;
t5 = -t53 * t12 + t74 * t9;
t52 = sin(qJ(6));
t41 = t46 ^ 2 / 0.2e1;
t37 = -pkin(8) * t62 + t45;
t32 = qJD(5) + t35;
t25 = t53 * t33 + t74 * t46;
t22 = qJD(6) + t23;
t20 = t33 * pkin(3) + t57;
t19 = t73 * t25 + t52 * t32;
t17 = t52 * t25 - t73 * t32;
t13 = -t46 * pkin(3) + t59;
t7 = t23 * pkin(5) - t25 * pkin(10) + t10;
t4 = t32 * pkin(10) + t6;
t3 = -t32 * pkin(5) - t5;
t2 = t73 * t4 + t52 * t7;
t1 = -t52 * t4 + t73 * t7;
t8 = [0, 0, 0, 0, 0, t56 / 0.2e1, 0, 0, 0, 0, t54 ^ 2 * t60, t54 * t64, t46 * t62, t55 ^ 2 * t60, t46 * t61, t41, pkin(1) * t64 + t37 * t46, -pkin(1) * t54 * t69 - t38 * t46 (-t37 * t54 + t38 * t55) * t68, t38 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t60, t65, -t72, t71, t66, -t70, t41, t15 * t46 + t39 * t33, -t16 * t46 + t39 * t35, -t15 * t35 - t16 * t33, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t41, -t71, t70, t65, -t72, t66, t13 * t35 + t14 * t33, t13 * t46 - t20 * t33, -t14 * t46 - t20 * t35, t20 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t32, t23 ^ 2 / 0.2e1, -t23 * t32, t32 ^ 2 / 0.2e1, t10 * t23 + t5 * t32, t10 * t25 - t6 * t32, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t22, t17 ^ 2 / 0.2e1, -t17 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t17, t3 * t19 - t2 * t22, -t1 * t19 - t2 * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
