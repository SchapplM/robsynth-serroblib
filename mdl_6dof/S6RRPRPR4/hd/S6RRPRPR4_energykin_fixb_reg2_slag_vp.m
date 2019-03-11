% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:16
% EndTime: 2019-03-09 10:26:16
% DurationCPUTime: 0.25s
% Computational Cost: add. (1397->73), mult. (3882->165), div. (0->0), fcn. (3100->12), ass. (0->57)
t64 = cos(qJ(2));
t71 = cos(pkin(6)) * qJD(1);
t69 = pkin(1) * t71;
t53 = t64 * t69;
t54 = qJD(2) + t71;
t63 = sin(qJ(2));
t58 = sin(pkin(6));
t72 = qJD(1) * t58;
t68 = t63 * t72;
t35 = t54 * pkin(2) + t53 + (-pkin(8) - qJ(3)) * t68;
t67 = t64 * t72;
t46 = pkin(8) * t67 + t63 * t69;
t38 = qJ(3) * t67 + t46;
t57 = sin(pkin(11));
t59 = cos(pkin(11));
t26 = t35 * t57 + t38 * t59;
t24 = pkin(9) * t54 + t26;
t42 = t57 * t68 - t59 * t67;
t44 = (t57 * t64 + t59 * t63) * t72;
t47 = qJD(3) + (-pkin(2) * t64 - pkin(1)) * t72;
t29 = t42 * pkin(3) - t44 * pkin(9) + t47;
t62 = sin(qJ(4));
t76 = cos(qJ(4));
t13 = t24 * t76 + t29 * t62;
t32 = t44 * t62 - t54 * t76;
t11 = -qJ(5) * t32 + t13;
t56 = sin(pkin(12));
t73 = cos(pkin(12));
t12 = -t24 * t62 + t29 * t76;
t34 = t44 * t76 + t54 * t62;
t41 = qJD(4) + t42;
t9 = pkin(4) * t41 - qJ(5) * t34 + t12;
t6 = t11 * t73 + t56 * t9;
t75 = cos(qJ(6));
t65 = qJD(1) ^ 2;
t74 = t58 ^ 2 * t65;
t70 = t64 * t74;
t66 = t74 / 0.2e1;
t20 = t32 * t73 + t34 * t56;
t25 = t35 * t59 - t38 * t57;
t23 = -pkin(3) * t54 - t25;
t5 = -t11 * t56 + t73 * t9;
t14 = pkin(4) * t32 + qJD(5) + t23;
t61 = sin(qJ(6));
t50 = t54 ^ 2 / 0.2e1;
t45 = -pkin(8) * t68 + t53;
t40 = t41 ^ 2 / 0.2e1;
t22 = -t32 * t56 + t34 * t73;
t18 = qJD(6) + t20;
t17 = t22 * t75 + t41 * t61;
t15 = t22 * t61 - t41 * t75;
t7 = pkin(5) * t20 - pkin(10) * t22 + t14;
t4 = pkin(10) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t75 + t61 * t7;
t1 = -t4 * t61 + t7 * t75;
t8 = [0, 0, 0, 0, 0, t65 / 0.2e1, 0, 0, 0, 0, t63 ^ 2 * t66, t63 * t70, t54 * t68, t64 ^ 2 * t66, t54 * t67, t50, pkin(1) * t70 + t45 * t54, -pkin(1) * t63 * t74 - t46 * t54 (-t45 * t63 + t46 * t64) * t72, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t66, t44 ^ 2 / 0.2e1, -t44 * t42, t44 * t54, t42 ^ 2 / 0.2e1, -t42 * t54, t50, t25 * t54 + t42 * t47, -t26 * t54 + t44 * t47, -t25 * t44 - t26 * t42, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t40, t12 * t41 + t23 * t32, -t13 * t41 + t23 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t41, t20 ^ 2 / 0.2e1, -t20 * t41, t40, t14 * t20 + t41 * t5, t14 * t22 - t41 * t6, -t20 * t6 - t22 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t18, t15 ^ 2 / 0.2e1, -t15 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t15 * t3, t17 * t3 - t18 * t2, -t1 * t17 - t15 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
