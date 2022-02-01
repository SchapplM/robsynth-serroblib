% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:24
% DurationCPUTime: 0.64s
% Computational Cost: add. (1968->116), mult. (4069->168), div. (0->0), fcn. (2412->8), ass. (0->92)
t60 = cos(pkin(9));
t63 = sin(qJ(2));
t106 = t60 * t63;
t59 = sin(pkin(9));
t66 = cos(qJ(2));
t74 = pkin(1) * (-t59 * t66 - t106);
t38 = qJD(1) * t74;
t107 = t59 * t63;
t73 = pkin(1) * (t60 * t66 - t107);
t40 = qJD(1) * t73;
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t117 = pkin(2) * t59;
t51 = t60 * pkin(2) + pkin(3);
t76 = -t62 * t117 + t65 * t51;
t102 = t76 * qJD(4) - t62 * t38 - t65 * t40;
t56 = qJD(1) + qJD(2);
t54 = qJD(4) + t56;
t120 = t102 * t54;
t100 = t65 * t117 + t62 * t51;
t37 = pkin(8) + t100;
t67 = qJD(5) ^ 2;
t101 = -t100 * qJD(4) - t65 * t38 + t62 * t40;
t91 = t101 * t54;
t119 = t37 * t67 - t91;
t97 = pkin(1) * qJD(1);
t93 = t66 * t97;
t46 = t56 * pkin(2) + t93;
t94 = t63 * t97;
t87 = t59 * t94;
t23 = t60 * t46 - t87;
t21 = t56 * pkin(3) + t23;
t24 = t59 * t46 + t60 * t94;
t16 = t62 * t21 + t65 * t24;
t14 = t54 * pkin(8) + t16;
t61 = sin(qJ(5));
t64 = cos(qJ(5));
t10 = t61 * qJD(3) + t64 * t14;
t105 = t62 * t24;
t39 = qJD(2) * t74;
t30 = qJD(1) * t39;
t31 = (t60 * t93 - t87) * qJD(2);
t5 = (qJD(4) * t21 + t31) * t65 - qJD(4) * t105 + t62 * t30;
t9 = t64 * qJD(3) - t61 * t14;
t2 = t9 * qJD(5) + t64 * t5;
t3 = -t10 * qJD(5) - t61 * t5;
t118 = t2 * t64 - t3 * t61 + (-t10 * t61 - t64 * t9) * qJD(5);
t52 = t66 * pkin(1) + pkin(2);
t85 = -pkin(1) * t107 + t60 * t52;
t35 = pkin(3) + t85;
t42 = pkin(1) * t106 + t59 * t52;
t103 = t62 * t35 + t65 * t42;
t81 = t10 * t64 - t61 * t9;
t6 = t16 * qJD(4) - t65 * t30 + t62 * t31;
t116 = t54 * pkin(4);
t41 = qJD(2) * t73;
t79 = t65 * t35 - t62 * t42;
t7 = t79 * qJD(4) + t62 * t39 + t65 * t41;
t114 = t7 * t54;
t8 = t103 * qJD(4) - t65 * t39 + t62 * t41;
t113 = t8 * t54;
t15 = t65 * t21 - t105;
t13 = -t15 - t116;
t96 = qJD(5) * t64;
t112 = t13 * t96 + t6 * t61;
t110 = t15 * t54;
t109 = t16 * t54;
t104 = t67 * t61;
t57 = t61 ^ 2;
t58 = t64 ^ 2;
t99 = t57 - t58;
t98 = t57 + t58;
t53 = t54 ^ 2;
t95 = t61 * t53 * t64;
t92 = -t13 * t54 - t5;
t86 = t61 * t54 * t96;
t84 = (-qJD(2) + t56) * t97;
t83 = pkin(1) * qJD(2) * (-qJD(1) - t56);
t82 = pkin(8) * t67 - t109;
t18 = pkin(8) + t103;
t80 = t18 * t67 + t113;
t78 = qJD(5) * (t15 - t116);
t17 = -pkin(4) - t79;
t77 = qJD(5) * (t17 * t54 - t7);
t36 = -pkin(4) - t76;
t75 = qJD(5) * (t36 * t54 - t102);
t55 = t67 * t64;
t45 = -0.2e1 * t86;
t44 = 0.2e1 * t86;
t26 = -0.2e1 * t99 * t54 * qJD(5);
t11 = t13 * qJD(5) * t61;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t83, t66 * t83, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t56 + t30, -t41 * t56 - t31, 0, t23 * t39 + t24 * t41 + t30 * t85 + t31 * t42, 0, 0, 0, 0, 0, 0, -t6 - t113, -t5 - t114, 0, t5 * t103 - t15 * t8 + t16 * t7 - t6 * t79, t44, t26, t55, t45, -t104, 0, t11 + t61 * t77 + (-t6 - t80) * t64, t80 * t61 + t64 * t77 + t112, t98 * t114 + t118, t118 * t18 + t13 * t8 + t6 * t17 + t81 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t84, t66 * t84, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t56 + t30, t40 * t56 - t31, 0, -t23 * t38 - t24 * t40 + (t30 * t60 + t31 * t59) * pkin(2), 0, 0, 0, 0, 0, 0, -t6 + t91, -t5 - t120, 0, t5 * t100 + t101 * t15 + t102 * t16 - t6 * t76, t44, t26, t55, t45, -t104, 0, t11 + t61 * t75 + (-t6 - t119) * t64, t119 * t61 + t64 * t75 + t112, t98 * t120 + t118, -t101 * t13 + t81 * t102 + t118 * t37 + t6 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t55, 0, t81 * qJD(5) + t2 * t61 + t3 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 + t109, -t5 + t110, 0, 0, t44, t26, t55, t45, -t104, 0, t11 + t61 * t78 + (-t6 - t82) * t64, t82 * t61 + t64 * t78 + t112, -t98 * t110 + t118, -t6 * pkin(4) + pkin(8) * t118 - t13 * t16 - t81 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t99 * t53, 0, t95, 0, 0, t92 * t61, t92 * t64, 0, 0;];
tauc_reg = t1;
