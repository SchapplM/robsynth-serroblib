% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:38
% EndTime: 2021-01-15 10:56:42
% DurationCPUTime: 0.76s
% Computational Cost: add. (872->160), mult. (2397->241), div. (0->0), fcn. (1655->6), ass. (0->99)
t102 = cos(pkin(7));
t68 = cos(qJ(2));
t87 = t102 * t68;
t54 = qJD(1) * t87;
t64 = sin(pkin(7));
t66 = sin(qJ(2));
t99 = qJD(1) * t66;
t37 = t64 * t99 - t54;
t35 = qJD(4) + t37;
t122 = qJD(4) - t35;
t95 = qJD(1) * qJD(2);
t121 = -0.2e1 * t95;
t88 = t102 * t66;
t47 = t64 * t68 + t88;
t101 = qJD(1) * t47;
t105 = -qJ(3) - pkin(5);
t86 = qJD(2) * t105;
t36 = t68 * qJD(3) + t66 * t86;
t30 = t36 * qJD(1);
t74 = t66 * qJD(3) - t68 * t86;
t71 = t102 * t74;
t4 = qJD(1) * t71 + t64 * t30;
t57 = t64 * pkin(2) + pkin(6);
t93 = pkin(2) * t99;
t120 = (pkin(3) * t101 + t37 * pkin(6) + qJD(4) * t57 + t93) * t35 + t4;
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t26 = t65 * qJD(2) + t101 * t67;
t91 = t66 * t95;
t53 = t64 * t91;
t33 = qJD(2) * t54 - t53;
t7 = qJD(4) * t26 + t65 * t33;
t72 = t64 * t74;
t11 = t102 * t36 - t72;
t52 = t105 * t68;
t50 = qJD(1) * t52;
t111 = t64 * t50;
t103 = qJD(2) * pkin(2);
t92 = t105 * t66;
t49 = qJD(1) * t92;
t45 = t49 + t103;
t18 = t102 * t45 + t111;
t14 = -qJD(2) * pkin(3) - t18;
t60 = -t68 * pkin(2) - pkin(1);
t76 = -t64 * t66 + t87;
t17 = -pkin(3) * t76 - t47 * pkin(6) + t60;
t42 = t76 * qJD(2);
t89 = t102 * t30;
t5 = -qJD(1) * t72 + t89;
t23 = -t102 * t52 + t64 * t92;
t39 = t47 * qJD(2);
t32 = qJD(1) * t39;
t79 = -t23 * t32 + t4 * t47;
t100 = qJD(1) * t60;
t51 = qJD(3) + t100;
t9 = t37 * pkin(3) - pkin(6) * t101 + t51;
t119 = t14 * t42 - (qJD(4) * t17 + t11) * t35 + (qJD(4) * t9 + t5) * t76 + t79;
t118 = pkin(2) * t66;
t96 = t67 * qJD(2);
t97 = qJD(4) * t65;
t6 = qJD(4) * t96 - t101 * t97 + t67 * t33;
t117 = t6 * t65;
t116 = t17 * t32;
t24 = t101 * t65 - t96;
t115 = t24 * t35;
t114 = t26 * t35;
t113 = t26 * t101;
t112 = t101 * t24;
t110 = t65 * t32;
t28 = t67 * t32;
t70 = qJD(1) ^ 2;
t108 = t68 * t70;
t69 = qJD(2) ^ 2;
t107 = t69 * t66;
t106 = t69 * t68;
t43 = t102 * t50;
t19 = t64 * t45 - t43;
t104 = t66 ^ 2 - t68 ^ 2;
t98 = qJD(4) * t47;
t94 = t66 * t103;
t85 = t35 * t67;
t84 = 0.2e1 * t101;
t81 = pkin(1) * t121;
t15 = qJD(2) * pkin(6) + t19;
t2 = t67 * t15 + t65 * t9;
t80 = t65 * t15 - t67 * t9;
t78 = t28 + (-t37 * t65 - t97) * t35;
t77 = t67 * t42 - t47 * t97;
t21 = t102 * t49 + t111;
t73 = -t57 * t32 + (t14 + t21) * t35;
t58 = -t102 * pkin(2) - pkin(3);
t55 = pkin(2) * t91;
t22 = -t105 * t88 - t64 * t52;
t20 = t64 * t49 - t43;
t13 = t39 * pkin(3) - t42 * pkin(6) + t94;
t10 = t64 * t36 + t71;
t8 = t32 * pkin(3) - t33 * pkin(6) + t55;
t3 = t67 * t8;
t1 = [0, 0, 0, 0.2e1 * t68 * t91, t104 * t121, t106, -t107, 0, -pkin(5) * t106 + t66 * t81, pkin(5) * t107 + t68 * t81, t60 * t32 + t51 * t39 + (-t10 + (-qJD(1) * t76 + t37) * t118) * qJD(2), t60 * t33 + t51 * t42 + (t84 * t118 - t11) * qJD(2), t10 * t101 - t11 * t37 - t18 * t42 - t19 * t39 + t22 * t33 + t5 * t76 + t79, -t18 * t10 + t19 * t11 + t4 * t22 + t5 * t23 + (t51 + t100) * t94, t6 * t67 * t47 + t77 * t26, (-t24 * t67 - t26 * t65) * t42 + (-t117 - t67 * t7 + (t24 * t65 - t26 * t67) * qJD(4)) * t47, t26 * t39 + t47 * t28 + t77 * t35 - t6 * t76, -t47 * t110 - t24 * t39 + t7 * t76 + (-t65 * t42 - t67 * t98) * t35, -t32 * t76 + t35 * t39, -t80 * t39 + t10 * t24 + t22 * t7 - t3 * t76 + (t13 * t35 + t116 + (t14 * t47 + t15 * t76 - t23 * t35) * qJD(4)) * t67 + t119 * t65, t10 * t26 - t2 * t39 + t22 * t6 + (-(-qJD(4) * t23 + t13) * t35 - t116 + (-qJD(4) * t15 + t8) * t76 - t14 * t98) * t65 + t119 * t67; 0, 0, 0, -t66 * t108, t104 * t70, 0, 0, 0, t70 * pkin(1) * t66, pkin(1) * t108, t20 * qJD(2) - t101 * t51 - t37 * t93 - t4, -t89 + t21 * qJD(2) + t51 * t37 + (-t101 * t118 + t72) * qJD(1), (t19 - t20) * t101 + (-t18 + t21) * t37 + (-t102 * t33 - t32 * t64) * pkin(2), t18 * t20 - t19 * t21 + (-t102 * t4 + t5 * t64 - t51 * t99) * pkin(2), t26 * t85 + t117, (t6 - t115) * t67 + (-t7 - t114) * t65, t35 * t85 + t110 - t113, t78 + t112, -t35 * t101, t101 * t80 - t120 * t67 - t20 * t24 + t58 * t7 + t73 * t65, t101 * t2 + t120 * t65 - t20 * t26 + t58 * t6 + t73 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * qJD(2), -t53 + (t54 - t37) * qJD(2), -t101 ^ 2 - t37 ^ 2, t101 * t18 + t19 * t37 + t55, 0, 0, 0, 0, 0, t78 - t112, -t35 ^ 2 * t67 - t110 - t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t24, -t24 ^ 2 + t26 ^ 2, t6 + t115, t114 - t7, t32, -t122 * t2 - t14 * t26 - t65 * t5 + t3, t122 * t80 + t14 * t24 - t67 * t5 - t65 * t8;];
tauc_reg = t1;
