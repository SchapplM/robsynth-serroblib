% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:40
% DurationCPUTime: 0.90s
% Computational Cost: add. (651->164), mult. (1518->240), div. (0->0), fcn. (850->6), ass. (0->107)
t68 = sin(qJ(3));
t108 = t68 * qJD(1);
t52 = qJD(5) + t108;
t134 = qJD(5) - t52;
t53 = sin(pkin(8)) * pkin(1) + pkin(6);
t44 = t53 * qJD(1);
t70 = cos(qJ(3));
t117 = t70 * qJD(2) - t68 * t44;
t133 = -qJD(4) + t117;
t17 = -qJD(3) * pkin(3) - t133;
t69 = cos(qJ(5));
t106 = t69 * qJD(3);
t114 = qJD(1) * t70;
t67 = sin(qJ(5));
t42 = -t67 * t114 + t106;
t104 = qJD(3) * qJ(4);
t25 = t68 * qJD(2) + t70 * t44;
t19 = -t104 - t25;
t99 = pkin(4) * t108;
t105 = t99 - t133;
t109 = t67 * qJD(3);
t64 = t70 ^ 2;
t115 = qJD(1) * t64;
t123 = t52 * t68;
t110 = qJD(5) * t69;
t95 = t70 * t110;
t132 = -(t115 - t123) * t109 - t52 * t95;
t71 = -pkin(3) - pkin(7);
t102 = qJD(2) * qJD(3);
t113 = qJD(3) * t68;
t119 = -t70 * t102 + t44 * t113;
t3 = (qJD(4) - t99) * qJD(3) - t119;
t131 = t3 * t67;
t130 = t3 * t69;
t129 = pkin(4) + t53;
t112 = qJD(3) * t70;
t40 = t69 * t114 + t109;
t103 = qJD(1) * qJD(3);
t92 = t68 * t103;
t13 = -qJD(5) * t40 + t67 * t92;
t128 = t42 * t112 + t13 * t68;
t127 = t13 * t69;
t126 = t40 * t52;
t125 = t40 * t70;
t124 = t42 * t52;
t122 = t52 * t71;
t14 = t42 * qJD(5) - t69 * t92;
t121 = t68 * t14;
t72 = qJD(3) ^ 2;
t120 = t72 * t68;
t62 = t72 * t70;
t22 = t68 * t102 + t44 * t112;
t101 = t69 * t123;
t111 = qJD(5) * t67;
t97 = t52 * t111;
t118 = qJD(3) * t101 + t70 * t97;
t63 = t68 ^ 2;
t116 = t63 - t64;
t54 = -cos(pkin(8)) * pkin(1) - pkin(2);
t79 = -t68 * qJ(4) + t54;
t31 = -t70 * pkin(3) + t79;
t20 = qJD(1) * t31;
t45 = qJD(1) * t54;
t107 = t68 * qJD(4);
t73 = qJD(1) ^ 2;
t100 = t68 * t73 * t70;
t98 = t69 * t115;
t96 = t52 * t110;
t51 = pkin(3) * t92;
t86 = pkin(7) * t68 - qJ(4) * t70;
t75 = t86 * qJD(3) - t107;
t6 = t75 * qJD(1) + t51;
t91 = t70 * t103;
t7 = pkin(4) * t91 + t22;
t93 = -t67 * t6 + t69 * t7;
t33 = t129 * t70;
t88 = t25 * qJD(3) - t22;
t23 = t71 * t70 + t79;
t10 = t23 * qJD(1);
t5 = t71 * qJD(3) + t105;
t2 = t69 * t10 + t67 * t5;
t87 = t67 * t10 - t69 * t5;
t32 = t129 * t68;
t85 = t69 * t23 + t67 * t32;
t82 = -0.2e1 * qJD(3) * t20;
t81 = 0.2e1 * qJD(3) * t45;
t80 = t52 * t67;
t57 = pkin(4) * t114;
t9 = -t19 + t57;
t78 = t71 * t112 + t68 * t9;
t76 = -t70 * t104 - t107;
t21 = t76 * qJD(1) + t51;
t58 = pkin(3) * t113;
t30 = t58 + t76;
t77 = qJD(1) * t30 + t53 * t72 + t21;
t12 = -qJD(3) * qJD(4) + t119;
t74 = -t12 * t70 + t22 * t68 + (t17 * t70 + t19 * t68) * qJD(3);
t59 = pkin(3) * t108;
t47 = t69 * t91;
t43 = -qJ(4) * t114 + t59;
t29 = qJD(3) * t33;
t28 = t129 * t113;
t26 = t86 * qJD(1) + t59;
t18 = t58 + t75;
t16 = t57 + t25;
t11 = t20 * t108;
t1 = [0, 0, 0, 0, 0.2e1 * t68 * t91, -0.2e1 * t116 * t103, t62, -t120, 0, -t53 * t62 + t68 * t81, t53 * t120 + t70 * t81, t74, t68 * t82 + t77 * t70, -t77 * t68 + t70 * t82, t20 * t30 + t21 * t31 + t74 * t53, -t13 * t67 * t70 + (t68 * t109 - t95) * t42, (-t40 * t67 + t42 * t69) * t113 + (-t127 + t14 * t67 + (t40 * t69 + t42 * t67) * qJD(5)) * t70, t128 + t132, -t121 + (-t98 - t125) * qJD(3) + t118, (t52 + t108) * t112, (-t67 * t18 + t69 * t29) * t52 - t28 * t40 + t33 * t14 + (-t9 * t106 + t93) * t68 + (-t2 * t68 - t85 * t52) * qJD(5) + (-t9 * t111 + t130 + ((-t67 * t23 + t69 * t32) * qJD(1) - t87) * qJD(3)) * t70, t33 * t13 - t28 * t42 + (-(qJD(5) * t32 + t18) * t52 - (qJD(5) * t5 + t6) * t68) * t69 + (-(-qJD(5) * t23 + t29) * t52 + (t9 * qJD(3) + qJD(5) * t10 - t7) * t68) * t67 + (-t9 * t110 - t131 + (-t85 * qJD(1) - t2) * qJD(3)) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t62, 0, t120, t62, -t12 * t68 - t22 * t70 + (t17 * t68 - t19 * t70) * qJD(3), 0, 0, 0, 0, 0, t121 + (-t98 + t125) * qJD(3) + t118, t128 - t132; 0, 0, 0, 0, -t100, t116 * t73, 0, 0, 0, -t45 * t108 + t88, qJD(3) * t117 - t45 * t114 + t119, 0, -t43 * t114 + t11 - t88, (0.2e1 * qJD(4) - t117) * qJD(3) + (t20 * t70 + t43 * t68) * qJD(1) - t119, -t22 * pkin(3) - t12 * qJ(4) + t133 * t19 - t17 * t25 - t20 * t43, -t42 * t80 + t127, (-t14 - t124) * t69 + (-t13 + t126) * t67, -t97 + t47 + (-t67 * t123 - t42 * t70) * qJD(1), -t96 + (-t101 + (t40 - t109) * t70) * qJD(1), -t52 * t114, qJ(4) * t14 + t131 - (t69 * t16 - t67 * t26) * t52 + t105 * t40 + (-t67 * t122 + t9 * t69) * qJD(5) + (t78 * t69 + t70 * t87) * qJD(1), qJ(4) * t13 + t130 + (t67 * t16 + t69 * t26) * t52 + t105 * t42 + (-t69 * t122 - t9 * t67) * qJD(5) + (t2 * t70 - t78 * t67) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t63 * t73 - t72, t19 * qJD(3) + t11 + t22, 0, 0, 0, 0, 0, -qJD(3) * t40 - t52 * t80 + t47, -t96 - qJD(3) * t42 + (-t70 * t109 - t101) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t40, -t40 ^ 2 + t42 ^ 2, t13 + t126, -t14 + t124, t91, -t134 * t2 - t9 * t42 + t93, t134 * t87 + t9 * t40 - t69 * t6 - t67 * t7;];
tauc_reg = t1;
