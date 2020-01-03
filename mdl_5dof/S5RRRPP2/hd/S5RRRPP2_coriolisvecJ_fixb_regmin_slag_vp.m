% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:58
% EndTime: 2019-12-31 20:52:02
% DurationCPUTime: 0.85s
% Computational Cost: add. (1015->191), mult. (1704->236), div. (0->0), fcn. (777->4), ass. (0->130)
t74 = qJD(3) * qJ(4);
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t148 = t77 * qJD(4) + t79 * t74;
t144 = pkin(3) + pkin(4);
t101 = t144 * qJD(3);
t72 = qJD(1) + qJD(2);
t125 = qJ(5) * t72;
t128 = pkin(1) * qJD(1);
t78 = sin(qJ(2));
t108 = t78 * t128;
t47 = t72 * pkin(7) + t108;
t38 = t77 * t47;
t23 = t77 * t125 - t38;
t114 = qJD(4) - t23;
t12 = -t101 + t114;
t122 = qJD(3) * t77;
t127 = pkin(1) * qJD(2);
t103 = qJD(1) * t127;
t80 = cos(qJ(2));
t94 = t80 * t103;
t53 = t79 * t94;
t73 = qJD(3) * qJD(4);
t15 = -t47 * t122 + t53 + t73;
t121 = qJD(3) * t79;
t52 = t77 * t94;
t20 = t47 * t121 + t52;
t147 = t15 * t79 + t20 * t77;
t136 = t72 * t79;
t39 = t79 * t47;
t24 = -t79 * t125 + t39;
t19 = t24 + t74;
t67 = t77 * qJ(4);
t70 = t79 * pkin(3);
t146 = t67 + t70;
t100 = pkin(2) + t67;
t107 = t80 * t128;
t7 = t107 + qJD(5) + (t144 * t79 + t100) * t72;
t123 = qJD(5) + t7;
t145 = t123 * t77;
t137 = t72 * t77;
t95 = t78 * t103;
t85 = t148 * t72 - t95;
t5 = -t101 * t137 + t85;
t143 = t7 * t121 + t5 * t77;
t82 = qJD(3) ^ 2;
t142 = pkin(7) * t82;
t141 = t72 * pkin(2);
t140 = t80 * pkin(1);
t139 = t19 * t79;
t112 = pkin(2) + t146;
t69 = t79 * pkin(4);
t41 = t69 + t112;
t138 = t41 * t72;
t135 = t82 * t77;
t66 = t82 * t79;
t134 = pkin(7) - qJ(5);
t48 = -t107 - t141;
t133 = t48 * t121 + t77 * t95;
t120 = qJD(3) * t80;
t106 = t77 * t120;
t98 = t72 * t108;
t132 = t106 * t128 + t79 * t98;
t131 = t53 + 0.2e1 * t73;
t75 = t77 ^ 2;
t76 = t79 ^ 2;
t130 = t75 - t76;
t129 = t75 + t76;
t126 = qJ(4) * t79;
t62 = t78 * pkin(1) + pkin(7);
t124 = -qJ(5) + t62;
t119 = t19 * qJD(3);
t30 = t39 + t74;
t118 = t30 * qJD(3);
t117 = t77 * qJD(5);
t116 = t79 * qJD(5);
t115 = -qJD(1) - t72;
t113 = qJ(5) * qJD(3);
t71 = t72 ^ 2;
t111 = t77 * t71 * t79;
t33 = pkin(3) * t122 - t148;
t110 = t80 * t127;
t109 = t78 * t127;
t105 = t72 * t122;
t104 = t72 * t121;
t63 = -pkin(2) - t140;
t60 = t77 * t113;
t102 = t79 * t113;
t56 = t134 * t79;
t43 = t124 * t79;
t99 = 0.2e1 * t104;
t97 = -qJD(3) * pkin(3) + qJD(4);
t96 = t129 * t140;
t93 = (-qJD(2) + t72) * t128;
t92 = t115 * t127;
t91 = -qJD(5) + t110;
t27 = t38 + t97;
t90 = t27 * t77 + t30 * t79;
t40 = t63 - t146;
t8 = pkin(3) * t105 - t85;
t89 = -t33 * t72 - t142 - t8;
t26 = t109 + t33;
t88 = -t26 * t72 - t62 * t82 - t8;
t87 = t40 * t72 - t110;
t25 = -pkin(4) * t122 - t33;
t86 = -t77 * t118 + t27 * t121 + t147;
t57 = t72 * t60;
t1 = -t72 * t116 + t15 + t57;
t4 = (-t102 - t117) * t72 + t20;
t84 = (-qJD(3) * t12 - t1) * t79 + (-t4 + t119) * t77;
t83 = (t27 * t79 - t30 * t77) * qJD(3) + t147;
t55 = t134 * t77;
t50 = -t75 * t71 - t82;
t46 = t77 * t99;
t44 = t77 * t98;
t42 = t124 * t77;
t36 = t48 * t122;
t34 = (pkin(3) * t77 - t126) * t72;
t32 = qJD(3) * t56 - t117;
t31 = -pkin(7) * t122 - t116 + t60;
t29 = -t40 + t69;
t28 = -0.2e1 * t130 * t72 * qJD(3);
t22 = (-t144 * t77 + t126) * t72;
t18 = -t107 + (-t100 - t70) * t72;
t17 = t25 - t109;
t14 = qJD(3) * t43 + t91 * t77;
t13 = -t62 * t122 + t91 * t79 + t60;
t10 = t18 * t122;
t3 = t5 * t79;
t2 = [0, 0, 0, 0, t78 * t92, t80 * t92, t46, t28, t66, -t135, 0, t63 * t105 - t62 * t66 + t36 + (t115 * t79 * t78 - t106) * t127, t63 * t104 + t62 * t135 + (-t79 * t120 + t78 * t137) * t127 + t133, t87 * t122 + t88 * t79 + t10, t72 * qJD(2) * t96 + t86, t88 * t77 + (-t18 - t87) * t121, t90 * t110 + t18 * t26 + t8 * t40 + t83 * t62, t17 * t136 + t3 + (-t14 + (-t29 * t72 - t7) * t77) * qJD(3), t17 * t137 + (t29 * t136 + t13) * qJD(3) + t143, (-t13 * t79 - t14 * t77 + (-t42 * t79 + t43 * t77) * qJD(3)) * t72 + t84, t1 * t43 + t12 * t14 + t19 * t13 + t7 * t17 + t5 * t29 + t4 * t42; 0, 0, 0, 0, t78 * t93, t80 * t93, t46, t28, t66, -t135, 0, -pkin(2) * t105 + t36 + (-t95 - t142) * t79 + t132, pkin(7) * t135 - t44 + (t107 - t141) * t121 + t133, -t105 * t112 + t89 * t79 + t10 + t132, -t129 * t72 * t107 + t86, t44 + t89 * t77 + (t112 * t72 - t107 - t18) * t121, t18 * t33 - t8 * t112 + (-t18 * t78 - t90 * t80) * t128 + t83 * pkin(7), t25 * t136 + t3 + (-t32 + (-t7 - t138) * t77) * qJD(3) + t132, t25 * t137 + t44 + (t31 + (-t107 + t138) * t79) * qJD(3) + t143, (-t31 * t79 - t32 * t77 + (-t55 * t79 + t56 * t77) * qJD(3) + qJD(1) * t96) * t72 + t84, t1 * t56 + t12 * t32 + t19 * t31 + t7 * t25 + t4 * t55 + t5 * t41 + (t7 * t78 + (-t12 * t77 - t139) * t80) * t128; 0, 0, 0, 0, 0, 0, -t111, t130 * t71, 0, 0, 0, -t48 * t137 - t52, -t48 * t136 - t53, -t52 + (-t18 * t77 + t34 * t79) * t72, ((t30 - t74) * t77 + (-t27 + t97) * t79) * t72, (t18 * t79 + t34 * t77) * t72 + t131, -t27 * t39 - t20 * pkin(3) + t15 * qJ(4) - t18 * t34 + (qJD(4) + t38) * t30, t24 * qJD(3) + ((-t22 + t113) * t79 + t145) * t72 - t20, t57 + (-t23 - t38) * qJD(3) + (-t123 * t79 - t22 * t77) * t72 + t131, 0, t1 * qJ(4) + t114 * t19 - t12 * t24 - t144 * t4 - t7 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, 0, t50, t18 * t137 - t118 + t20, -t111, t50, 0, -t119 + (-t102 - t145) * t72 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t105, t99, -t129 * t71, (t139 + (t12 - t101) * t77) * t72 + t85;];
tauc_reg = t2;
