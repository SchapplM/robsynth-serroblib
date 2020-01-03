% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:59
% EndTime: 2019-12-31 20:52:02
% DurationCPUTime: 1.10s
% Computational Cost: add. (1089->210), mult. (1938->250), div. (0->0), fcn. (872->4), ass. (0->136)
t87 = sin(qJ(3));
t85 = t87 ^ 2;
t89 = cos(qJ(3));
t86 = t89 ^ 2;
t140 = t85 + t86;
t82 = qJD(1) + qJD(2);
t162 = t140 * t82;
t84 = qJD(3) * qJ(4);
t161 = t87 * qJD(4) + t89 * t84;
t156 = pkin(3) + pkin(4);
t113 = t156 * qJD(3);
t136 = qJ(5) * t82;
t139 = pkin(1) * qJD(1);
t88 = sin(qJ(2));
t120 = t88 * t139;
t53 = t82 * pkin(7) + t120;
t40 = t87 * t53;
t23 = t87 * t136 - t40;
t125 = qJD(4) - t23;
t12 = -t113 + t125;
t133 = qJD(3) * t87;
t138 = pkin(1) * qJD(2);
t115 = qJD(1) * t138;
t90 = cos(qJ(2));
t106 = t90 * t115;
t61 = t89 * t106;
t83 = qJD(3) * qJD(4);
t15 = -t53 * t133 + t61 + t83;
t132 = qJD(3) * t89;
t59 = t87 * t106;
t20 = t53 * t132 + t59;
t160 = t15 * t89 + t20 * t87;
t148 = t82 * t89;
t41 = t89 * t53;
t24 = -t89 * t136 + t41;
t19 = t24 + t84;
t77 = t87 * qJ(4);
t80 = t89 * pkin(3);
t159 = t77 + t80;
t112 = pkin(2) + t77;
t119 = t90 * t139;
t7 = t119 + qJD(5) + (t156 * t89 + t112) * t82;
t134 = qJD(5) + t7;
t158 = t134 * t87;
t141 = -t85 + t86;
t157 = 0.2e1 * t141 * t82 * qJD(3);
t149 = t82 * t87;
t105 = t88 * t115;
t95 = t161 * t82 - t105;
t5 = -t113 * t149 + t95;
t155 = t7 * t132 + t5 * t87;
t92 = qJD(3) ^ 2;
t154 = pkin(7) * t92;
t153 = t82 * pkin(2);
t152 = t90 * pkin(1);
t151 = t19 * t89;
t123 = pkin(2) + t159;
t79 = t89 * pkin(4);
t43 = t79 + t123;
t150 = t43 * t82;
t75 = t92 * t87;
t76 = t92 * t89;
t147 = pkin(7) - qJ(5);
t54 = -t119 - t153;
t146 = t87 * t105 + t54 * t132;
t122 = t90 * t138;
t145 = t122 * t162;
t110 = t82 * t120;
t131 = qJD(3) * t90;
t118 = t87 * t131;
t144 = t89 * t110 + t118 * t139;
t143 = t140 * t106;
t142 = t61 + 0.2e1 * t83;
t137 = qJ(4) * t89;
t71 = t88 * pkin(1) + pkin(7);
t135 = -qJ(5) + t71;
t130 = t19 * qJD(3);
t31 = t41 + t84;
t129 = t31 * qJD(3);
t128 = t87 * qJD(5);
t127 = t89 * qJD(5);
t126 = -qJD(1) - t82;
t124 = qJ(5) * qJD(3);
t35 = pkin(3) * t133 - t161;
t121 = t88 * t138;
t117 = t82 * t133;
t116 = t82 * t132;
t72 = -pkin(2) - t152;
t69 = t87 * t124;
t114 = t89 * t124;
t64 = t147 * t89;
t46 = t135 * t89;
t111 = qJD(1) * t140;
t108 = -qJD(3) * pkin(3) + qJD(4);
t107 = t87 * t116;
t104 = (-qJD(2) + t82) * t139;
t103 = t126 * t138;
t102 = -qJD(5) + t122;
t27 = t108 + t40;
t101 = t27 * t87 + t31 * t89;
t42 = t72 - t159;
t8 = pkin(3) * t117 - t95;
t100 = -t35 * t82 - t154 - t8;
t26 = t121 + t35;
t99 = -t26 * t82 - t71 * t92 - t8;
t98 = t42 * t82 - t122;
t97 = t119 * t162;
t25 = -pkin(4) * t133 - t35;
t96 = -t87 * t129 + t27 * t132 + t160;
t65 = t82 * t69;
t1 = -t82 * t127 + t15 + t65;
t4 = (-t114 - t128) * t82 + t20;
t94 = (-qJD(3) * t12 - t1) * t89 + (-t4 + t130) * t87;
t93 = (t27 * t89 - t31 * t87) * qJD(3) + t160;
t81 = t82 ^ 2;
t67 = t87 * t81 * t89;
t63 = t147 * t87;
t56 = -t85 * t81 - t92;
t52 = -0.2e1 * t107;
t51 = 0.2e1 * t107;
t48 = t87 * t110;
t45 = t141 * t81;
t44 = t135 * t87;
t38 = t54 * t133;
t36 = (pkin(3) * t87 - t137) * t82;
t34 = qJD(3) * t64 - t128;
t33 = -pkin(7) * t133 - t127 + t69;
t30 = -t42 + t79;
t22 = (-t156 * t87 + t137) * t82;
t18 = -t119 + (-t112 - t80) * t82;
t17 = t25 - t121;
t14 = qJD(3) * t46 + t102 * t87;
t13 = t102 * t89 - t71 * t133 + t69;
t10 = t18 * t133;
t3 = t5 * t89;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t103, t90 * t103, 0, 0, t51, t157, t76, t52, -t75, 0, t72 * t117 - t71 * t76 + t38 + (t126 * t89 * t88 - t118) * t138, t72 * t116 + t71 * t75 + (-t89 * t131 + t88 * t149) * t138 + t146, t143 + t145, ((qJD(1) * t72 + t54) * t88 + (t71 * t111 + t140 * t53) * t90) * t138, t51, t76, -t157, 0, t75, t52, t98 * t133 + t99 * t89 + t10, t96 + t145, t99 * t87 + (-t18 - t98) * t132, t101 * t122 + t18 * t26 + t8 * t42 + t93 * t71, t51, -t157, -t76, t52, -t75, 0, t17 * t148 + t3 + (-t14 + (-t30 * t82 - t7) * t87) * qJD(3), t17 * t149 + (t30 * t148 + t13) * qJD(3) + t155, (-t13 * t89 - t14 * t87 + (-t44 * t89 + t46 * t87) * qJD(3)) * t82 + t94, t1 * t46 + t12 * t14 + t19 * t13 + t7 * t17 + t5 * t30 + t4 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t104, t90 * t104, 0, 0, t51, t157, t76, t52, -t75, 0, -pkin(2) * t117 + t38 + (-t105 - t154) * t89 + t144, pkin(7) * t75 - t48 + (t119 - t153) * t132 + t146, -t97 + t143, ((-pkin(2) * qJD(2) - t54) * t88 + (pkin(7) * qJD(2) - t53) * t90 * t140) * t139, t51, t76, -t157, 0, t75, t52, t100 * t89 - t117 * t123 + t10 + t144, -t97 + t96, t48 + t100 * t87 + (t123 * t82 - t119 - t18) * t132, t18 * t35 - t8 * t123 + (-t101 * t90 - t18 * t88) * t139 + t93 * pkin(7), t51, -t157, -t76, t52, -t75, 0, t25 * t148 + t3 + (-t34 + (-t7 - t150) * t87) * qJD(3) + t144, t25 * t149 + t48 + (t33 + (-t119 + t150) * t89) * qJD(3) + t155, (-t33 * t89 - t34 * t87 + (-t63 * t89 + t64 * t87) * qJD(3) + t111 * t152) * t82 + t94, t1 * t64 + t12 * t34 + t19 * t33 + t7 * t25 + t4 * t63 + t5 * t43 + (t7 * t88 + (-t12 * t87 - t151) * t90) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t45, 0, t67, 0, 0, -t54 * t149 - t59, -t54 * t148 - t61, 0, 0, -t67, 0, t45, 0, 0, t67, -t59 + (-t18 * t87 + t36 * t89) * t82, ((t31 - t84) * t87 + (t108 - t27) * t89) * t82, (t18 * t89 + t36 * t87) * t82 + t142, -t27 * t41 - t20 * pkin(3) + t15 * qJ(4) - t18 * t36 + (qJD(4) + t40) * t31, -t67, t45, 0, t67, 0, 0, t24 * qJD(3) + ((-t22 + t124) * t89 + t158) * t82 - t20, t65 + (-t23 - t40) * qJD(3) + (-t134 * t89 - t22 * t87) * t82 + t142, 0, t1 * qJ(4) - t12 * t24 + t125 * t19 - t156 * t4 - t7 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, t56, t18 * t149 - t129 + t20, 0, 0, 0, 0, 0, 0, -t67, t56, 0, -t130 + (-t114 - t158) * t82 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t117, 0.2e1 * t116, -t140 * t81, (t151 + (t12 - t113) * t87) * t82 + t95;];
tauc_reg = t2;
