% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:52
% EndTime: 2019-03-09 01:35:56
% DurationCPUTime: 1.56s
% Computational Cost: add. (2432->256), mult. (4232->365), div. (0->0), fcn. (2277->6), ass. (0->150)
t118 = qJD(1) * qJD(2);
t61 = sin(pkin(9));
t105 = t61 * t118;
t65 = sin(qJ(5));
t119 = qJ(2) * qJD(1);
t68 = -pkin(1) - pkin(2);
t47 = qJD(1) * t68 + qJD(2);
t62 = cos(pkin(9));
t33 = -t61 * t119 + t47 * t62;
t28 = qJD(1) * pkin(3) + qJD(4) - t33;
t24 = qJD(1) * pkin(7) + t28;
t67 = cos(qJ(5));
t82 = qJD(3) * t65 - t24 * t67;
t11 = -qJD(5) * t82 + t105 * t65;
t57 = qJD(1) * qJD(4);
t121 = t62 * qJD(2);
t96 = -pkin(5) * t67 - pkin(8) * t65;
t74 = qJD(5) * t96 + t121;
t25 = qJD(1) * t74 - t57;
t64 = sin(qJ(6));
t66 = cos(qJ(6));
t19 = qJD(3) * t67 + t24 * t65;
t17 = qJD(5) * pkin(8) + t19;
t147 = t61 * t47;
t54 = t62 * qJ(2);
t58 = qJD(1) * qJ(4);
t95 = -pkin(5) * t65 + pkin(8) * t67;
t20 = t147 - t58 + (t95 + t54) * qJD(1);
t89 = t17 * t64 - t20 * t66;
t1 = -qJD(6) * t89 + t66 * t11 + t64 * t25;
t132 = qJD(1) * t65;
t49 = -qJD(6) + t132;
t173 = -t89 * t49 + t1;
t44 = t67 * t105;
t12 = t19 * qJD(5) - t44;
t10 = t12 * t67;
t91 = t11 * t65 - t10;
t60 = t67 ^ 2;
t172 = (qJD(1) * t60 - t49 * t65) * t64;
t6 = t17 * t66 + t20 * t64;
t2 = -qJD(6) * t6 - t64 * t11 + t66 * t25;
t171 = t6 * t49 - t2;
t92 = -t6 * t66 - t64 * t89;
t170 = qJD(5) * t92 + t12;
t169 = -qJD(5) * (-t19 * t67 - t65 * t82) + t91;
t166 = t12 * t64;
t165 = t12 * t66;
t16 = -qJD(5) * pkin(5) + t82;
t164 = t16 * t64;
t163 = t16 * t66;
t116 = qJD(5) * qJD(6);
t123 = qJD(6) * t67;
t110 = t64 * t123;
t128 = qJD(5) * t66;
t75 = t128 * t65 + t110;
t26 = qJD(1) * t75 + t116 * t66;
t162 = t26 * t64;
t161 = t26 * t65;
t117 = qJD(1) * qJD(5);
t106 = t65 * t117;
t109 = t66 * t123;
t27 = qJD(1) * t109 + (-t106 - t116) * t64;
t160 = t27 * t65;
t159 = t27 * t66;
t131 = qJD(1) * t67;
t40 = t131 * t64 + t128;
t113 = t66 * t131;
t120 = t64 * qJD(5);
t41 = t113 - t120;
t158 = t40 * t41;
t157 = t40 * t49;
t156 = t40 * t64;
t155 = t40 * t66;
t154 = t40 * t67;
t153 = t41 * t49;
t152 = t41 * t64;
t151 = t41 * t66;
t150 = t41 * t67;
t149 = t49 * t64;
t148 = t49 * t66;
t70 = qJD(1) ^ 2;
t146 = t61 * t70;
t145 = t62 * t70;
t144 = t64 * t65;
t143 = t65 * t66;
t142 = t65 * t70;
t141 = t67 * t26;
t140 = t67 * t27;
t69 = qJD(5) ^ 2;
t139 = t69 * t67;
t127 = qJD(5) * t67;
t111 = t62 * t127;
t36 = t144 * t62 + t61 * t66;
t138 = qJD(6) * t36 - t111 * t66 - (t143 * t61 + t62 * t64) * qJD(1);
t80 = t143 * t62 - t61 * t64;
t137 = qJD(6) * t80 + t111 * t64 - (-t144 * t61 + t62 * t66) * qJD(1);
t136 = t61 * t68 + t54;
t59 = t65 ^ 2;
t135 = t59 - t60;
t134 = t59 + t60;
t133 = t69 + t70;
t130 = qJD(2) * t61;
t129 = qJD(5) * t65;
t126 = qJD(6) * t64;
t125 = qJD(6) * t65;
t124 = qJD(6) * t66;
t34 = t119 * t62 + t147;
t29 = t34 - t58;
t122 = t29 * qJD(1);
t115 = t49 * t143;
t114 = t67 * t142;
t39 = -qJ(4) + t136;
t112 = t62 * t129;
t108 = 0.2e1 * t118;
t107 = t62 * t133;
t104 = t62 * t118;
t103 = t67 * t117;
t102 = -t61 * qJ(2) + t62 * t68;
t101 = t29 + t130;
t100 = t49 * t109;
t99 = -qJD(1) + t125;
t98 = pkin(3) - t102;
t97 = t65 * t103;
t38 = pkin(7) + t98;
t94 = -t125 * t38 - qJD(4) + t74;
t93 = t6 * t64 - t66 * t89;
t88 = t19 * t65 - t67 * t82;
t43 = -t57 + t104;
t48 = -qJD(4) + t121;
t86 = t29 * t48 + t43 * t39;
t85 = t33 * t61 - t34 * t62;
t84 = t67 * (t49 + t132);
t81 = t66 * t60 * t117 - t110 * t49;
t79 = pkin(8) * t127 - t16 * t65;
t78 = -qJD(1) * t48 - t38 * t69 - t43;
t77 = t19 * t127 + t129 * t82 + t91;
t76 = qJD(1) * t39 - t130 + t29;
t30 = t95 + t39;
t73 = qJD(6) * t30 + t127 * t38 + t130 * t65;
t72 = -qJD(6) * t93 + t1 * t66 - t2 * t64;
t71 = qJD(5) * t16 + t72;
t56 = t69 * t65;
t42 = t96 * qJD(1);
t14 = t143 * t38 + t30 * t64;
t13 = -t144 * t38 + t30 * t66;
t8 = t42 * t64 - t66 * t82;
t7 = t42 * t66 + t64 * t82;
t4 = -t64 * t73 + t66 * t94;
t3 = t64 * t94 + t66 * t73;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, qJ(2) * t108, 0, 0, 0, 0, 0, 0, 0.2e1 * t105, 0.2e1 * t104, 0 ((-t102 * t61 + t136 * t62) * qJD(1) - t85) * qJD(2), 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t105, t57 + (-t48 - t121) * qJD(1) (qJD(1) * t98 + t28) * t130 + t86, -0.2e1 * t97, 0.2e1 * t135 * t117, t56, 0.2e1 * t97, t139, 0, -t127 * t76 + t65 * t78, t129 * t76 + t67 * t78, t105 * t134 + t77, t88 * t130 + t169 * t38 + t86, -t141 * t66 - t41 * t75 (t152 + t155) * t129 + (t162 - t159 + (-t151 + t156) * qJD(6)) * t67, -t161 + (-t115 + t150) * qJD(5) + t81, t64 * t140 + (-t120 * t65 + t109) * t40, -t100 - t160 + (-t154 - t172) * qJD(5), qJD(5) * t84, -t4 * t49 + (-t2 + (-t38 * t40 + t164) * qJD(5)) * t65 + (t40 * t130 - t16 * t124 - t166 + t27 * t38 + (-qJD(1) * t13 + t89) * qJD(5)) * t67, t3 * t49 + (t1 + (-t38 * t41 + t163) * qJD(5)) * t65 + (t41 * t130 + t16 * t126 - t165 - t26 * t38 + (qJD(1) * t14 + t6) * qJD(5)) * t67, -t13 * t26 + t14 * t27 + t3 * t40 + t4 * t41 - t93 * t129 + (-qJD(6) * t92 + t1 * t64 + t2 * t66) * t67, -t38 * t10 + t1 * t14 + t2 * t13 + t6 * t3 - t89 * t4 + (t129 * t38 - t130 * t67) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t70 * qJ(2), 0, 0, 0, 0, 0, 0, -t146, -t145, 0, t85 * qJD(1), 0, 0, 0, 0, 0, 0, 0, t146, t145, t43 * t61 + (-t29 * t62 + (-t28 - t121) * t61) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t103 * t61 + t107 * t65, 0.2e1 * t106 * t61 + t107 * t67, -t134 * t146 (-qJD(1) * t88 + t43) * t61 + (-t122 - t169) * t62, 0, 0, 0, 0, 0, 0, t40 * t112 - t137 * t49 + (-t27 * t62 + (-qJD(5) * t36 - t40 * t61) * qJD(1)) * t67, t41 * t112 + t138 * t49 + (t26 * t62 + (-qJD(5) * t80 - t41 * t61) * qJD(1)) * t67, t137 * t41 + t138 * t40 - t36 * t26 - t27 * t80, t62 * t10 - t1 * t80 + t2 * t36 + t138 * t6 - t137 * t89 + (t131 * t61 - t112) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t56, 0, -qJD(5) * t88 + t11 * t67 + t12 * t65, 0, 0, 0, 0, 0, 0, t100 - t160 + (-t154 + t172) * qJD(5), t161 + (-t115 - t150) * qJD(5) + t81 (t152 - t155) * t129 + (t162 + t159 + (-t151 - t156) * qJD(6)) * t67, t170 * t65 + t71 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t101 * qJD(1), 0, 0, 0, 0, 0, 0, -t56 - t142, -t133 * t67, 0, t77 + t122, 0, 0, 0, 0, 0, 0, t140 + t99 * t148 + (-t40 * t65 + t64 * t84) * qJD(5), -t141 - t99 * t149 + (t67 * t148 + (-t41 + t113) * t65) * qJD(5) (t127 * t40 - t41 * t99 + t160) * t66 + (-t127 * t41 - t40 * t99 + t161) * t64, t93 * qJD(1) - t170 * t67 + t71 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t135 * t70, 0, -t114, 0, 0, t122 * t67 + t44, -t101 * t132, 0, 0, t148 * t41 + t162 (t26 - t157) * t66 + (t27 - t153) * t64, -t49 * t124 + (t115 + (-t41 - t120) * t67) * qJD(1), t149 * t40 + t159, t49 * t126 + (-t49 * t144 + (t40 - t128) * t67) * qJD(1), -t49 * t131, pkin(5) * t27 - t165 + t19 * t40 + t49 * t7 + (pkin(8) * t148 + t164) * qJD(6) + (t64 * t79 - t67 * t89) * qJD(1), -pkin(5) * t26 + t166 + t19 * t41 - t49 * t8 + (-pkin(8) * t149 + t163) * qJD(6) + (-t6 * t67 + t66 * t79) * qJD(1), -t40 * t8 - t41 * t7 + ((-qJD(6) * t41 + t27) * pkin(8) + t173) * t66 + ((-qJD(6) * t40 + t26) * pkin(8) + t171) * t64, -t12 * pkin(5) + pkin(8) * t72 - t16 * t19 - t6 * t8 + t7 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, -t40 ^ 2 + t41 ^ 2, t26 + t157, -t158, t27 + t153, -t103, t16 * t41 - t171, -t16 * t40 - t173, 0, 0;];
tauc_reg  = t5;
