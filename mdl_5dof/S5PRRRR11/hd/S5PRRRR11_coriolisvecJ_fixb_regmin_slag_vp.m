% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:59
% DurationCPUTime: 1.47s
% Computational Cost: add. (1762->200), mult. (5880->296), div. (0->0), fcn. (4564->8), ass. (0->151)
t104 = sin(qJ(5));
t155 = qJD(5) * t104;
t105 = sin(qJ(4));
t106 = sin(qJ(3));
t102 = sin(pkin(5));
t188 = pkin(7) + pkin(8);
t151 = t102 * t188;
t132 = t106 * t151;
t103 = cos(pkin(5));
t109 = cos(qJ(3));
t160 = qJD(2) * t109;
t152 = pkin(2) * t160;
t162 = qJD(1) * t102;
t95 = t109 * t162;
t175 = t103 * t152 + t95;
t58 = -qJD(2) * t132 + t175;
t154 = t103 * qJD(2);
t96 = qJD(3) + t154;
t49 = t96 * pkin(3) + t58;
t108 = cos(qJ(4));
t144 = t106 * t162;
t153 = t103 * t106 * pkin(2);
t168 = t102 * t109;
t73 = t188 * t168 + t153;
t59 = t73 * qJD(2) + t144;
t57 = t108 * t59;
t124 = -t105 * t49 - t57;
t147 = t102 * t160;
t161 = qJD(2) * t102;
t148 = t106 * t161;
t194 = -t105 * t148 + t108 * t147;
t186 = t194 * pkin(9);
t18 = -t124 + t186;
t16 = t18 * t155;
t107 = cos(qJ(5));
t70 = t107 * t194;
t76 = -t105 * t147 - t108 * t148;
t36 = t104 * t76 + t70;
t80 = (-pkin(3) * t109 - pkin(2)) * t102;
t98 = t103 * qJD(1);
t77 = qJD(2) * t80 + t98;
t45 = -pkin(4) * t194 + t77;
t196 = -t45 * t36 + t16;
t121 = qJD(3) * t132;
t158 = qJD(3) * t109;
t145 = t103 * t158;
t128 = qJD(2) * t145;
t176 = pkin(2) * t128 + qJD(3) * t95;
t52 = -qJD(2) * t121 + t176;
t113 = -t109 * t151 - t153;
t53 = (t113 * qJD(2) - t144) * qJD(3);
t140 = -t105 * t52 + t108 * t53;
t112 = t124 * qJD(4) + t140;
t190 = qJD(3) + qJD(4);
t38 = t194 * t190;
t3 = -t38 * pkin(9) + t112;
t92 = qJD(4) + t96;
t88 = qJD(5) + t92;
t195 = (-t18 * t88 - t3) * t104 + t196;
t125 = -t104 * t194 + t107 * t76;
t182 = t125 * t36;
t6 = t125 ^ 2 - t36 ^ 2;
t120 = t105 * t109 + t106 * t108;
t189 = t102 * t190;
t51 = t120 * t189;
t39 = qJD(2) * t51;
t10 = qJD(5) * t70 - t104 * t39 + t107 * t38 + t76 * t155;
t4 = -t36 * t88 + t10;
t111 = t125 * qJD(5) - t104 * t38 - t107 * t39;
t5 = -t125 * t88 + t111;
t156 = qJD(4) * t108;
t157 = qJD(4) * t105;
t137 = -t105 * t53 - t108 * t52 - t49 * t156 + t59 * t157;
t2 = -pkin(9) * t39 - t137;
t149 = -t104 * t2 + t107 * t3;
t118 = t125 * t45 + t149;
t55 = t105 * t59;
t141 = t108 * t49 - t55;
t72 = t76 * pkin(9);
t17 = t141 + t72;
t193 = t106 * t109;
t164 = qJD(3) - t96;
t192 = qJD(5) - t88;
t191 = -t105 * t106 + t108 * t109;
t187 = pkin(3) * t88;
t185 = t76 * pkin(4);
t184 = pkin(7) * t106;
t78 = t191 * t102;
t79 = t120 * t102;
t43 = t104 * t79 - t107 * t78;
t50 = t191 * t189;
t13 = -t43 * qJD(5) - t104 * t51 + t107 * t50;
t183 = t13 * t88;
t180 = t50 * t92;
t179 = t76 * t194;
t178 = t77 * t76;
t177 = t108 * t58 - t55;
t174 = t103 * t111;
t173 = t103 * t39;
t172 = t107 * t18;
t82 = -pkin(2) * t161 + t98;
t171 = t109 * t82;
t99 = t102 ^ 2;
t170 = qJD(2) ^ 2 * t99;
t169 = t102 * t106;
t166 = t105 * t107;
t163 = t106 ^ 2 - t109 ^ 2;
t159 = qJD(3) * t106;
t15 = t92 * pkin(4) + t17;
t150 = -pkin(4) * t88 - t15;
t146 = t102 * t159;
t143 = qJD(5) * t15 + t2;
t139 = -t105 * t58 - t57;
t94 = pkin(2) * t145;
t67 = t94 - t121;
t68 = t113 * qJD(3);
t138 = -t105 * t67 + t108 * t68;
t136 = t96 + t154;
t135 = 0.2e1 * qJD(2) * qJD(3) * t99;
t134 = pkin(3) * t146;
t133 = pkin(3) * t148;
t131 = t102 * t96 * t158;
t127 = t164 * t169;
t126 = -t104 * t15 - t172;
t44 = t104 * t78 + t107 * t79;
t66 = (pkin(2) * t109 + pkin(3)) * t103 - t132;
t123 = -t105 * t66 - t108 * t73;
t122 = t102 * t136;
t117 = -t194 * t77 + t137;
t115 = t105 * t68 + t108 * t67 + t66 * t156 - t73 * t157;
t97 = t108 * pkin(3) + pkin(4);
t81 = t102 * t128;
t61 = t133 - t185;
t60 = -t78 * pkin(4) + t80;
t29 = t51 * t92;
t28 = t38 * t103;
t27 = t51 * pkin(4) + t134;
t26 = t39 * pkin(4) + qJD(3) * t133;
t25 = -t194 ^ 2 + t76 ^ 2;
t24 = t78 * pkin(9) - t123;
t23 = -t190 * t161 * t120 - t76 * t92;
t22 = -t194 * t92 + t38;
t21 = t103 * pkin(4) - t79 * pkin(9) - t105 * t73 + t108 * t66;
t20 = t72 + t177;
t19 = t139 - t186;
t14 = t44 * qJD(5) + t104 * t50 + t107 * t51;
t12 = t14 * t88;
t9 = t10 * t103;
t8 = -t50 * pkin(9) + t123 * qJD(4) + t138;
t7 = -t51 * pkin(9) + t115;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, (-t96 + t154) * t146, t81 - t131, 0, 0, 0, 0, 0, -t29 + t173, t28 - t180, 0, 0, 0, 0, 0, -t12 - t174, t9 - t183; 0, 0, 0, 0, t135 * t193, -t163 * t135, t81 + t131, -t122 * t159, 0, (-t109 * pkin(7) * t122 + ((t82 - t98) * t102 + (-t103 * t96 + (-t103 ^ 2 - t99) * qJD(2)) * pkin(2)) * t106) * qJD(3), -t94 * t96 - t176 * t103 + (-t99 * t152 + (t136 * t184 + t171) * t102) * qJD(3), t38 * t79 - t50 * t76, t194 * t50 + t38 * t78 - t39 * t79 + t51 * t76, t28 + t180, -t29 - t173, 0, t138 * t92 + t140 * t103 + t80 * t39 + t77 * t51 + (t124 * t103 + t123 * t92) * qJD(4) + (-qJD(2) * t78 - t194) * t134, -t115 * t92 + t137 * t103 + t80 * t38 + t77 * t50 + (qJD(2) * t79 - t76) * t134, t10 * t44 - t125 * t13, -t10 * t43 + t111 * t44 + t125 * t14 + t13 * t36, t9 + t183, -t12 + t174, 0, (-t104 * t7 + t107 * t8) * t88 + t149 * t103 - t27 * t36 - t60 * t111 + t26 * t43 + t45 * t14 + ((-t104 * t21 - t107 * t24) * t88 + t126 * t103) * qJD(5), t60 * t10 + t16 * t103 + t45 * t13 + t26 * t44 - t27 * t125 + (-(-qJD(5) * t24 + t8) * t88 - t3 * t103) * t104 + (-(qJD(5) * t21 + t7) * t88 - t143 * t103) * t107; 0, 0, 0, 0, -t170 * t193, t163 * t170, t164 * t147, -qJD(2) * t127, 0, -qJD(1) * t127 + (-t82 * t169 + t164 * (-pkin(7) * t168 - t153)) * qJD(2), t175 * t96 + (t164 * t184 - t171) * t161 - t176, t179, t25, t22, t23, 0, -t139 * t92 + t194 * t133 + t178 + (-t57 + (-pkin(3) * t92 - t49) * t105) * qJD(4) + t140, t177 * t92 + (t76 * t148 - t92 * t156) * pkin(3) + t117, t182, t6, t4, t5, 0, -(-t104 * t20 + t107 * t19) * t88 + t61 * t36 + (-t104 * t108 - t166) * qJD(4) * t187 + ((-pkin(3) * t166 - t104 * t97) * t88 + t126) * qJD(5) + t118, t61 * t125 + (t19 * t88 - t3 - (-qJD(4) - qJD(5)) * t105 * t187) * t104 + ((-pkin(3) * t156 - qJD(5) * t97 + t20) * t88 - t143) * t107 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t25, t22, t23, 0, -t124 * t92 + t112 + t178, t141 * t92 + t117, t182, t6, t4, t5, 0, -(-t104 * t17 - t172) * t88 - t36 * t185 + (t150 * t104 - t172) * qJD(5) + t118, -t125 * t185 + (t150 * qJD(5) + t17 * t88 - t2) * t107 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t6, t4, t5, 0, t192 * t126 + t118, (-t192 * t15 - t2) * t107 + t195;];
tauc_reg = t1;
