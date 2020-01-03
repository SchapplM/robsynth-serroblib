% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:41
% EndTime: 2020-01-03 11:57:44
% DurationCPUTime: 1.49s
% Computational Cost: add. (2547->262), mult. (4039->352), div. (0->0), fcn. (2482->14), ass. (0->171)
t114 = qJ(1) + qJ(2);
t102 = pkin(8) + t114;
t92 = sin(t102);
t93 = cos(t102);
t202 = g(2) * t93 + g(3) * t92;
t108 = qJDD(1) + qJDD(2);
t116 = sin(pkin(8));
t118 = cos(pkin(8));
t123 = cos(qJ(2));
t199 = pkin(1) * qJD(2);
t164 = qJD(1) * t199;
t120 = sin(qJ(2));
t171 = qJDD(1) * t120;
t212 = pkin(1) * t171 + t123 * t164;
t209 = pkin(1) * t123;
t99 = qJDD(1) * t209;
t50 = pkin(2) * t108 - t120 * t164 + t99;
t168 = t212 * t116 - t118 * t50;
t152 = qJDD(4) + t168;
t22 = -t108 * pkin(3) + t152;
t137 = -t22 - t202;
t111 = qJD(1) + qJD(2);
t115 = sin(pkin(9));
t109 = t115 ^ 2;
t117 = cos(pkin(9));
t110 = t117 ^ 2;
t176 = t109 + t110;
t211 = t111 * t176;
t200 = pkin(1) * qJD(1);
t165 = t123 * t200;
t166 = t120 * t200;
t78 = t116 * t166;
t55 = t118 * t165 - t78;
t177 = qJD(4) - t55;
t210 = pkin(1) * t120;
t208 = pkin(2) * t116;
t207 = pkin(2) * t118;
t206 = pkin(4) * t117;
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t174 = qJD(4) * t117;
t180 = t117 * t122;
t181 = t117 * t119;
t136 = -pkin(7) * t115 - pkin(3) - t206;
t59 = t136 - t207;
t91 = qJ(4) + t208;
t34 = t122 * t59 - t91 * t181;
t189 = qJD(5) * t34;
t179 = t118 * t120;
t130 = pkin(1) * (t116 * t123 + t179);
t53 = qJD(1) * t130;
t205 = -t119 * t53 + t122 * t174 - t55 * t180 + t189;
t35 = t119 * t59 + t91 * t180;
t188 = qJD(5) * t35;
t204 = -t119 * t174 - t122 * t53 + t55 * t181 - t188;
t194 = t115 * t93;
t195 = t115 * t92;
t203 = g(2) * t194 + g(3) * t195;
t201 = -g(2) * t92 + g(3) * t93;
t98 = pkin(2) + t209;
t58 = pkin(1) * t179 + t116 * t98;
t62 = pkin(2) * t111 + t165;
t41 = t118 * t62 - t78;
t144 = qJD(4) - t41;
t23 = t136 * t111 + t144;
t42 = t116 * t62 + t118 * t166;
t37 = qJ(4) * t111 + t42;
t32 = qJD(3) * t115 + t117 * t37;
t8 = t119 * t23 + t122 * t32;
t198 = qJD(5) * t8;
t197 = t111 * t53;
t54 = qJD(2) * t130;
t196 = t111 * t54;
t30 = t116 * t50 + t212 * t118;
t20 = qJ(4) * t108 + qJD(4) * t111 + t30;
t14 = -t117 * qJDD(3) + t115 * t20;
t193 = t14 * t115;
t31 = -t117 * qJD(3) + t115 * t37;
t192 = t31 * t115;
t182 = t117 * t108;
t72 = -qJDD(5) + t182;
t191 = t72 * t117;
t89 = t116 * t210;
t57 = t118 * t98 - t89;
t38 = t136 - t57;
t51 = qJ(4) + t58;
t18 = t122 * t38 - t51 * t181;
t190 = qJD(5) * t18;
t187 = t108 * t119;
t186 = t108 * t122;
t107 = t111 ^ 2;
t185 = t109 * t107;
t184 = t111 * t117;
t183 = t115 * t108;
t178 = t119 * t122;
t112 = t119 ^ 2;
t113 = t122 ^ 2;
t175 = t112 - t113;
t173 = qJD(5) * t119;
t172 = qJD(5) * t122;
t161 = t115 * t173;
t7 = -t119 * t32 + t122 * t23;
t170 = t7 * t161 - t203;
t169 = t22 * t115 + t203;
t104 = cos(t114);
t97 = pkin(2) * t104;
t167 = t93 * pkin(3) + t92 * qJ(4) + t97;
t75 = -qJD(5) + t184;
t163 = t75 * t173;
t160 = -t14 * t117 - g(1);
t103 = sin(t114);
t159 = g(2) * t103 - g(3) * t104;
t158 = t72 - t182;
t157 = t72 + t182;
t156 = t108 * t176;
t155 = t111 * (-qJD(5) - t75);
t154 = qJD(1) * (-qJD(2) + t111);
t153 = qJD(2) * (-qJD(1) - t111);
t150 = t178 * t185;
t149 = t111 * t119 * t172;
t96 = pkin(2) * t103;
t148 = t92 * pkin(3) - qJ(4) * t93 + t96;
t147 = pkin(7) * t194 + t93 * t206 + t167;
t146 = t168 + t202;
t145 = -t30 - t201;
t143 = -g(2) * t104 - g(3) * t103;
t121 = sin(qJ(1));
t124 = cos(qJ(1));
t142 = -g(2) * t124 - g(3) * t121;
t141 = t119 * t7 - t122 * t8;
t52 = -pkin(3) - t57;
t140 = t108 * t52 + t196;
t94 = -pkin(3) - t207;
t139 = t108 * t94 - t197;
t15 = qJDD(3) * t115 + t117 * t20;
t138 = t15 * t117 + t193 + t201;
t56 = t118 * t123 * t199 - qJD(2) * t89;
t135 = t143 + t99;
t134 = g(1) * t115 - qJD(5) * t23 - t15;
t19 = t119 * t38 + t51 * t180;
t133 = pkin(7) * t195 + t92 * t206 + t148;
t132 = -qJD(5) * t32 - t111 * t192;
t49 = qJD(4) + t56;
t131 = (t14 * t51 + t31 * t49) * t115;
t129 = (t14 * t91 + t177 * t31) * t115;
t13 = t136 * t108 + t152;
t11 = t122 * t13;
t3 = -t119 * t15 + t11 - t198;
t44 = -t119 * t93 + t92 * t180;
t46 = t119 * t92 + t93 * t180;
t128 = -g(2) * t46 - g(3) * t44 - t117 * t3 + t119 * t193 + t172 * t192;
t127 = -t75 ^ 2 - t185;
t2 = t7 * qJD(5) + t119 * t13 + t122 * t15;
t43 = -t122 * t93 - t92 * t181;
t45 = -t92 * t122 + t93 * t181;
t126 = g(2) * t45 - g(3) * t43 + t2 * t117 + t122 * t193 - t31 * t161;
t106 = t124 * pkin(1);
t105 = t121 * pkin(1);
t88 = t110 * t108;
t87 = t109 * t108;
t63 = 0.2e1 * t115 * t182;
t61 = t161 * t184;
t40 = (t108 * t113 - 0.2e1 * t149) * t109;
t39 = (t108 * t112 + 0.2e1 * t149) * t109;
t36 = -pkin(3) * t111 + t144;
t33 = 0.2e1 * (t175 * t111 * qJD(5) - t108 * t178) * t109;
t17 = (t157 * t119 + (t75 + t184) * t172) * t115;
t16 = t61 + (-t157 * t122 + t163) * t115;
t5 = -qJD(5) * t19 + t122 * t54 - t49 * t181;
t4 = t119 * t54 + t49 * t180 + t190;
t1 = [0, 0, 0, 0, 0, qJDD(1), t142, g(2) * t121 - g(3) * t124, 0, 0, 0, 0, 0, 0, 0, t108, (t108 * t123 + t120 * t153) * pkin(1) + t135, ((-qJDD(1) - t108) * t120 + t123 * t153) * pkin(1) + t159, 0, (t142 + (t120 ^ 2 + t123 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t108, t108 * t57 - t146 - t196, -t108 * t58 - t111 * t56 + t145, 0, t30 * t58 + t42 * t56 - t168 * t57 - t41 * t54 - g(2) * (t97 + t106) - g(3) * (t96 + t105), t87, t63, 0, t88, 0, 0, (t137 - t140) * t117, t140 * t115 + t169, t51 * t156 + t49 * t211 + t138, t22 * t52 + t36 * t54 - g(2) * (t106 + t167) - g(3) * (t105 + t148) + (t15 * t51 + t32 * t49) * t117 + t131, t40, t33, t16, t39, t17, t191, -t18 * t72 - t5 * t75 + (t51 * t187 + (t119 * t49 + t51 * t172) * t111) * t109 + t128, t19 * t72 + t4 * t75 + (t51 * t186 + (t122 * t49 - t51 * t173) * t111) * t109 + t126, ((-t108 * t19 - t2 + (-t4 + t190) * t111) * t119 + (-t108 * t18 - t111 * t5 - t3 + (-t111 * t19 - t8) * qJD(5)) * t122) * t115 + t170, t2 * t19 + t8 * t4 + t3 * t18 + t7 * t5 - g(2) * (t106 + t147) - g(3) * (t105 + t133) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t154 * t210 + t135, (t123 * t154 - t171) * pkin(1) + t159, 0, 0, 0, 0, 0, 0, 0, t108, t108 * t207 - t146 + t197, -t108 * t208 + t111 * t55 + t145, 0, t41 * t53 - t42 * t55 + (t116 * t30 - t118 * t168 + t143) * pkin(2), t87, t63, 0, t88, 0, 0, (t137 - t139) * t117, t139 * t115 + t169, t91 * t156 + t177 * t211 + t138, t22 * t94 - t36 * t53 - g(2) * t167 - g(3) * t148 + (t15 * t91 + t177 * t32) * t117 + t129, t40, t33, t16, t39, t17, t191, -t34 * t72 - t204 * t75 + (t91 * t187 + (t177 * t119 + t91 * t172) * t111) * t109 + t128, t35 * t72 + t205 * t75 + (t91 * t186 + (t177 * t122 - t91 * t173) * t111) * t109 + t126, ((-t108 * t35 - t2) * t119 + (-t108 * t34 - t198 - t3) * t122 + ((-t188 - t204) * t122 + (t189 - t205) * t119) * t111) * t115 + t170, -g(2) * t147 - g(3) * t133 + t2 * t35 + t204 * t7 + t205 * t8 + t3 * t34 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t15 + t160, 0, 0, 0, 0, 0, 0, (t158 * t119 + (t75 - t184) * t172) * t115, t61 + (t122 * t158 - t163) * t115, 0, (-t119 * t3 + t122 * t2 + (-t119 * t8 - t122 * t7) * qJD(5)) * t115 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t183, -t176 * t107, (-t117 * t32 - t192) * t111 - t137, 0, 0, 0, 0, 0, 0, t119 * t127 - t122 * t72, t119 * t72 + t122 * t127, (-t112 - t113) * t183, t2 * t119 + t3 * t122 - t141 * qJD(5) + (t117 * t141 - t192) * t111 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t175 * t185, (t119 * t155 + t186) * t115, -t150, (t122 * t155 - t187) * t115, -t72, -g(2) * t43 - g(3) * t45 + t119 * t134 + t122 * t132 - t8 * t75 + t11, g(2) * t44 - g(3) * t46 - t7 * t75 + t134 * t122 + (-t13 - t132) * t119, 0, 0;];
tau_reg = t1;
