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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:57
% DurationCPUTime: 1.61s
% Computational Cost: add. (2547->266), mult. (4039->355), div. (0->0), fcn. (2482->14), ass. (0->175)
t121 = cos(qJ(2));
t196 = pkin(1) * qJD(2);
t160 = qJD(1) * t196;
t118 = sin(qJ(2));
t168 = qJDD(1) * t118;
t214 = pkin(1) * t168 + t121 * t160;
t112 = qJ(1) + qJ(2);
t102 = sin(t112);
t103 = cos(t112);
t213 = g(1) * t102 - g(2) * t103;
t109 = qJD(1) + qJD(2);
t113 = sin(pkin(9));
t107 = t113 ^ 2;
t115 = cos(pkin(9));
t108 = t115 ^ 2;
t173 = t107 + t108;
t212 = t109 * t173;
t116 = cos(pkin(8));
t197 = pkin(1) * qJD(1);
t161 = t121 * t197;
t114 = sin(pkin(8));
t162 = t118 * t197;
t76 = t114 * t162;
t55 = t116 * t161 - t76;
t174 = qJD(4) - t55;
t101 = pkin(8) + t112;
t90 = sin(t101);
t83 = g(1) * t90;
t91 = cos(t101);
t82 = g(2) * t91;
t211 = pkin(1) * t118;
t210 = pkin(2) * t102;
t208 = t114 * pkin(2);
t207 = t115 * pkin(4);
t206 = t116 * pkin(2);
t119 = sin(qJ(1));
t205 = t119 * pkin(1);
t204 = t121 * pkin(1);
t117 = sin(qJ(5));
t120 = cos(qJ(5));
t171 = qJD(4) * t115;
t177 = t115 * t120;
t178 = t115 * t117;
t132 = -t113 * pkin(7) - pkin(3) - t207;
t59 = t132 - t206;
t89 = qJ(4) + t208;
t34 = t120 * t59 - t89 * t178;
t186 = qJD(5) * t34;
t176 = t116 * t118;
t128 = pkin(1) * (t114 * t121 + t176);
t53 = qJD(1) * t128;
t203 = -t117 * t53 + t120 * t171 - t55 * t177 + t186;
t35 = t117 * t59 + t89 * t177;
t185 = qJD(5) * t35;
t202 = -t117 * t171 - t120 * t53 + t55 * t178 - t185;
t194 = t113 * t91;
t201 = g(2) * t194 - t113 * t83;
t200 = -g(1) * t91 - g(2) * t90;
t199 = t82 - t83;
t97 = pkin(2) + t204;
t58 = pkin(1) * t176 + t114 * t97;
t198 = g(1) * t103 + g(2) * t102;
t62 = t109 * pkin(2) + t161;
t41 = t116 * t62 - t76;
t140 = qJD(4) - t41;
t23 = t132 * t109 + t140;
t42 = t114 * t62 + t116 * t162;
t37 = t109 * qJ(4) + t42;
t32 = t113 * qJD(3) + t115 * t37;
t8 = t117 * t23 + t120 * t32;
t195 = qJD(5) * t8;
t193 = t115 * t90;
t106 = qJDD(1) + qJDD(2);
t98 = qJDD(1) * t204;
t50 = t106 * pkin(2) - t118 * t160 + t98;
t30 = t114 * t50 + t214 * t116;
t20 = t106 * qJ(4) + t109 * qJD(4) + t30;
t14 = -t115 * qJDD(3) + t113 * t20;
t192 = t14 * t113;
t31 = -t115 * qJD(3) + t113 * t37;
t191 = t31 * t113;
t190 = t53 * t109;
t54 = qJD(2) * t128;
t189 = t54 * t109;
t180 = t115 * t106;
t71 = -qJDD(5) + t180;
t188 = t71 * t115;
t87 = t114 * t211;
t57 = t116 * t97 - t87;
t38 = t132 - t57;
t51 = qJ(4) + t58;
t18 = t120 * t38 - t51 * t178;
t187 = qJD(5) * t18;
t184 = t106 * t117;
t183 = t106 * t120;
t105 = t109 ^ 2;
t182 = t107 * t105;
t181 = t113 * t106;
t179 = t115 * t109;
t175 = t117 * t120;
t110 = t117 ^ 2;
t111 = t120 ^ 2;
t172 = t110 - t111;
t170 = qJD(5) * t117;
t169 = qJD(5) * t120;
t167 = pkin(7) * t83;
t156 = t113 * t170;
t7 = -t117 * t32 + t120 * t23;
t166 = t7 * t156 - t201;
t164 = t214 * t114 - t116 * t50;
t148 = qJDD(4) + t164;
t22 = -t106 * pkin(3) + t148;
t165 = t22 * t113 + t201;
t96 = pkin(2) * t103;
t163 = t91 * pkin(3) + t90 * qJ(4) + t96;
t73 = -qJD(5) + t179;
t159 = t73 * t170;
t157 = -t22 - t82;
t155 = -t14 * t115 - g(3);
t154 = t71 - t180;
t153 = t71 + t180;
t152 = t106 * t173;
t151 = t109 * (-qJD(5) - t73);
t150 = qJD(1) * (-qJD(2) + t109);
t149 = qJD(2) * (-qJD(1) - t109);
t146 = t175 * t182;
t145 = t98 + t213;
t144 = t109 * t117 * t169;
t143 = pkin(7) * t194 + t91 * t207 + t163;
t142 = t164 + t199;
t141 = -t30 - t200;
t122 = cos(qJ(1));
t139 = g(1) * t119 - g(2) * t122;
t138 = t117 * t7 - t120 * t8;
t56 = t116 * t121 * t196 - qJD(2) * t87;
t49 = qJD(4) + t56;
t137 = t14 * t51 + t31 * t49;
t52 = -pkin(3) - t57;
t136 = t106 * t52 + t189;
t92 = -pkin(3) - t206;
t135 = t106 * t92 - t190;
t15 = t113 * qJDD(3) + t115 * t20;
t134 = t15 * t115 + t192 + t200;
t133 = -t90 * pkin(3) + t91 * qJ(4) - t210;
t131 = g(3) * t113 - qJD(5) * t23 - t15;
t19 = t117 * t38 + t51 * t177;
t130 = -qJD(5) * t32 - t109 * t191;
t129 = t14 * t89 + t174 * t31;
t13 = t132 * t106 + t148;
t11 = t120 * t13;
t3 = -t117 * t15 + t11 - t195;
t44 = t91 * t117 - t90 * t177;
t46 = t90 * t117 + t91 * t177;
t127 = -g(1) * t44 - g(2) * t46 - t3 * t115 + t117 * t192 + t169 * t191;
t126 = -t73 ^ 2 - t182;
t125 = -pkin(4) * t193 + t133;
t2 = t7 * qJD(5) + t117 * t13 + t120 * t15;
t43 = t91 * t120 + t90 * t178;
t45 = t90 * t120 - t91 * t178;
t124 = -g(1) * t43 - g(2) * t45 + t2 * t115 + t120 * t192 - t31 * t156;
t104 = t122 * pkin(1);
t86 = t108 * t106;
t85 = t107 * t106;
t68 = g(1) * t193;
t63 = 0.2e1 * t113 * t180;
t61 = t156 * t179;
t40 = (t106 * t111 - 0.2e1 * t144) * t107;
t39 = (t106 * t110 + 0.2e1 * t144) * t107;
t36 = -t109 * pkin(3) + t140;
t33 = 0.2e1 * (t172 * t109 * qJD(5) - t106 * t175) * t107;
t17 = (t153 * t117 + (t73 + t179) * t169) * t113;
t16 = t61 + (-t153 * t120 + t159) * t113;
t5 = -t19 * qJD(5) + t120 * t54 - t49 * t178;
t4 = t117 * t54 + t49 * t177 + t187;
t1 = [0, 0, 0, 0, 0, qJDD(1), t139, g(1) * t122 + g(2) * t119, 0, 0, 0, 0, 0, 0, 0, t106, (t106 * t121 + t118 * t149) * pkin(1) + t145, ((-qJDD(1) - t106) * t118 + t121 * t149) * pkin(1) + t198, 0, (t139 + (t118 ^ 2 + t121 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t106, t57 * t106 - t142 - t189, -t58 * t106 - t56 * t109 + t141, 0, t30 * t58 + t42 * t56 - t164 * t57 - t41 * t54 - g(1) * (-t205 - t210) - g(2) * (t96 + t104), t85, t63, 0, t86, 0, 0, t68 + (-t136 + t157) * t115, t136 * t113 + t165, t51 * t152 + t49 * t212 + t134, t22 * t52 + t36 * t54 - g(1) * (t133 - t205) - g(2) * (t104 + t163) + (t15 * t51 + t32 * t49) * t115 + t137 * t113, t40, t33, t16, t39, t17, t188, -t18 * t71 - t5 * t73 + (t51 * t184 + (t117 * t49 + t51 * t169) * t109) * t107 + t127, t19 * t71 + t4 * t73 + (t51 * t183 + (t120 * t49 - t51 * t170) * t109) * t107 + t124, ((-t106 * t19 - t2 + (-t4 + t187) * t109) * t117 + (-t106 * t18 - t109 * t5 - t3 + (-t109 * t19 - t8) * qJD(5)) * t120) * t113 + t166, t2 * t19 + t8 * t4 + t3 * t18 + t7 * t5 - g(1) * (t125 - t205) - g(2) * (t104 + t143) + (t137 + t167) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t150 * t211 + t145, (t121 * t150 - t168) * pkin(1) + t198, 0, 0, 0, 0, 0, 0, 0, t106, t106 * t206 - t142 + t190, -t106 * t208 + t55 * t109 + t141, 0, t41 * t53 - t42 * t55 + (t114 * t30 - t116 * t164 + t213) * pkin(2), t85, t63, 0, t86, 0, 0, t68 + (-t135 + t157) * t115, t135 * t113 + t165, t89 * t152 + t174 * t212 + t134, t22 * t92 - t36 * t53 - g(1) * t133 - g(2) * t163 + (t15 * t89 + t174 * t32) * t115 + t129 * t113, t40, t33, t16, t39, t17, t188, -t34 * t71 - t202 * t73 + (t89 * t184 + (t174 * t117 + t89 * t169) * t109) * t107 + t127, t35 * t71 + t203 * t73 + (t89 * t183 + (t174 * t120 - t89 * t170) * t109) * t107 + t124, ((-t106 * t35 - t2) * t117 + (-t106 * t34 - t195 - t3) * t120 + ((-t185 - t202) * t120 + (t186 - t203) * t117) * t109) * t113 + t166, t2 * t35 + t3 * t34 - g(1) * t125 - g(2) * t143 + t203 * t8 + t202 * t7 + (t129 + t167) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t113 + t155, 0, 0, 0, 0, 0, 0, (t154 * t117 + (t73 - t179) * t169) * t113, t61 + (t154 * t120 - t159) * t113, 0, (-t117 * t3 + t120 * t2 + (-t117 * t8 - t120 * t7) * qJD(5)) * t113 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t181, -t173 * t105, (-t32 * t115 - t191) * t109 + t22 + t199, 0, 0, 0, 0, 0, 0, t126 * t117 - t120 * t71, t117 * t71 + t126 * t120, (-t110 - t111) * t181, t2 * t117 + t3 * t120 - t138 * qJD(5) + (t138 * t115 - t191) * t109 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, -t172 * t182, (t117 * t151 + t183) * t113, -t146, (t120 * t151 - t184) * t113, -t71, -g(1) * t45 + g(2) * t43 + t131 * t117 + t130 * t120 - t8 * t73 + t11, g(1) * t46 - g(2) * t44 - t7 * t73 + t131 * t120 + (-t13 - t130) * t117, 0, 0;];
tau_reg = t1;
