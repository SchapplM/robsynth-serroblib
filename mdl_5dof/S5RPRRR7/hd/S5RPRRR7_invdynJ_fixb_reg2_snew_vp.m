% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:11
% EndTime: 2019-12-31 19:04:17
% DurationCPUTime: 2.38s
% Computational Cost: add. (9399->313), mult. (18340->434), div. (0->0), fcn. (12150->10), ass. (0->192)
t168 = sin(qJ(5));
t170 = sin(qJ(3));
t200 = qJD(1) * qJD(3);
t154 = t170 * t200;
t174 = cos(qJ(3));
t198 = t174 * qJDD(1);
t137 = -t154 + t198;
t128 = -qJDD(4) + t137;
t125 = -qJDD(5) + t128;
t169 = sin(qJ(4));
t173 = cos(qJ(4));
t204 = qJD(1) * t170;
t129 = -t173 * qJD(3) + t169 * t204;
t131 = t169 * qJD(3) + t173 * t204;
t172 = cos(qJ(5));
t110 = t172 * t129 + t168 * t131;
t112 = -t168 * t129 + t172 * t131;
t83 = t112 * t110;
t230 = -t125 - t83;
t233 = t168 * t230;
t232 = t172 * t230;
t176 = qJD(1) ^ 2;
t194 = t174 * t200;
t199 = t170 * qJDD(1);
t136 = t194 + t199;
t182 = -t169 * qJDD(3) - t173 * t136;
t106 = -t129 * qJD(4) - t182;
t189 = -t173 * qJDD(3) + t169 * t136;
t180 = t131 * qJD(4) + t189;
t63 = -t110 * qJD(5) + t172 * t106 - t168 * t180;
t150 = t174 * qJD(1) - qJD(4);
t144 = -qJD(5) + t150;
t96 = t110 * t144;
t231 = t96 + t63;
t213 = t131 * t129;
t177 = -t128 - t213;
t229 = t169 * t177;
t228 = t173 * t177;
t120 = t129 * t150;
t88 = t106 - t120;
t191 = t168 * t106 + t172 * t180;
t47 = (qJD(5) + t144) * t112 + t191;
t84 = (qJD(4) + t150) * t131 + t189;
t108 = t110 ^ 2;
t109 = t112 ^ 2;
t227 = t129 ^ 2;
t127 = t131 ^ 2;
t143 = t144 ^ 2;
t148 = t150 ^ 2;
t226 = qJD(3) ^ 2;
t171 = sin(qJ(1));
t175 = cos(qJ(1));
t193 = t171 * g(1) - t175 * g(2);
t132 = qJDD(1) * pkin(1) + t193;
t186 = t175 * g(1) + t171 * g(2);
t133 = -t176 * pkin(1) - t186;
t164 = sin(pkin(9));
t165 = cos(pkin(9));
t190 = t165 * t132 - t164 * t133;
t103 = -qJDD(1) * pkin(2) - t176 * pkin(6) - t190;
t184 = -t137 + t154;
t185 = t136 + t194;
t75 = pkin(3) * t184 - pkin(7) * t185 + t103;
t205 = t164 * t132 + t165 * t133;
t104 = -t176 * pkin(2) + qJDD(1) * pkin(6) + t205;
t187 = -t174 * pkin(3) - t170 * pkin(7);
t188 = t176 * t187 + t104;
t206 = -g(3) + qJDD(2);
t192 = t170 * t206;
t81 = -t226 * pkin(3) + qJDD(3) * pkin(7) + t174 * t188 + t192;
t44 = t169 * t81 - t173 * t75;
t33 = t177 * pkin(4) - t88 * pkin(8) - t44;
t117 = -t150 * pkin(4) - t131 * pkin(8);
t45 = t169 * t75 + t173 * t81;
t34 = -t227 * pkin(4) - pkin(8) * t180 + t150 * t117 + t45;
t13 = t168 * t34 - t172 * t33;
t14 = t168 * t33 + t172 * t34;
t7 = -t172 * t13 + t168 * t14;
t225 = pkin(4) * t7;
t50 = -t96 + t63;
t26 = -t168 * t47 - t172 * t50;
t224 = pkin(4) * t26;
t223 = t169 * t7;
t222 = t173 * t7;
t153 = t174 * t206;
t80 = -qJDD(3) * pkin(3) - t226 * pkin(7) + t170 * t188 - t153;
t43 = pkin(4) * t180 - t227 * pkin(8) + t131 * t117 + t80;
t221 = t168 * t43;
t70 = t125 - t83;
t220 = t168 * t70;
t219 = t169 * t80;
t99 = t128 - t213;
t218 = t169 * t99;
t217 = t172 * t43;
t216 = t172 * t70;
t215 = t173 * t80;
t214 = t173 * t99;
t212 = t144 * t168;
t211 = t144 * t172;
t210 = t150 * t169;
t149 = t174 * t176 * t170;
t141 = qJDD(3) + t149;
t209 = t170 * t141;
t208 = t173 * t150;
t142 = qJDD(3) - t149;
t207 = t174 * t142;
t203 = qJD(4) - t150;
t197 = t174 * t83;
t196 = t174 * t213;
t195 = pkin(1) * t164 + pkin(6);
t8 = t168 * t13 + t172 * t14;
t24 = t169 * t44 + t173 * t45;
t92 = t170 * t104 - t153;
t93 = t174 * t104 + t192;
t64 = t170 * t92 + t174 * t93;
t183 = t169 * t45 - t173 * t44;
t76 = -t143 - t108;
t36 = t168 * t76 + t232;
t181 = pkin(4) * t36 - t13;
t91 = -t109 - t143;
t53 = t172 * t91 + t220;
t179 = pkin(4) * t53 - t14;
t178 = -pkin(1) * t165 - pkin(2) + t187;
t161 = t174 ^ 2;
t160 = t170 ^ 2;
t158 = t161 * t176;
t156 = t160 * t176;
t147 = -t158 - t226;
t146 = -t156 - t226;
t140 = t156 + t158;
t139 = (t160 + t161) * qJDD(1);
t138 = -0.2e1 * t154 + t198;
t135 = 0.2e1 * t194 + t199;
t119 = -t127 + t148;
t118 = -t148 + t227;
t116 = -t170 * t146 - t207;
t115 = t174 * t147 - t209;
t114 = t127 - t227;
t113 = -t127 - t148;
t107 = -t148 - t227;
t98 = t127 + t227;
t95 = -t109 + t143;
t94 = t108 - t143;
t89 = t203 * t129 + t182;
t87 = t106 + t120;
t85 = -t203 * t131 - t189;
t82 = t109 - t108;
t78 = -t169 * t113 + t214;
t77 = t173 * t113 + t218;
t74 = t173 * t107 - t229;
t73 = t169 * t107 + t228;
t67 = (t110 * t172 - t112 * t168) * t144;
t66 = (t110 * t168 + t112 * t172) * t144;
t65 = -t108 - t109;
t62 = -t112 * qJD(5) - t191;
t61 = t169 * t88 - t173 * t84;
t59 = t172 * t94 + t220;
t58 = -t168 * t95 + t232;
t57 = t168 * t94 - t216;
t56 = t172 * t95 + t233;
t55 = -t170 * t89 + t174 * t78;
t54 = -t168 * t91 + t216;
t52 = -t170 * t85 + t174 * t74;
t46 = (qJD(5) - t144) * t112 + t191;
t41 = t112 * t212 + t172 * t63;
t40 = -t112 * t211 + t168 * t63;
t39 = -t110 * t211 - t168 * t62;
t38 = -t110 * t212 + t172 * t62;
t37 = t172 * t76 - t233;
t31 = -t169 * t53 + t173 * t54;
t30 = t169 * t54 + t173 * t53;
t29 = -pkin(8) * t53 + t217;
t28 = t168 * t50 - t172 * t47;
t27 = -t168 * t231 - t172 * t46;
t25 = -t168 * t46 + t172 * t231;
t22 = -pkin(8) * t36 + t221;
t21 = -t169 * t36 + t173 * t37;
t20 = t169 * t37 + t173 * t36;
t18 = t170 * t231 + t174 * t31;
t17 = -pkin(4) * t231 + pkin(8) * t54 + t221;
t16 = -pkin(4) * t46 + pkin(8) * t37 - t217;
t15 = t170 * t46 + t174 * t21;
t11 = -t169 * t26 + t173 * t28;
t10 = t169 * t28 + t173 * t26;
t9 = t174 * t11 + t170 * t65;
t6 = -pkin(4) * t43 + pkin(8) * t8;
t5 = -pkin(8) * t26 - t7;
t4 = -pkin(4) * t65 + pkin(8) * t28 + t8;
t3 = t173 * t8 - t223;
t2 = t169 * t8 + t222;
t1 = t170 * t43 + t174 * t3;
t12 = [0, 0, 0, 0, 0, qJDD(1), t193, t186, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t165 * qJDD(1) - t164 * t176) + t190, pkin(1) * (-t164 * qJDD(1) - t165 * t176) - t205, 0, pkin(1) * (t164 * t205 + t165 * t190), t185 * t170, t174 * t135 + t170 * t138, t209 + t174 * (-t156 + t226), -t184 * t174, t170 * (t158 - t226) + t207, 0, -t174 * t103 + pkin(2) * t138 + pkin(6) * t115 + pkin(1) * (t164 * t115 + t165 * t138), t170 * t103 - pkin(2) * t135 + pkin(6) * t116 + pkin(1) * (t164 * t116 - t165 * t135), pkin(2) * t140 + pkin(6) * t139 + pkin(1) * (t164 * t139 + t165 * t140) + t64, -pkin(2) * t103 + pkin(6) * t64 + pkin(1) * (-t165 * t103 + t164 * t64), t170 * (t173 * t106 + t131 * t210) - t196, t170 * (-t169 * t87 + t173 * t85) - t174 * t114, t170 * (-t169 * t119 + t228) - t174 * t88, t170 * (-t129 * t208 + t169 * t180) + t196, t170 * (t173 * t118 + t218) + t174 * t84, t174 * t128 + t170 * (t129 * t173 - t131 * t169) * t150, t170 * (-pkin(7) * t73 + t219) + t174 * (-pkin(3) * t73 + t44) - pkin(2) * t73 + pkin(6) * t52 + pkin(1) * (t164 * t52 - t165 * t73), t170 * (-pkin(7) * t77 + t215) + t174 * (-pkin(3) * t77 + t45) - pkin(2) * t77 + pkin(6) * t55 + pkin(1) * (t164 * t55 - t165 * t77), -t170 * t183 + t195 * (-t170 * t98 + t174 * t61) + t178 * (-t169 * t84 - t173 * t88), t195 * (t170 * t80 + t174 * t24) + t178 * t183, t170 * (-t169 * t40 + t173 * t41) - t197, t170 * (-t169 * t25 + t173 * t27) - t174 * t82, t170 * (-t169 * t56 + t173 * t58) - t174 * t50, t170 * (-t169 * t38 + t173 * t39) + t197, t170 * (-t169 * t57 + t173 * t59) + t174 * t47, t170 * (-t169 * t66 + t173 * t67) + t174 * t125, t170 * (-pkin(7) * t20 - t169 * t16 + t173 * t22) + t174 * (-pkin(3) * t20 - t181) - pkin(2) * t20 + pkin(6) * t15 + pkin(1) * (t164 * t15 - t165 * t20), t170 * (-pkin(7) * t30 - t169 * t17 + t173 * t29) + t174 * (-pkin(3) * t30 - t179) - pkin(2) * t30 + pkin(6) * t18 + pkin(1) * (t164 * t18 - t165 * t30), t170 * (-pkin(7) * t10 - t169 * t4 + t173 * t5) + t174 * (-pkin(3) * t10 - t224) - pkin(2) * t10 + pkin(6) * t9 + pkin(1) * (-t165 * t10 + t164 * t9), t170 * (-pkin(7) * t2 - pkin(8) * t222 - t169 * t6) + t174 * (-pkin(3) * t2 - t225) - pkin(2) * t2 + pkin(6) * t1 + pkin(1) * (t164 * t1 - t165 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, 0, 0, 0, 0, 0, t174 * t141 + t170 * t147, -t170 * t142 + t174 * t146, 0, t170 * t93 - t174 * t92, 0, 0, 0, 0, 0, 0, t170 * t74 + t174 * t85, t170 * t78 + t174 * t89, t170 * t61 + t174 * t98, t170 * t24 - t174 * t80, 0, 0, 0, 0, 0, 0, t170 * t21 - t174 * t46, t170 * t31 - t174 * t231, t170 * t11 - t174 * t65, t170 * t3 - t174 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t156 - t158, t199, t149, t198, qJDD(3), -t92, -t93, 0, 0, t169 * t106 - t131 * t208, t169 * t85 + t173 * t87, t173 * t119 + t229, -t129 * t210 - t173 * t180, t169 * t118 - t214, (t129 * t169 + t131 * t173) * t150, pkin(3) * t85 + pkin(7) * t74 - t215, pkin(3) * t89 + pkin(7) * t78 + t219, pkin(3) * t98 + pkin(7) * t61 + t24, -pkin(3) * t80 + pkin(7) * t24, t169 * t41 + t173 * t40, t169 * t27 + t173 * t25, t169 * t58 + t173 * t56, t169 * t39 + t173 * t38, t169 * t59 + t173 * t57, t169 * t67 + t173 * t66, -pkin(3) * t46 + pkin(7) * t21 + t173 * t16 + t169 * t22, -pkin(3) * t231 + pkin(7) * t31 + t169 * t29 + t173 * t17, -pkin(3) * t65 + pkin(7) * t11 + t169 * t5 + t173 * t4, -pkin(3) * t43 + pkin(7) * t3 - pkin(8) * t223 + t173 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t114, t88, -t213, -t84, -t128, -t44, -t45, 0, 0, t83, t82, t50, -t83, -t47, -t125, t181, t179, t224, t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t82, t50, -t83, -t47, -t125, -t13, -t14, 0, 0;];
tauJ_reg = t12;
