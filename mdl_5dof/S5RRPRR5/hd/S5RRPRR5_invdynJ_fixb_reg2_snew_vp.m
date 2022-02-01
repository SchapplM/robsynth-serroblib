% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:56
% EndTime: 2022-01-20 11:03:06
% DurationCPUTime: 3.09s
% Computational Cost: add. (18064->272), mult. (25253->401), div. (0->0), fcn. (17906->10), ass. (0->190)
t262 = sin(qJ(1));
t263 = cos(qJ(1));
t217 = g(1) * t263 + g(2) * t262;
t181 = -qJD(1) ^ 2 * pkin(1) - t217;
t202 = sin(qJ(2));
t205 = cos(qJ(2));
t216 = g(1) * t262 - g(2) * t263;
t214 = qJDD(1) * pkin(1) + t216;
t157 = t205 * t181 + t202 * t214;
t196 = qJDD(1) + qJDD(2);
t239 = qJD(1) + qJD(2);
t238 = t239 ^ 2;
t269 = -pkin(2) * t238 + qJ(3) * t196 + 0.2e1 * qJD(3) * t239 + t157;
t200 = sin(qJ(5));
t198 = sin(pkin(9));
t199 = cos(pkin(9));
t204 = cos(qJ(4));
t226 = t204 * t239;
t201 = sin(qJ(4));
t227 = t201 * t239;
t169 = t198 * t227 - t199 * t226;
t171 = t198 * t226 + t199 * t227;
t203 = cos(qJ(5));
t134 = t169 * t203 + t171 * t200;
t136 = -t169 * t200 + t171 * t203;
t110 = t136 * t134;
t195 = qJDD(4) + qJDD(5);
t266 = -t110 + t195;
t273 = t200 * t266;
t155 = t171 * t169;
t265 = qJDD(4) - t155;
t272 = t201 * t265;
t271 = t203 * t266;
t270 = t204 * t265;
t207 = t198 ^ 2;
t209 = t199 ^ 2;
t224 = t209 * t238;
t180 = t207 * t238 + t224;
t221 = -t198 * g(3) + t199 * t269;
t251 = t199 * t196;
t116 = -pkin(3) * t224 + pkin(7) * t251 + t221;
t261 = t199 * g(3);
t212 = -t261 + (pkin(3) * t199 * t238 - t196 * pkin(7) - t269) * t198;
t84 = t204 * t116 + t201 * t212;
t197 = qJD(4) + qJD(5);
t131 = t197 * t134;
t252 = t198 * t196;
t121 = -t201 * t252 + t204 * t251;
t240 = t171 * qJD(4);
t151 = t121 - t240;
t166 = (t198 * t204 + t199 * t201) * t196;
t241 = t169 * qJD(4);
t153 = t166 - t241;
t99 = -qJD(5) * t134 + t151 * t200 + t153 * t203;
t268 = -t131 + t99;
t101 = t199 * t221 + t198 * (t198 * t269 + t261);
t83 = t201 * t116 - t204 * t212;
t264 = -t83 + (-t153 - t241) * pkin(8);
t132 = t134 ^ 2;
t133 = t136 ^ 2;
t167 = t169 ^ 2;
t168 = t171 ^ 2;
t194 = t197 ^ 2;
t52 = t201 * t84 - t204 * t83;
t260 = t198 * t52;
t156 = -t202 * t181 + t205 * t214;
t145 = -pkin(2) * t196 - qJ(3) * t238 + qJDD(3) - t156;
t120 = -pkin(3) * t251 - pkin(7) * t180 + t145;
t158 = qJD(4) * pkin(4) - pkin(8) * t171;
t78 = -pkin(4) * t151 - pkin(8) * t167 + t158 * t171 + t120;
t259 = t200 * t78;
t211 = pkin(4) * t265 + t264;
t55 = -t167 * pkin(4) + t151 * pkin(8) - qJD(4) * t158 + t84;
t36 = t200 * t55 - t203 * t211;
t257 = t203 * t55;
t37 = t200 * t211 + t257;
t19 = t200 * t37 - t203 * t36;
t258 = t201 * t19;
t256 = t203 * t78;
t255 = t204 * t19;
t254 = t197 * t200;
t253 = t197 * t203;
t107 = t110 + t195;
t250 = t200 * t107;
t249 = t201 * t120;
t148 = qJDD(4) + t155;
t248 = t201 * t148;
t247 = t203 * t107;
t245 = t204 * t120;
t244 = t204 * t148;
t243 = t205 * t196;
t242 = -pkin(2) * t145 + qJ(3) * t101;
t175 = t180 * t199;
t237 = pkin(2) * t251 - qJ(3) * t175 - t145 * t199;
t20 = t200 * t36 + t203 * t37;
t17 = -pkin(4) * t78 + pkin(8) * t20;
t7 = t20 * t201 + t255;
t8 = t20 * t204 - t258;
t4 = -t198 * t7 + t199 * t8;
t236 = t198 * (-pkin(7) * t7 - pkin(8) * t255 - t17 * t201) + t199 * (-pkin(3) * t78 + pkin(7) * t8 - pkin(8) * t258 + t17 * t204) - pkin(2) * t78 + qJ(3) * t4;
t228 = -t151 * t203 + t200 * t153;
t215 = (-qJD(5) + t197) * t136 - t228;
t74 = t131 + t99;
t45 = t200 * t215 - t203 * t74;
t47 = t200 * t74 + t203 * t215;
t25 = t201 * t47 + t204 * t45;
t26 = -t201 * t45 + t204 * t47;
t13 = -t198 * t25 + t199 * t26;
t14 = -pkin(8) * t45 - t19;
t97 = -t132 - t133;
t9 = -pkin(4) * t97 + pkin(8) * t47 + t20;
t235 = t198 * (-pkin(7) * t25 + t14 * t204 - t201 * t9) + t199 * (-pkin(3) * t97 + pkin(7) * t26 + t14 * t201 + t204 * t9) - pkin(2) * t97 + qJ(3) * t13;
t53 = t201 * t83 + t204 * t84;
t105 = -t194 - t132;
t57 = t105 * t200 + t271;
t58 = t105 * t203 - t273;
t41 = t201 * t58 + t204 * t57;
t42 = -t201 * t57 + t204 * t58;
t22 = -t198 * t41 + t199 * t42;
t69 = (qJD(5) + t197) * t136 + t228;
t33 = -pkin(4) * t69 + pkin(8) * t58 - t256;
t48 = -pkin(8) * t57 + t259;
t233 = t198 * (-pkin(7) * t41 - t201 * t33 + t204 * t48) + t199 * (-pkin(3) * t69 + pkin(7) * t42 + t201 * t48 + t204 * t33) - pkin(2) * t69 + qJ(3) * t22;
t127 = -t133 - t194;
t81 = t127 * t203 - t250;
t82 = -t127 * t200 - t247;
t50 = t201 * t82 + t204 * t81;
t51 = -t201 * t81 + t204 * t82;
t28 = -t198 * t50 + t199 * t51;
t34 = -pkin(4) * t268 + pkin(8) * t82 + t259;
t49 = -pkin(8) * t81 + t256;
t232 = t198 * (-pkin(7) * t50 - t201 * t34 + t204 * t49) + t199 * (-pkin(3) * t268 + pkin(7) * t51 + t201 * t49 + t204 * t34) - pkin(2) * t268 + qJ(3) * t28;
t113 = t121 * t201 - t166 * t204;
t114 = t121 * t204 + t166 * t201;
t125 = -t167 - t168;
t76 = -t113 * t198 + t114 * t199;
t231 = t198 * (-pkin(7) * t113 - t52) + t199 * (-pkin(3) * t125 + pkin(7) * t114 + t53) - pkin(2) * t125 + qJ(3) * t76;
t206 = qJD(4) ^ 2;
t146 = -t206 - t167;
t111 = t146 * t201 + t270;
t112 = t146 * t204 - t272;
t150 = -t121 + 0.2e1 * t240;
t65 = -t111 * t198 + t112 * t199;
t230 = t198 * (-pkin(7) * t111 + t249) + t199 * (-pkin(3) * t150 + pkin(7) * t112 - t245) - pkin(2) * t150 + qJ(3) * t65;
t161 = -t168 - t206;
t117 = t161 * t204 - t248;
t118 = -t161 * t201 - t244;
t152 = t166 - 0.2e1 * t241;
t92 = -t117 * t198 + t118 * t199;
t229 = t198 * (-pkin(7) * t117 + t245) + t199 * (-pkin(3) * t152 + pkin(7) * t118 + t249) - pkin(2) * t152 + qJ(3) * t92;
t190 = t207 * t196;
t191 = t209 * t196;
t178 = t191 + t190;
t222 = pkin(2) * t180 + qJ(3) * t178 + t101;
t174 = t180 * t198;
t219 = -pkin(2) * t252 + qJ(3) * t174 + t145 * t198;
t30 = t199 * t53 - t260;
t218 = -pkin(7) * t260 + qJ(3) * t30 - pkin(2) * t120 + t199 * (-pkin(3) * t120 + pkin(7) * t53);
t182 = 0.2e1 * t198 * t251;
t160 = -t168 + t206;
t159 = t167 - t206;
t129 = -t133 + t194;
t128 = t132 - t194;
t109 = t133 - t132;
t104 = (t198 * (-t169 * t204 + t171 * t201) + t199 * (-t169 * t201 - t171 * t204)) * qJD(4);
t103 = (-t134 * t203 + t136 * t200) * t197;
t102 = (-t134 * t200 - t136 * t203) * t197;
t98 = -qJD(5) * t136 - t228;
t94 = t198 * (t153 * t204 - t201 * t240) + t199 * (t153 * t201 + t204 * t240);
t93 = t198 * (-t151 * t201 + t204 * t241) + t199 * (t151 * t204 + t201 * t241);
t91 = t198 * (-t160 * t201 + t270) + t199 * (t160 * t204 + t272);
t90 = t198 * (t159 * t204 - t248) + t199 * (t159 * t201 + t244);
t88 = t128 * t203 - t250;
t87 = -t129 * t200 + t271;
t86 = t128 * t200 + t247;
t85 = t129 * t203 + t273;
t75 = t198 * (-t150 * t204 - t152 * t201) + t199 * (-t150 * t201 + t152 * t204);
t62 = -t136 * t254 + t203 * t99;
t61 = t136 * t253 + t200 * t99;
t60 = t134 * t253 - t200 * t98;
t59 = t134 * t254 + t203 * t98;
t46 = -t200 * t268 - t203 * t69;
t44 = -t200 * t69 + t203 * t268;
t39 = t198 * (-t102 * t201 + t103 * t204) + t199 * (t102 * t204 + t103 * t201);
t32 = t198 * (-t201 * t86 + t204 * t88) + t199 * (t201 * t88 + t204 * t86);
t31 = t198 * (-t201 * t85 + t204 * t87) + t199 * (t201 * t87 + t204 * t85);
t24 = t198 * (-t201 * t61 + t204 * t62) + t199 * (t201 * t62 + t204 * t61);
t23 = t198 * (-t201 * t59 + t204 * t60) + t199 * (t201 * t60 + t204 * t59);
t12 = t198 * (-t201 * t44 + t204 * t46) + t199 * (t201 * t46 + t204 * t44);
t1 = [0, 0, 0, 0, 0, qJDD(1), t216, t217, 0, 0, 0, 0, 0, 0, 0, t196, pkin(1) * (-t202 * t238 + t243) + t156, pkin(1) * (-t202 * t196 - t205 * t238) - t157, 0, pkin(1) * (t156 * t205 + t157 * t202), t190, t182, 0, t191, 0, 0, pkin(1) * (-t175 * t202 + t199 * t243) + t237, pkin(1) * (t174 * t202 - t198 * t243) + t219, pkin(1) * (t178 * t202 + t180 * t205) + t222, pkin(1) * (t101 * t202 - t145 * t205) + t242, t94, t75, t91, t93, t90, t104, pkin(1) * (-t150 * t205 + t202 * t65) + t230, pkin(1) * (-t152 * t205 + t202 * t92) + t229, pkin(1) * (-t125 * t205 + t202 * t76) + t231, pkin(1) * (-t120 * t205 + t202 * t30) + t218, t24, t12, t31, t23, t32, t39, pkin(1) * (t202 * t22 - t205 * t69) + t233, pkin(1) * (t202 * t28 - t205 * t268) + t232, pkin(1) * (t13 * t202 - t205 * t97) + t235, pkin(1) * (t202 * t4 - t205 * t78) + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t156, -t157, 0, 0, t190, t182, 0, t191, 0, 0, t237, t219, t222, t242, t94, t75, t91, t93, t90, t104, t230, t229, t231, t218, t24, t12, t31, t23, t32, t39, t233, t232, t235, t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, t252, -t180, t145, 0, 0, 0, 0, 0, 0, t150, t152, t125, t120, 0, 0, 0, 0, 0, 0, t69, t268, t97, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t168 - t167, t166, -t155, t121, qJDD(4), -t83, -t84, 0, 0, t110, t109, t74, -t110, t215, t195, pkin(4) * t57 - t36, -t257 - t200 * t264 + (-t200 * t265 + t81) * pkin(4), pkin(4) * t45, pkin(4) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t109, t74, -t110, t215, t195, -t36, -t37, 0, 0;];
tauJ_reg = t1;
