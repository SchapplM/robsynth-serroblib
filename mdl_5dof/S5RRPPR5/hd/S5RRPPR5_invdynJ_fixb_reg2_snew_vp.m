% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:58
% EndTime: 2019-12-31 19:30:07
% DurationCPUTime: 3.99s
% Computational Cost: add. (8483->358), mult. (20193->448), div. (0->0), fcn. (13600->8), ass. (0->203)
t186 = sin(qJ(2));
t189 = cos(qJ(2));
t184 = cos(pkin(8));
t183 = sin(pkin(8));
t219 = qJD(1) * t189;
t220 = qJD(1) * t186;
t151 = t183 * t220 - t184 * t219;
t153 = t183 * t219 + t184 * t220;
t232 = t153 * t151;
t252 = qJDD(2) + t232;
t228 = t183 * t252;
t150 = t153 ^ 2;
t190 = qJD(2) ^ 2;
t256 = -t150 - t190;
t78 = -t184 * t256 + t228;
t225 = t184 * t252;
t80 = t183 * t256 + t225;
t293 = pkin(6) * (t186 * t78 - t189 * t80);
t292 = pkin(2) * t78;
t291 = qJ(3) * t78;
t290 = qJ(3) * t80;
t251 = t151 ^ 2;
t135 = t251 - t190;
t288 = t186 * (-t184 * t135 + t228) - t189 * (t183 * t135 + t225);
t168 = t186 * qJDD(1);
t211 = qJD(1) * qJD(2);
t208 = t189 * t211;
t159 = t168 + t208;
t169 = t189 * qJDD(1);
t209 = t186 * t211;
t160 = t169 - t209;
t128 = t159 * t184 + t160 * t183;
t218 = qJD(2) * t151;
t101 = t128 + t218;
t127 = t159 * t183 - t160 * t184;
t147 = qJD(2) * t153;
t97 = t127 - t147;
t274 = -t101 * t184 - t183 * t97;
t287 = pkin(2) * t274;
t286 = qJ(3) * t274;
t273 = t101 * t183 - t184 * t97;
t94 = -t251 - t150;
t283 = -pkin(2) * t94 + qJ(3) * t273;
t282 = pkin(6) * (-t186 * t274 + t189 * t273) - pkin(1) * t94;
t253 = qJDD(2) - t232;
t227 = t183 * t253;
t255 = -t190 - t251;
t262 = t184 * t255 - t227;
t279 = qJ(3) * t262;
t112 = t184 * t253;
t264 = t183 * t255 + t112;
t278 = qJ(3) * t264;
t214 = qJD(3) * t153;
t275 = pkin(2) * t264 - 0.2e1 * t214;
t257 = t128 - t218;
t96 = t127 + t147;
t272 = t186 * (t183 * t257 + t184 * t96) - t189 * (-t183 * t96 + t184 * t257);
t271 = pkin(6) * (-t186 * t264 + t189 * t262) - pkin(1) * t96;
t136 = -t150 + t190;
t270 = t186 * (-t136 * t183 + t112) + t189 * (t136 * t184 + t227);
t185 = sin(qJ(5));
t175 = qJDD(2) - qJDD(5);
t188 = cos(qJ(5));
t109 = -t151 * t188 + t153 * t185;
t111 = t151 * t185 + t153 * t188;
t76 = t111 * t109;
t259 = -t76 - t175;
t267 = t185 * t259;
t266 = t188 * t259;
t265 = t257 * qJ(4);
t176 = qJD(2) - qJD(5);
t106 = t109 * t176;
t64 = -qJD(5) * t109 + t127 * t185 + t128 * t188;
t260 = t106 + t64;
t181 = t189 ^ 2;
t191 = qJD(1) ^ 2;
t197 = qJD(2) * pkin(2) - qJ(3) * t220;
t187 = sin(qJ(1));
t248 = cos(qJ(1));
t206 = t187 * g(1) - g(2) * t248;
t198 = qJDD(1) * pkin(1) + t206;
t93 = t160 * pkin(2) - qJDD(3) - t197 * t220 + (qJ(3) * t181 + pkin(6)) * t191 + t198;
t258 = -pkin(3) * t147 + t93;
t222 = t186 * t191;
t199 = g(1) * t248 + g(2) * t187;
t233 = qJDD(1) * pkin(6);
t156 = -pkin(1) * t191 - t199 + t233;
t224 = t186 * t156;
t88 = qJDD(2) * pkin(2) - t159 * qJ(3) - t224 + (pkin(2) * t222 + qJ(3) * t211 - g(3)) * t189;
t131 = -t186 * g(3) + t189 * t156;
t171 = t181 * t191;
t89 = -pkin(2) * t171 + t160 * qJ(3) - qJD(2) * t197 + t131;
t56 = -0.2e1 * qJD(3) * t151 + t183 * t88 + t184 * t89;
t254 = -t251 + t150;
t204 = -t127 * t188 + t185 * t128;
t47 = (qJD(5) + t176) * t111 + t204;
t107 = t109 ^ 2;
t108 = t111 ^ 2;
t174 = t176 ^ 2;
t250 = 2 * qJD(4);
t249 = pkin(3) + pkin(4);
t247 = pkin(3) * t183;
t246 = pkin(3) * t184;
t245 = t127 * pkin(3);
t243 = t183 * t93;
t241 = t184 * t93;
t133 = -qJD(2) * pkin(4) - pkin(7) * t153;
t193 = t258 + t265;
t32 = -t245 - t127 * pkin(4) - t251 * pkin(7) + (t250 + t133) * t153 + t193;
t239 = t185 * t32;
t72 = -t76 + t175;
t238 = t185 * t72;
t205 = t183 * t89 - t184 * t88;
t55 = t205 + 0.2e1 * t214;
t29 = t183 * t56 - t184 * t55;
t237 = t186 * t29;
t236 = t188 * t32;
t235 = t188 * t72;
t234 = qJ(4) * t184;
t231 = t176 * t185;
t230 = t176 * t188;
t167 = t189 * t222;
t223 = t186 * (qJDD(2) + t167);
t221 = t189 * (qJDD(2) - t167);
t217 = qJD(2) * t183;
t216 = qJD(2) * t184;
t113 = pkin(3) * t151 - qJ(4) * t153;
t213 = 0.2e1 * qJD(3) + t113;
t210 = t153 * t250;
t207 = -qJ(4) * t183 - pkin(2);
t30 = t183 * t55 + t184 * t56;
t196 = -qJDD(2) * pkin(3) - qJ(4) * t190 + qJDD(4) + t205;
t28 = -qJDD(2) * pkin(4) - t101 * pkin(7) + (pkin(4) * t151 + t213) * t153 + t196;
t200 = qJDD(2) * qJ(4) + qJD(2) * t250 - t113 * t151 + t56;
t37 = -pkin(3) * t190 + t200;
t31 = -pkin(4) * t251 + pkin(7) * t127 + qJD(2) * t133 + t37;
t14 = t185 * t31 - t188 * t28;
t130 = t189 * g(3) + t224;
t203 = t186 * t130 + t131 * t189;
t15 = t185 * t28 + t188 * t31;
t7 = -t14 * t188 + t15 * t185;
t8 = t185 * t14 + t188 * t15;
t38 = t213 * t153 + t196;
t195 = t186 * (t127 * t183 + t151 * t216) + t189 * (-t127 * t184 + t151 * t217);
t134 = t153 * t217;
t194 = t186 * t134 + (-t186 * t151 * t184 + t189 * (-t151 * t183 - t153 * t184)) * qJD(2);
t192 = t193 + t210;
t180 = t186 ^ 2;
t170 = t180 * t191;
t161 = t169 - 0.2e1 * t209;
t158 = t168 + 0.2e1 * t208;
t155 = t191 * pkin(6) + t198;
t104 = -t108 + t174;
t103 = t107 - t174;
t102 = -t108 - t174;
t75 = t108 - t107;
t71 = -t174 - t107;
t70 = (t109 * t188 - t111 * t185) * t176;
t69 = (-t109 * t185 - t111 * t188) * t176;
t63 = -qJD(5) * t111 - t204;
t62 = -t107 - t108;
t61 = t186 * (t128 * t184 - t134) + t189 * (t128 * t183 + t153 * t216);
t60 = t103 * t188 + t238;
t59 = -t104 * t185 + t266;
t58 = -t103 * t185 + t235;
t57 = -t104 * t188 - t267;
t54 = -t102 * t185 + t235;
t53 = t102 * t188 + t238;
t51 = -t106 + t64;
t46 = (qJD(5) - t176) * t111 + t204;
t45 = t111 * t231 + t188 * t64;
t44 = t111 * t230 - t185 * t64;
t43 = -t109 * t230 - t185 * t63;
t42 = t109 * t231 - t188 * t63;
t41 = t192 - t245;
t40 = t188 * t71 - t267;
t39 = t185 * t71 + t266;
t36 = (-t127 - t96) * pkin(3) + t192;
t35 = t210 - t245 + t258 + 0.2e1 * t265;
t34 = -qJ(4) * t94 + t38;
t33 = (-t190 - t94) * pkin(3) + t200;
t27 = t183 * t53 + t184 * t54;
t26 = t183 * t54 - t184 * t53;
t24 = t185 * t51 - t188 * t47;
t23 = -t185 * t260 - t188 * t46;
t22 = -t185 * t47 - t188 * t51;
t21 = t185 * t46 - t188 * t260;
t20 = t183 * t39 + t184 * t40;
t19 = t183 * t40 - t184 * t39;
t17 = t183 * t37 - t184 * t38;
t16 = -pkin(7) * t53 + qJ(4) * t260 + t236;
t13 = -pkin(7) * t39 + qJ(4) * t46 + t239;
t12 = t183 * t22 + t184 * t24;
t11 = t183 * t24 - t184 * t22;
t10 = -pkin(7) * t54 + t249 * t260 - t239;
t9 = -pkin(7) * t40 + t249 * t46 + t236;
t6 = -pkin(7) * t7 + qJ(4) * t32;
t5 = -pkin(7) * t22 + qJ(4) * t62 - t7;
t4 = -pkin(7) * t8 + t249 * t32;
t3 = -pkin(7) * t24 + t249 * t62 - t8;
t2 = t183 * t7 + t184 * t8;
t1 = t183 * t8 - t184 * t7;
t18 = [0, 0, 0, 0, 0, qJDD(1), t206, t199, 0, 0, (t159 + t208) * t186, t158 * t189 + t161 * t186, t223 + t189 * (-t170 + t190), (t160 - t209) * t189, t186 * (t171 - t190) + t221, 0, t189 * t155 + pkin(1) * t161 + pkin(6) * (t189 * (-t171 - t190) - t223), -t186 * t155 - pkin(1) * t158 + pkin(6) * (-t221 - t186 * (-t170 - t190)), pkin(1) * (t170 + t171) + (t180 + t181) * t233 + t203, pkin(1) * t155 + pkin(6) * t203, t61, -t272, t270, t195, -t288, t194, t186 * (-t243 - t278) + t189 * (-pkin(2) * t96 + t241 + t279) + t271, t186 * (-t241 + t291) + t189 * (-pkin(2) * t257 - t243 - t290) - pkin(1) * t257 + t293, t186 * (-t29 - t286) + t189 * (t283 + t30) + t282, -qJ(3) * t237 + t189 * (pkin(2) * t93 + qJ(3) * t30) + pkin(1) * t93 + pkin(6) * (t189 * t30 - t237), t61, t270, t272, t194, t288, t195, t186 * (-t183 * t36 - t234 * t96 - t278) + t189 * (t184 * t36 + t207 * t96 + t279) + t271, t186 * (-t183 * t33 + t184 * t34 - t286) + t189 * (t183 * t34 + t184 * t33 + t283) + t282, t186 * (t184 * t35 - t291) + t189 * (t183 * t35 + t290) - t293 + (-t186 * t247 + t189 * (pkin(2) + t246) + pkin(1)) * t257, (t186 * (t234 - t247) + t189 * (-t207 + t246) + pkin(1)) * t41 + (pkin(6) + qJ(3)) * (-t186 * t17 + t189 * (t183 * t38 + t184 * t37)), t186 * (-t183 * t44 + t184 * t45) + t189 * (t183 * t45 + t184 * t44), t186 * (-t183 * t21 + t184 * t23) + t189 * (t183 * t23 + t184 * t21), t186 * (-t183 * t57 + t184 * t59) + t189 * (t183 * t59 + t184 * t57), t186 * (-t183 * t42 + t184 * t43) + t189 * (t183 * t43 + t184 * t42), t186 * (-t183 * t58 + t184 * t60) + t189 * (t183 * t60 + t184 * t58), t186 * (-t183 * t69 + t184 * t70) + t189 * (t183 * t70 + t184 * t69), t186 * (-qJ(3) * t19 + t13 * t184 - t183 * t9) + t189 * (pkin(2) * t46 + qJ(3) * t20 + t13 * t183 + t184 * t9) + pkin(1) * t46 + pkin(6) * (-t186 * t19 + t189 * t20), t186 * (-qJ(3) * t26 - t10 * t183 + t16 * t184) + t189 * (pkin(2) * t260 + qJ(3) * t27 + t10 * t184 + t16 * t183) + pkin(1) * t260 + pkin(6) * (-t186 * t26 + t189 * t27), t186 * (-qJ(3) * t11 - t183 * t3 + t184 * t5) + t189 * (pkin(2) * t62 + qJ(3) * t12 + t183 * t5 + t184 * t3) + pkin(1) * t62 + pkin(6) * (-t11 * t186 + t12 * t189), t186 * (-qJ(3) * t1 - t183 * t4 + t184 * t6) + t189 * (pkin(2) * t32 + qJ(3) * t2 + t183 * t6 + t184 * t4) + pkin(1) * t32 + pkin(6) * (-t1 * t186 + t189 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, -t171 + t170, t168, t167, t169, qJDD(2), -t130, -t131, 0, 0, t232, t254, t101, -t232, -t97, qJDD(2), -t205 + t275, -t56 - t292, t287, pkin(2) * t29, t232, t101, -t254, qJDD(2), t97, -t232, pkin(3) * t253 + qJ(4) * t255 - t153 * t113 - t196 + t275, -pkin(3) * t101 - qJ(4) * t97 + t287, t292 + qJ(4) * t252 + (-t256 - t190) * pkin(3) + t200, pkin(2) * t17 - pkin(3) * t38 + qJ(4) * t37, -t76, -t75, -t51, t76, t47, t175, pkin(2) * t19 + qJ(4) * t40 - t249 * t39 + t14, pkin(2) * t26 + qJ(4) * t54 - t249 * t53 + t15, pkin(2) * t11 + qJ(4) * t24 - t22 * t249, pkin(2) * t1 + qJ(4) * t8 - t249 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t257, t94, -t93, 0, 0, 0, 0, 0, 0, t96, t94, -t257, -t41, 0, 0, 0, 0, 0, 0, -t46, -t260, -t62, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t101, t256, t38, 0, 0, 0, 0, 0, 0, t39, t53, t22, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t75, t51, -t76, -t47, -t175, -t14, -t15, 0, 0;];
tauJ_reg = t18;
