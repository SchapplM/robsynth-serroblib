% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:44
% EndTime: 2021-01-15 23:10:58
% DurationCPUTime: 3.21s
% Computational Cost: add. (4948->357), mult. (11962->473), div. (0->0), fcn. (8877->14), ass. (0->202)
t173 = qJ(2) + qJ(3);
t164 = pkin(9) + t173;
t153 = sin(t164);
t178 = sin(qJ(1));
t182 = cos(qJ(1));
t216 = g(1) * t182 + g(2) * t178;
t289 = t216 * t153;
t165 = sin(t173);
t166 = cos(t173);
t288 = -g(3) * t166 + t216 * t165;
t154 = cos(t164);
t276 = g(3) * t154;
t168 = qJDD(2) + qJDD(3);
t174 = sin(pkin(9));
t176 = sin(qJ(3));
t180 = cos(qJ(3));
t181 = cos(qJ(2));
t243 = qJD(1) * t181;
t177 = sin(qJ(2));
t244 = qJD(1) * t177;
t106 = -t176 * t243 - t180 * t244;
t279 = pkin(6) + pkin(7);
t132 = t279 * t177;
t120 = qJD(1) * t132;
t261 = qJD(2) * pkin(2);
t111 = -t120 + t261;
t133 = t279 * t181;
t122 = qJD(1) * t133;
t251 = t180 * t122;
t210 = -t176 * t111 - t251;
t238 = qJD(1) * qJD(2);
t228 = t181 * t238;
t237 = t177 * qJDD(1);
t284 = t228 + t237;
t82 = qJDD(2) * pkin(2) - t279 * t284;
t229 = t177 * t238;
t236 = t181 * qJDD(1);
t84 = t279 * (-t229 + t236);
t189 = t210 * qJD(3) - t176 * t84 + t180 * t82;
t170 = qJD(2) + qJD(3);
t230 = t180 * t243;
t231 = t176 * t244;
t61 = qJD(3) * t230 - t170 * t231 + t176 * t236 + t284 * t180;
t19 = t168 * pkin(3) - t61 * qJ(4) + t106 * qJD(4) + t189;
t104 = -t230 + t231;
t242 = qJD(3) * t176;
t281 = (qJD(3) * t111 + t84) * t180 - t122 * t242 + t176 * t82;
t212 = t176 * t237 - t180 * t236;
t116 = t176 * t181 + t180 * t177;
t81 = t170 * t116;
t62 = qJD(1) * t81 + t212;
t23 = -t62 * qJ(4) - t104 * qJD(4) + t281;
t257 = cos(pkin(9));
t5 = -t174 * t23 + t257 * t19;
t3 = -t168 * pkin(4) - t5;
t287 = t3 + t276;
t167 = t181 * pkin(2);
t265 = pkin(1) + t167;
t155 = t174 * pkin(3) + pkin(8);
t202 = -t174 * t104 - t257 * t106;
t274 = t106 * pkin(3);
t73 = -t257 * t104 + t174 * t106;
t43 = pkin(4) * t202 - t73 * pkin(8) - t274;
t69 = qJD(5) - t73;
t283 = (qJD(5) * t155 + t43) * t69;
t162 = pkin(2) * t244;
t160 = t180 * pkin(2) + pkin(3);
t226 = t257 * t176;
t99 = pkin(2) * t226 + t174 * t160;
t94 = pkin(8) + t99;
t282 = (qJD(5) * t94 + t162 + t43) * t69;
t100 = t106 * qJ(4);
t107 = t176 * t122;
t221 = t180 * t111 - t107;
t59 = t100 + t221;
t247 = -t176 * t132 + t180 * t133;
t175 = sin(qJ(5));
t179 = cos(qJ(5));
t36 = -t174 * t62 + t257 * t61;
t65 = t175 * t170 + t179 * t202;
t22 = t65 * qJD(5) - t179 * t168 + t175 * t36;
t232 = qJD(2) * t279;
t121 = t177 * t232;
t123 = t181 * t232;
t188 = -t247 * qJD(3) + t176 * t121 - t180 * t123;
t115 = t176 * t177 - t180 * t181;
t80 = t170 * t115;
t185 = t80 * qJ(4) - t116 * qJD(4) + t188;
t241 = qJD(3) * t180;
t198 = -t180 * t121 - t176 * t123 - t132 * t241 - t133 * t242;
t39 = -t81 * qJ(4) - t115 * qJD(4) + t198;
t15 = t174 * t185 + t257 * t39;
t131 = t265 * qJD(1);
t83 = t104 * pkin(3) + qJD(4) - t131;
t38 = -pkin(4) * t73 - pkin(8) * t202 + t83;
t6 = t174 * t19 + t257 * t23;
t227 = t168 * pkin(8) + qJD(5) * t38 + t6;
t256 = t104 * qJ(4);
t60 = -t210 - t256;
t260 = t174 * t60;
t54 = t170 * pkin(3) + t59;
t28 = t257 * t54 - t260;
t26 = -t170 * pkin(4) - t28;
t35 = t174 * t61 + t257 * t62;
t32 = qJDD(5) + t35;
t77 = t257 * t115 + t174 * t116;
t78 = -t174 * t115 + t257 * t116;
t86 = t115 * pkin(3) - t265;
t44 = t77 * pkin(4) - t78 * pkin(8) + t86;
t208 = -t180 * t132 - t176 * t133;
t194 = -t116 * qJ(4) + t208;
t70 = -t115 * qJ(4) + t247;
t46 = t174 * t194 + t257 * t70;
t52 = -t174 * t81 - t257 * t80;
t280 = -(qJD(5) * t44 + t15) * t69 - t227 * t77 + t26 * t52 + t3 * t78 - t46 * t32;
t141 = g(3) * t153;
t273 = t26 * t73;
t272 = t26 * t78;
t271 = t44 * t32;
t63 = -t179 * t170 + t175 * t202;
t270 = t63 * t69;
t269 = t65 * t69;
t268 = t65 * t202;
t267 = t69 * t202;
t266 = t202 * t63;
t55 = t257 * t60;
t29 = t174 * t54 + t55;
t209 = t176 * t120 - t251;
t195 = t209 + t256;
t262 = pkin(2) * qJD(3);
t248 = -t180 * t120 - t107;
t66 = t100 + t248;
t264 = -t174 * t66 + t257 * t195 + (t174 * t180 + t226) * t262;
t254 = t174 * t176;
t263 = -t174 * t195 - t257 * t66 + (t257 * t180 - t254) * t262;
t259 = t175 * t32;
t239 = qJD(5) * t179;
t240 = qJD(5) * t175;
t21 = t175 * t168 + t170 * t239 + t179 * t36 - t202 * t240;
t258 = t21 * t175;
t255 = t106 * t104;
t253 = t178 * t175;
t252 = t178 * t179;
t250 = t182 * t175;
t249 = t182 * t179;
t246 = pkin(3) * t166 + t167;
t171 = t177 ^ 2;
t245 = -t181 ^ 2 + t171;
t163 = t177 * t261;
t235 = t78 * t240;
t75 = t81 * pkin(3) + t163;
t101 = pkin(2) * t229 - qJDD(1) * t265;
t48 = t62 * pkin(3) + qJDD(4) + t101;
t11 = t35 * pkin(4) - t36 * pkin(8) + t48;
t27 = t170 * pkin(8) + t29;
t225 = qJD(5) * t27 - t11;
t219 = t179 * t69;
t215 = g(1) * t178 - g(2) * t182;
t214 = t202 * t29 + t28 * t73;
t213 = t32 * t78 + t52 * t69;
t13 = t175 * t38 + t179 * t27;
t211 = t13 * t202 + t287 * t175 + t26 * t239;
t12 = -t175 * t27 + t179 * t38;
t207 = -t12 * t202 + t179 * t289 + t26 * t240;
t206 = t179 * t32 + (t175 * t73 - t240) * t69;
t205 = -t227 + t141;
t203 = -0.2e1 * pkin(1) * t238 - pkin(6) * qJDD(2);
t34 = t257 * t59 - t260;
t197 = -t155 * t32 + t34 * t69 - t273;
t196 = t216 * t154 - t83 * t73 + t141 - t6;
t98 = -pkin(2) * t254 + t257 * t160;
t193 = -t263 * t69 - t94 * t32 - t273;
t183 = qJD(2) ^ 2;
t192 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t183 + t215;
t184 = qJD(1) ^ 2;
t191 = pkin(1) * t184 - pkin(6) * qJDD(1) + t216;
t190 = -t202 * t83 - t276 + t289 + t5;
t187 = g(3) * t165 - t131 * t104 + t216 * t166 - t281;
t186 = -t131 * t106 + t189 + t288;
t169 = -qJ(4) - t279;
t156 = -t257 * pkin(3) - pkin(4);
t119 = pkin(1) + t246;
t93 = -pkin(4) - t98;
t92 = t154 * t249 + t253;
t91 = -t154 * t250 + t252;
t90 = -t154 * t252 + t250;
t89 = t154 * t253 + t249;
t85 = t162 - t274;
t67 = -t104 ^ 2 + t106 ^ 2;
t51 = -t174 * t80 + t257 * t81;
t50 = -t212 + (-qJD(1) * t116 - t106) * t170;
t49 = t104 * t170 + t61;
t45 = t174 * t70 - t257 * t194;
t33 = t174 * t59 + t55;
t18 = t51 * pkin(4) - t52 * pkin(8) + t75;
t14 = t174 * t39 - t257 * t185;
t10 = t179 * t11;
t9 = t65 * t219 + t258;
t8 = t69 * t219 + t259 - t268;
t7 = t206 + t266;
t1 = (t21 - t270) * t179 + (-t22 - t269) * t175;
t2 = [qJDD(1), t215, t216, t171 * qJDD(1) + 0.2e1 * t177 * t228, 0.2e1 * t177 * t236 - 0.2e1 * t245 * t238, qJDD(2) * t177 + t183 * t181, qJDD(2) * t181 - t183 * t177, 0, t177 * t203 + t181 * t192, -t177 * t192 + t181 * t203, t106 * t80 + t61 * t116, t80 * t104 + t106 * t81 - t61 * t115 - t116 * t62, t116 * t168 - t80 * t170, -t115 * t168 - t81 * t170, 0, t101 * t115 + t104 * t163 - t131 * t81 + t166 * t215 + t168 * t208 + t170 * t188 - t265 * t62, t101 * t116 - t106 * t163 + t131 * t80 - t165 * t215 - t168 * t247 - t170 * t198 - t265 * t61, -t14 * t170 + t154 * t215 - t45 * t168 + t86 * t35 + t48 * t77 + t83 * t51 - t73 * t75, -t15 * t170 - t153 * t215 - t46 * t168 + t202 * t75 + t86 * t36 + t48 * t78 + t83 * t52, t14 * t202 + t15 * t73 - t28 * t52 - t29 * t51 - t46 * t35 + t45 * t36 - t5 * t78 - t6 * t77 - t216, t6 * t46 + t29 * t15 - t5 * t45 - t28 * t14 + t48 * t86 + t83 * t75 - g(1) * (-t178 * t119 - t182 * t169) - g(2) * (t182 * t119 - t178 * t169), -t65 * t235 + (t21 * t78 + t52 * t65) * t179, (-t175 * t65 - t179 * t63) * t52 + (-t258 - t179 * t22 + (t175 * t63 - t179 * t65) * qJD(5)) * t78, t179 * t213 + t21 * t77 - t235 * t69 + t65 * t51, -t239 * t69 * t78 - t175 * t213 - t22 * t77 - t63 * t51, t32 * t77 + t69 * t51, -g(1) * t90 - g(2) * t92 + t10 * t77 + t12 * t51 + t14 * t63 + t45 * t22 + (t18 * t69 + t271 + (-t27 * t77 - t46 * t69 + t272) * qJD(5)) * t179 + t280 * t175, -g(1) * t89 - g(2) * t91 - t13 * t51 + t14 * t65 + t45 * t21 + (-(-qJD(5) * t46 + t18) * t69 - t271 + t225 * t77 - qJD(5) * t272) * t175 + t280 * t179; 0, 0, 0, -t177 * t184 * t181, t245 * t184, t237, t236, qJDD(2), -g(3) * t181 + t177 * t191, g(3) * t177 + t181 * t191, -t255, t67, t49, t50, t168, -t209 * t170 + (-t104 * t244 + t180 * t168 - t170 * t242) * pkin(2) + t186, t248 * t170 + (t106 * t244 - t176 * t168 - t170 * t241) * pkin(2) + t187, t98 * t168 - t264 * t170 + t73 * t85 + t190, -t99 * t168 - t263 * t170 - t202 * t85 + t196, t202 * t264 + t263 * t73 - t99 * t35 - t98 * t36 + t214, t6 * t99 + t5 * t98 - t83 * t85 - g(3) * t246 + t263 * t29 - t264 * t28 - t216 * (-t177 * pkin(2) - pkin(3) * t165), t9, t1, t8, t7, -t267, t93 * t22 + t264 * t63 + (-t287 - t282) * t179 + t193 * t175 + t207, t93 * t21 + t264 * t65 + t193 * t179 + (-t289 + t282) * t175 + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t67, t49, t50, t168, -t170 * t210 + t186, t170 * t221 + t187, t33 * t170 + (-t106 * t73 + t257 * t168) * pkin(3) + t190, t34 * t170 + (t106 * t202 - t168 * t174) * pkin(3) + t196, -t33 * t202 - t34 * t73 + (-t174 * t35 - t257 * t36) * pkin(3) + t214, t28 * t33 - t29 * t34 + (t106 * t83 + t174 * t6 + t257 * t5 + t288) * pkin(3), t9, t1, t8, t7, -t267, t156 * t22 - t33 * t63 + t197 * t175 + (-t287 - t283) * t179 + t207, t156 * t21 - t33 * t65 + t197 * t179 + (-t289 + t283) * t175 + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170 * t202 + t35, t73 * t170 + t36, -t202 ^ 2 - t73 ^ 2, t202 * t28 - t29 * t73 - t215 + t48, 0, 0, 0, 0, 0, t206 - t266, -t179 * t69 ^ 2 - t259 - t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t63 ^ 2 + t65 ^ 2, t21 + t270, -t22 + t269, t32, -g(1) * t91 + g(2) * t89 + t13 * t69 + t175 * t205 - t239 * t27 - t26 * t65 + t10, g(1) * t92 - g(2) * t90 + t12 * t69 + t175 * t225 + t179 * t205 + t26 * t63;];
tau_reg = t2;
