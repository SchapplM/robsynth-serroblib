% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:51
% EndTime: 2019-12-31 20:55:57
% DurationCPUTime: 2.98s
% Computational Cost: add. (5081->390), mult. (12310->470), div. (0->0), fcn. (8832->12), ass. (0->204)
t174 = qJD(2) + qJD(3);
t182 = cos(qJ(2));
t282 = cos(qJ(3));
t237 = t282 * t182;
t224 = qJD(1) * t237;
t179 = sin(qJ(3));
t180 = sin(qJ(2));
t247 = qJD(1) * t180;
t236 = t179 * t247;
t105 = -t224 + t236;
t116 = t179 * t182 + t282 * t180;
t107 = t116 * qJD(1);
t178 = sin(pkin(8));
t257 = cos(pkin(8));
t66 = t257 * t105 + t107 * t178;
t258 = t66 * t174;
t251 = t179 * t180;
t213 = t174 * t251;
t228 = qJDD(1) * t282;
t243 = t182 * qJDD(1);
t225 = -t174 * t224 - t179 * t243 - t180 * t228;
t59 = qJD(1) * t213 + t225;
t244 = t180 * qJDD(1);
t211 = t179 * t244 - t182 * t228;
t81 = t174 * t116;
t60 = qJD(1) * t81 + t211;
t31 = -t178 * t60 - t257 * t59;
t22 = t31 + t258;
t177 = qJ(2) + qJ(3);
t169 = sin(t177);
t170 = cos(t177);
t181 = sin(qJ(1));
t183 = cos(qJ(1));
t221 = g(1) * t183 + g(2) * t181;
t289 = -g(3) * t170 + t169 * t221;
t206 = -t178 * t105 + t107 * t257;
t288 = t66 * t206;
t271 = t206 ^ 2;
t283 = pkin(7) + pkin(6);
t255 = t105 * qJ(4);
t134 = t283 * t180;
t120 = qJD(1) * t134;
t135 = t283 * t182;
t122 = qJD(1) * t135;
t238 = t282 * t122;
t78 = t179 * t120 - t238;
t198 = t78 + t255;
t235 = qJD(3) * t282;
t242 = t178 * t179 * pkin(2);
t101 = t107 * qJ(4);
t108 = t179 * t122;
t79 = -t282 * t120 - t108;
t61 = -t101 + t79;
t264 = -qJD(3) * t242 - t178 * t198 + (pkin(2) * t235 - t61) * t257;
t263 = qJD(2) * pkin(2);
t113 = -t120 + t263;
t71 = t282 * t113 - t108;
t57 = -t101 + t71;
t87 = -t179 * t134 + t282 * t135;
t168 = pkin(8) + t177;
t155 = sin(t168);
t156 = cos(t168);
t212 = t156 * pkin(4) + t155 * qJ(5);
t286 = g(1) * t181 - g(2) * t183;
t172 = qJDD(2) + qJDD(3);
t285 = -t172 * pkin(4) + qJDD(5);
t72 = t179 * t113 + t238;
t245 = qJD(1) * qJD(2);
t232 = t182 * t245;
t83 = qJDD(2) * pkin(2) + t283 * (-t232 - t244);
t233 = t180 * t245;
t85 = t283 * (-t233 + t243);
t39 = -qJD(3) * t72 - t179 * t85 + t282 * t83;
t14 = t172 * pkin(3) + t59 * qJ(4) - t107 * qJD(4) + t39;
t246 = qJD(3) * t179;
t226 = -t113 * t235 + t122 * t246 - t179 * t83 - t282 * t85;
t17 = -qJ(4) * t60 - qJD(4) * t105 - t226;
t267 = -t257 * t14 + t178 * t17;
t208 = -g(3) * t156 + t221 * t155 - t267;
t273 = t182 * pkin(2);
t162 = pkin(1) + t273;
t133 = t162 * qJD(1);
t84 = pkin(3) * t105 + qJD(4) - t133;
t33 = pkin(4) * t66 - qJ(5) * t206 + t84;
t190 = -t206 * t33 + t208 - t285;
t86 = -t282 * t134 - t179 * t135;
t197 = -t116 * qJ(4) + t86;
t115 = -t237 + t251;
t64 = -qJ(4) * t115 + t87;
t43 = t178 * t197 + t257 * t64;
t240 = qJD(2) * t283;
t121 = t180 * t240;
t123 = t182 * t240;
t51 = -t87 * qJD(3) + t179 * t121 - t282 * t123;
t80 = -qJD(2) * t237 - t182 * t235 + t213;
t188 = t80 * qJ(4) - t116 * qJD(4) + t51;
t50 = -t282 * t121 - t179 * t123 - t134 * t235 - t135 * t246;
t34 = -qJ(4) * t81 - qJD(4) * t115 + t50;
t9 = t178 * t188 + t257 * t34;
t284 = t155 * t286 + t172 * t43 + t174 * t9;
t281 = pkin(3) * t107;
t280 = pkin(3) * t169;
t279 = pkin(4) * t155;
t274 = g(3) * t182;
t58 = t72 - t255;
t53 = t257 * t58;
t27 = t178 * t57 + t53;
t272 = t27 * t206;
t269 = t66 ^ 2;
t4 = t178 * t14 + t257 * t17;
t52 = pkin(3) * t174 + t57;
t26 = t178 * t52 + t53;
t266 = qJD(5) + t264;
t227 = t257 * t179;
t265 = -t178 * t61 + t198 * t257 + (t282 * t178 + t227) * qJD(3) * pkin(2);
t260 = t178 * t58;
t28 = t257 * t57 - t260;
t262 = t174 * t28;
t261 = t174 * t206;
t259 = t27 * t174;
t256 = pkin(6) * qJDD(1);
t254 = t107 * t105;
t253 = t156 * t181;
t252 = t156 * t183;
t250 = qJD(5) - t28;
t161 = t282 * pkin(2) + pkin(3);
t100 = pkin(2) * t227 + t178 * t161;
t175 = t180 ^ 2;
t176 = t182 ^ 2;
t249 = t175 - t176;
t248 = t175 + t176;
t165 = t180 * t263;
t186 = qJD(1) ^ 2;
t241 = t180 * t186 * t182;
t239 = t265 * t206;
t70 = pkin(3) * t81 + t165;
t234 = t265 * t174;
t159 = pkin(3) * t170;
t231 = t159 + t273;
t124 = -pkin(2) * t180 - t280;
t230 = t124 - t279;
t30 = -t178 * t59 + t257 * t60;
t223 = t180 * t232;
t222 = -t279 - t280;
t25 = t257 * t52 - t260;
t23 = -t174 * pkin(4) + qJD(5) - t25;
t24 = qJ(5) * t174 + t26;
t219 = t206 * t24 + t23 * t66;
t218 = t206 * t26 - t25 * t66;
t47 = -t178 * t80 + t257 * t81;
t74 = t115 * t257 + t116 * t178;
t217 = t30 * t74 + t47 * t66;
t216 = -t269 - t271;
t19 = -t269 + t271;
t215 = g(1) * t252 + g(2) * t253 + g(3) * t155 - t4;
t210 = t172 * t74 + t174 * t47;
t89 = pkin(3) * t115 - t162;
t209 = t30 + t261;
t21 = -t30 + t261;
t207 = -0.2e1 * pkin(1) * t245 - pkin(6) * qJDD(2);
t102 = pkin(2) * t233 - qJDD(1) * t162;
t40 = pkin(4) * t206 + qJ(5) * t66 + t281;
t203 = t66 * t84 + t215;
t99 = t161 * t257 - t242;
t163 = t172 * qJ(5);
t202 = -t33 * t66 + t163 - t215;
t201 = -t206 * t84 + t208;
t199 = -t31 + t258;
t42 = t178 * t64 - t197 * t257;
t8 = t178 * t34 - t188 * t257;
t196 = g(1) * t253 - g(2) * t252 - t172 * t42 - t174 * t8;
t48 = -t178 * t81 - t257 * t80;
t75 = -t178 * t115 + t116 * t257;
t195 = t206 * t47 + t30 * t75 + t31 * t74 + t48 * t66;
t194 = g(3) * t169 - t133 * t105 + t221 * t170 + t226;
t185 = qJD(2) ^ 2;
t193 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t185 + t286;
t192 = pkin(1) * t186 + t221 - t256;
t44 = pkin(3) * t60 + qJDD(4) + t102;
t191 = t206 * t8 - t43 * t30 + t31 * t42 - t9 * t66 - t221;
t5 = pkin(4) * t30 - qJ(5) * t31 - qJD(5) * t206 + t44;
t189 = t133 * t107 + t289 + t39;
t173 = -qJ(4) - t283;
t166 = t174 * qJD(5);
t164 = pkin(2) * t247;
t157 = -pkin(3) * t257 - pkin(4);
t154 = pkin(3) * t178 + qJ(5);
t126 = qJ(5) * t252;
t125 = qJ(5) * t253;
t119 = pkin(1) + t231;
t112 = t183 * t119;
t95 = -pkin(4) - t99;
t94 = qJ(5) + t100;
t88 = t164 + t281;
t62 = -t105 ^ 2 + t107 ^ 2;
t45 = -t225 + (t105 - t236) * t174;
t41 = pkin(4) * t74 - qJ(5) * t75 + t89;
t37 = t164 + t40;
t29 = t172 * t75 + t174 * t48;
t10 = pkin(4) * t47 - qJ(5) * t48 - qJD(5) * t75 + t70;
t7 = t206 * t48 + t31 * t75;
t2 = t267 + t285;
t1 = t163 + t166 + t4;
t3 = [0, 0, 0, 0, 0, qJDD(1), t286, t221, 0, 0, qJDD(1) * t175 + 0.2e1 * t223, 0.2e1 * t180 * t243 - 0.2e1 * t245 * t249, qJDD(2) * t180 + t182 * t185, qJDD(1) * t176 - 0.2e1 * t223, qJDD(2) * t182 - t180 * t185, 0, t180 * t207 + t182 * t193, -t180 * t193 + t182 * t207, 0.2e1 * t248 * t256 - t221, -g(1) * (-pkin(1) * t181 + pkin(6) * t183) - g(2) * (pkin(1) * t183 + pkin(6) * t181) + (pkin(6) ^ 2 * t248 + pkin(1) ^ 2) * qJDD(1), -t107 * t80 - t116 * t59, t105 * t80 - t107 * t81 + t115 * t59 - t116 * t60, t116 * t172 - t174 * t80, t105 * t81 + t115 * t60, -t115 * t172 - t174 * t81, 0, t102 * t115 + t105 * t165 - t133 * t81 - t162 * t60 + t170 * t286 + t172 * t86 + t174 * t51, t102 * t116 + t107 * t165 + t133 * t80 + t162 * t59 - t169 * t286 - t172 * t87 - t174 * t50, -t105 * t50 - t107 * t51 + t115 * t226 - t116 * t39 + t59 * t86 - t60 * t87 + t71 * t80 - t72 * t81 - t221, -t226 * t87 + t72 * t50 + t39 * t86 + t71 * t51 - t102 * t162 - t133 * t165 - g(1) * (-t162 * t181 + t183 * t283) - g(2) * (t162 * t183 + t181 * t283), t7, -t195, t29, t217, -t210, 0, t30 * t89 + t44 * t74 + t47 * t84 + t66 * t70 + t196, t206 * t70 + t31 * t89 + t44 * t75 + t48 * t84 - t284, -t25 * t48 - t26 * t47 + t267 * t75 - t4 * t74 + t191, t4 * t43 + t26 * t9 + t267 * t42 - t25 * t8 + t44 * t89 + t84 * t70 - g(1) * (-t119 * t181 - t173 * t183) - g(2) * (-t173 * t181 + t112), t7, t29, t195, 0, t210, t217, t10 * t66 + t30 * t41 + t33 * t47 + t5 * t74 + t196, -t1 * t74 + t2 * t75 + t23 * t48 - t24 * t47 + t191, -t10 * t206 - t31 * t41 - t33 * t48 - t5 * t75 + t284, -g(2) * t112 + t1 * t43 + t33 * t10 + t2 * t42 + t23 * t8 + t24 * t9 + t5 * t41 + (g(1) * t173 - g(2) * t212) * t183 + (-g(1) * (-t119 - t212) + g(2) * t173) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t249 * t186, t244, t241, t243, qJDD(2), t180 * t192 - t274, g(3) * t180 + t182 * t192, 0, 0, t254, t62, t45, -t254, -t211, t172, -t78 * t174 + (-t105 * t247 + t282 * t172 - t174 * t246) * pkin(2) + t189, t79 * t174 + (-t107 * t247 - t172 * t179 - t174 * t235) * pkin(2) + t194, (t72 + t78) * t107 + (-t71 + t79) * t105 + (t282 * t59 - t179 * t60 + (-t282 * t105 + t107 * t179) * qJD(3)) * pkin(2), -t71 * t78 - t72 * t79 + (t282 * t39 - t274 - t179 * t226 + (-t179 * t71 + t282 * t72) * qJD(3) + (qJD(1) * t133 + t221) * t180) * pkin(2), t288, t19, t22, -t288, t21, t172, t172 * t99 - t66 * t88 + t201 - t234, -t100 * t172 - t174 * t264 - t206 * t88 + t203, -t100 * t30 - t264 * t66 - t31 * t99 + t218 + t239, -g(3) * t231 + t4 * t100 - t124 * t221 - t25 * t265 + t26 * t264 - t267 * t99 - t84 * t88, t288, t22, -t19, t172, -t21, -t288, -t172 * t95 - t37 * t66 + t190 - t234, -t266 * t66 - t30 * t94 + t31 * t95 + t219 + t239, t172 * t94 + t174 * t266 + t206 * t37 + t166 + t202, t1 * t94 + t2 * t95 - t33 * t37 - g(1) * (t183 * t230 + t126) - g(2) * (t181 * t230 + t125) - g(3) * (t231 + t212) + t266 * t24 + t265 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, t62, t45, -t254, -t211, t172, t72 * t174 + t189, t174 * t71 + t194, 0, 0, t288, t19, t22, -t288, t21, t172, t259 + (-t107 * t66 + t172 * t257) * pkin(3) + t201, t262 + (-t107 * t206 - t172 * t178) * pkin(3) + t203, -t272 + t28 * t66 + (-t178 * t30 - t257 * t31) * pkin(3) + t218, t25 * t27 - t26 * t28 + (-t107 * t84 + t178 * t4 - t257 * t267 + t289) * pkin(3), t288, t22, -t19, t172, -t21, -t288, -t157 * t172 - t40 * t66 + t190 + t259, -t154 * t30 + t157 * t31 - t250 * t66 + t219 - t272, t154 * t172 + t206 * t40 + 0.2e1 * t166 + t202 - t262, t1 * t154 + t2 * t157 - t33 * t40 - t23 * t27 - g(1) * (t183 * t222 + t126) - g(2) * (t181 * t222 + t125) - g(3) * (t159 + t212) + t250 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -t199, t216, t206 * t25 + t26 * t66 - t286 + t44, 0, 0, 0, 0, 0, 0, t209, t216, t199, -t206 * t23 + t24 * t66 - t286 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172 + t288, t22, -t174 ^ 2 - t271, -t174 * t24 - t190;];
tau_reg = t3;
