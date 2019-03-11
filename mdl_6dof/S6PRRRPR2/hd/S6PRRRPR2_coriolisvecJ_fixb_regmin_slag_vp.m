% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:08:58
% EndTime: 2019-03-08 23:09:08
% DurationCPUTime: 3.86s
% Computational Cost: add. (5312->352), mult. (13391->526), div. (0->0), fcn. (10452->12), ass. (0->209)
t180 = qJD(3) + qJD(4);
t188 = sin(qJ(4));
t192 = cos(qJ(3));
t293 = cos(qJ(4));
t237 = qJD(2) * t293;
t189 = sin(qJ(3));
t256 = qJD(2) * t189;
t302 = -t188 * t256 + t192 * t237;
t107 = t302 * t180;
t190 = sin(qJ(2));
t184 = sin(pkin(6));
t259 = qJD(1) * t184;
t243 = t190 * t259;
t285 = qJD(3) * pkin(3);
t202 = t189 * t285 - t243;
t138 = qJD(6) - t302;
t294 = -pkin(9) - pkin(8);
t245 = qJD(3) * t294;
t154 = t189 * t245;
t155 = t192 * t245;
t209 = -t188 * t189 + t293 * t192;
t193 = cos(qJ(2));
t239 = t193 * t259;
t161 = t294 * t189;
t162 = t294 * t192;
t296 = t293 * t161 + t188 * t162;
t304 = -qJD(4) * t296 - t293 * t154 - t188 * t155 + t209 * t239;
t119 = t180 * t209;
t264 = t188 * t192;
t152 = t293 * t189 + t264;
t120 = t180 * t152;
t303 = pkin(4) * t120 - qJ(5) * t119 - qJD(5) * t152 + t202;
t185 = cos(pkin(12));
t191 = cos(qJ(6));
t263 = t191 * t185;
t183 = sin(pkin(12));
t187 = sin(qJ(6));
t268 = t183 * t187;
t149 = -t263 + t268;
t283 = t138 * t149;
t150 = t183 * t191 + t185 * t187;
t301 = t138 * t150;
t146 = -qJD(2) * t264 - t189 * t237;
t129 = t146 * t183 + t185 * t180;
t300 = t129 * t191;
t127 = t146 * t185 - t180 * t183;
t214 = t127 * t191 - t129 * t187;
t299 = t138 * t214;
t288 = t304 * t183 + t303 * t185;
t287 = t303 * t183 - t304 * t185;
t232 = -t294 * qJD(2) + t243;
t186 = cos(pkin(6));
t258 = qJD(1) * t186;
t124 = -t232 * t189 + t192 * t258;
t116 = t124 + t285;
t125 = t189 * t258 + t232 * t192;
t236 = qJD(4) * t293;
t254 = qJD(4) * t188;
t257 = qJD(2) * t184;
t235 = qJD(1) * t257;
t229 = t193 * t235;
t87 = qJD(3) * t124 + t192 * t229;
t88 = -qJD(3) * t125 - t189 * t229;
t197 = t116 * t236 - t125 * t254 + t188 * t88 + t293 * t87;
t20 = t180 * qJD(5) + t197;
t108 = t120 * qJD(2);
t251 = qJD(2) * qJD(3);
t234 = t189 * t251;
t139 = pkin(3) * t234 + t190 * t235;
t48 = pkin(4) * t108 - qJ(5) * t107 + qJD(5) * t146 + t139;
t10 = t183 * t48 + t185 * t20;
t9 = -t183 * t20 + t185 * t48;
t223 = t10 * t185 - t183 * t9;
t133 = t188 * t161 - t293 * t162;
t284 = -t133 * qJD(4) + t152 * t239 - t188 * t154 + t293 * t155;
t166 = pkin(3) * t236 + qJD(5);
t114 = t188 * t125;
t70 = t293 * t124 - t114;
t109 = -pkin(4) * t146 - qJ(5) * t302;
t97 = pkin(3) * t256 + t109;
t42 = -t183 * t70 + t185 * t97;
t298 = t166 * t183 + t42;
t43 = t183 * t97 + t185 * t70;
t297 = -t166 * t185 + t43;
t115 = t293 * t125;
t69 = t188 * t124 + t115;
t227 = pkin(3) * t254 - t69;
t295 = -qJD(6) + t138;
t30 = -t214 * qJD(6) + t150 * t107;
t291 = t185 * pkin(5);
t179 = t185 * pkin(10);
t290 = pkin(5) * t120 - t119 * t179 + t288;
t277 = t119 * t183;
t289 = pkin(10) * t277 - t287;
t68 = t188 * t116 + t115;
t62 = qJ(5) * t180 + t68;
t178 = -pkin(3) * t192 - pkin(2);
t136 = t178 * qJD(2) - t239;
t83 = -pkin(4) * t302 + qJ(5) * t146 + t136;
t34 = t183 * t83 + t185 * t62;
t286 = qJD(2) * pkin(2);
t281 = pkin(5) * t277 - t284;
t67 = t293 * t116 - t114;
t47 = t183 * t109 + t185 * t67;
t252 = qJD(6) * t191;
t280 = t107 * t263 + t129 * t252;
t279 = t107 * t183;
t278 = t107 * t185;
t276 = t138 * t146;
t275 = t302 * t183;
t274 = t302 * t185;
t273 = t146 * t302;
t272 = t152 * t183;
t271 = t152 * t185;
t267 = t184 * t190;
t266 = t184 * t193;
t195 = qJD(2) ^ 2;
t265 = t184 * t195;
t194 = qJD(3) ^ 2;
t262 = t194 * t189;
t261 = t194 * t192;
t111 = -pkin(4) * t209 - qJ(5) * t152 + t178;
t66 = t183 * t111 + t185 * t133;
t260 = t189 ^ 2 - t192 ^ 2;
t255 = qJD(2) * t190;
t253 = qJD(6) * t152;
t250 = pkin(10) * t275;
t247 = t190 * t265;
t3 = pkin(5) * t108 - pkin(10) * t278 + t9;
t7 = -pkin(10) * t279 + t10;
t246 = -t187 * t7 + t191 * t3;
t242 = t184 * t255;
t241 = t193 * t257;
t22 = t116 * t254 + t125 * t236 + t188 * t87 - t293 * t88;
t233 = -t146 * t34 + t22 * t183;
t33 = -t183 * t62 + t185 * t83;
t46 = t185 * t109 - t183 * t67;
t65 = t185 * t111 - t133 * t183;
t231 = t189 * t241;
t230 = t192 * t241;
t177 = -t293 * pkin(3) - pkin(4);
t135 = pkin(5) * t275;
t228 = -t135 + t227;
t226 = -t146 * pkin(5) - pkin(10) * t274;
t225 = t187 * t3 + t191 * t7;
t222 = t33 * t274 + t34 * t275 + t223;
t221 = t146 * t33 - t185 * t22;
t18 = -pkin(5) * t302 + pkin(10) * t127 + t33;
t23 = pkin(10) * t129 + t34;
t5 = t18 * t191 - t187 * t23;
t6 = t18 * t187 + t191 * t23;
t220 = -t183 * t33 + t185 * t34;
t52 = -pkin(5) * t209 - pkin(10) * t271 + t65;
t57 = -pkin(10) * t272 + t66;
t219 = -t187 * t57 + t191 * t52;
t218 = t187 * t52 + t191 * t57;
t142 = t186 * t192 - t189 * t267;
t143 = t186 * t189 + t192 * t267;
t99 = t188 * t142 + t293 * t143;
t81 = -t183 * t99 - t185 * t266;
t82 = -t183 * t266 + t185 * t99;
t217 = -t187 * t82 + t191 * t81;
t216 = t187 * t81 + t191 * t82;
t213 = qJD(6) * t127 - t279;
t174 = pkin(3) * t188 + qJ(5);
t147 = (-pkin(10) - t174) * t183;
t212 = -qJD(6) * t147 - t250 + t297;
t148 = t174 * t185 + t179;
t211 = qJD(6) * t148 + t226 + t298;
t210 = t293 * t142 - t188 * t143;
t159 = (-pkin(10) - qJ(5)) * t183;
t208 = -qJD(5) * t185 - qJD(6) * t159 - t250 + t47;
t160 = qJ(5) * t185 + t179;
t207 = qJD(5) * t183 + qJD(6) * t160 + t226 + t46;
t16 = pkin(5) * t279 + t22;
t61 = -t180 * pkin(4) + qJD(5) - t67;
t51 = -pkin(5) * t129 + t61;
t206 = t5 * t146 + t16 * t149 + t301 * t51;
t205 = -t6 * t146 + t16 * t150 - t283 * t51;
t204 = t136 * t146 - t22;
t203 = t286 * qJD(2);
t201 = -t107 * t296 + t119 * t61 + t152 * t22;
t200 = -0.2e1 * qJD(3) * t286;
t199 = -pkin(4) * t107 - qJ(5) * t108 - (-qJD(5) + t61) * t302;
t198 = t107 * t177 - t108 * t174 - (-t166 + t61) * t302;
t29 = t213 * t187 + t280;
t196 = -t136 * t302 - t197;
t175 = -pkin(4) - t291;
t158 = t177 - t291;
t123 = -t143 * qJD(3) - t231;
t122 = t142 * qJD(3) + t230;
t103 = t149 * t152;
t102 = t150 * t152;
t93 = pkin(5) * t272 - t296;
t84 = t146 ^ 2 - t302 ^ 2;
t78 = (-qJD(2) * t152 - t146) * t180;
t72 = -t127 * t187 - t300;
t56 = t135 + t68;
t50 = t99 * qJD(4) + t188 * t122 - t293 * t123;
t49 = t210 * qJD(4) + t293 * t122 + t188 * t123;
t41 = t150 * t119 + t252 * t271 - t253 * t268;
t40 = -t149 * t119 - t150 * t253;
t37 = t183 * t242 + t185 * t49;
t36 = -t183 * t49 + t185 * t242;
t12 = -t149 * t108 - t138 * t301 - t72 * t146;
t11 = t150 * t108 - t283 * t138 - t146 * t214;
t4 = t29 * t150 + t214 * t283;
t1 = -t149 * t29 - t150 * t30 + t214 * t301 + t283 * t72;
t2 = [0, 0, -t247, -t193 * t265, 0, 0, 0, 0, 0, -t192 * t247 + (t123 - t231) * qJD(3), t189 * t247 + (-t122 - t230) * qJD(3), 0, 0, 0, 0, 0, -t50 * t180 + (-t108 * t193 - t255 * t302) * t184, -t49 * t180 + (-t107 * t193 - t146 * t255) * t184, t81 * t108 - t129 * t50 - t210 * t279 - t302 * t36, -t82 * t108 - t127 * t50 - t210 * t278 + t302 * t37, t129 * t37 + t127 * t36 + (-t183 * t82 - t185 * t81) * t107, t10 * t82 - t210 * t22 + t33 * t36 + t34 * t37 + t50 * t61 + t81 * t9, 0, 0, 0, 0, 0 (-qJD(6) * t216 - t187 * t37 + t191 * t36) * t138 + t217 * t108 + t50 * t72 - t210 * t30 -(qJD(6) * t217 + t187 * t36 + t191 * t37) * t138 - t216 * t108 - t50 * t214 - t210 * t29; 0, 0, 0, 0, 0.2e1 * t192 * t234, -0.2e1 * t260 * t251, t261, -t262, 0, -pkin(8) * t261 + t189 * t200, pkin(8) * t262 + t192 * t200, t107 * t152 - t119 * t146, t107 * t209 - t108 * t152 + t119 * t302 + t120 * t146, t119 * t180, -t120 * t180, 0, t108 * t178 + t120 * t136 - t139 * t209 + t284 * t180 - t202 * t302, t107 * t178 + t119 * t136 + t139 * t152 - t202 * t146 + t304 * t180, t108 * t65 + t120 * t33 + t129 * t284 + t201 * t183 - t209 * t9 - t288 * t302, t10 * t209 - t108 * t66 - t120 * t34 + t127 * t284 + t185 * t201 + t287 * t302, t288 * t127 + t287 * t129 + (-t107 * t65 - t119 * t33 - t152 * t9) * t185 + (-t10 * t152 - t107 * t66 - t119 * t34) * t183, t10 * t66 - t22 * t296 - t284 * t61 + t287 * t34 + t288 * t33 + t65 * t9, -t103 * t29 - t214 * t40, -t102 * t29 + t103 * t30 + t214 * t41 - t40 * t72, -t103 * t108 - t120 * t214 + t138 * t40 - t209 * t29, -t102 * t108 - t120 * t72 - t138 * t41 + t209 * t30, -t108 * t209 + t120 * t138, t219 * t108 - t246 * t209 + t5 * t120 + t93 * t30 + t16 * t102 + t51 * t41 + t281 * t72 + (t187 * t289 + t191 * t290) * t138 + (-t138 * t218 + t209 * t6) * qJD(6), -t218 * t108 + t225 * t209 - t6 * t120 + t93 * t29 - t16 * t103 + t51 * t40 - t281 * t214 + (-t187 * t290 + t191 * t289) * t138 + (-t138 * t219 + t209 * t5) * qJD(6); 0, 0, 0, 0, -t189 * t195 * t192, t260 * t195, 0, 0, 0, t189 * t203, t192 * t203, t273, t84, 0, t78, 0, t180 * t69 + (-t180 * t254 + t256 * t302) * pkin(3) + t204, t70 * t180 + (t146 * t256 - t180 * t236) * pkin(3) + t196, -t129 * t227 + t183 * t198 + t302 * t42 + t221, -t127 * t227 + t185 * t198 - t302 * t43 + t233, -t127 * t298 - t129 * t297 + t222, t166 * t220 + t174 * t223 + t177 * t22 + t227 * t61 - t33 * t42 - t34 * t43, t4, t1, t11, t12, t276 (t147 * t191 - t148 * t187) * t108 + t158 * t30 + t228 * t72 + (t187 * t212 - t191 * t211) * t138 + t206 -(t147 * t187 + t148 * t191) * t108 + t158 * t29 - t228 * t214 + (t187 * t211 + t191 * t212) * t138 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, t84, 0, t78, 0, t180 * t68 + t204, t67 * t180 + t196, t129 * t68 + t183 * t199 + t302 * t46 + t221, t127 * t68 + t185 * t199 - t302 * t47 + t233, -t129 * t47 - t127 * t46 + (-t127 * t183 + t129 * t185) * qJD(5) + t222, -pkin(4) * t22 + qJ(5) * t223 + qJD(5) * t220 - t33 * t46 - t34 * t47 - t61 * t68, t4, t1, t11, t12, t276 (t159 * t191 - t160 * t187) * t108 + t175 * t30 - t56 * t72 + (t187 * t208 - t191 * t207) * t138 + t206 -(t159 * t187 + t160 * t191) * t108 + t175 * t29 + t56 * t214 + (t187 * t207 + t191 * t208) * t138 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * t302 + t279, -t129 * t302 + t278, -t127 ^ 2 - t129 ^ 2, -t127 * t33 - t129 * t34 + t22, 0, 0, 0, 0, 0, t30 - t299, t138 * t300 + (t127 * t138 + t213) * t187 + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214 * t72, t214 ^ 2 - t72 ^ 2, t138 * t72 + t29, -t30 - t299, t108, t214 * t51 + t295 * t6 + t246, t295 * t5 + t51 * t72 - t225;];
tauc_reg  = t2;
