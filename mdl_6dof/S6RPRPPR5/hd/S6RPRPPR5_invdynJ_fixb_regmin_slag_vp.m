% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:41
% EndTime: 2019-03-09 02:51:52
% DurationCPUTime: 4.88s
% Computational Cost: add. (4467->428), mult. (10453->527), div. (0->0), fcn. (8053->14), ass. (0->219)
t186 = sin(pkin(9));
t188 = cos(pkin(9));
t192 = sin(qJ(3));
t302 = cos(qJ(3));
t141 = t302 * t186 + t192 * t188;
t313 = t141 * qJD(1);
t122 = qJD(6) + t313;
t185 = sin(pkin(10));
t187 = cos(pkin(10));
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t138 = t185 * t194 + t187 * t191;
t286 = t122 * t138;
t314 = -t185 * t191 + t187 * t194;
t260 = qJD(3) * t192;
t248 = t186 * t260;
t250 = t302 * t188;
t238 = qJD(1) * t250;
t244 = qJDD(1) * t302;
t255 = t188 * qJDD(1);
t254 = qJD(3) * t238 + t186 * t244 + t192 * t255;
t90 = qJD(1) * t248 - t254;
t88 = -qJDD(6) + t90;
t321 = -t286 * t122 - t314 * t88;
t271 = t186 * t192;
t139 = -t250 + t271;
t82 = t314 * t139;
t258 = qJD(6) * t194;
t259 = qJD(6) * t191;
t285 = -t185 * t259 + t187 * t258 + t314 * t313;
t240 = -t122 * t285 + t138 * t88;
t249 = qJD(1) * t271;
t129 = -t238 + t249;
t100 = qJD(3) * t185 - t187 * t129;
t102 = qJD(3) * t187 + t129 * t185;
t225 = t100 * t191 - t102 * t194;
t320 = t122 * t225;
t282 = qJDD(1) * pkin(1);
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t316 = g(1) * t193 - g(2) * t195;
t219 = -qJDD(2) + t282 + t316;
t181 = pkin(9) + qJ(3);
t174 = sin(t181);
t217 = t316 * t174;
t294 = pkin(7) + qJ(2);
t148 = t294 * t186;
t149 = t294 * t188;
t247 = qJD(3) * t302;
t74 = (qJD(2) * t186 + qJD(3) * t149) * t192 - qJD(2) * t250 + t148 * t247;
t99 = -t192 * t148 + t302 * t149;
t319 = qJD(3) * t74 - qJDD(3) * t99 - t217;
t176 = cos(t181);
t299 = g(1) * t195;
t236 = g(2) * t193 + t299;
t209 = -g(3) * t174 - t176 * t236;
t257 = qJD(1) * qJD(2);
t307 = t294 * qJDD(1) + t257;
t113 = t307 * t186;
t114 = t307 * t188;
t142 = qJD(1) * t148;
t143 = qJD(1) * t149;
t242 = t192 * t113 - t302 * t114 + t142 * t247 + t143 * t260;
t199 = t209 - t242;
t264 = -t302 * t142 - t192 * t143;
t318 = qJD(3) * t264 - t199;
t312 = qJD(4) - t264;
t256 = t186 * qJDD(1);
t231 = -t188 * t244 + t192 * t256;
t317 = 0.2e1 * t313 * qJD(3) + t231;
t262 = t176 * pkin(3) + t174 * qJ(4);
t182 = qJDD(3) * qJ(4);
t183 = qJD(3) * qJD(4);
t315 = -t182 - t183;
t311 = -qJD(6) + t122;
t310 = qJ(2) * qJDD(1);
t136 = t141 * qJD(3);
t91 = qJD(1) * t136 + t231;
t309 = -pkin(4) * t91 + qJDD(5);
t305 = t313 ^ 2;
t308 = t185 * t90 - t187 * t305;
t78 = qJDD(3) * t185 - t187 * t91;
t79 = qJDD(3) * t187 + t185 * t91;
t18 = -qJD(6) * t225 + t191 * t79 + t194 * t78;
t306 = t129 ^ 2;
t304 = pkin(3) * t91;
t301 = pkin(8) * t187;
t169 = g(3) * t176;
t296 = pkin(3) + qJ(5);
t295 = pkin(4) + t294;
t293 = -pkin(8) - t296;
t171 = pkin(2) * t188 + pkin(1);
t146 = -qJDD(1) * t171 + qJDD(2);
t207 = qJ(4) * t90 + t146;
t200 = -qJD(4) * t313 + t207;
t14 = qJD(5) * t129 + t296 * t91 + t200;
t213 = t302 * t113 + t192 * t114 - t142 * t260 + t143 * t247;
t211 = qJDD(4) + t213;
t26 = -t90 * pkin(4) - qJD(3) * qJD(5) - t296 * qJDD(3) + t211;
t7 = t187 * t14 + t185 * t26;
t135 = -t188 * t247 + t248;
t223 = qJ(4) * t135 - qJD(4) * t141;
t39 = qJD(5) * t139 + t296 * t136 + t223;
t75 = qJD(2) * t141 + qJD(3) * t99;
t48 = -t135 * pkin(4) + t75;
t16 = t185 * t48 + t187 * t39;
t147 = -qJD(1) * t171 + qJD(2);
t208 = -qJ(4) * t313 + t147;
t46 = t296 * t129 + t208;
t266 = pkin(4) * t313 + t312;
t52 = -t296 * qJD(3) + t266;
t25 = t185 * t52 + t187 * t46;
t280 = t129 * qJ(4);
t61 = t296 * t313 + t280;
t95 = -t192 * t142 + t302 * t143;
t65 = -pkin(4) * t129 + t95;
t30 = t185 * t65 + t187 * t61;
t220 = -qJ(4) * t141 - t171;
t63 = t296 * t139 + t220;
t98 = t302 * t148 + t192 * t149;
t80 = t141 * pkin(4) + t98;
t33 = t185 * t80 + t187 * t63;
t53 = t194 * t100 + t102 * t191;
t291 = t122 * t53;
t290 = t129 * t53;
t288 = t296 * t90;
t287 = t225 * t129;
t284 = qJ(5) * t176;
t281 = qJDD(3) * pkin(3);
t279 = t313 * t129;
t278 = t174 * t193;
t277 = t174 * t195;
t180 = pkin(10) + qJ(6);
t175 = cos(t180);
t276 = t175 * t193;
t275 = t175 * t195;
t274 = t176 * t193;
t273 = t176 * t195;
t269 = t95 * qJD(3);
t252 = -pkin(5) * t187 - pkin(4);
t267 = -t252 * t313 + t312;
t184 = qJD(3) * qJ(4);
t57 = qJD(5) + t184 + t65;
t265 = qJD(5) - t57;
t261 = t186 ^ 2 + t188 ^ 2;
t253 = -g(1) * t277 - g(2) * t278 + t169;
t6 = -t14 * t185 + t187 * t26;
t2 = -pkin(5) * t90 - pkin(8) * t79 + t6;
t5 = -pkin(8) * t78 + t7;
t251 = -t191 * t5 + t194 * t2;
t24 = -t185 * t46 + t187 * t52;
t243 = t261 * qJD(1) ^ 2;
t241 = 0.2e1 * t261;
t239 = g(2) * (pkin(3) * t273 + qJ(4) * t277 + t195 * t171);
t237 = g(1) * t274 - g(2) * t273;
t234 = t185 * t7 + t187 * t6;
t233 = -t185 * t6 + t187 * t7;
t232 = t191 * t2 + t194 * t5;
t10 = pkin(5) * t313 - pkin(8) * t102 + t24;
t13 = -pkin(8) * t100 + t25;
t3 = t10 * t194 - t13 * t191;
t4 = t10 * t191 + t13 * t194;
t230 = -t185 * t25 - t187 * t24;
t229 = -t185 * t24 + t187 * t25;
t73 = t187 * t80;
t20 = t141 * pkin(5) + t73 + (-pkin(8) * t139 - t63) * t185;
t28 = t139 * t301 + t33;
t228 = t191 * t28 - t194 * t20;
t227 = t191 * t20 + t194 * t28;
t226 = -t305 * t185 - t187 * t90;
t218 = -t171 - t262;
t34 = t242 + t315;
t144 = t293 * t185;
t60 = t187 * t65;
t216 = qJD(5) * t187 + qJD(6) * t144 - t129 * pkin(5) + t60 + (-pkin(8) * t313 - t61) * t185;
t145 = t293 * t187;
t215 = qJD(5) * t185 - qJD(6) * t145 + t301 * t313 + t30;
t17 = -t100 * t258 - t102 * t259 - t191 * t78 + t194 * t79;
t83 = t138 * t139;
t210 = t219 + t282;
t27 = -t34 + t309;
t206 = t136 * t57 + t139 * t27 + t236;
t204 = -qJD(3) * t75 - qJDD(3) * t98 + t237;
t203 = -t213 - t253;
t202 = t27 + t209;
t201 = t241 * t257 - t236;
t71 = pkin(3) * t129 + t208;
t198 = t313 * t71 + qJDD(4) - t203;
t173 = sin(t180);
t167 = pkin(5) * t185 + qJ(4);
t152 = qJ(4) * t273;
t150 = qJ(4) * t274;
t117 = qJD(3) * t129;
t112 = -t173 * t278 + t275;
t111 = t173 * t195 + t174 * t276;
t110 = t173 * t277 + t276;
t109 = -t173 * t193 + t174 * t275;
t89 = pkin(3) * t139 + t220;
t87 = pkin(3) * t313 + t280;
t86 = -t184 - t95;
t85 = -qJD(3) * pkin(3) + t312;
t81 = -t139 * pkin(4) + t99;
t62 = pkin(3) * t136 + t223;
t51 = t139 * t252 + t99;
t47 = -pkin(4) * t136 - t74;
t45 = t187 * t48;
t41 = pkin(5) * t100 + t57;
t40 = t136 * t252 - t74;
t38 = t211 - t281;
t37 = qJD(6) * t83 - t314 * t136;
t36 = qJD(6) * t82 + t136 * t138;
t32 = -t185 * t63 + t73;
t31 = t200 + t304;
t29 = -t185 * t61 + t60;
t15 = -t185 * t39 + t45;
t11 = pkin(5) * t78 + t27;
t9 = t136 * t301 + t16;
t8 = -t135 * pkin(5) + t45 + (-pkin(8) * t136 - t39) * t185;
t1 = [qJDD(1), t316, t236, t210 * t188, -t210 * t186, t241 * t310 + t201, pkin(1) * t219 + (t261 * t310 + t201) * qJ(2), -t135 * t313 - t141 * t90, t129 * t135 - t136 * t313 + t139 * t90 - t141 * t91, -qJD(3) * t135 + qJDD(3) * t141, -qJD(3) * t136 - qJDD(3) * t139, 0, t136 * t147 + t139 * t146 - t171 * t91 + t204, -t135 * t147 + t141 * t146 + t171 * t90 + t319, t129 * t74 - t135 * t85 + t136 * t86 + t139 * t34 + t141 * t38 + t313 * t75 - t90 * t98 - t91 * t99 - t236, -t129 * t62 - t136 * t71 - t139 * t31 - t89 * t91 - t204, t135 * t71 - t141 * t31 - t313 * t62 + t89 * t90 - t319, t31 * t89 + t71 * t62 - t34 * t99 + t86 * t74 + t38 * t98 + t85 * t75 - t294 * t299 - t239 + (-g(1) * t218 - g(2) * t294) * t193, t47 * t100 - t24 * t135 + t6 * t141 + t15 * t313 + t185 * t217 - t187 * t206 - t32 * t90 + t81 * t78, t47 * t102 + t25 * t135 - t7 * t141 - t16 * t313 + t185 * t206 + t187 * t217 + t33 * t90 + t81 * t79, -t100 * t16 - t102 * t15 + t136 * t229 + t139 * t233 - t32 * t79 - t33 * t78 + t237, t7 * t33 + t25 * t16 + t6 * t32 + t24 * t15 + t27 * t81 + t57 * t47 - t239 + (-g(1) * t295 - g(2) * t284) * t195 + (-g(1) * (t218 - t284) - g(2) * t295) * t193, t17 * t83 - t225 * t36, t17 * t82 - t18 * t83 + t225 * t37 - t36 * t53, t122 * t36 + t135 * t225 + t141 * t17 - t83 * t88, -t122 * t37 + t135 * t53 - t141 * t18 - t82 * t88, -t122 * t135 - t141 * t88 (-t191 * t9 + t194 * t8) * t122 + t228 * t88 + t251 * t141 - t3 * t135 + t40 * t53 + t51 * t18 - t11 * t82 + t41 * t37 - g(1) * t112 - g(2) * t110 + (-t122 * t227 - t141 * t4) * qJD(6) -(t191 * t8 + t194 * t9) * t122 + t227 * t88 - t232 * t141 + t4 * t135 - t40 * t225 + t51 * t17 + t11 * t83 + t41 * t36 + g(1) * t111 - g(2) * t109 + (t122 * t228 - t141 * t3) * qJD(6); 0, 0, 0, -t255, t256, -t243, -qJ(2) * t243 - t219, 0, 0, 0, 0, 0, t317 (-t129 - t249) * qJD(3) + t254, -t305 - t306, -t317, t117 + t90, t304 - t129 * t86 + (-qJD(4) - t85) * t313 + t207 - t316, t100 * t129 + t308, t102 * t129 - t226, t185 * t79 - t187 * t78 + (t100 * t185 + t102 * t187) * t313, t129 * t57 + t230 * t313 + t233 - t316, 0, 0, 0, 0, 0, t240 + t290, -t287 - t321; 0, 0, 0, 0, 0, 0, 0, t279, t305 - t306 (t129 - t249) * qJD(3) + t254, -t231, qJDD(3), -t147 * t313 + t203 + t269, t129 * t147 + t318, pkin(3) * t90 - qJ(4) * t91 + (-t86 - t95) * t313 + (t85 - t312) * t129, t87 * t129 + t198 - t269 - 0.2e1 * t281, -t129 * t71 + t313 * t87 + 0.2e1 * t182 + 0.2e1 * t183 - t318, -t34 * qJ(4) - t38 * pkin(3) - t71 * t87 - t85 * t95 - g(1) * (-pkin(3) * t277 + t152) - g(2) * (-pkin(3) * t278 + t150) - g(3) * t262 - t312 * t86, t187 * t288 + qJ(4) * t78 + t129 * t24 + t266 * t100 + (-t187 * t265 - t29) * t313 + t202 * t185, -t185 * t288 + qJ(4) * t79 - t129 * t25 + t266 * t102 + (t185 * t265 + t30) * t313 + t202 * t187, t100 * t30 + t102 * t29 + (qJD(5) * t102 - t25 * t313 + t296 * t79 - t6) * t187 + (qJD(5) * t100 + t24 * t313 + t296 * t78 - t7) * t185 - t253, t27 * qJ(4) - t25 * t30 - t24 * t29 - g(1) * t152 - g(2) * t150 - g(3) * (t262 + t284) + t266 * t57 + t230 * qJD(5) + (t174 * t236 - t234) * t296, t17 * t314 + t225 * t286, -t138 * t17 - t18 * t314 + t225 * t285 + t286 * t53, -t287 + t321, t240 - t290, t122 * t129 -(-t144 * t191 + t145 * t194) * t88 + t167 * t18 + t11 * t138 + t3 * t129 + t267 * t53 + t285 * t41 + (t191 * t215 - t194 * t216) * t122 + t209 * t173 (t144 * t194 + t145 * t191) * t88 + t167 * t17 + t11 * t314 - t4 * t129 - t267 * t225 - t286 * t41 + (t191 * t216 + t194 * t215) * t122 + t209 * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 - t90, qJDD(3) - t279, -qJD(3) ^ 2 - t305, t86 * qJD(3) + t198 - t281, -qJD(3) * t100 + t226, -qJD(3) * t102 + t308, -t185 * t78 - t187 * t79 + (-t100 * t187 + t102 * t185) * t313, -qJD(3) * t57 + t229 * t313 + t234 + t253, 0, 0, 0, 0, 0, -qJD(3) * t53 + t321, qJD(3) * t225 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t313 + t78, -t100 * t313 + t79, -t100 ^ 2 - t102 ^ 2, t100 * t25 + t102 * t24 + t199 + t309 - t315, 0, 0, 0, 0, 0, t18 - t320, t17 - t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225 * t53, t225 ^ 2 - t53 ^ 2, t17 + t291, -t18 - t320, -t88, -g(1) * t109 - g(2) * t111 + t175 * t169 + t225 * t41 + t311 * t4 + t251, g(1) * t110 - g(2) * t112 - t173 * t169 + t311 * t3 + t41 * t53 - t232;];
tau_reg  = t1;
