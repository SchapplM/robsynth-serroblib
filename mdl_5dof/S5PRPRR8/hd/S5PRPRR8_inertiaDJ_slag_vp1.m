% Calculate time derivative of joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:28
% EndTime: 2019-12-05 16:02:48
% DurationCPUTime: 9.69s
% Computational Cost: add. (25917->762), mult. (73768->1136), div. (0->0), fcn. (83660->10), ass. (0->294)
t259 = cos(pkin(5));
t258 = cos(pkin(9));
t263 = cos(qJ(2));
t285 = qJD(2) * t263;
t276 = t258 * t285;
t256 = sin(pkin(9));
t261 = sin(qJ(2));
t286 = qJD(2) * t261;
t279 = t256 * t286;
t229 = -t259 * t276 + t279;
t297 = t259 * t263;
t243 = t256 * t297 + t258 * t261;
t231 = t243 * qJD(2);
t298 = t259 * t261;
t242 = t256 * t263 + t258 * t298;
t244 = -t256 * t298 + t258 * t263;
t257 = sin(pkin(5));
t311 = cos(qJ(4));
t283 = t257 * t311;
t310 = sin(qJ(4));
t209 = t243 * t310 + t256 * t283;
t260 = sin(qJ(5));
t262 = cos(qJ(5));
t163 = -t209 * t260 + t244 * t262;
t164 = t209 * t262 + t244 * t260;
t282 = t257 * t310;
t208 = -t243 * t311 + t256 * t282;
t107 = Icges(6,5) * t164 + Icges(6,6) * t163 + Icges(6,3) * t208;
t109 = Icges(6,4) * t164 + Icges(6,2) * t163 + Icges(6,6) * t208;
t111 = Icges(6,1) * t164 + Icges(6,4) * t163 + Icges(6,5) * t208;
t232 = -t259 * t279 + t276;
t159 = -t208 * qJD(4) + t232 * t310;
t117 = -t164 * qJD(5) - t159 * t260 - t231 * t262;
t118 = t163 * qJD(5) + t159 * t262 - t231 * t260;
t160 = t209 * qJD(4) - t232 * t311;
t73 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t160;
t75 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t160;
t77 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t160;
t16 = t107 * t160 + t109 * t117 + t111 * t118 + t163 * t75 + t164 * t77 + t208 * t73;
t241 = t256 * t261 - t258 * t297;
t326 = -t241 * t310 + t258 * t283;
t165 = t242 * t262 + t260 * t326;
t166 = t242 * t260 - t262 * t326;
t210 = t241 * t311 + t258 * t282;
t108 = Icges(6,5) * t166 + Icges(6,6) * t165 - Icges(6,3) * t210;
t110 = Icges(6,4) * t166 + Icges(6,2) * t165 - Icges(6,6) * t210;
t112 = Icges(6,1) * t166 + Icges(6,4) * t165 - Icges(6,5) * t210;
t230 = t242 * qJD(2);
t161 = t210 * qJD(4) + t230 * t310;
t119 = -t166 * qJD(5) - t161 * t260 - t229 * t262;
t120 = t165 * qJD(5) + t161 * t262 - t229 * t260;
t162 = t326 * qJD(4) + t230 * t311;
t74 = Icges(6,5) * t120 + Icges(6,6) * t119 - Icges(6,3) * t162;
t76 = Icges(6,4) * t120 + Icges(6,2) * t119 - Icges(6,6) * t162;
t78 = Icges(6,1) * t120 + Icges(6,4) * t119 - Icges(6,5) * t162;
t17 = t108 * t160 + t110 * t117 + t112 * t118 + t163 * t76 + t164 * t78 + t208 * t74;
t266 = -t259 * t310 - t263 * t283;
t278 = t257 * t286;
t212 = t266 * qJD(4) + t310 * t278;
t246 = t259 * t311 - t263 * t282;
t299 = t257 * t261;
t215 = t246 * t262 + t260 * t299;
t277 = t257 * t285;
t145 = -t215 * qJD(5) - t212 * t260 + t262 * t277;
t214 = -t246 * t260 + t262 * t299;
t146 = t214 * qJD(5) + t212 * t262 + t260 * t277;
t213 = t246 * qJD(4) - t311 * t278;
t100 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t213;
t141 = Icges(6,5) * t215 + Icges(6,6) * t214 - Icges(6,3) * t266;
t142 = Icges(6,4) * t215 + Icges(6,2) * t214 - Icges(6,6) * t266;
t143 = Icges(6,1) * t215 + Icges(6,4) * t214 - Icges(6,5) * t266;
t98 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t213;
t99 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t213;
t32 = t100 * t164 + t117 * t142 + t118 * t143 + t141 * t160 + t163 * t99 + t208 * t98;
t50 = t107 * t208 + t109 * t163 + t111 * t164;
t51 = t108 * t208 + t110 * t163 + t112 * t164;
t66 = t141 * t208 + t142 * t163 + t143 * t164;
t3 = t16 * t244 + t17 * t242 - t51 * t229 - t50 * t231 + (t261 * t32 + t66 * t285) * t257;
t121 = Icges(5,5) * t159 - Icges(5,6) * t160 - Icges(5,3) * t231;
t123 = Icges(5,4) * t159 - Icges(5,2) * t160 - Icges(5,6) * t231;
t125 = Icges(5,1) * t159 - Icges(5,4) * t160 - Icges(5,5) * t231;
t133 = Icges(5,5) * t209 - Icges(5,6) * t208 + Icges(5,3) * t244;
t135 = Icges(5,4) * t209 - Icges(5,2) * t208 + Icges(5,6) * t244;
t137 = Icges(5,1) * t209 - Icges(5,4) * t208 + Icges(5,5) * t244;
t41 = t121 * t244 - t123 * t208 + t125 * t209 - t133 * t231 - t135 * t160 + t137 * t159;
t122 = Icges(5,5) * t161 + Icges(5,6) * t162 - Icges(5,3) * t229;
t124 = Icges(5,4) * t161 + Icges(5,2) * t162 - Icges(5,6) * t229;
t126 = Icges(5,1) * t161 + Icges(5,4) * t162 - Icges(5,5) * t229;
t134 = -Icges(5,5) * t326 + Icges(5,6) * t210 + Icges(5,3) * t242;
t136 = -Icges(5,4) * t326 + Icges(5,2) * t210 + Icges(5,6) * t242;
t138 = -Icges(5,1) * t326 + Icges(5,4) * t210 + Icges(5,5) * t242;
t42 = t122 * t244 - t124 * t208 + t126 * t209 - t134 * t231 - t136 * t160 + t138 * t159;
t147 = Icges(5,5) * t212 - Icges(5,6) * t213 + Icges(5,3) * t277;
t148 = Icges(5,4) * t212 - Icges(5,2) * t213 + Icges(5,6) * t277;
t149 = Icges(5,1) * t212 - Icges(5,4) * t213 + Icges(5,5) * t277;
t179 = Icges(5,5) * t246 + Icges(5,6) * t266 + Icges(5,3) * t299;
t180 = Icges(5,4) * t246 + Icges(5,2) * t266 + Icges(5,6) * t299;
t181 = Icges(5,1) * t246 + Icges(5,4) * t266 + Icges(5,5) * t299;
t54 = t147 * t244 - t148 * t208 + t149 * t209 + t159 * t181 - t160 * t180 - t179 * t231;
t82 = t133 * t244 - t135 * t208 + t137 * t209;
t83 = t134 * t244 - t136 * t208 + t138 * t209;
t94 = t179 * t244 - t180 * t208 + t181 * t209;
t328 = -t83 * t229 - t82 * t231 + t42 * t242 + t41 * t244 + (t261 * t54 + t94 * t285) * t257 + t3;
t18 = -t107 * t162 + t109 * t119 + t111 * t120 + t165 * t75 + t166 * t77 - t210 * t73;
t19 = -t108 * t162 + t110 * t119 + t112 * t120 + t165 * t76 + t166 * t78 - t210 * t74;
t33 = t100 * t166 + t119 * t142 + t120 * t143 - t141 * t162 + t165 * t99 - t210 * t98;
t52 = -t107 * t210 + t109 * t165 + t111 * t166;
t53 = -t108 * t210 + t110 * t165 + t112 * t166;
t67 = -t141 * t210 + t142 * t165 + t143 * t166;
t4 = t18 * t244 + t19 * t242 - t53 * t229 - t52 * t231 + (t261 * t33 + t67 * t285) * t257;
t43 = t121 * t242 + t123 * t210 - t125 * t326 - t133 * t229 + t135 * t162 + t137 * t161;
t44 = t122 * t242 + t124 * t210 - t126 * t326 - t134 * t229 + t136 * t162 + t138 * t161;
t55 = t147 * t242 + t148 * t210 - t149 * t326 + t161 * t181 + t162 * t180 - t179 * t229;
t84 = t133 * t242 + t135 * t210 - t137 * t326;
t85 = t134 * t242 + t136 * t210 - t138 * t326;
t95 = t179 * t242 + t180 * t210 - t181 * t326;
t327 = -t85 * t229 - t84 * t231 + t44 * t242 + t43 * t244 + (t261 * t55 + t95 * t285) * t257 + t4;
t325 = 2 * m(5);
t324 = 2 * m(6);
t323 = t259 ^ 2;
t322 = t160 / 0.2e1;
t321 = -t162 / 0.2e1;
t320 = t208 / 0.2e1;
t319 = -t210 / 0.2e1;
t318 = t213 / 0.2e1;
t317 = -t229 / 0.2e1;
t316 = -t231 / 0.2e1;
t315 = -t266 / 0.2e1;
t314 = t256 / 0.2e1;
t313 = -t258 / 0.2e1;
t312 = t259 / 0.2e1;
t309 = pkin(7) * t229;
t308 = pkin(7) * t231;
t79 = rSges(6,1) * t118 + rSges(6,2) * t117 + rSges(6,3) * t160;
t307 = pkin(4) * t159 + pkin(8) * t160 + t79;
t80 = rSges(6,1) * t120 + rSges(6,2) * t119 - rSges(6,3) * t162;
t306 = pkin(4) * t161 - pkin(8) * t162 + t80;
t305 = Icges(3,4) * t261;
t304 = Icges(3,4) * t263;
t303 = Icges(4,6) * t261;
t302 = Icges(4,6) * t263;
t301 = t256 * t257;
t300 = t257 * t258;
t101 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t213;
t296 = pkin(4) * t212 + pkin(8) * t213 + t101;
t113 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t208;
t295 = pkin(4) * t209 + pkin(8) * t208 + t113;
t114 = rSges(6,1) * t166 + rSges(6,2) * t165 - rSges(6,3) * t210;
t294 = -pkin(4) * t326 - pkin(8) * t210 + t114;
t144 = rSges(6,1) * t215 + rSges(6,2) * t214 - rSges(6,3) * t266;
t293 = pkin(4) * t246 - pkin(8) * t266 + t144;
t157 = -pkin(2) * t229 + qJ(3) * t230 + qJD(3) * t241;
t158 = -pkin(2) * t231 + qJ(3) * t232 + qJD(3) * t243;
t292 = t157 * t301 + t158 * t300;
t202 = pkin(2) * t242 + qJ(3) * t241;
t203 = pkin(2) * t244 + qJ(3) * t243;
t291 = t202 * t301 + t203 * t300;
t201 = t259 * t203;
t216 = pkin(3) * t301 + pkin(7) * t244;
t290 = t259 * t216 + t201;
t217 = -pkin(3) * t300 + pkin(7) * t242;
t289 = -t202 - t217;
t247 = (pkin(2) * t261 - qJ(3) * t263) * t257;
t288 = -pkin(3) * t259 - pkin(7) * t299 - t247;
t287 = qJD(2) * t257;
t280 = t257 ^ 2 * t285;
t275 = t257 * (-t259 * rSges(4,1) - (-rSges(4,2) * t261 - rSges(4,3) * t263) * t257 - t247);
t218 = (-qJD(3) * t263 + (pkin(2) * t263 + qJ(3) * t261) * qJD(2)) * t257;
t274 = (-t218 - (-rSges(4,2) * t263 + rSges(4,3) * t261) * t287) * t257;
t273 = pkin(7) * t280;
t272 = t216 * t300 + t217 * t301 + t291;
t269 = t307 - t308;
t182 = rSges(5,1) * t246 + rSges(5,2) * t266 + rSges(5,3) * t299;
t268 = (-t182 + t288) * t257;
t267 = (t288 - t293) * t257;
t150 = t212 * rSges(5,1) - t213 * rSges(5,2) + rSges(5,3) * t277;
t265 = -t273 + (-t150 - t218) * t257;
t264 = -t273 + (-t218 - t296) * t257;
t240 = (rSges(3,1) * t263 - rSges(3,2) * t261) * t287;
t238 = (Icges(3,1) * t263 - t305) * t287;
t237 = (-Icges(3,2) * t261 + t304) * t287;
t236 = (-Icges(4,4) * t263 + Icges(4,5) * t261) * t287;
t235 = (Icges(3,5) * t263 - Icges(3,6) * t261) * t287;
t234 = (-Icges(4,2) * t263 + t303) * t287;
t233 = (Icges(4,3) * t261 - t302) * t287;
t226 = t259 * rSges(3,3) + (rSges(3,1) * t261 + rSges(3,2) * t263) * t257;
t225 = Icges(4,4) * t259 + (-Icges(4,2) * t261 - t302) * t257;
t224 = Icges(4,5) * t259 + (-Icges(4,3) * t263 - t303) * t257;
t223 = Icges(3,5) * t259 + (Icges(3,1) * t261 + t304) * t257;
t222 = Icges(3,6) * t259 + (Icges(3,2) * t263 + t305) * t257;
t219 = t259 * t309;
t200 = -rSges(3,1) * t231 - rSges(3,2) * t232;
t199 = rSges(4,2) * t231 + rSges(4,3) * t232;
t198 = -rSges(3,1) * t229 - rSges(3,2) * t230;
t197 = rSges(4,2) * t229 + rSges(4,3) * t230;
t196 = -Icges(3,1) * t231 - Icges(3,4) * t232;
t195 = -Icges(3,1) * t229 - Icges(3,4) * t230;
t194 = -Icges(3,4) * t231 - Icges(3,2) * t232;
t193 = -Icges(3,4) * t229 - Icges(3,2) * t230;
t192 = Icges(4,4) * t231 + Icges(4,5) * t232;
t191 = Icges(4,4) * t229 + Icges(4,5) * t230;
t190 = -Icges(3,5) * t231 - Icges(3,6) * t232;
t189 = -Icges(3,5) * t229 - Icges(3,6) * t230;
t188 = Icges(4,2) * t231 + Icges(4,6) * t232;
t187 = Icges(4,2) * t229 + Icges(4,6) * t230;
t186 = Icges(4,6) * t231 + Icges(4,3) * t232;
t185 = Icges(4,6) * t229 + Icges(4,3) * t230;
t178 = rSges(3,1) * t244 - rSges(3,2) * t243 + rSges(3,3) * t301;
t177 = rSges(3,1) * t242 - rSges(3,2) * t241 - rSges(3,3) * t300;
t176 = -rSges(4,1) * t300 - rSges(4,2) * t242 + rSges(4,3) * t241;
t175 = rSges(4,1) * t301 - rSges(4,2) * t244 + rSges(4,3) * t243;
t174 = Icges(3,1) * t244 - Icges(3,4) * t243 + Icges(3,5) * t301;
t173 = Icges(3,1) * t242 - Icges(3,4) * t241 - Icges(3,5) * t300;
t172 = Icges(3,4) * t244 - Icges(3,2) * t243 + Icges(3,6) * t301;
t171 = Icges(3,4) * t242 - Icges(3,2) * t241 - Icges(3,6) * t300;
t170 = -Icges(4,4) * t300 - Icges(4,2) * t242 + Icges(4,6) * t241;
t169 = Icges(4,4) * t301 - Icges(4,2) * t244 + Icges(4,6) * t243;
t168 = -Icges(4,5) * t300 - Icges(4,6) * t242 + Icges(4,3) * t241;
t167 = Icges(4,5) * t301 - Icges(4,6) * t244 + Icges(4,3) * t243;
t153 = t259 * t158;
t140 = -rSges(5,1) * t326 + rSges(5,2) * t210 + rSges(5,3) * t242;
t139 = rSges(5,1) * t209 - rSges(5,2) * t208 + rSges(5,3) * t244;
t131 = (t198 * t256 + t200 * t258) * t257;
t128 = rSges(5,1) * t161 + rSges(5,2) * t162 - rSges(5,3) * t229;
t127 = rSges(5,1) * t159 - rSges(5,2) * t160 - rSges(5,3) * t231;
t116 = (-t176 - t202) * t259 + t258 * t275;
t115 = t175 * t259 + t256 * t275 + t201;
t106 = t139 * t299 - t182 * t244;
t105 = -t140 * t299 + t182 * t242;
t104 = (-t157 - t197) * t259 + t258 * t274;
t103 = t199 * t259 + t256 * t274 + t153;
t102 = t179 * t299 + t180 * t266 + t181 * t246;
t97 = (t175 * t258 + t176 * t256) * t257 + t291;
t96 = -t139 * t242 + t140 * t244;
t93 = (t197 * t256 + t199 * t258) * t257 + t292;
t92 = (-t140 + t289) * t259 + t258 * t268;
t91 = t139 * t259 + t256 * t268 + t290;
t90 = t114 * t266 - t144 * t210;
t89 = -t113 * t266 - t144 * t208;
t88 = t134 * t299 + t136 * t266 + t138 * t246;
t87 = t133 * t299 + t135 * t266 + t137 * t246;
t86 = -t141 * t266 + t142 * t214 + t143 * t215;
t81 = (t139 * t258 + t140 * t256) * t257 + t272;
t72 = t113 * t210 + t114 * t208;
t71 = t219 + (-t128 - t157) * t259 + t265 * t258;
t70 = t153 + (t127 - t308) * t259 + t265 * t256;
t69 = -t293 * t244 + t295 * t299;
t68 = t293 * t242 - t294 * t299;
t65 = -t244 * t150 + t231 * t182 + (t127 * t261 + t139 * t285) * t257;
t64 = t242 * t150 - t229 * t182 + (-t128 * t261 - t140 * t285) * t257;
t63 = (t289 - t294) * t259 + t258 * t267;
t62 = t256 * t267 + t295 * t259 + t290;
t61 = t266 * t148 + t246 * t149 - t213 * t180 + t212 * t181 + (t147 * t261 + t179 * t285) * t257;
t60 = (t127 * t258 + t128 * t256 + (-t229 * t256 - t231 * t258) * pkin(7)) * t257 + t292;
t59 = -t295 * t242 + t294 * t244;
t58 = -t108 * t266 + t110 * t214 + t112 * t215;
t57 = -t107 * t266 + t109 * t214 + t111 * t215;
t56 = -t127 * t242 + t128 * t244 + t139 * t229 - t140 * t231;
t49 = (t294 * t256 + t295 * t258) * t257 + t272;
t48 = t266 * t124 + t246 * t126 - t213 * t136 + t212 * t138 + (t122 * t261 + t134 * t285) * t257;
t47 = t266 * t123 + t246 * t125 - t213 * t135 + t212 * t137 + (t121 * t261 + t133 * t285) * t257;
t46 = t219 + (-t157 - t306) * t259 + t264 * t258;
t45 = t264 * t256 + t269 * t259 + t153;
t40 = -t101 * t210 - t114 * t213 - t144 * t162 + t266 * t80;
t39 = -t101 * t208 + t113 * t213 - t144 * t160 - t266 * t79;
t38 = t100 * t215 + t141 * t213 + t142 * t145 + t143 * t146 + t214 * t99 - t266 * t98;
t37 = (t269 * t258 + (t306 - t309) * t256) * t257 + t292;
t36 = -t296 * t244 + t293 * t231 + (t307 * t261 + t295 * t285) * t257;
t35 = t296 * t242 - t293 * t229 + (-t306 * t261 - t294 * t285) * t257;
t34 = t113 * t162 + t114 * t160 + t208 * t80 + t210 * t79;
t31 = t259 * t86 + (t256 * t57 - t258 * t58) * t257;
t30 = t242 * t58 + t244 * t57 + t86 * t299;
t29 = t208 * t57 - t210 * t58 - t266 * t86;
t28 = t295 * t229 - t294 * t231 - t307 * t242 + t306 * t244;
t27 = t259 * t67 + (t256 * t52 - t258 * t53) * t257;
t26 = t259 * t66 + (t256 * t50 - t258 * t51) * t257;
t25 = t242 * t53 + t244 * t52 + t67 * t299;
t24 = t242 * t51 + t244 * t50 + t66 * t299;
t23 = t208 * t52 - t210 * t53 - t266 * t67;
t22 = t208 * t50 - t210 * t51 - t266 * t66;
t21 = t108 * t213 + t110 * t145 + t112 * t146 + t214 * t76 + t215 * t78 - t266 * t74;
t20 = t107 * t213 + t109 * t145 + t111 * t146 + t214 * t75 + t215 * t77 - t266 * t73;
t15 = t259 * t61 + (t256 * t47 - t258 * t48) * t257;
t14 = t259 * t55 + (t256 * t43 - t258 * t44) * t257;
t13 = t259 * t54 + (t256 * t41 - t258 * t42) * t257;
t12 = -t88 * t229 - t87 * t231 + t48 * t242 + t47 * t244 + (t102 * t285 + t261 * t61) * t257;
t9 = t259 * t38 + (t20 * t256 - t21 * t258) * t257;
t8 = t259 * t33 + (t18 * t256 - t19 * t258) * t257;
t7 = t259 * t32 + (t16 * t256 - t17 * t258) * t257;
t6 = t20 * t244 + t21 * t242 - t58 * t229 - t57 * t231 + (t261 * t38 + t86 * t285) * t257;
t5 = t160 * t57 - t162 * t58 + t20 * t208 - t21 * t210 + t213 * t86 - t266 * t38;
t2 = t160 * t52 - t162 * t53 + t18 * t208 - t19 * t210 + t213 * t67 - t266 * t33;
t1 = t16 * t208 + t160 * t50 - t162 * t51 - t17 * t210 + t213 * t66 - t266 * t32;
t10 = [0; m(3) * t131 + m(4) * t93 + m(5) * t60 + m(6) * t37; -t8 * t300 + t7 * t301 - t14 * t300 + t13 * t301 - ((t167 * t230 + t169 * t229 + t186 * t241 - t188 * t242 - t192 * t300) * t301 - (t168 * t230 + t170 * t229 + t185 * t241 - t187 * t242 - t191 * t300) * t300 + (t224 * t230 + t225 * t229 + t233 * t241 - t234 * t242 - t236 * t300) * t259) * t300 + ((-t172 * t232 - t174 * t231 + t190 * t301 - t194 * t243 + t244 * t196) * t301 - (-t171 * t232 - t173 * t231 + t189 * t301 - t193 * t243 + t244 * t195) * t300 + (-t222 * t232 - t223 * t231 + t235 * t301 - t237 * t243 + t238 * t244) * t259) * t301 + ((t167 * t232 + t169 * t231 + t186 * t243 - t188 * t244 + t192 * t301) * t301 - (t168 * t232 + t170 * t231 + t185 * t243 - t187 * t244 + t191 * t301) * t300 + (t224 * t232 + t225 * t231 + t233 * t243 - t234 * t244 + t236 * t301) * t259) * t301 - ((-t172 * t230 - t174 * t229 - t190 * t300 - t194 * t241 + t196 * t242) * t301 - (-t171 * t230 - t173 * t229 - t189 * t300 - t193 * t241 + t195 * t242) * t300 + (-t222 * t230 - t223 * t229 - t235 * t300 - t237 * t241 + t238 * t242) * t259) * t300 + (t37 * t49 + t45 * t62 + t46 * t63) * t324 + t259 * t9 + (t60 * t81 + t70 * t91 + t71 * t92) * t325 + t259 * t15 + 0.2e1 * m(4) * (t103 * t115 + t104 * t116 + t93 * t97) + t259 * (t323 * t235 + (((t194 * t263 + t196 * t261) * t256 - (t193 * t263 + t195 * t261) * t258 + ((-t172 * t261 + t174 * t263) * t256 - (-t171 * t261 + t173 * t263) * t258) * qJD(2)) * t257 + (-t189 * t258 + t190 * t256 + t237 * t263 + t238 * t261 + (-t222 * t261 + t223 * t263) * qJD(2)) * t259) * t257) + t259 * (t323 * t236 + (((-t186 * t263 - t188 * t261) * t256 - (-t185 * t263 - t187 * t261) * t258 + ((t167 * t261 - t169 * t263) * t256 - (t168 * t261 - t170 * t263) * t258) * qJD(2)) * t257 + (-t191 * t258 + t192 * t256 - t233 * t263 - t234 * t261 + (t224 * t261 - t225 * t263) * qJD(2)) * t259) * t257) + 0.2e1 * m(3) * ((-t177 * t259 - t226 * t300) * (-t198 * t259 - t240 * t300) + (t178 * t259 - t226 * t301) * (t200 * t259 - t240 * t301) + (t177 * t256 + t178 * t258) * t257 * t131); (m(4) + m(5) + m(6)) * t278; m(6) * (t230 * t62 + t232 * t63 + t241 * t45 + t243 * t46 + (-t263 * t37 + t49 * t286) * t257) + m(5) * (t230 * t91 + t232 * t92 + t241 * t70 + t243 * t71 + (-t263 * t60 + t81 * t286) * t257) + m(4) * (t241 * t103 + t243 * t104 + t230 * t115 + t232 * t116 + (-t263 * t93 + t97 * t286) * t257); 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t241 * t230 + t243 * t232 - t261 * t280); m(5) * t56 + m(6) * t28; t27 * t317 + t26 * t316 + (t7 / 0.2e1 + t13 / 0.2e1) * t244 + (t8 / 0.2e1 + t14 / 0.2e1) * t242 + m(6) * (t28 * t49 + t35 * t63 + t36 * t62 + t37 * t59 + t45 * t69 + t46 * t68) + m(5) * (t105 * t71 + t106 * t70 + t56 * t81 + t60 * t96 + t64 * t92 + t65 * t91) + (t6 / 0.2e1 + t12 / 0.2e1 + t95 * t317 + t94 * t316) * t259 + ((t256 * t84 - t258 * t85) * t317 + (t256 * t82 - t258 * t83) * t316 + (t9 / 0.2e1 + t15 / 0.2e1) * t261 + (t31 / 0.2e1 + t102 * t312 + (t256 * t87 - t258 * t88) * t257 / 0.2e1) * t285 + t328 * t314 + t327 * t313) * t257; m(5) * (t105 * t232 + t106 * t230 + t65 * t241 + t64 * t243 + (-t263 * t56 + t96 * t286) * t257) + m(6) * (t69 * t230 + t68 * t232 + t36 * t241 + t35 * t243 + (-t263 * t28 + t59 * t286) * t257); (t12 + t6) * t299 + t328 * t244 + t327 * t242 + (-t242 * t83 - t244 * t82 - t94 * t299 - t24) * t231 + (-t242 * t85 - t244 * t84 - t95 * t299 - t25) * t229 + (t102 * t299 + t242 * t88 + t244 * t87 + t30) * t277 + (t28 * t59 + t35 * t68 + t36 * t69) * t324 + (t105 * t64 + t106 * t65 + t56 * t96) * t325; m(6) * t34; m(6) * (t34 * t49 + t37 * t72 + t39 * t62 + t40 * t63 + t45 * t89 + t46 * t90) + t31 * t318 + t9 * t315 + t5 * t312 + t27 * t321 + t8 * t319 + t26 * t322 + t7 * t320 + (t1 * t314 + t2 * t313) * t257; m(6) * (t89 * t230 + t90 * t232 + t39 * t241 + t40 * t243 + (-t263 * t34 + t72 * t286) * t257); m(6) * (t28 * t72 + t34 * t59 + t35 * t90 + t36 * t89 + t39 * t69 + t40 * t68) + t25 * t321 + t4 * t319 + t24 * t322 + t3 * t320 + t22 * t316 + t244 * t1 / 0.2e1 + t23 * t317 + t242 * t2 / 0.2e1 + t30 * t318 + t6 * t315 + (t29 * t285 / 0.2e1 + t261 * t5 / 0.2e1) * t257; (t34 * t72 + t39 * t89 + t40 * t90) * t324 + t160 * t22 + t208 * t1 - t162 * t23 - t210 * t2 + t213 * t29 - t266 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
