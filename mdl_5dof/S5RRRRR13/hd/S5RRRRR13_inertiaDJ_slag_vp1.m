% Calculate time derivative of joint inertia matrix for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR13_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR13_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:24
% EndTime: 2024-09-27 17:30:30
% DurationCPUTime: 4.89s
% Computational Cost: add. (34660->521), mult. (20140->722), div. (0->0), fcn. (16146->16), ass. (0->266)
t250 = qJD(4) + qJD(5);
t320 = t250 / 0.2e1;
t252 = qJ(4) + qJ(5);
t241 = pkin(5) + t252;
t232 = sin(t241);
t326 = t232 * t320;
t251 = qJD(1) + qJD(2);
t240 = qJD(3) + t251;
t255 = cos(pkin(5));
t258 = cos(qJ(4));
t294 = qJD(4) * t258;
t289 = pkin(4) * t294;
t275 = t255 * t289;
t239 = t258 * pkin(4) + pkin(3);
t317 = pkin(3) - t239;
t325 = t317 * t240 - t275;
t324 = 2 * m(3);
t323 = 2 * m(4);
t322 = 2 * m(5);
t321 = 2 * m(6);
t254 = sin(pkin(5));
t249 = t254 ^ 2;
t253 = qJ(1) + qJ(2);
t244 = sin(t253);
t319 = pkin(2) * t244;
t247 = qJ(3) + t253;
t238 = cos(t247);
t230 = t238 * pkin(3);
t257 = sin(qJ(1));
t318 = t257 * pkin(1);
t316 = pkin(1) * qJD(1);
t256 = sin(qJ(4));
t315 = Icges(5,4) * t256;
t314 = Icges(5,4) * t258;
t234 = cos(t241);
t217 = t234 * t320;
t242 = pkin(5) - t252;
t235 = cos(t242);
t303 = t250 * t235;
t179 = t217 - t303 / 0.2e1;
t233 = sin(t242);
t218 = t233 * t320;
t181 = t218 + t326;
t152 = rSges(6,1) * t181 + rSges(6,2) * t179;
t237 = sin(t247);
t313 = t152 * t237;
t312 = t237 * t240;
t311 = t237 * t254;
t260 = pkin(9) + pkin(10);
t301 = t255 * t256;
t212 = pkin(4) * t301 - t254 * t260;
t310 = t238 * t212;
t309 = t238 * t240;
t308 = t238 * t254;
t307 = t240 * t254;
t306 = t244 * t251;
t246 = cos(t253);
t305 = t246 * t251;
t302 = t254 * t256;
t300 = t255 * t258;
t228 = t235 / 0.2e1;
t210 = t228 + t234 / 0.2e1;
t243 = sin(t252);
t163 = t210 * t238 - t237 * t243;
t227 = t232 / 0.2e1;
t209 = t227 - t233 / 0.2e1;
t245 = cos(t252);
t164 = t238 * t209 + t237 * t245;
t123 = t164 * rSges(6,1) + t163 * rSges(6,2) - rSges(6,3) * t308;
t165 = -t210 * t237 - t238 * t243;
t166 = -t237 * t209 + t238 * t245;
t124 = t166 * rSges(6,1) + t165 * rSges(6,2) + rSges(6,3) * t311;
t71 = t123 * t311 + t124 * t308;
t208 = t227 + t233 / 0.2e1;
t211 = t228 - t234 / 0.2e1;
t158 = t211 * rSges(6,1) + t208 * rSges(6,2) + t255 * rSges(6,3);
t196 = pkin(4) * t302 + (-pkin(9) + t260) * t255;
t299 = -t158 - t196;
t295 = qJD(4) * t256;
t290 = pkin(4) * t295;
t298 = -t212 * t312 - t237 * t290;
t297 = pkin(9) * t311 + t230;
t296 = qJD(4) * t254;
t293 = pkin(2) * t306;
t292 = pkin(2) * t305;
t284 = t238 * t307;
t269 = -t210 * t240 - t245 * t250;
t278 = t240 * t243 - t218 + t326;
t110 = t278 * t237 + t269 * t238;
t270 = -t209 * t240 - t243 * t250;
t279 = t240 * t245 + t217 + t303 / 0.2e1;
t111 = -t279 * t237 + t270 * t238;
t67 = t111 * rSges(6,1) + t110 * rSges(6,2) + rSges(6,3) * t284;
t112 = t269 * t237 - t278 * t238;
t113 = t270 * t237 + t279 * t238;
t286 = t237 * t307;
t68 = t113 * rSges(6,1) + t112 * rSges(6,2) + rSges(6,3) * t286;
t291 = t123 * t284 + t67 * t308 + t68 * t311;
t288 = t257 * t316;
t259 = cos(qJ(1));
t287 = t259 * t316;
t285 = t158 * t311;
t186 = -t237 * t256 + t238 * t300;
t267 = t237 * t301 - t238 * t258;
t141 = t267 * qJD(4) - t186 * t240;
t187 = t237 * t258 + t238 * t301;
t188 = -t237 * t300 - t238 * t256;
t264 = t188 * qJD(4);
t142 = -t187 * t240 + t264;
t89 = t142 * rSges(5,1) + t141 * rSges(5,2) + rSges(5,3) * t284;
t138 = -rSges(5,1) * t267 + t188 * rSges(5,2) + rSges(5,3) * t311;
t283 = t311 / 0.2e1;
t282 = -t308 / 0.2e1;
t281 = t307 / 0.2e1;
t214 = t246 * rSges(3,1) - t244 * rSges(3,2);
t204 = t238 * rSges(4,1) - t237 * rSges(4,2);
t280 = t299 * t254;
t277 = -t237 * t212 + t238 * t239;
t276 = t249 * t289;
t274 = t237 * t281;
t273 = t238 * t281;
t194 = -rSges(3,1) * t305 + rSges(3,2) * t306;
t173 = -rSges(4,1) * t309 + rSges(4,2) * t312;
t236 = pkin(2) * t246;
t185 = t204 + t236;
t213 = -t244 * rSges(3,1) - t246 * rSges(3,2);
t203 = -t237 * rSges(4,1) - t238 * rSges(4,2);
t143 = -t187 * qJD(4) + t188 * t240;
t144 = t186 * qJD(4) - t267 * t240;
t272 = -t144 * rSges(5,1) - t143 * rSges(5,2);
t149 = Icges(6,5) * t181 + Icges(6,6) * t179;
t150 = Icges(6,4) * t181 + Icges(6,2) * t179;
t151 = Icges(6,1) * t181 + Icges(6,4) * t179;
t156 = Icges(6,4) * t211 + Icges(6,2) * t208 + Icges(6,6) * t255;
t157 = Icges(6,1) * t211 + Icges(6,4) * t208 + Icges(6,5) * t255;
t271 = t255 * t149 + t208 * t150 + t211 * t151 + t179 * t156 + t181 * t157;
t126 = t138 + t297;
t268 = -t237 * t239 - t310;
t137 = t187 * rSges(5,1) + t186 * rSges(5,2) - rSges(5,3) * t308;
t116 = t236 + t126;
t120 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t311;
t122 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t311;
t61 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t284;
t63 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t284;
t65 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t284;
t10 = t179 * t120 + t181 * t122 + t208 * t63 + t211 * t65 + t255 * t61;
t117 = Icges(6,5) * t164 + Icges(6,6) * t163 - Icges(6,3) * t308;
t118 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t311;
t119 = Icges(6,4) * t164 + Icges(6,2) * t163 - Icges(6,6) * t308;
t121 = Icges(6,1) * t164 + Icges(6,4) * t163 - Icges(6,5) * t308;
t155 = Icges(6,5) * t211 + Icges(6,6) * t208 + Icges(6,3) * t255;
t18 = t110 * t156 + t111 * t157 + t165 * t150 + t166 * t151 + (t149 * t237 + t155 * t309) * t254;
t19 = t112 * t156 + t113 * t157 + t163 * t150 + t164 * t151 + (-t149 * t238 + t155 * t312) * t254;
t22 = -t117 * t308 + t163 * t119 + t164 * t121;
t23 = -t118 * t308 + t163 * t120 + t164 * t122;
t24 = t117 * t311 + t165 * t119 + t166 * t121;
t25 = t118 * t311 + t165 * t120 + t166 * t122;
t26 = t271 * t255;
t36 = t255 * t117 + t208 * t119 + t211 * t121;
t37 = t255 * t118 + t208 * t120 + t211 * t122;
t62 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t286;
t64 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t286;
t66 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t286;
t69 = -t155 * t308 + t163 * t156 + t164 * t157;
t70 = t155 * t311 + t165 * t156 + t166 * t157;
t9 = t179 * t119 + t181 * t121 + t208 * t64 + t211 * t66 + t255 * t62;
t266 = -(t19 * t255 + ((t112 * t120 + t113 * t122 + t163 * t63 + t164 * t65) * t237 + t23 * t309 - (t112 * t119 + t113 * t121 + t163 * t64 + t164 * t66) * t238 + t22 * t312 + ((t118 * t312 - t238 * t61) * t237 - (t117 * t312 - t62 * t238) * t238) * t254) * t254) * t308 + (t18 * t255 + ((t110 * t120 + t111 * t122 + t165 * t63 + t166 * t65) * t237 + t25 * t309 - (t110 * t119 + t111 * t121 + t165 * t64 + t166 * t66) * t238 + t24 * t312 + ((t118 * t309 + t237 * t61) * t237 - (t117 * t309 + t237 * t62) * t238) * t254) * t254) * t311 + t255 * (t26 + ((t240 * t37 - t9) * t238 + (t240 * t36 + t10) * t237) * t254) + (t69 * t255 + (-t22 * t238 + t23 * t237) * t254) * t286 + (t70 * t255 + (t237 * t25 - t238 * t24) * t254) * t284;
t193 = t213 * t251;
t172 = t203 * t240;
t98 = t277 + t124;
t184 = t203 - t319;
t207 = pkin(9) * t284;
t76 = -pkin(3) * t312 + t207 + t89;
t265 = t26 + (t10 + t18) * t283 + (t19 + t9) * t282 + (t36 + t69) * t274 + (t37 + t70) * t273;
t94 = t236 + t98;
t160 = t173 - t292;
t191 = Icges(5,6) * t255 + (Icges(5,2) * t258 + t315) * t254;
t192 = Icges(5,5) * t255 + (Icges(5,1) * t256 + t314) * t254;
t197 = (Icges(5,5) * t258 - Icges(5,6) * t256) * t296;
t198 = (-Icges(5,2) * t256 + t314) * t296;
t199 = (Icges(5,1) * t258 - t315) * t296;
t263 = t255 * t197 + t199 * t302 + (-t191 * t295 + t192 * t294 + t198 * t258) * t254;
t225 = pkin(9) * t308;
t125 = -t237 * pkin(3) - t137 + t225;
t159 = t172 - t293;
t74 = t76 - t293;
t97 = -t123 + t268;
t115 = t125 - t319;
t93 = t97 - t319;
t262 = t263 + t271;
t77 = (-t230 + (-rSges(5,3) - pkin(9)) * t311) * t240 + t272;
t133 = Icges(5,4) * t187 + Icges(5,2) * t186 - Icges(5,6) * t308;
t135 = Icges(5,1) * t187 + Icges(5,4) * t186 - Icges(5,5) * t308;
t82 = Icges(5,5) * t144 + Icges(5,6) * t143 + Icges(5,3) * t286;
t84 = Icges(5,4) * t144 + Icges(5,2) * t143 + Icges(5,6) * t286;
t86 = Icges(5,1) * t144 + Icges(5,4) * t143 + Icges(5,5) * t286;
t14 = t255 * t82 + (t256 * t86 + t258 * t84 + (-t133 * t256 + t135 * t258) * qJD(4)) * t254;
t134 = -Icges(5,4) * t267 + Icges(5,2) * t188 + Icges(5,6) * t311;
t136 = -Icges(5,1) * t267 + Icges(5,4) * t188 + Icges(5,5) * t311;
t81 = Icges(5,5) * t142 + Icges(5,6) * t141 + Icges(5,3) * t284;
t83 = Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t284;
t85 = Icges(5,1) * t142 + Icges(5,4) * t141 + Icges(5,5) * t284;
t15 = t255 * t81 + (t256 * t85 + t258 * t83 + (-t134 * t256 + t136 * t258) * qJD(4)) * t254;
t190 = Icges(5,3) * t255 + (Icges(5,5) * t256 + Icges(5,6) * t258) * t254;
t32 = t141 * t191 + t142 * t192 + t188 * t198 - t267 * t199 + (t190 * t309 + t197 * t237) * t254;
t33 = t143 * t191 + t144 * t192 + t186 * t198 + t187 * t199 + (t190 * t312 - t197 * t238) * t254;
t131 = Icges(5,5) * t187 + Icges(5,6) * t186 - Icges(5,3) * t308;
t59 = t255 * t131 + (t133 * t258 + t135 * t256) * t254;
t132 = -Icges(5,5) * t267 + Icges(5,6) * t188 + Icges(5,3) * t311;
t60 = t255 * t132 + (t134 * t258 + t136 * t256) * t254;
t80 = t263 * t255;
t95 = t186 * t191 + t187 * t192 - t190 * t308;
t96 = t188 * t191 + t190 * t311 - t192 * t267;
t261 = t265 + t80 + (t15 + t32) * t283 + (t14 + t33) * t282 + (t59 + t95) * t274 + (t60 + t96) * t273;
t75 = t77 - t292;
t42 = pkin(4) * t264 + t268 * t240 + t67;
t43 = (-t239 * t240 - t275) * t238 - t68 - t298;
t40 = t42 - t293;
t41 = t43 - t292;
t248 = t259 * pkin(1);
t202 = t214 + t248;
t201 = t213 - t318;
t200 = (rSges(5,1) * t258 - rSges(5,2) * t256) * t296;
t195 = t255 * rSges(5,3) + (rSges(5,1) * t256 + rSges(5,2) * t258) * t254;
t171 = t194 - t287;
t170 = t193 - t288;
t169 = t185 + t248;
t168 = t184 - t318;
t154 = t160 - t287;
t153 = t159 - t288;
t148 = t277 - t297;
t147 = -t317 * t237 + t225 + t310;
t146 = t240 * t285;
t114 = t255 * t124;
t109 = t248 + t116;
t108 = t115 - t318;
t103 = -t255 * t137 - t195 * t308;
t102 = t255 * t138 - t195 * t311;
t100 = -pkin(9) * t286 - t325 * t238 + t298;
t99 = -t207 + (-t212 * t240 - t290) * t238 + t325 * t237;
t92 = t248 + t94;
t91 = t93 - t318;
t90 = rSges(5,3) * t286 - t272;
t88 = -t255 * t123 - t158 * t308;
t87 = t114 - t285;
t73 = t75 - t287;
t72 = t74 - t288;
t58 = t255 * t67;
t57 = t255 * t89 + (-t195 * t309 - t200 * t237) * t254;
t56 = -t255 * t90 + (t195 * t312 - t200 * t238) * t254;
t49 = (-t123 - t147) * t255 + t238 * t280;
t48 = t255 * t148 + t237 * t280 + t114;
t47 = t132 * t311 + t188 * t134 - t136 * t267;
t46 = t131 * t311 + t188 * t133 - t135 * t267;
t45 = -t132 * t308 + t186 * t134 + t187 * t136;
t44 = -t131 * t308 + t186 * t133 + t187 * t135;
t39 = t41 - t287;
t38 = t40 - t288;
t31 = (t147 * t237 + t148 * t238) * t254 + t71;
t30 = t58 + (-t158 * t309 - t313) * t254;
t29 = -t152 * t308 - t255 * t68 + t146;
t21 = -t237 * t276 + t255 * t99 + t58 + (t299 * t309 - t313) * t254;
t20 = t196 * t286 + t146 + (-t100 - t68) * t255 + (-t152 * t254 - t276) * t238;
t11 = -t124 * t286 + t291;
t4 = ((t147 * t240 + t99) * t238 + (t100 + (-t124 - t148) * t240) * t237) * t254 + t291;
t1 = [(t38 * t92 + t39 * t91) * t321 + (t108 * t73 + t109 * t72) * t322 + (t153 * t169 + t154 * t168) * t323 + (t170 * t202 + t171 * t201) * t324 + t262; m(6) * (t38 * t94 + t39 * t93 + t40 * t92 + t41 * t91) + m(5) * (t108 * t75 + t109 * t74 + t115 * t73 + t116 * t72) + m(4) * (t153 * t185 + t154 * t184 + t159 * t169 + t160 * t168) + m(3) * (t170 * t214 + t171 * t213 + t193 * t202 + t194 * t201) + t262; (t40 * t94 + t41 * t93) * t321 + (t115 * t75 + t116 * t74) * t322 + (t159 * t185 + t160 * t184) * t323 + (t193 * t214 + t194 * t213) * t324 + t262; m(6) * (t38 * t98 + t39 * t97 + t42 * t92 + t43 * t91) + m(5) * (t108 * t77 + t109 * t76 + t125 * t73 + t126 * t72) + m(4) * (t153 * t204 + t154 * t203 + t168 * t173 + t169 * t172) + t262; m(6) * (t40 * t98 + t41 * t97 + t42 * t94 + t43 * t93) + m(5) * (t115 * t77 + t116 * t76 + t125 * t75 + t126 * t74) + m(4) * (t159 * t204 + t160 * t203 + t172 * t185 + t173 * t184) + t262; (t42 * t98 + t43 * t97) * t321 + (t125 * t77 + t126 * t76) * t322 + (t172 * t204 + t173 * t203) * t323 + t262; m(6) * (t20 * t91 + t21 * t92 + t38 * t48 + t39 * t49) + m(5) * (t102 * t72 + t103 * t73 + t108 * t56 + t109 * t57) + t261; m(6) * (t20 * t93 + t21 * t94 + t40 * t48 + t41 * t49) + m(5) * (t102 * t74 + t103 * t75 + t115 * t56 + t116 * t57) + t261; m(6) * (t20 * t97 + t21 * t98 + t42 * t48 + t43 * t49) + m(5) * (t102 * t76 + t103 * t77 + t125 * t56 + t126 * t57) + t261; (t102 * t57 + t103 * t56 + (t137 * t237 + t138 * t238) * ((t137 * t240 + t89) * t238 + (-t138 * t240 + t90) * t237) * t249) * t322 + (t96 * t255 + (t237 * t47 - t238 * t46) * t254) * t284 + (t32 * t255 + ((t141 * t134 + t142 * t136 + t188 * t83 - t267 * t85) * t237 + t47 * t309 - (t141 * t133 + t142 * t135 + t188 * t84 - t267 * t86) * t238 + t46 * t312 + ((t132 * t309 + t237 * t81) * t237 - (t131 * t309 + t237 * t82) * t238) * t254) * t254) * t311 + (t95 * t255 + (t237 * t45 - t238 * t44) * t254) * t286 - (t33 * t255 + ((t143 * t134 + t144 * t136 + t186 * t83 + t187 * t85) * t237 + t45 * t309 - (t143 * t133 + t144 * t135 + t186 * t84 + t187 * t86) * t238 + t44 * t312 + ((t132 * t312 - t238 * t81) * t237 - (t131 * t312 - t82 * t238) * t238) * t254) * t254) * t308 + t255 * (t80 + ((t240 * t60 - t14) * t238 + (t240 * t59 + t15) * t237) * t254) + (t20 * t49 + t21 * t48 + t31 * t4) * t321 + t266; m(6) * (t29 * t91 + t30 * t92 + t38 * t87 + t39 * t88) + t265; m(6) * (t29 * t93 + t30 * t94 + t40 * t87 + t41 * t88) + t265; m(6) * (t29 * t97 + t30 * t98 + t42 * t87 + t43 * t88) + t265; m(6) * (t11 * t31 + t20 * t88 + t21 * t87 + t29 * t49 + t30 * t48 + t4 * t71) + t266; (t11 * t71 + t29 * t88 + t30 * t87) * t321 + t266;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
