% Calculate time derivative of joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:28
% DurationCPUTime: 5.88s
% Computational Cost: add. (9492->462), mult. (9772->657), div. (0->0), fcn. (7678->8), ass. (0->274)
t198 = sin(qJ(1));
t200 = cos(qJ(1));
t196 = qJ(2) + qJ(3);
t186 = qJ(4) + t196;
t177 = sin(t186);
t328 = rSges(5,2) * t177;
t178 = cos(t186);
t332 = rSges(5,1) * t178;
t247 = -t328 + t332;
t120 = t198 * rSges(5,3) + t200 * t247;
t188 = t198 * rSges(4,3);
t184 = sin(t196);
t185 = cos(t196);
t333 = rSges(4,1) * t185;
t248 = -rSges(4,2) * t184 + t333;
t132 = t200 * t248 + t188;
t201 = -pkin(6) - pkin(5);
t197 = sin(qJ(2));
t276 = qJD(2) * t197;
t269 = pkin(2) * t276;
t359 = qJD(1) * t201 + t269;
t327 = rSges(5,2) * t178;
t151 = rSges(5,1) * t177 + t327;
t192 = qJD(2) + qJD(3);
t183 = qJD(4) + t192;
t220 = t151 * t183;
t199 = cos(qJ(2));
t169 = rSges(3,1) * t197 + rSges(3,2) * t199;
t218 = qJD(2) * t169;
t358 = t198 * t218;
t322 = Icges(3,4) * t199;
t239 = -Icges(3,2) * t197 + t322;
t141 = Icges(3,6) * t198 + t200 * t239;
t323 = Icges(3,4) * t197;
t243 = Icges(3,1) * t199 - t323;
t143 = Icges(3,5) * t198 + t200 * t243;
t227 = t141 * t197 - t143 * t199;
t357 = t198 * t227;
t320 = Icges(4,4) * t185;
t237 = -Icges(4,2) * t184 + t320;
t128 = Icges(4,6) * t198 + t200 * t237;
t321 = Icges(4,4) * t184;
t241 = Icges(4,1) * t185 - t321;
t130 = Icges(4,5) * t198 + t200 * t241;
t229 = t128 * t184 - t130 * t185;
t356 = t198 * t229;
t318 = Icges(5,4) * t178;
t236 = -Icges(5,2) * t177 + t318;
t112 = Icges(5,6) * t198 + t200 * t236;
t319 = Icges(5,4) * t177;
t240 = Icges(5,1) * t178 - t319;
t114 = Icges(5,5) * t198 + t200 * t240;
t231 = t112 * t177 - t114 * t178;
t355 = t198 * t231;
t179 = t199 * pkin(2) + pkin(1);
t336 = pkin(1) - t179;
t354 = t198 * t336;
t140 = -Icges(3,6) * t200 + t198 * t239;
t142 = -Icges(3,5) * t200 + t198 * t243;
t228 = t140 * t197 - t142 * t199;
t353 = t200 * t228;
t127 = -Icges(4,6) * t200 + t198 * t237;
t129 = -Icges(4,5) * t200 + t198 * t241;
t230 = t127 * t184 - t129 * t185;
t352 = t200 * t230;
t111 = -Icges(5,6) * t200 + t198 * t236;
t113 = -Icges(5,5) * t200 + t198 * t240;
t232 = t111 * t177 - t113 * t178;
t351 = t200 * t232;
t233 = Icges(5,5) * t178 - Icges(5,6) * t177;
t109 = -Icges(5,3) * t200 + t198 * t233;
t350 = qJD(1) * t109;
t234 = Icges(4,5) * t185 - Icges(4,6) * t184;
t125 = -Icges(4,3) * t200 + t198 * t234;
t349 = qJD(1) * t125;
t235 = Icges(3,5) * t199 - Icges(3,6) * t197;
t138 = -Icges(3,3) * t200 + t198 * t235;
t148 = Icges(5,2) * t178 + t319;
t149 = Icges(5,1) * t177 + t318;
t226 = t148 * t177 - t149 * t178;
t348 = qJD(1) * t226 + t233 * t183;
t154 = Icges(4,2) * t185 + t321;
t155 = Icges(4,1) * t184 + t320;
t225 = t154 * t184 - t155 * t185;
t347 = qJD(1) * t225 + t234 * t192;
t193 = -pkin(7) + t201;
t284 = t193 - t201;
t161 = pkin(3) * t185 + t179;
t289 = t161 - t179;
t92 = t198 * t289 + t200 * t284;
t346 = 2 * m(3);
t345 = 2 * m(4);
t344 = 2 * m(5);
t194 = t198 ^ 2;
t195 = t200 ^ 2;
t343 = t198 / 0.2e1;
t342 = -t200 / 0.2e1;
t341 = m(3) * t169;
t329 = rSges(4,2) * t185;
t156 = rSges(4,1) * t184 + t329;
t340 = m(4) * t156;
t339 = pkin(2) * t197;
t338 = t198 * pkin(5);
t191 = t200 * pkin(5);
t303 = t148 * t183;
t60 = qJD(1) * t112 - t198 * t303;
t259 = t113 * t183 + t60;
t302 = t149 * t183;
t62 = qJD(1) * t114 - t198 * t302;
t261 = t111 * t183 - t62;
t27 = -t109 * t200 - t198 * t232;
t110 = Icges(5,3) * t198 + t200 * t233;
t28 = -t110 * t200 - t355;
t282 = qJD(1) * t110;
t298 = t178 * t183;
t299 = t177 * t183;
t147 = Icges(5,5) * t177 + Icges(5,6) * t178;
t214 = t183 * t147;
t57 = -t200 * t214 - t350;
t58 = -t198 * t214 + t282;
t59 = -t111 * qJD(1) - t200 * t303;
t61 = -t113 * qJD(1) - t200 * t302;
t2 = (t200 * t58 + (t28 + t351) * qJD(1)) * t200 + (t27 * qJD(1) + (-t112 * t298 - t114 * t299 - t177 * t59 + t178 * t61 + t282) * t198 + (-t57 + t261 * t178 + t259 * t177 + (-t109 - t231) * qJD(1)) * t200) * t198;
t337 = t200 * t2;
t335 = -pkin(5) - t201;
t119 = -rSges(5,3) * t200 + t198 * t247;
t65 = t198 * t119 + t200 * t120;
t334 = rSges(3,1) * t199;
t331 = rSges(3,2) * t197;
t326 = rSges(3,3) * t200;
t189 = t198 * rSges(3,3);
t325 = rSges(5,3) - t193;
t152 = t200 * t161;
t170 = t200 * t179;
t93 = -t198 * t284 + t152 - t170;
t324 = -t120 - t93;
t137 = t248 * t192;
t308 = t137 * t198;
t307 = t140 * t199;
t306 = t141 * t199;
t305 = t142 * t197;
t304 = t143 * t197;
t301 = t154 * t192;
t300 = t155 * t192;
t297 = t183 * t200;
t296 = t184 * t192;
t295 = t185 * t192;
t294 = t192 * t198;
t293 = t192 * t200;
t292 = t198 * t193;
t123 = t200 * t201 + t191 - t354;
t124 = -t200 * pkin(1) + t198 * t335 + t170;
t291 = t198 * t123 + t200 * t124;
t131 = -t200 * rSges(4,3) + t198 * t248;
t68 = t198 * t131 + t200 * t132;
t279 = qJD(1) * t198;
t267 = t184 * t279;
t290 = pkin(3) * t267 + t151 * t279;
t278 = qJD(1) * t200;
t288 = rSges(5,3) * t278 + t279 * t328;
t287 = rSges(4,2) * t267 + rSges(4,3) * t278;
t286 = t359 * t198;
t285 = t200 * t334 + t189;
t283 = t194 + t195;
t126 = Icges(4,3) * t198 + t200 * t234;
t281 = qJD(1) * t126;
t139 = Icges(3,3) * t198 + t200 * t235;
t280 = qJD(1) * t139;
t275 = qJD(2) * t199;
t258 = t114 * t183 + t59;
t260 = -t112 * t183 + t61;
t29 = t109 * t198 - t351;
t30 = t110 * t198 - t200 * t231;
t274 = t198 * ((t198 * t57 + (t29 + t355) * qJD(1)) * t198 + (t30 * qJD(1) + (t111 * t298 + t113 * t299 + t177 * t60 - t178 * t62 - t350) * t200 + (-t58 + t260 * t178 - t258 * t177 + (t110 - t232) * qJD(1)) * t198) * t200) + (t198 * t28 - t200 * t27) * t279 + (t198 * t30 - t200 * t29) * t278;
t273 = t198 * (t120 * qJD(1) - t198 * t220) + t200 * (-t297 * t327 + (-t177 * t297 - t178 * t279) * rSges(5,1) + t288) + t119 * t278;
t253 = t200 * t269;
t272 = t198 * ((-t200 * t336 - t338) * qJD(1) - t286) + t200 * (-t253 + (t200 * t335 + t354) * qJD(1)) + t123 * t278;
t219 = t156 * t192;
t271 = t131 * t278 + t198 * (t132 * qJD(1) - t198 * t219) + t200 * (-t293 * t329 + (-t184 * t293 - t185 * t279) * rSges(4,1) + t287);
t270 = t200 * t331;
t268 = pkin(2) * t275;
t266 = t197 * t279;
t265 = t279 / 0.2e1;
t264 = t278 / 0.2e1;
t263 = -t156 - t339;
t262 = -pkin(3) * t184 - t151;
t23 = t198 * t92 + t200 * t93 + t65;
t74 = qJD(1) * t130 - t198 * t300;
t257 = t127 * t192 - t74;
t73 = -t129 * qJD(1) - t200 * t300;
t256 = -t128 * t192 + t73;
t72 = qJD(1) * t128 - t198 * t301;
t255 = t129 * t192 + t72;
t71 = -t127 * qJD(1) - t200 * t301;
t254 = t130 * t192 + t71;
t118 = t247 * t183;
t252 = -pkin(3) * t295 - t118;
t251 = t274 - t337;
t31 = -t125 * t200 - t198 * t230;
t32 = -t126 * t200 - t356;
t33 = t125 * t198 - t352;
t34 = t126 * t198 - t200 * t229;
t153 = Icges(4,5) * t184 + Icges(4,6) * t185;
t211 = t192 * t153;
t69 = -t200 * t211 - t349;
t70 = -t198 * t211 + t281;
t250 = (t198 * t32 - t200 * t31) * t279 + (t198 * t34 - t200 * t33) * t278 + t198 * ((t198 * t69 + (t33 + t356) * qJD(1)) * t198 + (t34 * qJD(1) + (t127 * t295 + t129 * t296 + t184 * t72 - t185 * t74 - t349) * t200 + (-t70 + t256 * t185 - t254 * t184 + (t126 - t230) * qJD(1)) * t198) * t200) + t274;
t249 = -t331 + t334;
t150 = -pkin(3) * t296 - t269;
t144 = t200 * t150;
t244 = t198 * (t150 * t198 + (t200 * t289 - t292) * qJD(1) + t286) + t200 * (-t92 * qJD(1) + t144 + t253) + t92 * t278 + t273;
t242 = Icges(3,1) * t197 + t322;
t238 = Icges(3,2) * t199 + t323;
t224 = -pkin(1) - t249;
t223 = t262 - t339;
t222 = -t179 - t248;
t221 = -t161 - t247;
t116 = t236 * t183;
t117 = t240 * t183;
t203 = qJD(1) * t147 + (t117 - t303) * t178 + (-t116 - t302) * t177;
t217 = (t177 * t260 + t178 * t258 + t348 * t198 + t203 * t200) * t343 + (-t177 * t261 + t178 * t259 + t203 * t198 - t348 * t200) * t342 + (t111 * t178 + t113 * t177 - t147 * t200 - t198 * t226) * t265 + (t112 * t178 + t114 * t177 + t147 * t198 - t200 * t226) * t264;
t210 = qJD(2) * t242;
t209 = qJD(2) * t238;
t208 = qJD(2) * (-Icges(3,5) * t197 - Icges(3,6) * t199);
t95 = t223 * t200;
t207 = t252 - t268;
t4 = (t200 * t70 + (t32 + t352) * qJD(1)) * t200 + (t31 * qJD(1) + (-t128 * t295 - t130 * t296 - t184 * t71 + t185 * t73 + t281) * t198 + (-t69 + t257 * t185 + t255 * t184 + (-t125 - t229) * qJD(1)) * t200) * t198;
t206 = (-t4 - t2) * t200 + t250;
t205 = rSges(3,2) * t266 + rSges(3,3) * t278 - t200 * t218;
t134 = t237 * t192;
t135 = t241 * t192;
t202 = qJD(1) * t153 + (t135 - t301) * t185 + (-t134 - t300) * t184;
t204 = t217 + (t184 * t256 + t185 * t254 + t347 * t198 + t202 * t200) * t343 + (-t184 * t257 + t185 * t255 + t202 * t198 - t347 * t200) * t342 + (t127 * t185 + t129 * t184 - t153 * t200 - t198 * t225) * t265 + (t128 * t185 + t130 * t184 + t153 * t198 - t200 * t225) * t264;
t174 = pkin(2) * t266;
t160 = t249 * qJD(2);
t146 = -t270 + t285;
t145 = t198 * t249 - t326;
t122 = t263 * t200;
t121 = t263 * t198;
t104 = t338 + (pkin(1) - t331) * t200 + t285;
t103 = t198 * t224 + t191 + t326;
t101 = t262 * t200;
t100 = t262 * t198;
t94 = t223 * t198;
t91 = -t198 * t201 + t132 + t170;
t90 = (rSges(4,3) - t201) * t200 + t222 * t198;
t83 = t198 * t208 + t280;
t82 = -qJD(1) * t138 + t200 * t208;
t78 = t120 + t152 - t292;
t77 = t198 * t221 + t200 * t325;
t76 = t358 + ((-rSges(3,3) - pkin(5)) * t198 + t224 * t200) * qJD(1);
t75 = (t191 + (-pkin(1) - t334) * t198) * qJD(1) + t205;
t64 = -t156 * t278 - t308 + (-t197 * t278 - t198 * t275) * pkin(2);
t63 = t156 * t279 + t174 + (-t137 - t268) * t200;
t48 = -t151 * t278 - t118 * t198 + (-t184 * t278 - t185 * t294) * pkin(3);
t47 = t200 * t252 + t290;
t44 = t139 * t198 - t200 * t227;
t43 = t138 * t198 - t353;
t42 = -t139 * t200 - t357;
t41 = -t138 * t200 - t198 * t228;
t38 = t156 * t294 + (t200 * t222 - t188) * qJD(1) + t286;
t37 = (-t179 - t333) * t279 + (-t219 - t359) * t200 + t287;
t36 = qJD(1) * t95 + t198 * t207;
t35 = t200 * t207 + t174 + t290;
t26 = (-t150 + t220) * t198 + (-t198 * t325 + t200 * t221) * qJD(1);
t25 = t144 - t200 * t220 + (-t200 * t193 + (-t161 - t332) * t198) * qJD(1) + t288;
t24 = t68 + t291;
t22 = -t132 * t279 + t271;
t19 = -t120 * t279 + t273;
t16 = t23 + t291;
t7 = (-t124 - t132) * t279 + t271 + t272;
t6 = t279 * t324 + t244;
t5 = (-t124 + t324) * t279 + t244 + t272;
t1 = [(t103 * t76 + t104 * t75) * t346 + t155 * t295 + t184 * t135 - t154 * t296 + t185 * t134 + (t37 * t91 + t38 * t90) * t345 + t149 * t298 + t177 * t117 - t148 * t299 + t178 * t116 + (t25 * t78 + t26 * t77) * t344 + (t243 - t238) * t276 + (t242 + t239) * t275; (-qJD(2) * t227 + t197 * (-t142 * qJD(1) - t200 * t210) + t199 * (-t140 * qJD(1) - t200 * t209)) * t343 + (-qJD(2) * t228 + t197 * (qJD(1) * t143 - t198 * t210) + t199 * (qJD(1) * t141 - t198 * t209)) * t342 + ((-t104 * t341 + t306 / 0.2e1 + t304 / 0.2e1) * t200 + (t307 / 0.2e1 + t305 / 0.2e1 + t103 * t341) * t198) * qJD(1) + (t194 / 0.2e1 + t195 / 0.2e1) * t235 * qJD(2) + t204 + m(5) * (t25 * t94 + t26 * t95 + t35 * t77 + t36 * t78) + m(4) * (t121 * t37 + t122 * t38 + t63 * t90 + t64 * t91) + m(3) * ((-t198 * t75 - t200 * t76) * t169 + (-t103 * t200 - t104 * t198) * t160); (t16 * t5 + t35 * t95 + t36 * t94) * t344 - t337 + (t121 * t64 + t122 * t63 + t24 * t7) * t345 - t200 * t4 + (t198 * t44 - t200 * t43) * t278 + t198 * ((t198 * t82 + (t43 + t357) * qJD(1)) * t198 + (t44 * qJD(1) + (t140 * t275 + t142 * t276) * t200 + (-t83 + (-t304 - t306) * qJD(2) + (t139 - t228) * qJD(1)) * t198) * t200) + ((t145 * t198 + t146 * t200) * ((qJD(1) * t145 + t205) * t200 + (-t358 + (-t146 - t270 + t189) * qJD(1)) * t198) + t283 * t169 * t160) * t346 + (t198 * t42 - t200 * t41) * t279 - t200 * ((t200 * t83 + (t42 + t353) * qJD(1)) * t200 + (t41 * qJD(1) + (-t141 * t275 - t143 * t276 + t280) * t198 + (-t82 + (t305 + t307) * qJD(2) - t227 * qJD(1)) * t200) * t198) + t250; m(5) * (t100 * t25 + t101 * t26 + t47 * t77 + t48 * t78) + (-t198 * t37 - t200 * t38 + (t198 * t90 - t200 * t91) * qJD(1)) * t340 + m(4) * (-t198 * t91 - t200 * t90) * t137 + t204; m(5) * (t100 * t36 + t101 * t35 + t6 * t16 + t23 * t5 + t47 * t95 + t48 * t94) + m(4) * (-t122 * t137 * t200 - t121 * t308 + t22 * t24 + t68 * t7) + (-t198 * t64 - t200 * t63 + (-t121 * t200 + t122 * t198) * qJD(1)) * t340 + t206; (t137 * t156 * t283 + t22 * t68) * t345 + (t100 * t48 + t101 * t47 + t23 * t6) * t344 + t206; m(5) * ((-t198 * t78 - t200 * t77) * t118 + (-t198 * t25 - t200 * t26 + (t198 * t77 - t200 * t78) * qJD(1)) * t151) + t217; m(5) * (t19 * t16 + t5 * t65 + (-t198 * t94 - t200 * t95) * t118 + (-t198 * t36 - t200 * t35 + (t198 * t95 - t200 * t94) * qJD(1)) * t151) + t251; m(5) * (t19 * t23 + t6 * t65 + (-t100 * t198 - t101 * t200) * t118 + (-t198 * t48 - t200 * t47 + (-t100 * t200 + t101 * t198) * qJD(1)) * t151) + t251; (t118 * t151 * t283 + t19 * t65) * t344 + t251;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
