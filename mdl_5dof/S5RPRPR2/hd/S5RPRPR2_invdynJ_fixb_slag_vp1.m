% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:48
% EndTime: 2022-01-23 09:18:54
% DurationCPUTime: 4.32s
% Computational Cost: add. (9387->383), mult. (5723->480), div. (0->0), fcn. (4240->10), ass. (0->230)
t192 = qJ(1) + pkin(8);
t187 = qJ(3) + t192;
t178 = sin(t187);
t179 = cos(t187);
t191 = qJD(1) + qJD(3);
t259 = qJD(5) * t191;
t110 = -qJDD(5) * t179 + t178 * t259;
t190 = pkin(9) + qJ(5);
t183 = sin(t190);
t185 = cos(t190);
t137 = rSges(6,1) * t185 - rSges(6,2) * t183;
t117 = t137 * qJD(5);
t294 = rSges(6,2) * t185;
t135 = rSges(6,1) * t183 + t294;
t195 = -pkin(7) - qJ(4);
t276 = t178 * t195;
t139 = t191 * t276;
t194 = cos(pkin(9));
t180 = pkin(4) * t194 + pkin(3);
t144 = t179 * t180;
t189 = qJDD(1) + qJDD(3);
t184 = sin(t192);
t186 = cos(t192);
t198 = qJD(1) ^ 2;
t196 = sin(qJ(1));
t197 = cos(qJ(1));
t217 = (-qJDD(1) * t196 - t197 * t198) * pkin(1);
t203 = (-qJDD(1) * t184 - t186 * t198) * pkin(2) + t217;
t263 = qJD(4) * t191;
t201 = qJDD(4) * t178 + t179 * t263 + t203;
t165 = t179 * qJ(4);
t304 = pkin(3) * t178;
t118 = -t165 + t304;
t318 = -t178 * t180 - t179 * t195;
t74 = t118 + t318;
t279 = t178 * t183;
t145 = rSges(6,2) * t279;
t278 = t178 * t185;
t254 = rSges(6,1) * t278;
t87 = -t179 * rSges(6,3) - t145 + t254;
t257 = -t118 + t74 - t87;
t261 = qJD(5) * t179;
t316 = t179 * pkin(3) + t178 * qJ(4);
t274 = t179 * t185;
t224 = rSges(6,1) * t274 + t178 * rSges(6,3);
t251 = qJD(5) * t294;
t275 = t179 * t183;
t252 = rSges(6,2) * t275;
t260 = qJD(5) * t183;
t250 = -t191 * t252 + (-rSges(6,1) * t260 - t251) * t178;
t58 = t224 * t191 + t250;
t162 = qJD(4) * t179;
t90 = t191 * t316 - t162;
t7 = -t117 * t261 + t110 * t135 + t257 * t189 + (-t90 + t139 - t58 + (t316 - t144) * t191) * t191 + t201;
t325 = t7 - g(1);
t109 = qJDD(5) * t178 + t179 * t259;
t273 = t179 * t191;
t151 = qJ(4) * t273;
t177 = pkin(2) * t186;
t188 = t197 * pkin(1);
t181 = qJDD(1) * t188;
t305 = pkin(1) * t196;
t239 = -pkin(2) * t184 - t305;
t213 = qJDD(1) * t177 + t239 * t198 + t181;
t161 = qJD(4) * t178;
t267 = t151 + t161;
t277 = t178 * t191;
t204 = -qJDD(4) * t179 + t189 * t316 + t178 * t263 + t213 + t191 * (-pkin(3) * t277 + t267);
t262 = qJD(5) * t178;
t88 = -t252 + t224;
t64 = t144 - t276 + t88;
t301 = -t316 + t64;
t221 = rSges(6,3) * t273 + t191 * t145 - t179 * t251;
t248 = t179 * t260;
t57 = (-t185 * t277 - t248) * rSges(6,1) + t221;
t8 = -t117 * t262 - t109 * t135 + t301 * t189 + (-t151 + t57 + (t318 + t304) * t191) * t191 + t204;
t324 = t8 - g(2);
t173 = t179 * rSges(4,1);
t121 = -rSges(4,2) * t178 + t173;
t119 = rSges(4,1) * t178 + rSges(4,2) * t179;
t272 = t191 * t119;
t323 = t121 * t189 - t191 * t272 - g(2) + t213;
t193 = sin(pkin(9));
t295 = rSges(5,2) * t193;
t253 = t179 * t295;
t132 = t191 * t253;
t297 = rSges(5,1) * t194;
t255 = t178 * t297;
t155 = t178 * t295;
t266 = t179 * rSges(5,3) + t155;
t91 = t255 - t266;
t286 = -t118 - t91;
t317 = -t178 * rSges(5,3) - t179 * t297;
t23 = t286 * t189 + (t191 * t317 + t132 - t90) * t191 + t201;
t322 = t23 - g(1);
t269 = rSges(5,3) * t273 + t191 * t155;
t92 = -t253 - t317;
t24 = t189 * t92 + t191 * (-t191 * t255 + t269) + t204;
t321 = t24 - g(2);
t152 = rSges(4,2) * t277;
t104 = rSges(4,1) * t273 - t152;
t320 = t104 * t191 + t119 * t189 + g(1) - t203;
t76 = t191 * t87;
t319 = t191 * t74 - t76;
t138 = t186 * rSges(3,1) - rSges(3,2) * t184;
t113 = t138 + t188;
t175 = Icges(6,4) * t185;
t228 = -Icges(6,2) * t183 + t175;
t129 = Icges(6,1) * t183 + t175;
t315 = t177 + t188;
t126 = Icges(6,5) * t185 - Icges(6,6) * t183;
t125 = Icges(6,5) * t183 + Icges(6,6) * t185;
t209 = Icges(6,3) * t191 - t125 * qJD(5);
t219 = t228 * t179;
t83 = Icges(6,6) * t178 + t219;
t290 = t183 * t83;
t284 = Icges(6,4) * t183;
t130 = Icges(6,1) * t185 - t284;
t220 = t130 * t179;
t85 = Icges(6,5) * t178 + t220;
t230 = -t185 * t85 + t290;
t313 = -t126 * t277 + t209 * t179 + t230 * t191;
t218 = t126 * t179;
t82 = Icges(6,4) * t278 - Icges(6,2) * t279 - Icges(6,6) * t179;
t291 = t183 * t82;
t142 = Icges(6,4) * t279;
t84 = Icges(6,1) * t278 - Icges(6,5) * t179 - t142;
t231 = -t185 * t84 + t291;
t312 = t209 * t178 + (t218 + t231) * t191;
t127 = Icges(6,2) * t185 + t284;
t226 = t183 * t127 - t185 * t129;
t311 = t126 * qJD(5) + t226 * t191;
t80 = Icges(6,5) * t278 - Icges(6,6) * t279 - Icges(6,3) * t179;
t30 = -t231 * t178 - t179 * t80;
t299 = -Icges(6,2) * t278 - t142 + t84;
t300 = t129 * t178 + t82;
t310 = -t299 * t183 - t300 * t185;
t309 = t109 / 0.2e1;
t308 = t110 / 0.2e1;
t307 = t178 / 0.2e1;
t306 = -t179 / 0.2e1;
t303 = -t178 * t80 - t84 * t274;
t81 = Icges(6,3) * t178 + t218;
t302 = t178 * t81 + t85 * t274;
t298 = -t127 * t179 + t85;
t222 = t315 * qJD(1);
t79 = t121 * t191 + t222;
t293 = t119 * t79;
t223 = t239 * qJD(1);
t216 = t161 + t223;
t249 = t135 * t261;
t208 = t216 - t249;
t28 = t257 * t191 + t208;
t289 = t191 * t28;
t281 = t125 * t179;
t49 = -t226 * t178 - t281;
t288 = t49 * t191;
t287 = -t129 * t179 - t83;
t69 = t316 + t92;
t282 = t125 * t178;
t280 = t126 * t191;
t271 = -t127 + t130;
t270 = t129 + t228;
t111 = t191 * t118;
t265 = t161 - t111;
t258 = m(3) + m(4) + m(5);
t256 = t316 + t301;
t247 = -pkin(3) - t297;
t246 = -t262 / 0.2e1;
t245 = t262 / 0.2e1;
t244 = -t261 / 0.2e1;
t243 = t261 / 0.2e1;
t65 = t85 * t278;
t242 = t179 * t81 - t65;
t241 = -t80 + t290;
t159 = rSges(2,1) * t197 - rSges(2,2) * t196;
t158 = rSges(2,1) * t196 + rSges(2,2) * t197;
t136 = rSges(3,1) * t184 + rSges(3,2) * t186;
t107 = t135 * t262;
t214 = t222 - t162;
t29 = t256 * t191 - t107 + t214;
t235 = -t178 * t29 - t179 * t28;
t31 = -t83 * t279 - t242;
t234 = t31 * t178 - t30 * t179;
t32 = -t82 * t275 - t303;
t33 = -t83 * t275 + t302;
t233 = t33 * t178 - t32 * t179;
t232 = t178 * t87 + t179 * t88;
t42 = t183 * t84 + t185 * t82;
t43 = t183 * t85 + t185 * t83;
t227 = t127 * t185 + t129 * t183;
t63 = -t87 + t318;
t215 = -t298 * t183 + t287 * t185;
t68 = t247 * t178 + t165 + t266;
t78 = t223 - t272;
t212 = (-t270 * t183 + t271 * t185) * t191;
t211 = Icges(6,5) * t191 - qJD(5) * t129;
t210 = Icges(6,6) * t191 - t127 * qJD(5);
t11 = t234 * qJD(5) + t288;
t115 = t228 * qJD(5);
t116 = t130 * qJD(5);
t50 = -t226 * t179 + t282;
t44 = t50 * t191;
t12 = t233 * qJD(5) + t44;
t54 = t210 * t178 + t191 * t219;
t56 = t211 * t178 + t191 * t220;
t16 = -t231 * qJD(5) + t183 * t56 + t185 * t54;
t53 = t210 * t179 - t228 * t277;
t55 = -t130 * t277 + t211 * t179;
t17 = -t230 * qJD(5) + t183 * t55 + t185 * t53;
t202 = -t227 * qJD(5) - t115 * t183 + t116 * t185 + t125 * t191;
t20 = t178 * t311 + t202 * t179;
t21 = t202 * t178 - t179 * t311;
t207 = (t44 + ((t31 - t65 + (t81 + t291) * t179 + t303) * t179 + t302 * t178) * qJD(5)) * t243 + (-t226 * qJD(5) + t115 * t185 + t116 * t183) * t191 + (t43 + t50) * t309 + (t42 + t49) * t308 + (-t288 + ((t241 * t179 - t302 + t33) * t179 + (t241 * t178 + t242 + t32) * t178) * qJD(5) + t11) * t246 + (t17 + t20) * t245 + (t12 + t16 + t21) * t244 + (Icges(4,3) + t227 + Icges(5,2) * t194 ^ 2 + (Icges(5,1) * t193 + 0.2e1 * Icges(5,4) * t194) * t193) * t189;
t206 = -t43 * qJD(5) - t183 * t53 + t185 * t55 + t191 * t81;
t205 = -t42 * qJD(5) - t183 * t54 + t185 * t56 + t191 * t80;
t45 = t286 * t191 + t216;
t46 = t69 * t191 + t214;
t200 = t45 * (t132 + t162) + t46 * (t267 + t269) + (t45 * t247 * t179 + (t45 * (-rSges(5,3) - qJ(4)) + t46 * t247) * t178) * t191;
t199 = t28 * (t139 + t162 - t250) + t29 * (-rSges(6,1) * t248 + t161 + t221) + (t28 * (-t224 - t144) + t29 * (t318 - t254)) * t191;
t102 = t135 * t179;
t101 = t135 * t178;
t86 = t191 * t91;
t41 = t232 * qJD(5) + qJD(2);
t13 = t109 * t87 - t110 * t88 + qJDD(2) + (t178 * t58 + t179 * t57) * qJD(5);
t6 = t206 * t178 - t179 * t313;
t5 = t205 * t178 - t179 * t312;
t4 = t178 * t313 + t206 * t179;
t3 = t178 * t312 + t205 * t179;
t1 = [t207 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t78 * t152 + (-t78 * t173 - t293) * t191 + (t79 * t239 - t315 * t78) * qJD(1) + t323 * (t121 + t315) - t320 * (-t119 + t239)) * m(4) + ((qJDD(1) * t138 - g(2) + t181) * t113 + (-g(1) - qJDD(1) * t136 + t217 + (-0.2e1 * t138 + 0.2e1 * t113 - t188) * t198) * (-t136 - t305)) * m(3) + (g(1) * t158 - g(2) * t159 + (t158 ^ 2 + t159 ^ 2) * qJDD(1)) * m(2) + ((t29 * t239 - t28 * t315) * qJD(1) + t199 - (-t111 - t28 + t208 + t319) * t29 + t324 * (t64 + t315) + t325 * (t239 + t63)) * m(6) + ((t46 * t239 - t315 * t45) * qJD(1) + t200 - (-t111 - t45 - t86 + t216) * t46 + t321 * (t69 + t315) + t322 * (t239 + t68)) * m(5); m(6) * t13 + t258 * qJDD(2) + (-m(6) - t258) * g(3); t207 + (t199 - t28 * (t107 + t162) - t29 * (-t249 + t265 + t319) + t256 * t289 + t324 * t64 + t325 * t63) * m(6) + (t200 - t45 * t162 - t46 * (-t86 + t265) + (t45 * t191 + t321) * t69 + t322 * t68) * m(5) + (-t272 * t79 - t104 * t78 + t293 * t191 + (t78 * t191 + t323) * t121 + t320 * t119) * m(4); (-m(5) - m(6)) * (g(1) * t178 - g(2) * t179) + 0.2e1 * (t8 * t306 + t7 * t307) * m(6) + 0.2e1 * (t23 * t307 + t24 * t306) * m(5); t12 * t273 / 0.2e1 + (t109 * t33 + t110 * t32 + t50 * t189 + t20 * t191 + (t178 * t4 - t179 * t3) * qJD(5)) * t307 + t233 * t309 + ((t191 * t33 - t3) * t179 + (t191 * t32 + t4) * t178) * t245 + t11 * t277 / 0.2e1 + (t31 * t109 + t30 * t110 + t49 * t189 + t21 * t191 + (t178 * t6 - t179 * t5) * qJD(5)) * t306 + t234 * t308 + ((t191 * t31 - t5) * t179 + (t191 * t30 + t6) * t178) * t244 + t189 * (t43 * t178 - t42 * t179) / 0.2e1 + t191 * ((t191 * t43 - t16) * t179 + (t191 * t42 + t17) * t178) / 0.2e1 + ((-t262 * t281 + t280) * t178 + (t212 + (-t310 * t179 + (t282 + t215) * t178) * qJD(5)) * t179) * t246 + ((-t261 * t282 - t280) * t179 + (t212 + (t215 * t178 + (-t310 + t281) * t179) * qJD(5)) * t178) * t243 - t191 * ((t271 * t183 + t270 * t185) * t191 + ((t298 * t178 - t299 * t179) * t185 + (t287 * t178 + t300 * t179) * t183) * qJD(5)) / 0.2e1 + (t13 * t232 + t41 * ((t57 + t76) * t179 + (-t191 * t88 + t58) * t178) + t235 * t117 + ((-t191 * t29 - t7) * t179 + (-t8 + t289) * t178) * t135 - (t101 * t28 - t102 * t29) * t191 - (t41 * (-t101 * t178 - t102 * t179) + t235 * t137) * qJD(5) + g(1) * t102 + g(2) * t101 - g(3) * t137) * m(6);];
tau = t1;
