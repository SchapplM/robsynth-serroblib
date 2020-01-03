% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:29
% DurationCPUTime: 5.71s
% Computational Cost: add. (12170->459), mult. (10102->588), div. (0->0), fcn. (9456->10), ass. (0->232)
t223 = sin(pkin(9));
t224 = cos(pkin(9));
t225 = sin(qJ(5));
t227 = cos(qJ(5));
t222 = qJ(1) + qJ(2);
t213 = pkin(8) + t222;
t203 = sin(t213);
t204 = cos(t213);
t300 = t224 * t225;
t130 = -t203 * t227 + t204 * t300;
t299 = t224 * t227;
t131 = t203 * t225 + t204 * t299;
t304 = t204 * t223;
t72 = Icges(6,5) * t131 - Icges(6,6) * t130 + Icges(6,3) * t304;
t119 = Icges(6,4) * t131;
t76 = Icges(6,2) * t130 - Icges(6,6) * t304 - t119;
t118 = Icges(6,4) * t130;
t78 = Icges(6,1) * t131 + Icges(6,5) * t304 - t118;
t31 = -(t225 * t76 + t227 * t78) * t223 + t224 * t72;
t221 = qJD(1) + qJD(2);
t305 = t204 * t221;
t169 = qJ(4) * t305;
t193 = qJD(4) * t203;
t292 = t169 + t193;
t308 = t203 * t221;
t108 = pkin(3) * t308 - t292;
t288 = qJD(5) * t221;
t115 = (-qJDD(5) * t204 + t203 * t288) * t223;
t142 = -rSges(6,3) * t224 + (rSges(6,1) * t227 - rSges(6,2) * t225) * t223;
t220 = qJDD(1) + qJDD(2);
t191 = -qJDD(5) * t224 + t220;
t192 = -qJD(5) * t224 + t221;
t214 = sin(t222);
t215 = cos(t222);
t226 = sin(qJ(1));
t217 = t226 * pkin(1);
t228 = cos(qJ(1));
t229 = qJD(1) ^ 2;
t309 = pkin(1) * qJDD(1);
t269 = -t217 * t229 + t228 * t309;
t326 = pkin(2) * t220;
t219 = t221 ^ 2;
t327 = pkin(2) * t219;
t239 = -t214 * t327 + t215 * t326 + t269;
t237 = t221 * t193 + t239;
t157 = (-rSges(6,1) * t225 - rSges(6,2) * t227) * t223;
t146 = qJD(5) * t157;
t287 = qJD(5) * t223;
t254 = -t146 * t287 - qJDD(4);
t303 = t204 * t224;
t133 = pkin(4) * t303 + pkin(7) * t304;
t151 = t204 * pkin(3) + t203 * qJ(4);
t295 = t151 + t133;
t325 = pkin(4) * t224;
t128 = -t203 * t300 - t204 * t227;
t89 = qJD(5) * t131 + t128 * t221;
t129 = t203 * t299 - t204 * t225;
t90 = qJD(5) * t130 + t129 * t221;
t263 = rSges(6,1) * t90 + rSges(6,2) * t89;
t301 = t221 * t223;
t285 = t203 * t301;
t49 = rSges(6,3) * t285 + t263;
t81 = t131 * rSges(6,1) - t130 * rSges(6,2) + rSges(6,3) * t304;
t17 = t115 * t142 + t191 * t81 - t192 * t49 + t295 * t220 + t254 * t204 + (-t108 + (-pkin(7) * t223 - t325) * t308) * t221 + t237;
t361 = t17 - g(2);
t259 = t130 * t76 + t131 * t78;
t27 = t304 * t72 + t259;
t360 = t204 * t27;
t205 = pkin(2) * t214;
t260 = -t203 * rSges(4,1) - t204 * rSges(4,2);
t359 = t260 - t205;
t336 = t214 * rSges(3,1) + t215 * rSges(3,2);
t140 = t336 * t221;
t319 = pkin(1) * qJD(1);
t286 = t226 * t319;
t124 = t286 + t140;
t302 = t215 * t221;
t186 = pkin(2) * t302;
t281 = t142 * t287;
t358 = t192 * t81 + t186 + (-qJD(4) - t281) * t204;
t357 = -t128 * t76 + t129 * t78;
t135 = -Icges(6,3) * t224 + (Icges(6,5) * t227 - Icges(6,6) * t225) * t223;
t310 = Icges(6,4) * t227;
t136 = -Icges(6,6) * t224 + (-Icges(6,2) * t225 + t310) * t223;
t311 = Icges(6,4) * t225;
t137 = -Icges(6,5) * t224 + (Icges(6,1) * t227 - t311) * t223;
t238 = -t130 * t136 + t131 * t137 + t135 * t304;
t353 = t238 * t192;
t289 = qJD(4) * t204;
t262 = pkin(3) * t305 + qJ(4) * t308 - t289;
t351 = -t221 * t151 + t186 + t262;
t114 = (qJDD(5) * t203 + t204 * t288) * t223;
t306 = t203 * t224;
t307 = t203 * t223;
t132 = pkin(4) * t306 + pkin(7) * t307;
t198 = t203 * pkin(3);
t149 = -qJ(4) * t204 + t198;
t218 = t228 * pkin(1);
t290 = t229 * t218 + t226 * t309;
t274 = t214 * t326 + t215 * t327 + t290;
t251 = t220 * t149 + t221 * t262 + t274;
t283 = t221 * t303;
t284 = t204 * t301;
t339 = pkin(4) * t283 + pkin(7) * t284;
t91 = -qJD(5) * t129 - t130 * t221;
t92 = qJD(5) * t128 + t131 * t221;
t50 = t92 * rSges(6,1) + t91 * rSges(6,2) + rSges(6,3) * t284;
t80 = t129 * rSges(6,1) + t128 * rSges(6,2) + rSges(6,3) * t307;
t16 = -t114 * t142 + t220 * t132 + t191 * t80 + t192 * t50 + (-t289 + t339) * t221 + t254 * t203 + t251;
t350 = -g(3) + t16;
t175 = rSges(5,1) * t306;
t266 = -rSges(5,2) * t307 + t175;
t110 = -rSges(5,3) * t204 + t266;
t293 = rSges(5,1) * t283 + rSges(5,3) * t308;
t32 = -qJDD(4) * t203 + t220 * t110 + ((-rSges(5,2) * t301 - qJD(4)) * t204 + t293) * t221 + t251;
t349 = t32 - g(3);
t294 = rSges(5,2) * t285 + rSges(5,3) * t305;
t177 = rSges(5,2) * t304;
t111 = rSges(5,1) * t303 + rSges(5,3) * t203 - t177;
t296 = t151 + t111;
t33 = -qJDD(4) * t204 + (-t175 * t221 - t108 + t294) * t221 + t296 * t220 + t237;
t348 = t33 - g(2);
t172 = rSges(4,1) * t305;
t347 = -t220 * t260 + t221 * (-rSges(4,2) * t308 + t172) + t274 - g(3);
t195 = t203 * rSges(4,2);
t152 = rSges(4,1) * t204 - t195;
t346 = t220 * t152 + t219 * t260 - g(2) + t239;
t200 = t214 * rSges(3,2);
t141 = rSges(3,1) * t302 - t200 * t221;
t345 = t141 * t221 + t220 * t336 - g(3) + t290;
t164 = rSges(3,1) * t215 - t200;
t344 = -t140 * t221 + t164 * t220 - g(2) + t269;
t212 = t228 * t319;
t36 = t221 * t295 + t212 + t358;
t316 = t221 * t36;
t343 = (-t205 + (-t325 - pkin(3) + (-rSges(6,3) - pkin(7)) * t223) * t203) * t316;
t277 = t149 + t205;
t342 = t221 * (t110 + t277);
t341 = t221 * t359;
t267 = t132 + t277;
t139 = t221 * t152;
t338 = t172 - t139;
t337 = -t198 - t205;
t334 = -t221 * t133 + t339 + t351 - t358 + t50;
t275 = t186 - t289;
t333 = -t221 * t111 - t275 + t293 + t351;
t206 = pkin(2) * t215;
t282 = t206 + t151;
t332 = t282 + t133 + t81;
t312 = Icges(6,4) * t129;
t74 = Icges(6,2) * t128 + Icges(6,6) * t307 + t312;
t331 = t203 * (-Icges(6,1) * t128 + t312 + t74) + t204 * (Icges(6,1) * t130 + t119 - t76);
t117 = Icges(6,4) * t128;
t77 = Icges(6,1) * t129 + Icges(6,5) * t307 + t117;
t244 = t203 * (-Icges(6,2) * t129 + t117 + t77) - t204 * (Icges(6,2) * t131 + t118 - t78);
t330 = -m(5) - m(6);
t329 = t114 / 0.2e1;
t328 = t115 / 0.2e1;
t154 = (-Icges(6,5) * t225 - Icges(6,6) * t227) * t223;
t143 = qJD(5) * t154;
t155 = (-Icges(6,2) * t227 - t311) * t223;
t144 = qJD(5) * t155;
t156 = (-Icges(6,1) * t225 - t310) * t223;
t145 = qJD(5) * t156;
t41 = -t143 * t224 + (-t144 * t225 + t145 * t227 + (-t136 * t227 - t137 * t225) * qJD(5)) * t223;
t67 = -t135 * t224 + (-t136 * t225 + t137 * t227) * t223;
t324 = t67 * t191 + t41 * t192;
t323 = -t130 * t74 + t131 * t77;
t318 = t203 * t72;
t71 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t307;
t317 = t204 * t71;
t30 = -t224 * t71 + (-t225 * t74 + t227 * t77) * t223;
t315 = t30 * t114;
t314 = t31 * t115;
t298 = t136 - t156;
t297 = t137 + t155;
t24 = t128 * t74 + t129 * t77 + t71 * t307;
t25 = -t307 * t72 - t357;
t279 = -t287 / 0.2e1;
t278 = t287 / 0.2e1;
t125 = t164 * t221 + t212;
t273 = t203 * t279;
t272 = t203 * t278;
t271 = t204 * t279;
t270 = t204 * t278;
t265 = -t193 + t286;
t127 = t152 + t206;
t264 = t192 * t80 - t203 * t281 - t193;
t188 = rSges(2,1) * t228 - t226 * rSges(2,2);
t187 = rSges(2,1) * t226 + rSges(2,2) * t228;
t35 = t221 * t267 + t264 + t286;
t258 = -t203 * t35 - t204 * t36;
t257 = t203 * t49 + t204 * t50;
t256 = -t203 * t81 + t204 * t80;
t255 = t203 * (Icges(6,5) * t128 - Icges(6,6) * t129) - t204 * (Icges(6,5) * t130 + Icges(6,6) * t131);
t247 = (t203 * t24 - t204 * t25) * t223;
t26 = -t304 * t71 - t323;
t246 = (t203 * t26 - t360) * t223;
t245 = -t263 + t292;
t54 = t80 + t267;
t94 = t111 + t282;
t93 = (-rSges(5,3) - qJ(4)) * t204 + t266 - t337;
t103 = t286 - t341;
t104 = t212 + t186 + t139;
t232 = (-t103 * t195 + t104 * t359) * t221;
t10 = qJD(5) * t246 - t353;
t44 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t284;
t46 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t284;
t48 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t284;
t13 = -t224 * t44 + (-t225 * t46 + t227 * t48 + (-t225 * t77 - t227 * t74) * qJD(5)) * t223;
t43 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t285;
t45 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t285;
t47 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t285;
t14 = -t224 * t43 + (-t225 * t45 + t227 * t47 + (t225 * t78 - t227 * t76) * qJD(5)) * t223;
t21 = t130 * t144 - t131 * t145 + t136 * t89 + t137 * t90 + (t135 * t308 - t143 * t204) * t223;
t22 = t128 * t144 + t129 * t145 + t136 * t91 + t137 * t92 + (t135 * t305 + t143 * t203) * t223;
t51 = t128 * t136 + t129 * t137 + t135 * t307;
t42 = t51 * t192;
t9 = qJD(5) * t247 + t42;
t231 = (t42 + ((t24 - t259 + t27) * t203 + (t26 + (t317 - t318) * t223 - t25 + t323) * t204) * t287) * t270 + t315 / 0.2e1 + t314 / 0.2e1 + t51 * t329 - t238 * t328 + t324 + (t13 + t22) * t272 + (t353 + (t360 + (t323 + t25 + (t317 + t318) * t223 + t357) * t203) * t287 + t10) * t273 + (t14 + t21 + t9) * t271 + (Icges(5,2) * t224 ^ 2 + (Icges(5,1) * t223 + 0.2e1 * Icges(5,4) * t224) * t223 + Icges(4,3) + Icges(3,3)) * t220;
t65 = t265 + t342;
t66 = t221 * t296 + t212 + t275;
t230 = (t66 * (-t175 + t337) - t65 * t177) * t221;
t102 = rSges(6,1) * t130 + rSges(6,2) * t131;
t101 = rSges(6,1) * t128 - rSges(6,2) * t129;
t39 = t256 * t287 + qJD(3);
t15 = -t114 * t81 - t115 * t80 + t257 * t287 + qJDD(3);
t6 = t128 * t45 + t129 * t47 + t76 * t91 - t78 * t92 + (t203 * t43 - t305 * t72) * t223;
t5 = t128 * t46 + t129 * t48 + t74 * t91 + t77 * t92 + (t203 * t44 + t305 * t71) * t223;
t4 = t130 * t45 - t131 * t47 + t76 * t89 - t78 * t90 + (-t204 * t43 - t308 * t72) * t223;
t3 = t130 * t46 - t131 * t48 + t74 * t89 + t77 * t90 + (-t204 * t44 + t308 * t71) * t223;
t1 = [Icges(2,3) * qJDD(1) + t231 + (t344 * (t164 + t218) + t345 * (t217 + t336) + (-t125 + t141 + t212) * t124) * m(3) + ((t187 ^ 2 + t188 ^ 2) * qJDD(1) - g(2) * t188 - g(3) * t187) * m(2) + (t36 * (t245 - t286) + t350 * (t217 + t54) + (t36 + t334) * t35 + t343 + t361 * (t218 + t332)) * m(6) + (t66 * (t169 - t265 + t294) + t230 + t348 * (t218 + t94) + t349 * (t217 + t93) + (t66 + t333) * t65) * m(5) + (-t104 * t286 + t232 + t346 * (t127 + t218) + t347 * (t217 - t359) + (t104 + t338) * t103) * m(4); t231 + (t267 * t316 + t350 * t54 + (t245 + t264) * t36 + t334 * t35 + t343 + t361 * t332) * m(6) + (t230 + t348 * t94 + t349 * t93 + t333 * t65 + (-t193 + t292 + t294 + t342) * t66) * m(5) + (t338 * t103 - t104 * t341 + t346 * t127 - t347 * t359 + t232) * m(4) + (t124 * t141 - t125 * t140 + (-t124 * t221 + t344) * t164 + (t125 * t221 + t345) * t336) * m(3); m(6) * t15 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t330) * g(1); t330 * (-g(2) * t204 - g(3) * t203) + m(5) * (-t203 * t32 - t204 * t33) + m(6) * (-t16 * t203 - t17 * t204); -t224 * (t315 + t314 + (t13 * t203 - t14 * t204) * t287 + t324) / 0.2e1 + t191 * (-t224 * t67 + (t203 * t30 - t204 * t31) * t223) / 0.2e1 + t192 * (-t224 * t41 + ((t221 * t30 - t14) * t204 + (t221 * t31 + t13) * t203) * t223) / 0.2e1 + (t114 * t24 + t115 * t25 + t191 * t51 + t192 * t22 + (t203 * t5 - t204 * t6) * t287) * t307 / 0.2e1 + (-t224 * t51 + t247) * t329 + (-t22 * t224 + ((t221 * t24 - t6) * t204 + (t221 * t25 + t5) * t203) * t223) * t272 - (t114 * t26 + t115 * t27 - t191 * t238 + t192 * t21 + (t203 * t3 - t204 * t4) * t287) * t304 / 0.2e1 + (t224 * t238 + t246) * t328 + (-t21 * t224 + ((t221 * t26 - t4) * t204 + (t221 * t27 + t3) * t203) * t223) * t271 - t192 * (-t224 * t154 * t192 + ((-t225 * t297 - t227 * t298) * t192 + ((-t244 * t225 - t227 * t331) * t223 - t255 * t224) * qJD(5)) * t223) / 0.2e1 + ((t128 * t297 - t129 * t298 + t154 * t307) * t192 + (t128 * t244 - t129 * t331 + t255 * t307) * t287) * t273 + ((t130 * t297 + t131 * t298 - t154 * t304) * t192 + (t244 * t130 + t331 * t131 - t255 * t304) * t287) * t270 + (t203 * t10 + t204 * t9) * t301 / 0.2e1 + ((-t16 * t80 - t17 * t81 - t35 * t50 + t36 * t49) * t224 + (t15 * t256 + t39 * (-t305 * t81 - t308 * t80 + t257) + t258 * t146 + ((-t221 * t35 - t17) * t204 + (-t16 + t316) * t203) * t142) * t223 - (t101 * t35 - t102 * t36) * t192 - (t39 * (t101 * t204 + t102 * t203) + t258 * t157) * t287 - g(1) * t157 - g(2) * t101 - g(3) * t102) * m(6);];
tau = t1;
