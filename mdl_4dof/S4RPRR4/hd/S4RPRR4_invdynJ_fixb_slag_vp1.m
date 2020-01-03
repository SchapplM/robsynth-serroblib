% Calculate vector of inverse dynamics joint torques for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:34
% DurationCPUTime: 17.24s
% Computational Cost: add. (13435->710), mult. (17768->978), div. (0->0), fcn. (16764->8), ass. (0->340)
t258 = sin(qJ(1));
t439 = pkin(1) * t258;
t255 = qJ(1) + pkin(7);
t248 = sin(t255);
t475 = t248 * pkin(5);
t261 = cos(qJ(1));
t254 = t261 * pkin(1);
t257 = sin(qJ(3));
t260 = cos(qJ(3));
t256 = sin(qJ(4));
t431 = rSges(5,2) * t256;
t259 = cos(qJ(4));
t433 = rSges(5,1) * t259;
t321 = -t431 + t433;
t168 = -rSges(5,3) * t260 + t257 * t321;
t360 = qJD(4) * t257;
t249 = cos(t255);
t363 = qJD(3) * t249;
t178 = -t248 * t360 + t363;
t438 = t257 * pkin(3);
t224 = -pkin(6) * t260 + t438;
t359 = qJD(4) * t260;
t237 = qJD(1) - t359;
t390 = t256 * t260;
t392 = t249 * t259;
t150 = t248 * t390 + t392;
t389 = t259 * t260;
t394 = t249 * t256;
t151 = t248 * t389 - t394;
t322 = rSges(5,1) * t151 - rSges(5,2) * t150;
t398 = t248 * t257;
t90 = -rSges(5,3) * t398 - t322;
t474 = t168 * t178 + t224 * t363 - t237 * t90;
t364 = qJD(3) * t248;
t188 = rSges(3,1) * t248 + rSges(3,2) * t249;
t175 = -t188 - t439;
t144 = Icges(5,4) * t151;
t83 = -Icges(5,2) * t150 + Icges(5,6) * t398 + t144;
t143 = Icges(5,4) * t150;
t87 = -Icges(5,1) * t151 - Icges(5,5) * t398 + t143;
t469 = t256 * t83 + t259 * t87;
t80 = Icges(5,5) * t151 - Icges(5,6) * t150 + Icges(5,3) * t398;
t31 = -t257 * t469 - t260 * t80;
t262 = qJD(1) ^ 2;
t356 = t262 * t254;
t177 = t249 * t360 + t364;
t25 = -t150 * t83 - t151 * t87 + t398 * t80;
t397 = t248 * t259;
t152 = -t249 * t390 + t397;
t399 = t248 * t256;
t153 = t249 * t389 + t399;
t393 = t249 * t257;
t82 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t393;
t413 = Icges(5,4) * t153;
t85 = Icges(5,2) * t152 + Icges(5,6) * t393 + t413;
t145 = Icges(5,4) * t152;
t88 = Icges(5,1) * t153 + Icges(5,5) * t393 + t145;
t26 = -t150 * t85 + t151 * t88 + t398 * t82;
t302 = Icges(5,5) * t259 - Icges(5,6) * t256;
t154 = -Icges(5,3) * t260 + t257 * t302;
t411 = Icges(5,4) * t259;
t303 = -Icges(5,2) * t256 + t411;
t156 = -Icges(5,6) * t260 + t257 * t303;
t412 = Icges(5,4) * t256;
t305 = Icges(5,1) * t259 - t412;
t158 = -Icges(5,5) * t260 + t257 * t305;
t50 = -t150 * t156 + t151 * t158 + t154 * t398;
t10 = t177 * t26 - t178 * t25 + t237 * t50;
t27 = t152 * t83 - t153 * t87 + t393 * t80;
t28 = t152 * t85 + t153 * t88 + t393 * t82;
t51 = t152 * t156 + t153 * t158 + t154 * t393;
t11 = t177 * t28 - t178 * t27 + t237 * t51;
t241 = t248 * rSges(4,3);
t391 = t249 * t260;
t140 = rSges(4,1) * t391 - rSges(4,2) * t393 + t241;
t191 = pkin(2) * t249 + t475;
t462 = t254 + t191;
t113 = t140 + t462;
t465 = t10 * t248 + t11 * t249;
t189 = rSges(3,1) * t249 - rSges(3,2) * t248;
t176 = t189 + t254;
t250 = Icges(4,4) * t260;
t304 = -Icges(4,2) * t257 + t250;
t209 = Icges(4,1) * t257 + t250;
t253 = t260 * pkin(3);
t463 = pkin(6) * t257 + t253;
t206 = Icges(4,5) * t260 - Icges(4,6) * t257;
t205 = Icges(4,5) * t257 + Icges(4,6) * t260;
t283 = qJD(3) * t205;
t414 = Icges(4,4) * t257;
t210 = Icges(4,1) * t260 - t414;
t138 = Icges(4,5) * t248 + t210 * t249;
t136 = Icges(4,6) * t248 + t249 * t304;
t402 = t136 * t257;
t299 = -t138 * t260 + t402;
t405 = Icges(4,3) * t249;
t459 = -t249 * t283 + (-t206 * t248 + t299 + t405) * qJD(1);
t218 = Icges(4,4) * t398;
t396 = t248 * t260;
t410 = Icges(4,5) * t249;
t137 = Icges(4,1) * t396 - t218 - t410;
t407 = Icges(4,6) * t249;
t135 = Icges(4,4) * t396 - Icges(4,2) * t398 - t407;
t403 = t135 * t257;
t300 = -t137 * t260 + t403;
t134 = Icges(4,3) * t248 + t206 * t249;
t370 = qJD(1) * t134;
t458 = qJD(1) * t300 - t248 * t283 + t370;
t133 = Icges(4,5) * t396 - Icges(4,6) * t398 - t405;
t52 = -t133 * t249 - t248 * t300;
t207 = Icges(4,2) * t260 + t414;
t295 = t207 * t257 - t209 * t260;
t457 = t295 * qJD(1) + qJD(3) * t206;
t381 = -Icges(4,2) * t396 + t137 - t218;
t383 = t209 * t248 + t135;
t456 = -t257 * t381 - t260 * t383;
t155 = Icges(5,3) * t257 + t260 * t302;
t297 = -t156 * t256 + t158 * t259;
t307 = -t256 * t85 + t259 * t88;
t455 = t177 * (-t154 * t249 - t307) - t178 * (-t154 * t248 + t469) + t237 * (t155 - t297);
t185 = (-Icges(5,2) * t259 - t412) * t257;
t454 = t177 * (-Icges(5,2) * t153 + t145 + t88) - t178 * (-Icges(5,2) * t151 - t143 - t87) + t237 * (t158 + t185);
t358 = qJD(1) * qJD(3);
t180 = qJDD(3) * t248 + t249 * t358;
t361 = qJD(3) * t260;
t345 = t249 * t361;
t366 = qJD(1) * t257;
t350 = t248 * t366;
t280 = t345 - t350;
t357 = qJDD(4) * t257;
t102 = qJD(4) * t280 + t249 * t357 + t180;
t453 = t102 / 0.2e1;
t181 = -qJDD(3) * t249 + t248 * t358;
t349 = t249 * t366;
t282 = t248 * t361 + t349;
t103 = qJD(4) * t282 + t248 * t357 + t181;
t452 = t103 / 0.2e1;
t451 = -t177 / 0.2e1;
t450 = t177 / 0.2e1;
t449 = -t178 / 0.2e1;
t448 = t178 / 0.2e1;
t447 = t180 / 0.2e1;
t446 = t181 / 0.2e1;
t192 = qJD(3) * t360 - qJDD(4) * t260 + qJDD(1);
t445 = t192 / 0.2e1;
t444 = -t237 / 0.2e1;
t443 = t237 / 0.2e1;
t440 = -rSges(5,3) - pkin(6);
t362 = qJD(3) * t257;
t275 = t237 * t259 + t256 * t362;
t365 = qJD(1) * t260;
t329 = -qJD(4) + t365;
t77 = t248 * t275 - t329 * t394;
t274 = t237 * t256 - t259 * t362;
t78 = t248 * t274 + t329 * t392;
t42 = Icges(5,5) * t78 + Icges(5,6) * t77 + Icges(5,3) * t282;
t44 = Icges(5,4) * t78 + Icges(5,2) * t77 + Icges(5,6) * t282;
t46 = Icges(5,1) * t78 + Icges(5,4) * t77 + Icges(5,5) * t282;
t7 = (-qJD(3) * t469 - t42) * t260 + (qJD(3) * t80 - t256 * t44 + t259 * t46 + (t256 * t87 - t259 * t83) * qJD(4)) * t257;
t437 = t7 * t178;
t75 = t249 * t275 + t329 * t399;
t76 = t249 * t274 - t329 * t397;
t41 = Icges(5,5) * t76 + Icges(5,6) * t75 + Icges(5,3) * t280;
t43 = Icges(5,4) * t76 + Icges(5,2) * t75 + Icges(5,6) * t280;
t45 = Icges(5,1) * t76 + Icges(5,4) * t75 + Icges(5,5) * t280;
t8 = (qJD(3) * t307 - t41) * t260 + (qJD(3) * t82 - t256 * t43 + t259 * t45 + (-t256 * t88 - t259 * t85) * qJD(4)) * t257;
t436 = t8 * t177;
t184 = (-Icges(5,5) * t256 - Icges(5,6) * t259) * t257;
t116 = qJD(3) * t155 + qJD(4) * t184;
t157 = Icges(5,6) * t257 + t260 * t303;
t117 = qJD(3) * t157 + qJD(4) * t185;
t159 = Icges(5,5) * t257 + t260 * t305;
t186 = (-Icges(5,1) * t256 - t411) * t257;
t118 = qJD(3) * t159 + qJD(4) * t186;
t24 = (qJD(3) * t297 - t116) * t260 + (qJD(3) * t154 - t117 * t256 + t118 * t259 + (-t156 * t259 - t158 * t256) * qJD(4)) * t257;
t58 = -t154 * t260 + t257 * t297;
t435 = t192 * t58 + t237 * t24;
t434 = rSges(4,1) * t260;
t212 = rSges(4,1) * t257 + rSges(4,2) * t260;
t167 = t212 * t249;
t67 = qJD(1) * t113 - t212 * t364;
t428 = t167 * t67;
t371 = rSges(4,2) * t398 + rSges(4,3) * t249;
t139 = rSges(4,1) * t396 - t371;
t245 = t249 * pkin(5);
t190 = pkin(2) * t248 - t245;
t334 = -t190 - t439;
t327 = -t139 + t334;
t348 = t212 * t363;
t66 = qJD(1) * t327 - t348;
t424 = t248 * t66;
t171 = t463 * t248;
t326 = -t171 + t334;
t37 = qJD(1) * t326 - t474;
t423 = t249 * t37;
t422 = t249 * t66;
t251 = t257 * rSges(5,3);
t421 = t31 * t103;
t32 = t257 * t307 - t260 * t82;
t420 = t32 * t102;
t417 = t171 - t90;
t173 = pkin(3) * t391 + pkin(6) * t393;
t91 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t393;
t416 = t173 + t91;
t401 = t205 * t248;
t400 = t205 * t249;
t98 = -t248 * t295 - t400;
t387 = t98 * qJD(1);
t187 = (-rSges(5,1) * t256 - rSges(5,2) * t259) * t257;
t119 = qJD(4) * t187 + (t260 * t321 + t251) * qJD(3);
t199 = t463 * qJD(3);
t386 = -t119 - t199;
t385 = -t133 * t248 - t137 * t391;
t384 = t134 * t248 + t138 * t391;
t382 = -t209 * t249 - t136;
t380 = -t207 * t249 + t138;
t377 = t168 + t224;
t353 = t257 * t431;
t376 = rSges(5,3) * t396 + t248 * t353;
t375 = rSges(5,3) * t391 + t249 * t353;
t367 = qJD(1) * t249;
t374 = rSges(4,2) * t350 + rSges(4,3) * t367;
t373 = -t207 + t210;
t372 = t209 + t304;
t369 = qJD(1) * t206;
t368 = qJD(1) * t248;
t355 = rSges(5,1) * t76 + rSges(5,2) * t75 + rSges(5,3) * t345;
t354 = t257 * t433;
t346 = t249 * t362;
t344 = -pkin(2) - t434;
t341 = t367 / 0.2e1;
t340 = -t364 / 0.2e1;
t339 = t364 / 0.2e1;
t338 = -t363 / 0.2e1;
t337 = t363 / 0.2e1;
t332 = t245 - t439;
t120 = t138 * t396;
t331 = t134 * t249 - t120;
t330 = -t133 + t402;
t328 = qJDD(1) * t254 - t262 * t439;
t325 = -rSges(5,1) * t78 - rSges(5,2) * t77;
t215 = rSges(2,1) * t261 - rSges(2,2) * t258;
t213 = rSges(2,1) * t258 + rSges(2,2) * t261;
t214 = -rSges(4,2) * t257 + t434;
t72 = t136 * t260 + t138 * t257;
t284 = qJD(3) * t207;
t94 = -t249 * t284 + (-t248 * t304 + t407) * qJD(1);
t285 = qJD(3) * t209;
t96 = -t249 * t285 + (-t210 * t248 + t410) * qJD(1);
t265 = -qJD(3) * t72 - t257 * t94 + t260 * t96 + t370;
t71 = t135 * t260 + t137 * t257;
t95 = qJD(1) * t136 - t248 * t284;
t97 = qJD(1) * t138 - t248 * t285;
t266 = qJD(1) * t133 - qJD(3) * t71 - t257 * t95 + t260 * t97;
t320 = -(t248 * t458 + t249 * t266) * t249 + (t248 * t459 + t249 * t265) * t248;
t319 = -(t248 * t266 - t249 * t458) * t249 + (t248 * t265 - t249 * t459) * t248;
t318 = t248 * t26 - t249 * t25;
t317 = t248 * t25 + t249 * t26;
t316 = t248 * t28 - t249 * t27;
t315 = t248 * t27 + t249 * t28;
t314 = t248 * t32 - t249 * t31;
t313 = t248 * t31 + t249 * t32;
t53 = -t136 * t398 - t331;
t312 = t248 * t53 - t249 * t52;
t54 = -t135 * t393 - t385;
t55 = -t136 * t393 + t384;
t311 = t248 * t55 - t249 * t54;
t310 = -t248 * t67 - t422;
t309 = -t248 * t91 - t249 * t90;
t281 = -t248 * t365 - t346;
t100 = rSges(4,1) * t281 - rSges(4,2) * t345 + t374;
t166 = t212 * t248;
t101 = -qJD(3) * t166 + (t214 * t249 + t241) * qJD(1);
t301 = t100 * t249 + t101 * t248;
t298 = t139 * t248 + t140 * t249;
t296 = t207 * t260 + t209 * t257;
t169 = rSges(5,1) * t389 - rSges(5,2) * t390 + t251;
t292 = t257 * t41 + t361 * t82;
t291 = t257 * t42 + t361 * t80;
t234 = pkin(5) * t367;
t287 = qJD(1) * (-pkin(2) * t368 + t234) + qJDD(1) * t191 + t328;
t286 = t257 * t440 - pkin(2) - t253;
t279 = t154 * t237 + t177 * t82 - t178 * t80;
t278 = (-Icges(5,5) * t150 - Icges(5,6) * t151) * t178 - (Icges(5,5) * t152 - Icges(5,6) * t153) * t177 - t184 * t237;
t277 = -t257 * t380 + t260 * t382;
t276 = t257 * t278;
t270 = (-t257 * t372 + t260 * t373) * qJD(1);
t268 = (Icges(5,1) * t152 - t413 - t85) * t177 - (-Icges(5,1) * t150 - t144 - t83) * t178 + (-t156 + t186) * t237;
t30 = -t177 * t90 + t178 * t91 + qJD(2) + (t171 * t248 + t173 * t249) * qJD(3);
t38 = -t224 * t364 - t168 * t177 + t237 * t91 + (t173 + t462) * qJD(1);
t267 = t30 * t309 + (t248 * t37 - t249 * t38) * t168;
t194 = t304 * qJD(3);
t195 = t210 * qJD(3);
t264 = qJD(1) * t205 - qJD(3) * t296 - t194 * t257 + t195 * t260;
t263 = t455 * t257;
t228 = pkin(6) * t391;
t226 = pkin(6) * t396;
t204 = pkin(6) * t345;
t196 = t214 * qJD(3);
t183 = qJD(1) * t190;
t179 = t191 * qJD(1);
t172 = -pkin(3) * t393 + t228;
t170 = -pkin(3) * t398 + t226;
t130 = -t249 * t354 + t375;
t129 = -t248 * t354 + t376;
t128 = t158 * t249;
t127 = t158 * t248;
t126 = t156 * t249;
t125 = t156 * t248;
t115 = t282 * pkin(6) + (-t248 * t362 + t249 * t365) * pkin(3);
t114 = pkin(3) * t281 - pkin(6) * t350 + t204;
t111 = rSges(5,1) * t152 - rSges(5,2) * t153;
t110 = -rSges(5,1) * t150 - rSges(5,2) * t151;
t99 = -t249 * t295 + t401;
t79 = t99 * qJD(1);
t65 = qJD(3) * t298 + qJD(2);
t48 = rSges(5,3) * t282 - t325;
t47 = -rSges(5,3) * t350 + t355;
t40 = t248 * t264 - t249 * t457;
t39 = t248 * t457 + t249 * t264;
t36 = qJD(1) * t100 + qJDD(1) * t140 - t180 * t212 - t196 * t364 + t287;
t35 = -t356 - t196 * t363 + t181 * t212 + (-t101 - t179) * qJD(1) + t327 * qJDD(1);
t34 = -qJD(3) * t299 + t257 * t96 + t260 * t94;
t33 = -qJD(3) * t300 + t257 * t97 + t260 * t95;
t29 = qJD(3) * t301 + t139 * t180 - t140 * t181 + qJDD(2);
t23 = qJD(3) * t311 + t79;
t22 = qJD(3) * t312 + t387;
t16 = t116 * t398 - t117 * t150 + t118 * t151 + t154 * t282 + t156 * t77 + t158 * t78;
t15 = t116 * t393 + t117 * t152 + t118 * t153 + t154 * t280 + t156 * t75 + t158 * t76;
t14 = t177 * t32 - t178 * t31 + t237 * t58;
t13 = qJD(1) * t114 + qJDD(1) * t173 - t102 * t168 - t119 * t177 - t180 * t224 + t192 * t91 - t199 * t364 + t237 * t47 + t287;
t12 = -t356 - t199 * t363 + t103 * t168 - t119 * t178 + t181 * t224 + t192 * t90 - t237 * t48 + (-t115 - t179) * qJD(1) + t326 * qJDD(1);
t9 = -t102 * t90 - t103 * t91 + t171 * t180 - t173 * t181 + t177 * t48 + t178 * t47 + qJDD(2) + (t114 * t249 + t115 * t248) * qJD(3);
t6 = -t150 * t43 + t151 * t45 + t248 * t292 + t349 * t82 + t77 * t85 + t78 * t88;
t5 = -t150 * t44 + t151 * t46 + t248 * t291 + t349 * t80 + t77 * t83 - t78 * t87;
t4 = t152 * t43 + t153 * t45 + t249 * t292 - t350 * t82 + t75 * t85 + t76 * t88;
t3 = t152 * t44 + t153 * t46 + t249 * t291 - t350 * t80 + t75 * t83 - t76 * t87;
t2 = t102 * t26 + t103 * t25 + t16 * t237 + t177 * t6 - t178 * t5 + t192 * t50;
t1 = t102 * t28 + t103 * t27 + t15 * t237 + t177 * t4 - t178 * t3 + t192 * t51;
t17 = [t15 * t450 + t50 * t452 + t51 * t453 + t16 * t449 - t437 / 0.2e1 + t436 / 0.2e1 + t435 + t420 / 0.2e1 + t421 / 0.2e1 + (t79 + ((t53 - t120 + (t134 + t403) * t249 + t385) * t249 + t384 * t248) * qJD(3)) * t337 + (-qJD(3) * t295 + t194 * t260 + t195 * t257) * qJD(1) - m(2) * (-g(1) * t213 + g(2) * t215) + (t72 + t99) * t447 + (t71 + t98) * t446 + (-t387 + ((t249 * t330 - t384 + t55) * t249 + (t248 * t330 + t331 + t54) * t248) * qJD(3) + t22) * t340 + (t34 + t39) * t339 + (t449 + t448) * t11 + ((-t188 * t262 - g(2) + t328) * t176 + (-t356 + (-0.2e1 * t189 - t254 + t176) * t262 - g(1)) * t175) * m(3) + (t33 + t40 + t23) * t338 + (m(2) * (t213 ^ 2 + t215 ^ 2) + m(3) * (t175 ^ 2 + t189 * t176) + t296 + Icges(2,3) + Icges(3,3)) * qJDD(1) + (t286 * t423 * qJD(1) + (t13 - g(2)) * (t462 + t416) + (t325 + (t260 * t440 + t438) * t364 + (-t254 - t475) * qJD(1)) * t37 + (-pkin(3) * t346 + t183 + t204 + t234 + t355 + t37 + ((-pkin(2) - t463 - t251) * t248 + t171) * qJD(1) + t474) * t38 + (-g(1) + t12) * (t248 * t286 - t322 + t332)) * m(5) + (t67 * (t234 + t374) + (t212 * t424 - t428) * qJD(3) + ((-t258 * t67 - t261 * t66) * pkin(1) + (-pkin(2) - t214) * t422 + (t66 * (-rSges(4,3) - pkin(5)) + t67 * t344) * t248) * qJD(1) - (-t348 - t183 - t66 + (-t139 - t439) * qJD(1)) * t67 + (t36 - g(2)) * t113 + (t35 - g(1)) * (t248 * t344 + t332 + t371)) * m(4); m(3) * qJDD(2) + m(4) * t29 + m(5) * t9 + (-m(3) - m(4) - m(5)) * g(3); ((-t126 * t152 - t128 * t153) * t177 - (-t125 * t152 - t127 * t153) * t178 + (t152 * t157 + t153 * t159) * t237 + (t257 * t51 + t27 * t396) * qJD(4) + ((qJD(4) * t28 + t279) * t260 + t263) * t249) * t451 + ((t126 * t150 - t128 * t151) * t177 - (t125 * t150 - t127 * t151) * t178 + (-t150 * t157 + t151 * t159) * t237 + (t257 * t50 + t26 * t391) * qJD(4) + ((qJD(4) * t25 + t279) * t260 + t263) * t248) * t448 + (((t126 * t256 - t128 * t259 + t82) * t177 - (t125 * t256 - t127 * t259 + t80) * t178 + (-t157 * t256 + t159 * t259 + t154) * t237 + t58 * qJD(4)) * t257 + (qJD(4) * t313 - t455) * t260) * t444 + (qJD(1) * t315 + t248 * t4 - t249 * t3) * t450 + t318 * t452 + t316 * t453 + (qJD(1) * t313 + t248 * t8 - t249 * t7) * t443 + t314 * t445 + t312 * t446 + t311 * t447 + (qJD(1) * t317 + t248 * t6 - t249 * t5) * t449 - qJD(1) * ((t257 * t373 + t260 * t372) * qJD(1) + ((t248 * t380 - t249 * t381) * t260 + (t248 * t382 + t249 * t383) * t257) * qJD(3)) / 0.2e1 + (t23 + t11) * t341 - t14 * t360 / 0.2e1 - (qJD(1) * t40 + qJD(3) * t319 + qJDD(1) * t98 + t180 * t53 + t181 * t52 + t2) * t249 / 0.2e1 + (qJD(1) * t39 + qJD(3) * t320 + qJDD(1) * t99 + t180 * t55 + t181 * t54 + t1) * t248 / 0.2e1 + ((-t363 * t401 - t369) * t249 + (t270 + (t277 * t248 + (t400 - t456) * t249) * qJD(3)) * t248) * t337 + ((-t364 * t400 + t369) * t248 + (t270 + (-t456 * t249 + (t401 + t277) * t248) * qJD(3)) * t249) * t340 + (t22 + t10) * t368 / 0.2e1 + qJD(1) * (t248 * t34 - t249 * t33 + (t71 * t248 + t249 * t72) * qJD(1)) / 0.2e1 + (t29 * t298 + t65 * ((t139 * t249 - t140 * t248) * qJD(1) + t301) + t310 * t196 + (-t36 * t248 - t35 * t249 + (-t249 * t67 + t424) * qJD(1)) * t212 - (t166 * t66 - t428) * qJD(1) - (t65 * (-t166 * t248 - t167 * t249) + t310 * t214) * qJD(3) + g(1) * t167 + g(2) * t166 - g(3) * t214) * m(4) + ((t248 * t52 + t53 * t249) * qJD(1) + t319) * t338 - t465 * t359 / 0.2e1 + ((-t12 * t377 + t37 * t386 + t9 * t416 + t30 * (t114 + t47) + (t30 * t417 - t377 * t38) * qJD(1)) * t249 + (-t13 * t377 + t38 * t386 + t9 * t417 + t30 * (t115 + t48) + (-t30 * t416 + t37 * t377) * qJD(1)) * t248 - t37 * (-qJD(1) * t170 - t129 * t237 - t169 * t178 - t363 * t463) - t38 * (qJD(1) * t172 + t130 * t237 - t169 * t177 - t364 * t463) - t30 * (t129 * t177 + t130 * t178 + t170 * t364 + t172 * t363) - ((t37 * t90 + t38 * t91) * t257 + t267 * t260) * qJD(4) - g(1) * (t228 + t375) - g(2) * (t226 + t376) - g(3) * (t169 + t463) - (g(1) * t249 + g(2) * t248) * t257 * (-pkin(3) - t433)) * m(5) + ((t54 * t248 + t249 * t55) * qJD(1) + t320) * t339 + qJDD(1) * (t248 * t72 - t249 * t71) / 0.2e1; -t11 * t350 / 0.2e1 + t1 * t393 / 0.2e1 + (t257 * t315 - t260 * t51) * t453 + ((qJD(3) * t315 - t15) * t260 + (-qJD(1) * t316 + qJD(3) * t51 + t248 * t3 + t249 * t4) * t257) * t450 + t257 * t10 * t341 + t2 * t398 / 0.2e1 + (t257 * t317 - t260 * t50) * t452 + ((qJD(3) * t317 - t16) * t260 + (-qJD(1) * t318 + qJD(3) * t50 + t248 * t5 + t249 * t6) * t257) * t449 + t14 * t362 / 0.2e1 - t260 * (t420 + t421 + t435 + t436 - t437) / 0.2e1 + (t257 * t313 - t260 * t58) * t445 + ((qJD(3) * t313 - t24) * t260 + (-qJD(1) * t314 + qJD(3) * t58 + t248 * t7 + t249 * t8) * t257) * t443 + (t152 * t454 + t268 * t153 - t249 * t276) * t451 + (-t150 * t454 + t151 * t268 - t248 * t276) * t448 + (t278 * t260 + (-t256 * t454 + t259 * t268) * t257) * t444 + t465 * t361 / 0.2e1 + ((qJD(3) * t267 - t12 * t90 - t13 * t91 + t37 * t48 - t38 * t47) * t260 + (t37 * (qJD(3) * t90 + t119 * t248) + t38 * (qJD(3) * t91 - t119 * t249) + t9 * t309 + t30 * (-t248 * t47 + t249 * t48 - t367 * t91 + t368 * t90) + (t12 * t248 - t13 * t249 + (t248 * t38 + t423) * qJD(1)) * t168) * t257 - t37 * (-t110 * t237 - t178 * t187) - t38 * (t111 * t237 - t177 * t187) - t30 * (t110 * t177 + t111 * t178) - g(1) * t111 - g(2) * t110 - g(3) * t187) * m(5);];
tau = t17;
