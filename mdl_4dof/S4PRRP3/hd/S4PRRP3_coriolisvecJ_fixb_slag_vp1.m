% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:53
% DurationCPUTime: 7.87s
% Computational Cost: add. (5417->349), mult. (6457->448), div. (0->0), fcn. (5050->4), ass. (0->209)
t410 = Icges(4,3) + Icges(5,3);
t187 = pkin(6) + qJ(2);
t183 = sin(t187);
t190 = cos(qJ(3));
t282 = t183 * t190;
t189 = sin(qJ(3));
t283 = t183 * t189;
t184 = cos(t187);
t292 = Icges(5,6) * t184;
t87 = Icges(5,4) * t282 - Icges(5,2) * t283 - t292;
t293 = Icges(4,6) * t184;
t89 = Icges(4,4) * t282 - Icges(4,2) * t283 - t293;
t401 = t87 + t89;
t138 = Icges(5,5) * t190 - Icges(5,6) * t189;
t140 = Icges(4,5) * t190 - Icges(4,6) * t189;
t391 = t138 + t140;
t157 = Icges(5,4) * t283;
t296 = Icges(5,5) * t184;
t91 = Icges(5,1) * t282 - t157 - t296;
t158 = Icges(4,4) * t283;
t297 = Icges(4,5) * t184;
t93 = Icges(4,1) * t282 - t158 - t297;
t406 = t91 + t93;
t298 = Icges(5,4) * t189;
t146 = Icges(5,1) * t190 - t298;
t92 = Icges(5,5) * t183 + t146 * t184;
t299 = Icges(4,4) * t189;
t148 = Icges(4,1) * t190 - t299;
t94 = Icges(4,5) * t183 + t148 * t184;
t399 = t92 + t94;
t141 = Icges(5,2) * t190 + t298;
t143 = Icges(4,2) * t190 + t299;
t409 = t141 + t143;
t185 = Icges(5,4) * t190;
t145 = Icges(5,1) * t189 + t185;
t186 = Icges(4,4) * t190;
t147 = Icges(4,1) * t189 + t186;
t408 = t145 + t147;
t407 = t410 * t184;
t379 = t407 + (Icges(4,6) + Icges(5,6)) * t283 + (-Icges(4,5) - Icges(5,5)) * t282;
t394 = t410 * t183 + t391 * t184;
t225 = -Icges(5,2) * t189 + t185;
t88 = Icges(5,6) * t183 + t184 * t225;
t226 = -Icges(4,2) * t189 + t186;
t90 = Icges(4,6) * t183 + t184 * t226;
t400 = t88 + t90;
t405 = t146 + t148;
t404 = t225 + t226;
t403 = t401 * t189;
t375 = -t406 * t190 + t403;
t372 = -t183 * t375 + t379 * t184;
t402 = t399 * t282;
t137 = Icges(5,5) * t189 + Icges(5,6) * t190;
t139 = Icges(4,5) * t189 + Icges(4,6) * t190;
t397 = t139 + t137;
t396 = t409 * qJD(3);
t395 = t408 * qJD(3);
t390 = t409 * t189 - t408 * t190;
t393 = t404 * qJD(3);
t392 = t405 * qJD(3);
t389 = t394 * t184 - t402;
t388 = t400 * t189;
t280 = t184 * t190;
t347 = -t394 * t183 - t399 * t280;
t386 = t379 * t183 - t406 * t280;
t371 = -t400 * t283 - t389;
t281 = t184 * t189;
t370 = -t401 * t281 - t386;
t369 = -t400 * t281 - t347;
t353 = t406 * t189 + t401 * t190;
t352 = t399 * t189 + t400 * t190;
t385 = -t396 * t184 + (-t183 * t404 + t292 + t293) * qJD(2);
t384 = t400 * qJD(2) - t396 * t183;
t383 = -t395 * t184 + (-t183 * t405 + t296 + t297) * qJD(2);
t382 = -t399 * qJD(2) + t395 * t183;
t284 = t139 * t184;
t286 = t137 * t184;
t381 = -t390 * t183 - t284 - t286;
t285 = t139 * t183;
t287 = t137 * t183;
t380 = -t184 * t390 + t285 + t287;
t378 = t397 * qJD(3);
t377 = t399 * t190 - t388;
t376 = t392 * t190 - t393 * t189 + (-t408 * t189 - t409 * t190) * qJD(3) + t397 * qJD(2);
t374 = t394 * qJD(2);
t373 = qJD(2) * t390 + qJD(3) * t391;
t368 = t380 * qJD(2);
t367 = -t352 * qJD(3) - t385 * t189 + t383 * t190 + t374;
t366 = t379 * qJD(2) + t353 * qJD(3) + t384 * t189 + t382 * t190;
t365 = (t369 * t183 - t370 * t184) * qJD(3);
t364 = (t371 * t183 - t372 * t184) * qJD(3);
t363 = t381 * qJD(2);
t362 = t375 * qJD(2) - t378 * t183 + t374;
t361 = -t378 * t184 + (-t391 * t183 - t377 + t407) * qJD(2);
t360 = 0.2e1 * qJD(3);
t359 = t363 + t364;
t358 = t365 + t368;
t357 = t375 * qJD(3) + t382 * t189 - t384 * t190;
t356 = t377 * qJD(3) + t383 * t189 + t385 * t190;
t355 = t373 * t183 + t184 * t376;
t354 = t183 * t376 - t373 * t184;
t175 = t183 * rSges(5,3);
t336 = pkin(3) * t190;
t182 = pkin(2) + t336;
t188 = -qJ(4) - pkin(5);
t351 = rSges(5,1) * t280 - rSges(5,2) * t281 + t184 * t182 - t183 * t188 + t175;
t163 = t184 * t188;
t350 = -rSges(5,1) * t282 + rSges(5,2) * t283 + t184 * rSges(5,3) - t183 * t182 - t163;
t167 = qJD(4) * t183;
t322 = rSges(5,2) * t190;
t149 = rSges(5,1) * t189 + t322;
t248 = pkin(3) * t189 + t149;
t263 = qJD(3) * t184;
t208 = -t248 * t263 + t167;
t266 = qJD(2) * t183;
t260 = t189 * t266;
t265 = qJD(2) * t184;
t349 = rSges(5,2) * t260 + rSges(5,3) * t265 + t167;
t122 = t184 * pkin(2) + t183 * pkin(5);
t348 = t379 + t388;
t304 = -t147 * t183 - t89;
t308 = -Icges(4,2) * t282 - t158 + t93;
t340 = -t189 * t308 + t190 * t304;
t306 = -t145 * t183 - t87;
t310 = -Icges(5,2) * t282 - t157 + t91;
t339 = -t189 * t310 + t190 * t306;
t338 = t183 / 0.2e1;
t337 = -t184 / 0.2e1;
t334 = qJD(2) / 0.2e1;
t333 = pkin(2) - t182;
t166 = pkin(5) * t265;
t257 = t189 * t263;
t209 = -t190 * t266 - t257;
t262 = qJD(3) * t190;
t256 = t184 * t262;
t332 = -pkin(3) * t257 - t166 + (t183 * t333 - t163) * qJD(2) + rSges(5,1) * t209 - rSges(5,2) * t256 + t349;
t111 = t149 * t183;
t323 = rSges(5,1) * t190;
t151 = -rSges(5,2) * t189 + t323;
t168 = qJD(4) * t184;
t261 = pkin(3) * t283;
t274 = qJD(3) * t261 + t168;
t331 = -qJD(3) * t111 - t274 + ((-pkin(5) - t188) * t183 + t175 + (-t333 + t151) * t184) * qJD(2);
t180 = t184 * pkin(5);
t121 = pkin(2) * t183 - t180;
t326 = t121 + t350;
t325 = -t122 + t351;
t324 = rSges(4,1) * t190;
t150 = rSges(4,1) * t189 + rSges(4,2) * t190;
t114 = t150 * t184;
t264 = qJD(3) * t183;
t259 = t150 * t264;
t176 = t183 * rSges(4,3);
t98 = rSges(4,1) * t280 - rSges(4,2) * t281 + t176;
t302 = t122 + t98;
t42 = qJD(2) * t302 - t259;
t321 = t114 * t42;
t258 = t150 * t263;
t269 = rSges(4,2) * t283 + t184 * rSges(4,3);
t96 = rSges(4,1) * t282 - t269;
t41 = -t258 + (-t121 - t96) * qJD(2);
t320 = t183 * t41;
t319 = t184 * t41;
t28 = (-t121 + t326) * qJD(2) + t208;
t311 = t28 * t149;
t309 = -t141 * t184 + t92;
t307 = -t143 * t184 + t94;
t305 = -t145 * t184 - t88;
t303 = -t147 * t184 - t90;
t275 = rSges(4,2) * t260 + rSges(4,3) * t265;
t273 = -t141 + t146;
t272 = t145 + t225;
t271 = -t143 + t148;
t270 = t147 + t226;
t268 = qJD(2) * t138;
t267 = qJD(2) * t140;
t255 = -pkin(2) - t324;
t252 = -t264 / 0.2e1;
t249 = t263 / 0.2e1;
t247 = -t151 - t336;
t240 = -pkin(3) * t281 - t149 * t184;
t129 = t151 * qJD(3);
t237 = -pkin(3) * t262 - t129;
t235 = -rSges(4,2) * t189 + t324;
t234 = -t183 * t42 - t319;
t233 = t183 * t96 + t184 * t98;
t221 = (-t336 * qJD(3) - t129) * qJD(3);
t112 = t150 * t183;
t207 = -t189 * t309 + t190 * t305;
t206 = -t189 * t307 + t190 * t303;
t205 = (-t189 * t272 + t190 * t273) * qJD(2);
t204 = (-t189 * t270 + t190 * t271) * qJD(2);
t68 = rSges(4,1) * t209 - rSges(4,2) * t256 + t275;
t70 = -qJD(3) * t112 + (t184 * t235 + t176) * qJD(2);
t203 = t183 * t70 + t184 * t68 + (-t183 * t98 + t184 * t96) * qJD(2);
t130 = t235 * qJD(3);
t118 = qJD(2) * t121;
t117 = t122 * qJD(2);
t115 = qJD(2) * (-pkin(2) * t266 + t166);
t40 = qJD(3) * t233 + qJD(1);
t31 = -t130 * t263 + (-t117 - t70 + t259) * qJD(2);
t30 = -t130 * t264 + t115 + (t68 - t258) * qJD(2);
t29 = -t149 * t264 + (t122 + t325) * qJD(2) - t274;
t25 = qJD(1) + (-t183 * t326 + t184 * t325) * qJD(3);
t16 = t203 * qJD(3);
t15 = t221 * t184 + (t248 * t264 - t117 + t168 - t331) * qJD(2);
t14 = t115 + t221 * t183 + (t332 + t208) * qJD(2);
t1 = (t332 * t184 + t331 * t183 + (-t183 * t325 - t184 * t326) * qJD(2)) * qJD(3);
t2 = [m(4) * t16 + m(5) * t1; ((((t394 + t403) * t184 + t371 + t386 - t402) * t184 - t347 * t183) * qJD(3) + t368) * t249 + (-t390 * qJD(3) + t392 * t189 + t190 * t393) * qJD(2) + (t15 * t350 + t28 * t274 + t14 * t351 + t29 * t349 + (t183 * t311 + t29 * (-t322 + (-rSges(5,1) - pkin(3)) * t189) * t184) * qJD(3) + ((t28 * (-t151 - t182) - t29 * t188) * t184 + (t28 * (-rSges(5,3) + t188) + t29 * (-t182 - t323)) * t183) * qJD(2) - (qJD(2) * t326 - t118 + t208 - t28) * t29) * m(5) + (t31 * (t183 * t255 + t180 + t269) + t30 * t302 + t42 * (t166 + t275) + (t150 * t320 - t321) * qJD(3) + ((-pkin(2) - t235) * t319 + (t41 * (-rSges(4,3) - pkin(5)) + t42 * t255) * t183) * qJD(2) - (-qJD(2) * t96 - t118 - t258 - t41) * t42) * m(4) + (((t184 * t348 + t347 + t369) * t184 + (t183 * t348 + t370 + t389) * t183) * qJD(3) + t359 - t363) * t252 + (t355 + t356) * t264 / 0.2e1 - (t354 - t357 + t358) * t263 / 0.2e1 + ((t353 + t381) * t183 + (t352 + t380) * t184) * qJD(3) * t334; -((((-t308 - t310) * t184 + (t307 + t309) * t183) * t190 + ((-t304 - t306) * t184 + (t303 + t305) * t183) * t189) * qJD(3) + ((t270 + t272) * t190 + (t271 + t273) * t189) * qJD(2)) * qJD(2) / 0.2e1 + (t357 * t184 + t356 * t183 + (t183 * t353 + t184 * t352) * qJD(2)) * t334 + ((-t264 * t286 + t268) * t183 + (t205 + (-t339 * t184 + (t287 + t207) * t183) * qJD(3)) * t184 + (-t264 * t284 + t267) * t183 + (t204 + (-t340 * t184 + (t285 + t206) * t183) * qJD(3)) * t184) * t252 + ((-t263 * t287 - t268) * t184 + (t205 + (t207 * t183 + (t286 - t339) * t184) * qJD(3)) * t183 + (-t263 * t285 - t267) * t184 + (t204 + (t206 * t183 + (t284 - t340) * t184) * qJD(3)) * t183) * t249 + (t16 * t233 + t40 * t203 + t234 * t130 + (-t30 * t183 - t31 * t184 + (-t184 * t42 + t320) * qJD(2)) * t150 - (t112 * t41 - t321) * qJD(2) - (t40 * (-t112 * t183 - t114 * t184) + t234 * t235) * qJD(3)) * m(4) + (t355 * qJD(2) + ((t369 * qJD(2) + t366 * t184) * t184 + (t361 * t183 + t370 * qJD(2) + (-t362 + t367) * t184) * t183) * t360) * t338 + (t354 * qJD(2) + ((t371 * qJD(2) + t362 * t184) * t184 + (t367 * t183 + t372 * qJD(2) + (-t361 + t366) * t184) * t183) * t360) * t337 + (-(t28 * t111 + t240 * t29) * qJD(2) - ((t240 * t25 + t247 * t28) * t184 + (t29 * t247 + (-t111 - t261) * t25) * t183) * qJD(3) + (-t14 * t248 + t29 * t237 - t1 * t326 + t25 * t331 + (-t25 * t325 + t311) * qJD(2)) * t183 + (-t15 * t248 + t28 * t237 + t1 * t325 + t25 * t332 + (-t248 * t29 - t25 * t326) * qJD(2)) * t184) * m(5) + (t359 + t364) * t266 / 0.2e1 + (t358 + t365) * t265 / 0.2e1; 0.2e1 * (t14 * t337 + t15 * t338) * m(5);];
tauc = t2(:);
