% Calculate time derivative of joint inertia matrix for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:52
% DurationCPUTime: 8.73s
% Computational Cost: add. (15912->354), mult. (21983->553), div. (0->0), fcn. (20896->8), ass. (0->207)
t429 = Icges(5,4) - Icges(6,5);
t430 = Icges(5,1) + Icges(6,1);
t419 = Icges(6,4) + Icges(5,5);
t428 = Icges(5,2) + Icges(6,3);
t435 = Icges(6,2) + Icges(5,3);
t418 = Icges(6,6) - Icges(5,6);
t229 = cos(qJ(4));
t434 = t429 * t229;
t227 = sin(qJ(4));
t433 = t429 * t227;
t224 = qJ(2) + qJ(3);
t219 = sin(t224);
t225 = sin(pkin(8));
t226 = cos(pkin(8));
t220 = cos(t224);
t308 = qJD(4) * t229;
t299 = t220 * t308;
t223 = qJD(2) + qJD(3);
t325 = t223 * t225;
t157 = -t225 * t299 + (qJD(4) * t226 + t219 * t325) * t227;
t320 = t226 * t229;
t323 = t225 * t227;
t206 = t220 * t323 + t320;
t331 = t219 * t223;
t305 = t229 * t331;
t158 = -qJD(4) * t206 - t225 * t305;
t304 = t220 * t325;
t432 = -t418 * t157 + t419 * t158 + t435 * t304;
t309 = qJD(4) * t227;
t321 = t226 * t227;
t159 = -t225 * t309 - t226 * t299 + t321 * t331;
t322 = t225 * t229;
t208 = t220 * t321 - t322;
t160 = -qJD(4) * t208 - t226 * t305;
t324 = t223 * t226;
t303 = t220 * t324;
t431 = -t418 * t159 + t419 * t160 + t435 * t303;
t207 = t220 * t322 - t321;
t330 = t219 * t225;
t415 = t428 * t206 - t429 * t207 + t418 * t330;
t209 = t220 * t320 + t323;
t329 = t219 * t226;
t414 = t428 * t208 - t429 * t209 + t418 * t329;
t413 = t418 * t206 + t419 * t207 + t435 * t330;
t405 = t418 * t208 + t419 * t209 + t435 * t329;
t412 = -t429 * t206 + t430 * t207 + t419 * t330;
t411 = -t429 * t208 + t430 * t209 + t419 * t329;
t426 = t428 * t227 - t434;
t425 = -t418 * t227 - t419 * t229;
t424 = t430 * t229 - t433;
t410 = -t425 * t219 - t220 * t435;
t389 = t410 * t220;
t423 = t414 * t206 + t411 * t207 + t405 * t330;
t422 = t415 * t208 + t412 * t209 + t413 * t329;
t421 = t428 * t157 + t429 * t158 - t418 * t304;
t420 = -t428 * t159 - t429 * t160 + t418 * t303;
t417 = t429 * t157 + t430 * t158 + t419 * t304;
t416 = t429 * t159 + t430 * t160 + t419 * t303;
t392 = t426 * t219 - t418 * t220;
t409 = t424 * t219 - t419 * t220;
t326 = t220 * t223;
t408 = t431 * t219 + t405 * t326;
t407 = t432 * t219 + t413 * t326;
t402 = rSges(6,1) + pkin(4);
t395 = rSges(6,3) + qJ(5);
t404 = t414 * t227 + t411 * t229;
t403 = t415 * t227 + t412 * t229;
t401 = t415 * t157 - t412 * t158 + t421 * t206 - t417 * t207 - t407 * t225;
t400 = -t414 * t157 + t411 * t158 + t420 * t206 + t416 * t207 + t408 * t225;
t399 = -t415 * t159 + t412 * t160 - t421 * t208 + t417 * t209 + t407 * t226;
t398 = -t414 * t159 + t411 * t160 + t420 * t208 + t416 * t209 + t408 * t226;
t374 = t415 * t206 + t412 * t207 + t413 * t330;
t373 = t414 * t208 + t411 * t209 + t405 * t329;
t397 = t392 * t206 + t409 * t207 + t410 * t330;
t396 = t392 * t208 + t409 * t209 + t410 * t329;
t221 = t225 ^ 2;
t222 = t226 ^ 2;
t383 = t221 + t222;
t369 = t223 * t383;
t394 = t426 * t326 + (t418 * t223 + (t428 * t229 + t433) * qJD(4)) * t219;
t393 = t424 * t326 + (t419 * t223 + (-t430 * t227 - t434) * qJD(4)) * t219;
t390 = (t425 * t326 + (-t435 * t223 + (t419 * t227 - t418 * t229) * qJD(4)) * t219) * t220;
t388 = t423 * t226;
t387 = t422 * t225;
t386 = t403 * t219 - t220 * t413;
t385 = t404 * t219 - t405 * t220;
t344 = rSges(6,2) * t304 + qJD(5) * t206 - t395 * t157 + t402 * t158;
t384 = rSges(6,2) * t303 + qJD(5) * t208 - t395 * t159 + t402 * t160;
t382 = -t392 * t227 - t229 * t409;
t381 = -Icges(4,5) * t219 - Icges(4,6) * t220;
t380 = (((-t389 + t374) * t225 + t388) * t223 - t393 * t207 - t394 * t206 - t409 * t158 + t392 * t157) * t220 + (t400 * t226 + (t390 - t401) * t225 + t397 * t223) * t219;
t379 = (((-t389 + t373) * t226 + t387) * t223 - t393 * t209 - t394 * t208 - t409 * t160 + t392 * t159) * t220 + ((t390 + t398) * t226 + t399 * t225 + t396 * t223) * t219;
t378 = t400 * t225 + t401 * t226;
t377 = t398 * t225 - t399 * t226;
t376 = (-t403 * t223 + t432) * t220 + (-t417 * t229 + t421 * t227 - t413 * t223 + (t412 * t227 - t415 * t229) * qJD(4)) * t219;
t375 = (t404 * t223 - t431) * t220 + (t416 * t229 + t420 * t227 + t405 * t223 + (-t411 * t227 + t414 * t229) * qJD(4)) * t219;
t317 = rSges(6,2) * t330 + t395 * t206 + t402 * t207;
t316 = rSges(6,2) * t329 + t395 * t208 + t402 * t209;
t372 = t395 * t227 + t402 * t229;
t241 = t223 * t381;
t171 = t225 * t241;
t172 = t226 * t241;
t359 = t381 * t369;
t371 = -t222 * t171 - (-t226 * t172 + t359) * t225 - t378;
t368 = t219 * t382 + t389;
t107 = rSges(5,1) * t158 + rSges(5,2) * t157 + rSges(5,3) * t304;
t109 = rSges(5,1) * t160 + rSges(5,2) * t159 + rSges(5,3) * t303;
t367 = t225 * t107 + t226 * t109;
t364 = t225 * t386 + t226 * t385;
t228 = sin(qJ(2));
t230 = cos(qJ(2));
t363 = qJD(2) * (rSges(3,1) * t228 + rSges(3,2) * t230);
t362 = t225 * t344 + t226 * t384;
t358 = 2 * m(4);
t357 = 2 * m(5);
t356 = 2 * m(6);
t352 = pkin(2) * t228;
t349 = pkin(2) * qJD(2);
t319 = t372 * t326 + (rSges(6,2) * t223 + qJD(5) * t227 + (-t402 * t227 + t395 * t229) * qJD(4)) * t219;
t280 = rSges(5,1) * t229 - rSges(5,2) * t227;
t125 = t280 * t326 + (rSges(5,3) * t223 + (-rSges(5,1) * t227 - rSges(5,2) * t229) * qJD(4)) * t219;
t283 = pkin(3) * t220 + pkin(7) * t219;
t198 = t283 * t223;
t318 = -t125 - t198;
t210 = rSges(4,1) * t219 + rSges(4,2) * t220;
t117 = t210 * t369;
t211 = pkin(3) * t219 - pkin(7) * t220;
t315 = t211 * t369;
t314 = t383 * pkin(2) * t230;
t281 = rSges(4,1) * t220 - rSges(4,2) * t219;
t126 = t383 * t281;
t313 = -rSges(6,2) * t220 + t372 * t219;
t184 = -rSges(5,3) * t220 + t219 * t280;
t312 = -t184 - t211;
t311 = t383 * t283;
t307 = (t221 * t172 + (-t225 * t171 + t359) * t226 + t377) * t225;
t306 = t230 * t349;
t302 = t227 * t326;
t301 = -t198 - t319;
t300 = -t211 - t313;
t294 = -t210 - t352;
t293 = -t211 - t352;
t292 = t225 * t313;
t291 = t313 * t226;
t147 = rSges(5,1) * t207 - rSges(5,2) * t206 + rSges(5,3) * t330;
t149 = rSges(5,1) * t209 - rSges(5,2) * t208 + rSges(5,3) * t329;
t71 = t225 * t147 + t226 * t149 + t311;
t286 = -t184 + t293;
t197 = t281 * t223;
t285 = -t197 - t306;
t284 = -t198 - t306;
t259 = t147 * t226 - t149 * t225;
t253 = t293 - t313;
t252 = t383 * t228 * t349;
t251 = -t125 + t284;
t52 = t317 * t225 + t316 * t226 + t311;
t244 = t284 - t319;
t237 = qJD(2) * (-Icges(3,5) * t228 - Icges(3,6) * t230);
t236 = t371 * t226 + t307;
t235 = t219 * t308 + t302;
t234 = -t252 - t315;
t233 = -t225 * t316 + t317 * t226;
t232 = -(t375 * t225 + t376 * t226) * t220 / 0.2e1 + t379 * t225 / 0.2e1 - t380 * t226 / 0.2e1 + (t225 * t385 - t226 * t386) * t331 / 0.2e1 + t378 * t330 / 0.2e1 + t377 * t329 / 0.2e1 + ((t423 * t225 - t374 * t226) * t225 + (t373 * t225 - t422 * t226) * t226) * t326 / 0.2e1;
t200 = t226 * t237;
t199 = t225 * t237;
t186 = t294 * t226;
t185 = t294 * t225;
t162 = t285 * t226;
t161 = t285 * t225;
t152 = t383 * t363;
t151 = t312 * t226;
t150 = t312 * t225;
t129 = t286 * t226;
t128 = t286 * t225;
t116 = t300 * t226;
t115 = t300 * t225;
t114 = t253 * t226;
t113 = t253 * t225;
t112 = -t252 - t117;
t111 = t318 * t226;
t110 = t318 * t225;
t93 = -t149 * t220 - t184 * t329;
t92 = t147 * t220 + t184 * t330;
t87 = t251 * t226;
t86 = t251 * t225;
t79 = t126 + t314;
t78 = t259 * t219;
t73 = t301 * t226;
t72 = t301 * t225;
t70 = t244 * t226;
t69 = t244 * t225;
t68 = -t219 * t291 - t220 * t316;
t67 = t219 * t292 + t220 * t317;
t54 = t71 + t314;
t53 = t233 * t219;
t51 = -t315 + t367;
t50 = (-t184 * t324 - t109) * t220 + (-t125 * t226 + t149 * t223) * t219;
t49 = (t184 * t325 + t107) * t220 + (t125 * t225 - t147 * t223) * t219;
t48 = t234 + t367;
t47 = t52 + t314;
t46 = t259 * t326 + (t107 * t226 - t109 * t225) * t219;
t37 = -t315 + t362;
t36 = t234 + t362;
t35 = (-t223 * t291 - t384) * t220 + (t223 * t316 - t319 * t226) * t219;
t34 = (t223 * t292 + t344) * t220 + (-t223 * t317 + t225 * t319) * t219;
t21 = t233 * t326 + (-t225 * t384 + t344 * t226) * t219;
t1 = [0; -m(3) * t152 + m(4) * t112 + m(5) * t48 + m(6) * t36; (t113 * t69 + t114 * t70 + t36 * t47) * t356 + (t128 * t86 + t129 * t87 + t48 * t54) * t357 + (t112 * t79 + t161 * t185 + t162 * t186) * t358 + t225 * t221 * t200 + t307 + 0.2e1 * m(3) * (-t152 + t363) * t383 * (rSges(3,1) * t230 - rSges(3,2) * t228) + (-t222 * t199 + (-t225 * t199 + t226 * t200) * t225 + t371) * t226; -m(4) * t117 + m(5) * t51 + m(6) * t37; m(6) * (t113 * t72 + t114 * t73 + t115 * t69 + t116 * t70 + t36 * t52 + t37 * t47) + m(5) * (t110 * t128 + t111 * t129 + t150 * t86 + t151 * t87 + t48 * t71 + t51 * t54) + m(4) * (t112 * t126 - t117 * t79 + (-t161 * t225 - t162 * t226) * t210 + (-t185 * t225 - t186 * t226) * t197) + t236; (t115 * t72 + t116 * t73 + t37 * t52) * t356 + (t110 * t150 + t111 * t151 + t51 * t71) * t357 + (t197 * t210 * t383 - t117 * t126) * t358 + t236; m(5) * t46 + m(6) * t21; t232 + m(6) * (t113 * t35 + t114 * t34 + t21 * t47 + t36 * t53 + t67 * t70 + t68 * t69) + m(5) * (t128 * t50 + t129 * t49 + t46 * t54 + t48 * t78 + t86 * t93 + t87 * t92); t232 + m(6) * (t115 * t35 + t116 * t34 + t21 * t52 + t37 * t53 + t67 * t73 + t68 * t72) + m(5) * (t110 * t93 + t111 * t92 + t150 * t50 + t151 * t49 + t46 * t71 + t51 * t78); (t21 * t53 + t34 * t67 + t35 * t68) * t356 + (t46 * t78 + t49 * t92 + t50 * t93) * t357 + t380 * t330 + t379 * t329 + (t364 * t331 + (t374 * t225 + t388) * t304 + (t373 * t226 + t387) * t303) * t219 + (t390 * t220 + t368 * t331 - t397 * t304 - t396 * t303 + (-t220 * t382 - t364) * t326 + ((t394 * t227 + t393 * t229 + t392 * t308 - t309 * t409) * t220 - t375 * t226 + t376 * t225 + (t368 + t389) * t223) * t219) * t220; t235 * m(6); m(6) * (t47 * t302 - t113 * t157 - t114 * t159 + t206 * t69 + t208 * t70 + (t227 * t36 + t308 * t47) * t219); m(6) * (t52 * t302 - t115 * t157 - t116 * t159 + t206 * t72 + t208 * t73 + (t227 * t37 + t308 * t52) * t219); m(6) * (t53 * t302 - t157 * t68 - t159 * t67 + t206 * t35 + t208 * t34 + (t21 * t227 + t308 * t53) * t219); (t219 * t227 * t235 - t157 * t206 - t159 * t208) * t356;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
