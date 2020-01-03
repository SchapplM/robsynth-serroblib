% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:43
% EndTime: 2019-12-31 16:53:56
% DurationCPUTime: 9.33s
% Computational Cost: add. (26882->463), mult. (35787->697), div. (0->0), fcn. (39295->7), ass. (0->278)
t283 = pkin(7) + qJ(3);
t277 = sin(t283);
t288 = sin(qJ(1));
t387 = t277 * t288;
t278 = cos(t283);
t287 = sin(qJ(4));
t289 = cos(qJ(4));
t318 = Icges(5,5) * t289 - Icges(5,6) * t287;
t216 = -Icges(5,3) * t278 + t277 * t318;
t399 = Icges(5,4) * t289;
t320 = -Icges(5,2) * t287 + t399;
t218 = -Icges(5,6) * t278 + t277 * t320;
t400 = Icges(5,4) * t287;
t321 = Icges(5,1) * t289 - t400;
t220 = -Icges(5,5) * t278 + t277 * t321;
t290 = cos(qJ(1));
t376 = t289 * t290;
t380 = t288 * t287;
t249 = t278 * t380 + t376;
t375 = t290 * t287;
t379 = t288 * t289;
t250 = t278 * t379 - t375;
t116 = t216 * t387 - t218 * t249 + t220 * t250;
t169 = Icges(5,5) * t250 - Icges(5,6) * t249 + Icges(5,3) * t387;
t241 = Icges(5,4) * t250;
t172 = -Icges(5,2) * t249 + Icges(5,6) * t387 + t241;
t240 = Icges(5,4) * t249;
t176 = -Icges(5,1) * t250 - Icges(5,5) * t387 + t240;
t97 = t169 * t387 - t172 * t249 - t176 * t250;
t251 = -t278 * t375 + t379;
t252 = t278 * t376 + t380;
t386 = t277 * t290;
t171 = Icges(5,5) * t252 + Icges(5,6) * t251 + Icges(5,3) * t386;
t401 = Icges(5,4) * t252;
t174 = Icges(5,2) * t251 + Icges(5,6) * t386 + t401;
t242 = Icges(5,4) * t251;
t177 = Icges(5,1) * t252 + Icges(5,5) * t386 + t242;
t98 = t171 * t387 - t249 * t174 + t250 * t177;
t325 = t97 * t288 + t290 * t98;
t12 = t116 * t278 - t277 * t325;
t118 = t216 * t386 + t251 * t218 + t252 * t220;
t100 = t171 * t386 + t251 * t174 + t252 * t177;
t99 = t169 * t386 + t251 * t172 - t252 * t176;
t324 = t100 * t290 + t288 * t99;
t477 = -t118 * t278 + t277 * t324;
t49 = t98 * t288 - t290 * t97;
t327 = t252 * rSges(5,1) + t251 * rSges(5,2);
t180 = rSges(5,3) * t386 + t327;
t410 = rSges(5,1) * t289;
t326 = -rSges(5,2) * t287 + t410;
t222 = -rSges(5,3) * t278 + t277 * t326;
t141 = t278 * t180 + t222 * t386;
t455 = -t250 * rSges(5,1) + t249 * rSges(5,2);
t178 = rSges(5,3) * t387 - t455;
t473 = t178 * t278 + t222 * t387;
t102 = -t141 * t288 + t290 * t473;
t393 = t169 * t278;
t471 = t172 * t287 + t176 * t289;
t107 = t471 * t277 + t393;
t260 = pkin(3) * t277 - pkin(6) * t278;
t358 = t222 + t260;
t181 = t358 * t288;
t183 = t358 * t290;
t453 = t181 * t288 + t183 * t290;
t258 = rSges(4,1) * t277 + rSges(4,2) * t278;
t284 = t288 ^ 2;
t285 = t290 ^ 2;
t472 = t284 + t285;
t457 = t472 * t258;
t464 = -m(5) / 0.2e1;
t370 = t453 * t464 - m(4) * t457 / 0.2e1;
t430 = rSges(5,3) + pkin(6);
t165 = (-t430 * t278 + (pkin(3) + t326) * t277) * t288;
t383 = t278 * t290;
t271 = pkin(6) * t383;
t353 = t277 * rSges(5,2) * t375 + rSges(5,3) * t383;
t166 = t271 + (-pkin(3) - t410) * t386 + t353;
t244 = t258 * t288;
t245 = t258 * t290;
t185 = t288 * t244 + t245 * t290;
t371 = m(5) * (t288 * t165 - t166 * t290) / 0.2e1 + m(4) * t185 / 0.2e1;
t43 = t371 - t370;
t475 = t43 * qJD(1);
t474 = t100 * t288 - t290 * t99;
t416 = -pkin(5) - qJ(2);
t275 = t288 * t416;
t274 = cos(pkin(7)) * pkin(2) + pkin(1);
t418 = pkin(3) * t278;
t451 = t430 * t277 + t274 + t418;
t149 = t451 * t290 - t275 + t327;
t338 = t290 * t416;
t465 = -t451 * t288 - t338 + t455;
t105 = t149 * t288 + t290 * t465;
t384 = t278 * t288;
t226 = Icges(4,4) * t384 - Icges(4,2) * t387 - Icges(4,6) * t290;
t273 = Icges(4,4) * t278;
t397 = Icges(4,2) * t277;
t227 = Icges(4,6) * t288 + (t273 - t397) * t290;
t402 = Icges(4,4) * t277;
t257 = Icges(4,1) * t278 - t402;
t229 = Icges(4,5) * t288 + t257 * t290;
t208 = t229 * t384;
t253 = Icges(4,5) * t278 - Icges(4,6) * t277;
t225 = Icges(4,3) * t288 + t253 * t290;
t331 = t225 * t290 - t208;
t224 = Icges(4,5) * t384 - Icges(4,6) * t387 - Icges(4,3) * t290;
t267 = Icges(4,4) * t387;
t228 = Icges(4,1) * t384 - Icges(4,5) * t290 - t267;
t362 = -t288 * t224 - t228 * t383;
t466 = -t226 * t386 - t227 * t387 - t331 - t362;
t435 = t288 / 0.2e1;
t433 = -t290 / 0.2e1;
t431 = t290 / 0.2e1;
t459 = t222 * t288;
t217 = Icges(5,3) * t277 + t278 * t318;
t378 = t289 * t220;
t382 = t287 * t218;
t403 = Icges(4,1) * t277;
t458 = t278 * (-t378 / 0.2e1 + t382 / 0.2e1 - t273 - t403 / 0.2e1 + t397 / 0.2e1 + t217 / 0.2e1);
t234 = (-Icges(5,2) * t289 - t400) * t277;
t237 = (-Icges(5,1) * t287 - t399) * t277;
t450 = -t287 * (t220 / 0.2e1 + t234 / 0.2e1) + (t237 / 0.2e1 - t218 / 0.2e1) * t289;
t235 = -Icges(4,2) * t384 - t267;
t254 = Icges(4,2) * t278 + t402;
t236 = t254 * t290;
t322 = -t273 - t403;
t238 = t322 * t288;
t239 = t322 * t290;
t449 = (t290 * (t226 - t238) + (-t227 + t239) * t288) * t278 + ((t228 + t235) * t290 + (-t229 + t236) * t288) * t277;
t448 = 0.4e1 * qJD(1);
t447 = -t477 / 0.2e1;
t446 = t49 / 0.2e1;
t445 = t474 / 0.2e1;
t188 = -Icges(5,5) * t249 - Icges(5,6) * t250;
t365 = -Icges(5,2) * t250 - t176 - t240;
t367 = -Icges(5,1) * t249 - t172 - t241;
t81 = -t188 * t278 + (-t365 * t287 + t367 * t289) * t277;
t444 = t81 / 0.2e1;
t223 = rSges(5,3) * t277 + t278 * t326;
t109 = (t223 * t288 - t178) * t277;
t207 = -rSges(5,1) * t277 * t376 + t353;
t110 = (-t222 * t290 - t207) * t278 + (-t223 * t290 + t180) * t277;
t442 = m(5) * (t109 * t465 + t110 * t149 - t141 * t166 + t165 * t473);
t196 = -rSges(5,1) * t249 - rSges(5,2) * t250;
t197 = rSges(5,1) * t251 - rSges(5,2) * t252;
t243 = (-rSges(5,1) * t287 - rSges(5,2) * t289) * t277;
t440 = m(5) * (-t105 * t243 - t181 * t197 + t183 * t196);
t438 = m(5) * (t149 * t166 + t165 * t465);
t437 = t277 / 0.2e1;
t436 = -t278 / 0.2e1;
t434 = t288 / 0.4e1;
t432 = -t290 / 0.4e1;
t429 = m(3) * t472 * (rSges(3,3) + qJ(2));
t411 = rSges(4,1) * t278;
t335 = t274 + t411;
t352 = rSges(4,2) * t387 + t290 * rSges(4,3);
t194 = -t288 * t335 - t338 + t352;
t334 = -rSges(4,2) * t386 + t288 * rSges(4,3);
t195 = t290 * t335 - t275 + t334;
t428 = m(4) * (t194 * t244 - t195 * t245);
t427 = m(4) * (t194 * t290 + t195 * t288);
t424 = m(5) * t102;
t423 = m(5) * t105;
t142 = t288 * t196 + t197 * t290;
t419 = m(5) * t142;
t417 = qJD(3) / 0.2e1;
t407 = t288 * t12;
t406 = t290 * t477;
t392 = t171 * t278;
t390 = t216 * t278;
t388 = t226 * t277;
t231 = (-Icges(5,5) * t287 - Icges(5,6) * t289) * t277;
t385 = t278 * t231;
t219 = Icges(5,6) * t277 + t278 * t320;
t381 = t287 * t219;
t221 = Icges(5,5) * t277 + t278 * t321;
t377 = t289 * t221;
t65 = t109 * t288 - t110 * t290;
t374 = t65 * qJD(2);
t373 = t65 * qJD(4);
t366 = Icges(5,1) * t251 - t174 - t401;
t364 = -Icges(5,2) * t252 + t177 + t242;
t361 = t288 * t225 + t229 * t383;
t360 = -t218 + t237;
t359 = t220 + t234;
t261 = pkin(6) * t277 + t418;
t357 = -t223 - t261;
t350 = qJD(1) * t277;
t349 = qJD(4) * t277;
t346 = t447 + t477 / 0.2e1;
t345 = t387 / 0.4e1;
t330 = t227 * t277 - t224;
t189 = Icges(5,5) * t251 - Icges(5,6) * t252;
t82 = -t189 * t278 + (-t364 * t287 + t366 * t289) * t277;
t93 = t231 * t387 - t359 * t249 + t360 * t250;
t94 = t231 * t386 + t359 * t251 + t360 * t252;
t329 = t440 / 0.2e1 + (t82 + t94) * t434 + (t81 + t93) * t432;
t319 = -Icges(4,5) * t277 - Icges(4,6) * t278;
t315 = -t174 * t287 + t177 * t289;
t108 = t277 * t315 - t392;
t317 = -t107 * t288 + t108 * t290;
t314 = t178 * t290 - t180 * t288;
t313 = t378 - t382;
t311 = -m(5) * (t149 * t197 - t196 * t465) + t385 / 0.2e1;
t135 = -t227 * t386 + t361;
t310 = (t290 * t330 + t135 - t361) * t431 + (-t288 * (-t228 * t278 + t388) - t224 * t290) * t433 - t49 / 0.2e1 + t446 + (t288 * t330 + t331 + t466) * t435;
t309 = t135 * t435 - t361 * t288 / 0.2e1 + t445 - t474 / 0.2e1 + (-t208 + (t225 + t388) * t290 + t362 + t466) * t433;
t72 = t188 * t387 - t365 * t249 + t367 * t250;
t73 = t189 * t387 - t364 * t249 + t366 * t250;
t32 = t73 * t288 - t290 * t72;
t74 = t188 * t386 + t365 * t251 + t367 * t252;
t75 = t189 * t386 + t364 * t251 + t366 * t252;
t33 = t75 * t288 - t290 * t74;
t308 = t32 * t433 + t33 * t435;
t305 = -t216 * t288 + t471;
t304 = -t216 * t290 - t315;
t303 = t217 - t313;
t301 = t12 * t434 + t477 * t432 - t407 / 0.4e1 + t406 / 0.4e1 + (t345 - t387 / 0.4e1) * t474;
t128 = t277 * t313 - t390;
t201 = t218 * t288;
t203 = t220 * t288;
t76 = -t305 * t278 + (t201 * t287 - t203 * t289 + t169) * t277;
t202 = t218 * t290;
t204 = t220 * t290;
t77 = -t304 * t278 + (t202 * t287 - t204 * t289 + t171) * t277;
t294 = t277 * t303 + t390;
t86 = -t219 * t249 + t221 * t250 + t288 * t294;
t87 = t251 * t219 + t252 * t221 + t290 * t294;
t89 = -t303 * t278 + (t216 + t377 - t381) * t277;
t297 = t128 * t437 + t442 / 0.2e1 + t89 * t436 + (t76 + t86) * t345 + (t77 + t87) * t386 / 0.4e1 + (-t107 + t116) * t384 / 0.4e1 + (t108 + t118) * t383 / 0.4e1;
t296 = t277 * t305 + t393;
t295 = t277 * t304 + t392;
t292 = -t377 / 0.2e1 + t381 / 0.2e1 - t257 / 0.2e1 + t254 / 0.2e1 - t216 / 0.2e1;
t259 = -rSges(4,2) * t277 + t411;
t233 = t319 * t290;
t232 = t319 * t288;
t184 = t357 * t290;
t182 = t357 * t288;
t151 = -t278 * t197 - t243 * t386;
t150 = t196 * t278 + t243 * t387;
t138 = -t419 / 0.2e1;
t131 = (t196 * t290 - t197 * t288) * t277;
t123 = t314 * t277;
t119 = (-pkin(3) * t386 + t207 + t271) * t290 + (-t260 * t288 - t459) * t288;
t112 = (t261 * t288 + t178) * t288 + (t261 * t290 + t180) * t290;
t111 = -t385 + (-t359 * t287 + t360 * t289) * t277;
t101 = t424 / 0.2e1;
t90 = t314 * t278 + (-t207 * t288 - t290 * t459) * t277;
t69 = -t251 * t202 - t252 * t204 + t290 * t295;
t68 = -t251 * t201 - t252 * t203 + t290 * t296;
t67 = t202 * t249 - t204 * t250 + t288 * t295;
t66 = t201 * t249 - t203 * t250 + t288 * t296;
t64 = t423 + t427 + t429;
t63 = m(5) * t65 * t417;
t57 = t112 * t142 + t453 * t243;
t46 = t450 * t277 - t311;
t44 = t370 + t371;
t42 = -t128 * t278 + t277 * t317;
t31 = t101 + t419 / 0.2e1;
t30 = t138 + t101;
t29 = t138 - t424 / 0.2e1;
t28 = t69 * t288 - t290 * t68;
t27 = t67 * t288 - t290 * t66;
t26 = -t277 * t292 + t428 + t438 - t458;
t23 = -t94 * t278 + (t288 * t74 + t290 * t75) * t277;
t22 = -t93 * t278 + (t288 * t72 + t290 * t73) * t277;
t21 = t109 * t473 - t110 * t141 + t123 * t90;
t14 = (t317 - t89) * t278 + (t76 * t288 + t77 * t290 + t128) * t277;
t9 = (t324 - t87) * t278 + (t288 * t68 + t290 * t69 + t118) * t277;
t8 = (t325 - t86) * t278 + (t288 * t66 + t290 * t67 + t116) * t277;
t7 = m(5) * t57 + t308;
t6 = t346 * t387;
t5 = t288 * t310 + t290 * t309;
t4 = m(5) * t21 + (t406 / 0.2e1 - t407 / 0.2e1 - t14 / 0.2e1) * t278 + (t9 * t431 + t8 * t435 + t42 / 0.2e1) * t277;
t3 = t301 + (-t94 / 0.4e1 - t82 / 0.4e1) * t288 + (t93 / 0.4e1 + t81 / 0.4e1) * t290 - t440 / 0.2e1 + t297;
t2 = t301 - t442 / 0.2e1 + (t89 / 0.2e1 + (-t118 / 0.4e1 - t108 / 0.4e1) * t290 + (-t116 / 0.4e1 + t107 / 0.4e1) * t288) * t278 + (-t128 / 0.2e1 + (-t87 / 0.4e1 - t77 / 0.4e1) * t290 + (-t86 / 0.4e1 - t76 / 0.4e1) * t288) * t277 + t329;
t1 = t297 + t329;
t10 = [t64 * qJD(2) + t26 * qJD(3) + t46 * qJD(4), qJD(1) * t64 + qJD(3) * t44 + qJD(4) * t30, t26 * qJD(1) + t44 * qJD(2) + t1 * qJD(4) + (m(5) * (t149 * t182 - t165 * t183 - t166 * t181 + t184 * t465) + (m(4) * (-t194 * t259 - t244 * t258) - t76 / 0.2e1 - t86 / 0.2e1 + t253 * t431 + (-t228 / 0.2e1 - t235 / 0.2e1) * t278 + (t226 / 0.2e1 - t238 / 0.2e1) * t277 - t309) * t290 + (m(4) * (-t195 * t259 + t245 * t258) + t77 / 0.2e1 + t87 / 0.2e1 + t253 * t435 + (t229 / 0.2e1 - t236 / 0.2e1) * t278 + (-t227 / 0.2e1 + t239 / 0.2e1) * t277 - t310) * t288) * qJD(3), t46 * qJD(1) + t30 * qJD(2) + t1 * qJD(3) + (-t111 * t278 + m(5) * (-t141 * t197 + t149 * t151 + t150 * t465 - t196 * t473)) * qJD(4) + ((t82 / 0.2e1 + t94 / 0.2e1) * t290 + (t444 + t93 / 0.2e1 - t346) * t288) * t349; t43 * qJD(3) + t29 * qJD(4) + (-t429 / 0.4e1 - t427 / 0.4e1 - t423 / 0.4e1) * t448, 0, t475 + 0.2e1 * ((-t182 * t290 + t184 * t288) * t417 + t373 / 0.4e1) * m(5), t29 * qJD(1) + t63 + m(5) * (t150 * t288 - t151 * t290) * qJD(4); -t43 * qJD(2) + t5 * qJD(3) + t2 * qJD(4) + (-t428 / 0.4e1 - t438 / 0.4e1) * t448 + t292 * t350 + t458 * qJD(1), t373 * t464 - t475, t5 * qJD(1) + (m(4) * (t259 * t457 - (t288 * (rSges(4,1) * t384 - t352) + t290 * (rSges(4,1) * t383 + t334)) * t185) + m(5) * (t112 * t119 - t181 * t182 - t183 * t184) + (t284 * t233 + (-t288 * t232 + t449) * t290 + t28) * t435 + (t285 * t232 + (-t290 * t233 + t449) * t288 + t27) * t433) * qJD(3) + t7 * qJD(4), t2 * qJD(1) + t7 * qJD(3) + t374 * t464 + (-t42 / 0.2e1 + (t33 / 0.2e1 - t9 / 0.2e1) * t290 + (t32 / 0.2e1 - t8 / 0.2e1) * t288) * t349 + (t22 * t433 + t23 * t435 + (t14 / 0.2e1 + (t444 + t447) * t290 + (-t82 / 0.2e1 + t12 / 0.2e1) * t288) * t278 + (-t102 * t243 + t131 * t112 + t123 * t142 - t150 * t183 - t151 * t181 - t21) * m(5)) * qJD(4); t311 * qJD(1) + t31 * qJD(2) + t3 * qJD(3) + t6 * qJD(4) - t450 * t350, qJD(1) * t31 + t63, t3 * qJD(1) + (t9 * t435 + t383 * t445 + t28 * t386 / 0.2e1 + t8 * t433 + t384 * t446 + t27 * t387 / 0.2e1 + (t107 * t290 + t108 * t288) * t437 + (t77 * t288 - t76 * t290) * t436 - t308) * qJD(3) + t4 * qJD(4) + (t374 / 0.2e1 + (-t109 * t183 - t110 * t181 + t112 * t90 + t119 * t123 - t141 * t182 + t184 * t473 - t57) * qJD(3)) * m(5), t6 * qJD(1) + t4 * qJD(3) + (m(5) * (t123 * t131 - t141 * t151 + t150 * t473) + t278 ^ 2 * t111 / 0.2e1 + (t23 * t431 + t22 * t435 + (t288 * t81 + t290 * t82) * t436) * t277) * qJD(4);];
Cq = t10;
