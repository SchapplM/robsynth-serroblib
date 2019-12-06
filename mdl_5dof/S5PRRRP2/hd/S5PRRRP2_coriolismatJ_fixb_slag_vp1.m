% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP2
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:38
% EndTime: 2019-12-05 16:41:48
% DurationCPUTime: 4.33s
% Computational Cost: add. (21318->291), mult. (14532->376), div. (0->0), fcn. (12957->6), ass. (0->200)
t288 = pkin(8) + qJ(2);
t285 = qJ(3) + t288;
t281 = cos(t285);
t289 = sin(qJ(4));
t359 = t281 * t289;
t280 = sin(t285);
t290 = cos(qJ(4));
t391 = rSges(6,1) + pkin(4);
t317 = t391 * t289;
t378 = rSges(6,3) + qJ(5);
t452 = -t378 * t290 + t317;
t457 = t452 * t280;
t458 = t359 * t457;
t286 = Icges(6,5) * t289;
t307 = Icges(6,3) * t290 - t286;
t375 = Icges(5,4) * t289;
t455 = Icges(5,2) * t290 + t307 + t375;
t287 = Icges(5,4) * t290;
t260 = Icges(5,1) * t289 + t287;
t374 = Icges(6,5) * t290;
t456 = Icges(6,1) * t289 + t260 - t374;
t257 = -Icges(5,2) * t289 + t287;
t454 = t257 + t456;
t261 = Icges(5,1) * t290 - t375;
t429 = Icges(6,1) * t290 + t286;
t453 = t261 + t429;
t278 = t280 ^ 2;
t279 = t281 ^ 2;
t321 = t278 + t279;
t430 = t378 * t289 + t391 * t290;
t277 = t281 * pkin(7);
t425 = pkin(3) + t430;
t145 = t281 * rSges(6,2) - t280 * t425 + t277;
t387 = pkin(2) * sin(t288);
t140 = t145 - t387;
t450 = t140 - t145;
t448 = (-Icges(5,6) + Icges(6,6)) * t290 + (-Icges(6,4) - Icges(5,5)) * t289;
t196 = Icges(6,4) * t280 + t281 * t429;
t198 = Icges(5,5) * t280 + t261 * t281;
t447 = -t455 * t281 + t196 + t198;
t195 = -Icges(6,4) * t281 + t280 * t429;
t365 = t280 * t289;
t248 = Icges(5,4) * t365;
t364 = t280 * t290;
t197 = Icges(5,1) * t364 - Icges(5,5) * t281 - t248;
t446 = -Icges(5,2) * t364 - t307 * t280 + t195 + t197 - t248;
t358 = t281 * t290;
t247 = Icges(6,5) * t358;
t188 = Icges(6,6) * t280 + Icges(6,3) * t359 + t247;
t194 = Icges(5,6) * t280 + t257 * t281;
t445 = -Icges(6,1) * t359 - t260 * t281 + t188 - t194 + t247;
t253 = Icges(6,3) * t289 + t374;
t187 = -Icges(6,6) * t281 + t253 * t280;
t193 = Icges(5,4) * t364 - Icges(5,2) * t365 - Icges(5,6) * t281;
t444 = t456 * t280 - t187 + t193;
t146 = (rSges(6,2) + pkin(7)) * t280 + t425 * t281;
t161 = -t281 * t317 + t378 * t358;
t67 = t145 * t457 + t146 * t161;
t379 = rSges(5,1) * t290;
t316 = pkin(3) + t379;
t324 = rSges(5,2) * t365 + t281 * rSges(5,3);
t158 = -t280 * t316 + t277 + t324;
t250 = rSges(5,2) * t359;
t159 = -t250 + t316 * t281 + (rSges(5,3) + pkin(7)) * t280;
t264 = t289 * rSges(5,1) + rSges(5,2) * t290;
t218 = t264 * t280;
t220 = t264 * t281;
t84 = t158 * t218 - t159 * t220;
t443 = -m(5) * t84 - m(6) * t67;
t420 = m(5) / 0.2e1;
t419 = m(6) / 0.2e1;
t442 = -t289 / 0.2e1;
t386 = pkin(2) * cos(t288);
t389 = m(4) * (t386 * (-rSges(4,1) * t280 - rSges(4,2) * t281) + (rSges(4,1) * t281 - t280 * rSges(4,2)) * t387);
t318 = t146 * t359;
t80 = -t145 * t365 + t318;
t439 = m(6) * qJD(3) * t80;
t255 = Icges(6,4) * t290 + Icges(6,6) * t289;
t366 = t280 * t255;
t191 = -Icges(6,2) * t281 + t366;
t178 = t280 * t191;
t357 = t289 * t187;
t100 = t195 * t358 + t281 * t357 + t178;
t438 = t100 * t281;
t437 = (t453 - t455) * t290 + (t253 - t454) * t289;
t434 = t280 * (t195 * t290 + t357);
t156 = t158 - t387;
t157 = t159 + t386;
t75 = -t159 * t156 + t157 * t158;
t432 = t448 * t280;
t431 = t448 * t281;
t428 = -t447 * t289 + t445 * t290;
t427 = t446 * t289 + t444 * t290;
t141 = t146 + t386;
t184 = t452 * t281;
t385 = ((-t141 + t146) * t184 + t450 * t457) * t419 + ((-t157 + t159) * t281 + (t156 - t158) * t280) * t264 * t420;
t63 = t140 * t457 + t141 * t161;
t81 = t156 * t218 - t157 * t220;
t426 = (t67 + t63) * t419 + (t84 + t81) * t420;
t300 = t455 * t442 + t453 * t289 / 0.2e1 + (-t253 / 0.2e1 + t454 / 0.2e1) * t290;
t424 = 4 * qJD(2);
t422 = 2 * qJD(4);
t416 = m(5) * t75;
t414 = m(5) * t81;
t123 = t141 * t359;
t408 = m(6) * (-t450 * t365 + t123 - t318);
t407 = m(6) * (t123 + t318 + (-t140 - t145) * t365);
t42 = -t146 * t140 + t141 * t145;
t406 = m(6) * t42;
t337 = t161 * t365 + t458;
t340 = t140 * t358 + t141 * t364;
t405 = m(6) * (t337 + t340);
t338 = t145 * t358 + t146 * t364;
t403 = m(6) * (t337 + t338);
t313 = t184 * t365 - t458;
t402 = m(6) * (t313 + t340);
t401 = m(6) * (t313 + t338);
t400 = m(6) * t63;
t221 = t321 * t289;
t396 = t221 / 0.2e1;
t395 = -t280 / 0.2e1;
t394 = t280 / 0.2e1;
t393 = -t281 / 0.2e1;
t360 = t281 * t255;
t192 = Icges(6,2) * t280 + t360;
t356 = t289 * t188;
t101 = t280 * t192 + t196 * t358 + t281 * t356;
t362 = t281 * t191;
t302 = t101 + t362;
t312 = t281 * t192 - t196 * t364 - t280 * t356;
t16 = (t101 - t302) * t281 + (t100 + t312 - t178) * t280;
t189 = Icges(5,5) * t364 - Icges(5,6) * t365 - Icges(5,3) * t281;
t335 = -t280 * t189 - t197 * t358;
t355 = t289 * t193;
t102 = -t281 * t355 - t335;
t254 = Icges(5,5) * t290 - Icges(5,6) * t289;
t361 = t281 * t254;
t190 = Icges(5,3) * t280 + t361;
t334 = t280 * t190 + t198 * t358;
t354 = t289 * t194;
t103 = -t281 * t354 + t334;
t314 = -t189 + t354;
t169 = t198 * t364;
t315 = t281 * t190 - t169;
t17 = (t281 * t314 + t103 - t334) * t281 + (t280 * t314 + t102 + t315) * t280;
t96 = -t362 + t434;
t18 = -t438 + (t302 + t96 - t434) * t280;
t99 = -t280 * t354 - t315;
t19 = (t99 - t169 + (t190 + t355) * t281 + t335) * t281 + t334 * t280;
t56 = -t280 * t312 - t281 * t96;
t57 = t280 * t99 - t281 * (-t280 * (-t197 * t290 + t355) - t281 * t189);
t58 = t101 * t280 - t438;
t59 = -t102 * t281 + t103 * t280;
t2 = (t59 / 0.2e1 - t19 / 0.2e1 + t58 / 0.2e1 - t18 / 0.2e1) * t281 + (t17 / 0.2e1 + t57 / 0.2e1 + t16 / 0.2e1 + t56 / 0.2e1) * t280;
t390 = t2 * qJD(4);
t381 = m(6) * qJD(4);
t380 = m(6) * qJD(5);
t349 = t289 * t290;
t339 = (-t161 - t184) * t457;
t336 = -t184 * t358 - t364 * t457;
t325 = t321 * t349;
t202 = (t396 + t442) * m(6);
t320 = t202 * qJD(1);
t77 = -t140 * t365 + t123;
t319 = m(6) * t77 * qJD(2);
t301 = (-t218 * t281 + t220 * t280) * t264;
t293 = t300 + t426;
t292 = -t300 + ((t187 + t193) * t290 + (-t195 + t197) * t289) * (t394 + t395);
t291 = ((t18 + t19) * t281 / 0.2e1 + (t16 + t17 + t56 + t57) * t395 + (t280 * t254 + t437 * t281 + t445 * t289 + t447 * t290 + t366) * t394 + (t437 * t280 - t444 * t289 + t446 * t290 - t360 - t361 + t58 + t59) * t393) * qJD(4);
t267 = -t289 * rSges(5,2) + t379;
t201 = m(6) * t396 + t289 * t419;
t186 = t325 - t349;
t185 = t430 * t281;
t183 = t430 * t280;
t151 = -t218 * t280 - t220 * t281;
t106 = t161 * t281 - t278 * t452;
t88 = t321 * t430;
t62 = t221 * t88 + t336;
t60 = t401 / 0.2e1;
t54 = t402 / 0.2e1;
t47 = t403 / 0.2e1;
t43 = t405 / 0.2e1;
t37 = t407 / 0.2e1;
t36 = t408 / 0.2e1;
t27 = t300 - t443;
t26 = t300 + t400 + t414;
t23 = t389 + t406 + t416;
t22 = t60 - t403 / 0.2e1;
t21 = t60 + t47;
t20 = t47 - t401 / 0.2e1;
t11 = t54 - t405 / 0.2e1;
t10 = t54 + t43;
t9 = t43 - t402 / 0.2e1;
t8 = t37 - t408 / 0.2e1;
t7 = t37 + t36;
t6 = t36 - t407 / 0.2e1;
t5 = t293 - t385;
t4 = t293 + t385;
t3 = t292 + t385 - t426;
t1 = [0, 0, 0, (t106 * t419 + t151 * t420) * t422 + t201 * qJD(5), t201 * qJD(4); 0, t23 * qJD(3) + t26 * qJD(4) + t77 * t380, t23 * qJD(2) + t4 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t389 / 0.2e1 + t42 * t419 + t75 * t420) * qJD(3), t26 * qJD(2) + t4 * qJD(3) + t10 * qJD(5) + ((-t140 * t185 - t141 * t183 + t339) * t419 + ((-t156 * t281 - t157 * t280) * t267 + t301) * t420) * t422 + t291, t7 * qJD(3) + t10 * qJD(4) + t319; 0, t5 * qJD(4) + t8 * qJD(5) + (-t406 / 0.4e1 - t416 / 0.4e1 - t389 / 0.4e1) * t424, t27 * qJD(4) + t80 * t380, t5 * qJD(2) + t27 * qJD(3) + t21 * qJD(5) + ((-t145 * t185 - t146 * t183 + t339) * t419 + ((-t158 * t281 - t159 * t280) * t267 + t301) * t420) * t422 + t291, t8 * qJD(2) + t21 * qJD(4) + t439; qJD(5) * t202, t292 * qJD(2) + t3 * qJD(3) + t11 * qJD(5) + (-t400 / 0.4e1 - t414 / 0.4e1) * t424 + t390, t3 * qJD(2) + t22 * qJD(5) + t390 + (t292 + t443) * qJD(3), (m(5) * (t264 * t267 * t321 + (t280 * (rSges(5,1) * t364 - t324) + t281 * (rSges(5,1) * t358 + t280 * rSges(5,3) - t250)) * t151) + m(6) * (t106 * t88 + t183 * t457 + t184 * t185) + (t431 * t278 + (t427 * t281 + (t428 - t432) * t280) * t281) * t394 + (t432 * t279 + (t428 * t280 + (t427 - t431) * t281) * t280) * t393) * qJD(4) + t62 * t380 + (qJD(2) + qJD(3)) * t2, t320 + t11 * qJD(2) + t22 * qJD(3) + t62 * t381 + (-t221 * t290 - t186 + t325) * t380; -t202 * qJD(4), t6 * qJD(3) + t9 * qJD(4) - t319, t6 * qJD(2) + t20 * qJD(4) - t439, -t320 + t9 * qJD(2) + t20 * qJD(3) + (-t106 * t290 + (-t183 * t280 - t185 * t281 + t88) * t289 - t62 + t336) * t381 + t186 * t380, t186 * t381;];
Cq = t1;
