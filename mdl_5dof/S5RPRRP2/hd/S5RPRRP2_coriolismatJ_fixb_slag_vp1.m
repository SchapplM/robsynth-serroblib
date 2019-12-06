% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:24
% EndTime: 2019-12-05 18:01:34
% DurationCPUTime: 4.19s
% Computational Cost: add. (18685->282), mult. (12347->351), div. (0->0), fcn. (10902->8), ass. (0->185)
t438 = Icges(5,1) + Icges(6,1);
t437 = Icges(5,2) + Icges(6,2);
t275 = cos(qJ(4));
t269 = Icges(6,4) * t275;
t270 = Icges(5,4) * t275;
t273 = sin(qJ(4));
t425 = t438 * t273 + t269 + t270;
t271 = qJ(1) + pkin(8);
t268 = qJ(3) + t271;
t263 = sin(t268);
t264 = cos(t268);
t334 = t264 * t275;
t335 = t264 * t273;
t176 = Icges(6,4) * t334 - Icges(6,2) * t335 + Icges(6,6) * t263;
t178 = Icges(5,4) * t334 - Icges(5,2) * t335 + Icges(5,6) * t263;
t436 = t176 + t178;
t231 = Icges(6,4) * t335;
t180 = Icges(6,1) * t334 + Icges(6,5) * t263 - t231;
t233 = Icges(5,4) * t335;
t182 = Icges(5,1) * t334 + Icges(5,5) * t263 - t233;
t435 = t180 + t182;
t434 = (-Icges(5,4) - Icges(6,4)) * t273;
t433 = t437 * t275 - t434;
t432 = t438 * t275 + t434;
t431 = (Icges(5,6) + Icges(6,6)) * t275 + (Icges(5,5) + Icges(6,5)) * t273;
t430 = -t437 * t334 - t231 - t233 + t435;
t340 = t263 * t273;
t230 = Icges(6,4) * t340;
t339 = t263 * t275;
t179 = -Icges(6,1) * t339 + Icges(6,5) * t264 + t230;
t232 = Icges(5,4) * t340;
t181 = -Icges(5,1) * t339 + Icges(5,5) * t264 + t232;
t429 = t437 * t339 + t179 + t181 + t230 + t232;
t428 = t425 * t264 + t436;
t404 = Icges(6,2) * t273 - t269;
t175 = Icges(6,6) * t264 + t404 * t263;
t403 = Icges(5,2) * t273 - t270;
t177 = Icges(5,6) * t264 + t403 * t263;
t427 = t425 * t263 - t175 - t177;
t426 = t404 + t403;
t364 = pkin(4) * t275;
t265 = pkin(3) + t364;
t281 = rSges(6,1) * t334 - rSges(6,2) * t335 + t264 * t265;
t410 = -rSges(6,3) - qJ(5) - pkin(7);
t147 = t410 * t263 - t281;
t293 = -cos(qJ(1)) * pkin(1) - pkin(2) * cos(t271);
t135 = t147 + t293;
t424 = t135 - t147;
t291 = -rSges(6,2) * t340 + t410 * t264;
t357 = rSges(6,1) * t275;
t297 = t265 + t357;
t146 = t297 * t263 + t291;
t240 = pkin(4) * t340;
t307 = rSges(6,1) * t340 + rSges(6,2) * t339;
t169 = t240 + t307;
t356 = rSges(6,2) * t275;
t170 = (t356 + (rSges(6,1) + pkin(4)) * t273) * t264;
t65 = -t146 * t169 + t147 * t170;
t239 = rSges(5,2) * t335;
t358 = rSges(5,1) * t275;
t299 = pkin(3) + t358;
t149 = t239 - t299 * t264 + (-rSges(5,3) - pkin(7)) * t263;
t252 = rSges(5,1) * t273 + rSges(5,2) * t275;
t204 = t252 * t264;
t260 = t264 * pkin(7);
t308 = -rSges(5,2) * t340 - t264 * rSges(5,3);
t148 = t299 * t263 - t260 + t308;
t203 = t252 * t263;
t347 = t148 * t203;
t70 = t149 * t204 - t347;
t423 = -m(5) * t70 - m(6) * t65;
t391 = m(5) / 0.2e1;
t390 = m(6) / 0.2e1;
t336 = t264 * t146;
t79 = -t147 * t263 - t336;
t418 = t79 * m(6) * qJD(3);
t417 = t431 * t263;
t416 = t431 * t264;
t414 = -t429 * t273 + t427 * t275;
t413 = t430 * t273 + t428 * t275;
t412 = (-t432 + t433) * t275 + (t425 - t426) * t273;
t301 = qJD(1) + qJD(3);
t294 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t271);
t370 = m(4) * (t293 * (rSges(4,1) * t263 + rSges(4,2) * t264) - (-t264 * rSges(4,1) + t263 * rSges(4,2)) * t294);
t134 = t146 + t294;
t124 = t264 * t134;
t69 = -t135 * t263 - t124;
t411 = t69 * m(6) * qJD(1);
t255 = -rSges(5,2) * t273 + t358;
t409 = t255 * t391;
t408 = (t178 * t273 - t182 * t275) * t264;
t407 = (t176 * t273 - t180 * t275) * t264;
t145 = t149 + t293;
t251 = rSges(6,1) * t273 + t356;
t185 = t251 * t263 + t240;
t365 = pkin(4) * t273;
t187 = (t251 + t365) * t264;
t144 = t148 + t294;
t349 = t144 * t203;
t363 = (t424 * t187 + (-t134 + t146) * t185) * t390 + (-t349 + t347 + (t145 - t149) * t204) * t391;
t63 = -t134 * t169 + t135 * t170;
t68 = t145 * t204 - t349;
t402 = (t65 + t63) * t390 + (t70 + t68) * t391;
t280 = (-t426 / 0.2e1 + t425 / 0.2e1) * t275 + (-t433 / 0.2e1 + t432 / 0.2e1) * t273;
t261 = t263 ^ 2;
t262 = t264 ^ 2;
t395 = 0.4e1 * qJD(1);
t392 = 2 * qJD(4);
t62 = -t149 * t144 + t145 * t148;
t387 = m(5) * t62;
t385 = m(5) * t68;
t381 = m(6) * (-t263 * t424 - t124 + t336);
t380 = m(6) * (-t124 - t336 + (-t135 - t147) * t263);
t57 = -t147 * t134 + t135 * t146;
t379 = m(6) * t57;
t377 = m(6) * t63;
t373 = t263 / 0.2e1;
t372 = -t264 / 0.2e1;
t371 = t264 / 0.2e1;
t368 = m(6) * (t169 * t263 + t170 * t264);
t366 = m(6) * (-t185 * t263 - t264 * t187);
t359 = m(6) * qJD(5);
t345 = t175 * t273;
t344 = t177 * t273;
t241 = Icges(6,5) * t275 - Icges(6,6) * t273;
t342 = t241 * t263;
t242 = Icges(5,5) * t275 - Icges(5,6) * t273;
t341 = t242 * t263;
t321 = -t187 * t169 + t185 * t170;
t171 = Icges(6,3) * t264 - t342;
t320 = t264 * t171 + t175 * t340;
t173 = Icges(5,3) * t264 - t341;
t319 = t264 * t173 + t177 * t340;
t318 = t263 * t171 + t179 * t334;
t317 = t263 * t173 + t181 * t334;
t302 = t261 + t262;
t172 = Icges(6,5) * t334 - Icges(6,6) * t335 + Icges(6,3) * t263;
t296 = -t179 * t275 - t172;
t81 = t264 * t172 + t176 * t340 - t180 * t339;
t84 = -t175 * t335 + t318;
t85 = t172 * t263 - t407;
t10 = (t320 + t85 + t407) * t264 + (-t84 + (t296 - t345) * t264 + t81 + t318) * t263;
t174 = Icges(5,5) * t334 - Icges(5,6) * t335 + Icges(5,3) * t263;
t295 = -t181 * t275 - t174;
t83 = t264 * t174 + t178 * t340 - t182 * t339;
t86 = -t177 * t335 + t317;
t87 = t174 * t263 - t408;
t11 = (t319 + t87 + t408) * t264 + (-t86 + (t295 - t344) * t264 + t83 + t317) * t263;
t80 = -t179 * t339 + t320;
t12 = (t81 + (-t172 + t345) * t264 - t318) * t264 + (t263 * t296 + t320 - t80) * t263;
t82 = -t181 * t339 + t319;
t13 = (t83 + (-t174 + t344) * t264 - t317) * t264 + (t263 * t295 + t319 - t82) * t263;
t48 = t263 * t81 + t264 * t80;
t49 = t263 * t83 + t264 * t82;
t50 = t263 * t85 + t264 * t84;
t51 = t263 * t87 + t264 * t86;
t2 = (t13 / 0.2e1 + t51 / 0.2e1 + t12 / 0.2e1 + t50 / 0.2e1) * t264 + (-t49 / 0.2e1 + t11 / 0.2e1 - t48 / 0.2e1 + t10 / 0.2e1) * t263;
t120 = t366 / 0.2e1;
t60 = t120 - t368 / 0.2e1;
t300 = t2 * qJD(4) + t60 * qJD(5);
t298 = -rSges(6,2) * t273 + t357 + t364;
t279 = t280 + t402;
t278 = -t280 + (t435 * t273 + t436 * t275) * (t371 + t372);
t107 = t368 / 0.2e1;
t59 = t120 + t107;
t277 = t59 * qJD(5) + (-(t10 + t11) * t263 / 0.2e1 + (t12 + t13 + t50 + t51) * t372 + (t429 * t275 + t427 * t273 + (t241 + t242) * t264 + t412 * t263) * t371 + (-t412 * t264 - t428 * t273 + t430 * t275 + t341 + t342 + t48 + t49) * t373) * qJD(4);
t188 = t298 * t264;
t186 = t298 * t263;
t141 = -t203 * t263 - t204 * t264;
t110 = -t251 * t262 - t263 * t307 - t302 * t365;
t58 = t107 - t366 / 0.2e1;
t55 = t59 * qJD(4);
t53 = t58 * qJD(4);
t35 = t380 / 0.2e1;
t34 = t381 / 0.2e1;
t23 = t280 - t423;
t20 = t280 + t377 + t385;
t17 = t370 + t379 + t387;
t16 = t35 - t381 / 0.2e1;
t15 = t35 + t34;
t14 = t34 - t380 / 0.2e1;
t5 = t279 + t363;
t4 = t279 - t363;
t3 = t278 + t363 - t402;
t1 = [t17 * qJD(3) + t20 * qJD(4) + t69 * t359, 0, t17 * qJD(1) + t5 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t370 / 0.2e1 + t62 * t391 + t57 * t390) * qJD(3), t20 * qJD(1) + t5 * qJD(3) + ((t144 * t264 + t145 * t263) * t409 + (t134 * t188 + t135 * t186 + t321) * t390) * t392 + t277, t15 * qJD(3) + t411 + t55; 0, 0, 0, (t110 * t390 + t141 * t391) * t392, 0; t4 * qJD(4) + t16 * qJD(5) + (-t370 / 0.4e1 - t387 / 0.4e1 - t379 / 0.4e1) * t395, 0, t23 * qJD(4) + t79 * t359, t4 * qJD(1) + t23 * qJD(3) + ((t146 * t188 + t147 * t186 + t321) * t390 + (t148 * t264 + t149 * t263) * t409) * t392 + t277, t16 * qJD(1) + t418 + t55; t278 * qJD(1) + t3 * qJD(3) + (-t385 / 0.4e1 - t377 / 0.4e1) * t395 + t300, 0, t3 * qJD(1) + t300 + (t278 + t423) * qJD(3), (m(5) * (t252 * t255 * t302 + (t264 * (rSges(5,1) * t334 + t263 * rSges(5,3) - t239) - t263 * (-rSges(5,1) * t339 - t308)) * t141) + m(6) * (t110 * ((-t264 * pkin(3) + t281) * t264 + ((-pkin(7) - t410) * t264 + t260 + t291 + (-pkin(3) + t297) * t263) * t263) + t185 * t186 + t187 * t188) + (-t416 * t261 + (t414 * t264 + (-t413 + t417) * t263) * t264) * t373 + (t417 * t262 + (t413 * t263 + (-t414 - t416) * t264) * t263) * t371) * qJD(4) + t301 * t2, t301 * t60; t14 * qJD(3) - t411 + t53, 0, t14 * qJD(1) - t418 + t53, m(6) * (t186 * t264 - t188 * t263) * qJD(4) + t301 * t58, 0;];
Cq = t1;
