% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP1
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:49
% EndTime: 2019-12-05 16:39:59
% DurationCPUTime: 4.34s
% Computational Cost: add. (18525->259), mult. (12179->338), div. (0->0), fcn. (10740->6), ass. (0->173)
t269 = pkin(8) + qJ(2);
t266 = qJ(3) + t269;
t261 = sin(t266);
t262 = cos(t266);
t270 = sin(qJ(4));
t271 = cos(qJ(4));
t419 = (rSges(6,1) + pkin(4)) * t270 + rSges(6,2) * t271;
t426 = t419 * t262;
t427 = t419 * t261;
t430 = m(6) * (t261 * t427 + t262 * t426);
t118 = -t430 / 0.2e1;
t106 = t430 / 0.2e1;
t429 = Icges(5,2) + Icges(6,2);
t349 = Icges(6,4) * t270;
t350 = Icges(5,4) * t270;
t422 = t429 * t271 + t349 + t350;
t267 = Icges(6,4) * t271;
t268 = Icges(5,4) * t271;
t425 = t267 + t268 + (Icges(5,1) + Icges(6,1)) * t270;
t339 = t261 * t271;
t340 = t261 * t270;
t172 = Icges(6,4) * t339 - Icges(6,2) * t340 - Icges(6,6) * t262;
t174 = Icges(5,4) * t339 - Icges(5,2) * t340 - Icges(5,6) * t262;
t424 = t172 + t174;
t227 = Icges(6,4) * t340;
t176 = Icges(6,1) * t339 - Icges(6,5) * t262 - t227;
t228 = Icges(5,4) * t340;
t178 = Icges(5,1) * t339 - Icges(5,5) * t262 - t228;
t423 = t176 + t178;
t240 = Icges(6,1) * t271 - t349;
t242 = Icges(5,1) * t271 - t350;
t421 = t240 + t242;
t236 = -Icges(6,2) * t270 + t267;
t238 = -Icges(5,2) * t270 + t268;
t420 = -t238 - t236 - t425;
t416 = (-Icges(5,6) - Icges(6,6)) * t271 + (-Icges(5,5) - Icges(6,5)) * t270;
t177 = Icges(6,5) * t261 + t240 * t262;
t179 = Icges(5,5) * t261 + t242 * t262;
t415 = -t422 * t262 + t177 + t179;
t414 = -t429 * t339 - t227 - t228 + t423;
t173 = Icges(6,6) * t261 + t236 * t262;
t175 = Icges(5,6) * t261 + t238 * t262;
t413 = -t425 * t262 - t173 - t175;
t412 = t425 * t261 + t424;
t363 = pkin(4) * t271;
t263 = pkin(3) + t363;
t362 = -qJ(5) - pkin(7);
t144 = -rSges(6,1) * t339 + rSges(6,2) * t340 - t261 * t263 + (rSges(6,3) - t362) * t262;
t354 = rSges(6,1) * t271;
t295 = t263 + t354;
t334 = t262 * t270;
t302 = -rSges(6,2) * t334 - t261 * t362;
t145 = t261 * rSges(6,3) + t295 * t262 + t302;
t65 = t144 * t427 - t145 * t426;
t258 = t262 * pkin(7);
t355 = rSges(5,1) * t271;
t298 = pkin(3) + t355;
t303 = rSges(5,2) * t340 + t262 * rSges(5,3);
t148 = -t298 * t261 + t258 + t303;
t232 = rSges(5,2) * t334;
t149 = -t232 + t298 * t262 + (rSges(5,3) + pkin(7)) * t261;
t246 = t270 * rSges(5,1) + rSges(5,2) * t271;
t200 = t246 * t261;
t201 = t246 * t262;
t69 = t148 * t200 - t149 * t201;
t411 = -m(5) * t69 - m(6) * t65;
t394 = m(5) / 0.2e1;
t393 = m(6) / 0.2e1;
t79 = t144 * t262 + t145 * t261;
t407 = m(6) * qJD(3) * t79;
t406 = (t421 - t422) * t271 + t420 * t270;
t365 = pkin(2) * cos(t269);
t366 = pkin(2) * sin(t269);
t371 = m(4) * (t365 * (-rSges(4,1) * t261 - rSges(4,2) * t262) + (rSges(4,1) * t262 - t261 * rSges(4,2)) * t366);
t139 = t144 - t366;
t140 = t145 + t365;
t70 = t139 * t262 + t140 * t261;
t404 = t70 * m(6) * qJD(2);
t403 = t416 * t261;
t402 = t416 * t262;
t259 = t261 ^ 2;
t260 = t262 ^ 2;
t301 = t259 + t260;
t401 = -t415 * t270 + t413 * t271;
t400 = t414 * t270 + t412 * t271;
t300 = qJD(2) + qJD(3);
t146 = t148 - t366;
t147 = t149 + t365;
t361 = ((-t140 + t145) * t426 + (t139 - t144) * t427) * t393 + ((-t147 + t149) * t262 + (t146 - t148) * t261) * t246 * t394;
t63 = t139 * t427 - t140 * t426;
t68 = t146 * t200 - t147 * t201;
t399 = (t65 + t63) * t393 + (t69 + t68) * t394;
t281 = -t420 * t271 / 0.2e1 + (-t422 / 0.2e1 + t421 / 0.2e1) * t270;
t397 = 0.4e1 * qJD(2);
t395 = 2 * qJD(4);
t62 = -t149 * t146 + t147 * t148;
t390 = m(5) * t62;
t388 = m(5) * t68;
t383 = m(6) * (t70 - t79);
t382 = m(6) * (t79 + t70);
t60 = -t145 * t139 + t140 * t144;
t381 = m(6) * t60;
t379 = m(6) * t63;
t376 = -t261 / 0.2e1;
t375 = t261 / 0.2e1;
t374 = -t262 / 0.2e1;
t356 = m(6) * qJD(5);
t233 = Icges(6,5) * t271 - Icges(6,6) * t270;
t336 = t262 * t233;
t234 = Icges(5,5) * t271 - Icges(5,6) * t270;
t335 = t262 * t234;
t333 = t262 * t271;
t332 = t270 * t172;
t331 = t270 * t173;
t330 = t270 * t174;
t329 = t270 * t175;
t168 = Icges(6,5) * t339 - Icges(6,6) * t340 - Icges(6,3) * t262;
t315 = -t261 * t168 - t176 * t333;
t169 = Icges(6,3) * t261 + t336;
t314 = t261 * t169 + t177 * t333;
t170 = Icges(5,5) * t339 - Icges(5,6) * t340 - Icges(5,3) * t262;
t313 = -t261 * t170 - t178 * t333;
t171 = Icges(5,3) * t261 + t335;
t312 = t261 * t171 + t179 * t333;
t292 = -t168 + t331;
t152 = t177 * t339;
t294 = t262 * t169 - t152;
t84 = -t262 * t332 - t315;
t85 = -t262 * t331 + t314;
t10 = (t292 * t262 - t314 + t85) * t262 + (t292 * t261 + t294 + t84) * t261;
t291 = -t170 + t329;
t153 = t179 * t339;
t293 = t262 * t171 - t153;
t86 = -t262 * t330 - t313;
t87 = -t262 * t329 + t312;
t11 = (t291 * t262 - t312 + t87) * t262 + (t291 * t261 + t293 + t86) * t261;
t81 = -t261 * t331 - t294;
t12 = (t81 - t152 + (t169 + t332) * t262 + t315) * t262 + t314 * t261;
t83 = -t261 * t329 - t293;
t13 = (t83 - t153 + (t171 + t330) * t262 + t313) * t262 + t312 * t261;
t48 = t261 * t81 - t262 * (-(-t176 * t271 + t332) * t261 - t262 * t168);
t49 = t261 * t83 - t262 * (-(-t178 * t271 + t330) * t261 - t262 * t170);
t50 = t261 * t85 - t262 * t84;
t51 = t261 * t87 - t262 * t86;
t2 = (t51 / 0.2e1 - t13 / 0.2e1 + t50 / 0.2e1 - t12 / 0.2e1) * t262 + (t11 / 0.2e1 + t49 / 0.2e1 + t10 / 0.2e1 + t48 / 0.2e1) * t261;
t59 = 0.2e1 * t118;
t299 = t2 * qJD(4) + t59 * qJD(5);
t296 = t270 * rSges(6,2) - t354 - t363;
t282 = (-t200 * t262 + t201 * t261) * t246;
t274 = t281 + t399;
t273 = -t281 + (t423 * t270 + t424 * t271) * (t375 + t376);
t58 = t118 + t106;
t272 = t58 * qJD(5) + ((t12 + t13) * t262 / 0.2e1 + (t10 + t11 + t48 + t49) * t376 + (t415 * t271 + t413 * t270 + t406 * t262 + (t233 + t234) * t261) * t375 + (t406 * t261 - t412 * t270 + t414 * t271 - t335 - t336 + t50 + t51) * t374) * qJD(4);
t248 = -t270 * rSges(5,2) + t355;
t185 = t296 * t262;
t183 = t296 * t261;
t141 = -t200 * t261 - t201 * t262;
t110 = t419 * t301;
t57 = 0.2e1 * t106;
t55 = t58 * qJD(4);
t53 = t57 * qJD(4);
t39 = t382 / 0.2e1;
t38 = t383 / 0.2e1;
t23 = t281 - t411;
t20 = t281 + t379 + t388;
t17 = t371 + t381 + t390;
t16 = t39 - t383 / 0.2e1;
t15 = t39 + t38;
t14 = t38 - t382 / 0.2e1;
t5 = t274 - t361;
t4 = t274 + t361;
t3 = t273 + t361 - t399;
t1 = [0, 0, 0, (-t110 * t393 + t141 * t394) * t395, 0; 0, t17 * qJD(3) + t20 * qJD(4) + t70 * t356, t17 * qJD(2) + t4 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t60 * t393 + t62 * t394 + t371 / 0.2e1) * qJD(3), t20 * qJD(2) + t4 * qJD(3) + (((-t146 * t262 - t147 * t261) * t248 + t282) * t394 + (t139 * t185 + t140 * t183) * t393) * t395 + t272, t15 * qJD(3) + t404 + t55; 0, t5 * qJD(4) + t16 * qJD(5) + (-t381 / 0.4e1 - t390 / 0.4e1 - t371 / 0.4e1) * t397, t23 * qJD(4) + t79 * t356, t5 * qJD(2) + t23 * qJD(3) + (((-t148 * t262 - t149 * t261) * t248 + t282) * t394 + (t144 * t185 + t145 * t183) * t393) * t395 + t272, t16 * qJD(2) + t407 + t55; 0, t273 * qJD(2) + t3 * qJD(3) + (-t379 / 0.4e1 - t388 / 0.4e1) * t397 + t299, t3 * qJD(2) + t299 + (t273 + t411) * qJD(3), (m(5) * (t246 * t248 * t301 + (t261 * (rSges(5,1) * t339 - t303) + t262 * (rSges(5,1) * t333 + t261 * rSges(5,3) - t232)) * t141) + m(6) * (-t110 * ((-t261 * pkin(3) - t144 + t258) * t261 + ((-pkin(3) + t295) * t262 + (rSges(6,3) - pkin(7)) * t261 + t302) * t262) - t427 * t183 - t426 * t185) + (t402 * t259 + (t400 * t262 + (t401 - t403) * t261) * t262) * t375 + (t403 * t260 + (t401 * t261 + (t400 - t402) * t262) * t261) * t374) * qJD(4) + t300 * t2, t300 * t59; 0, t14 * qJD(3) - t404 + t53, t14 * qJD(2) - t407 + t53, m(6) * (-t183 * t262 + t185 * t261) * qJD(4) + t300 * t57, 0;];
Cq = t1;
