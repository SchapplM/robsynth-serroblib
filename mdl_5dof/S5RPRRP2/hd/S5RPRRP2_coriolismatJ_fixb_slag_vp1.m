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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:27:59
% DurationCPUTime: 4.17s
% Computational Cost: add. (18685->261), mult. (12347->340), div. (0->0), fcn. (10902->8), ass. (0->171)
t271 = qJ(1) + pkin(8);
t268 = qJ(3) + t271;
t263 = sin(t268);
t264 = cos(t268);
t272 = sin(qJ(4));
t274 = cos(qJ(4));
t422 = rSges(6,2) * t274 + (rSges(6,1) + pkin(4)) * t272;
t429 = t422 * t264;
t430 = t422 * t263;
t433 = m(6) * (t263 * t430 + t264 * t429);
t120 = -t433 / 0.2e1;
t107 = t433 / 0.2e1;
t432 = Icges(5,2) + Icges(6,2);
t354 = Icges(6,4) * t272;
t355 = Icges(5,4) * t272;
t425 = t274 * t432 + t354 + t355;
t269 = Icges(6,4) * t274;
t270 = Icges(5,4) * t274;
t428 = t269 + t270 + (Icges(5,1) + Icges(6,1)) * t272;
t338 = t263 * t274;
t339 = t263 * t272;
t172 = Icges(6,4) * t338 - Icges(6,2) * t339 - Icges(6,6) * t264;
t174 = Icges(5,4) * t338 - Icges(5,2) * t339 - Icges(5,6) * t264;
t427 = t172 + t174;
t227 = Icges(6,4) * t339;
t176 = Icges(6,1) * t338 - Icges(6,5) * t264 - t227;
t228 = Icges(5,4) * t339;
t178 = Icges(5,1) * t338 - Icges(5,5) * t264 - t228;
t426 = t176 + t178;
t240 = Icges(6,1) * t274 - t354;
t242 = Icges(5,1) * t274 - t355;
t424 = t240 + t242;
t236 = -Icges(6,2) * t272 + t269;
t238 = -Icges(5,2) * t272 + t270;
t423 = -t238 - t236 - t428;
t419 = (-Icges(5,6) - Icges(6,6)) * t274 + (-Icges(5,5) - Icges(6,5)) * t272;
t177 = Icges(6,5) * t263 + t240 * t264;
t179 = Icges(5,5) * t263 + t242 * t264;
t418 = -t264 * t425 + t177 + t179;
t417 = -t338 * t432 - t227 - t228 + t426;
t173 = Icges(6,6) * t263 + t236 * t264;
t175 = Icges(5,6) * t263 + t238 * t264;
t416 = -t264 * t428 - t173 - t175;
t415 = t263 * t428 + t427;
t368 = pkin(4) * t274;
t265 = pkin(3) + t368;
t359 = rSges(6,1) * t274;
t302 = t265 + t359;
t337 = t264 * t272;
t367 = -qJ(5) - pkin(7);
t309 = -rSges(6,2) * t337 - t263 * t367;
t147 = t263 * rSges(6,3) + t302 * t264 + t309;
t286 = rSges(6,1) * t338 - rSges(6,2) * t339 + t263 * t265 + (-rSges(6,3) + t367) * t264;
t65 = -t147 * t429 - t286 * t430;
t260 = t264 * pkin(7);
t360 = rSges(5,1) * t274;
t305 = pkin(3) + t360;
t310 = rSges(5,2) * t339 + rSges(5,3) * t264;
t148 = -t263 * t305 + t260 + t310;
t232 = rSges(5,2) * t337;
t149 = -t232 + t305 * t264 + (rSges(5,3) + pkin(7)) * t263;
t246 = rSges(5,1) * t272 + rSges(5,2) * t274;
t200 = t246 * t263;
t201 = t246 * t264;
t70 = t148 * t200 - t149 * t201;
t414 = -m(5) * t70 - m(6) * t65;
t397 = m(5) / 0.2e1;
t396 = m(6) / 0.2e1;
t79 = t147 * t263 - t264 * t286;
t410 = m(6) * qJD(3) * t79;
t409 = (t424 - t425) * t274 + t423 * t272;
t296 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t271);
t297 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t271);
t374 = m(4) * ((-rSges(4,1) * t263 - rSges(4,2) * t264) * t296 - (rSges(4,1) * t264 - rSges(4,2) * t263) * t297);
t134 = -t286 + t297;
t135 = t147 + t296;
t69 = t134 * t264 + t135 * t263;
t407 = m(6) * qJD(1) * t69;
t406 = t419 * t263;
t405 = t419 * t264;
t261 = t263 ^ 2;
t262 = t264 ^ 2;
t308 = t261 + t262;
t404 = -t272 * t418 + t274 * t416;
t403 = t272 * t417 + t274 * t415;
t307 = qJD(1) + qJD(3);
t144 = t148 + t297;
t145 = t149 + t296;
t366 = ((-t135 + t147) * t429 + (t134 + t286) * t430) * t396 + ((-t145 + t149) * t264 + (t144 - t148) * t263) * t246 * t397;
t63 = t134 * t430 - t135 * t429;
t68 = t144 * t200 - t145 * t201;
t402 = (t65 + t63) * t396 + (t70 + t68) * t397;
t285 = -t423 * t274 / 0.2e1 + (-t425 / 0.2e1 + t424 / 0.2e1) * t272;
t400 = 0.4e1 * qJD(1);
t398 = 2 * qJD(4);
t62 = -t144 * t149 + t145 * t148;
t393 = m(5) * t62;
t391 = m(5) * t68;
t386 = m(6) * (t69 - t79);
t385 = m(6) * (t79 + t69);
t57 = -t134 * t147 - t135 * t286;
t384 = m(6) * t57;
t382 = m(6) * t63;
t379 = -t263 / 0.2e1;
t378 = t263 / 0.2e1;
t377 = -t264 / 0.2e1;
t361 = m(6) * qJD(5);
t346 = t172 * t272;
t345 = t174 * t272;
t233 = Icges(6,5) * t274 - Icges(6,6) * t272;
t343 = t233 * t264;
t234 = Icges(5,5) * t274 - Icges(5,6) * t272;
t342 = t234 * t264;
t336 = t264 * t274;
t168 = Icges(6,5) * t338 - Icges(6,6) * t339 - Icges(6,3) * t264;
t322 = -t263 * t168 - t176 * t336;
t169 = Icges(6,3) * t263 + t343;
t321 = t169 * t263 + t177 * t336;
t170 = Icges(5,5) * t338 - Icges(5,6) * t339 - Icges(5,3) * t264;
t320 = -t263 * t170 - t178 * t336;
t171 = Icges(5,3) * t263 + t342;
t319 = t171 * t263 + t179 * t336;
t299 = t173 * t272 - t168;
t152 = t177 * t338;
t301 = t169 * t264 - t152;
t84 = -t172 * t337 - t322;
t85 = -t173 * t337 + t321;
t10 = (t264 * t299 - t321 + t85) * t264 + (t263 * t299 + t301 + t84) * t263;
t298 = t175 * t272 - t170;
t153 = t179 * t338;
t300 = t171 * t264 - t153;
t86 = -t174 * t337 - t320;
t87 = -t175 * t337 + t319;
t11 = (t264 * t298 - t319 + t87) * t264 + (t263 * t298 + t300 + t86) * t263;
t81 = -t173 * t339 - t301;
t12 = (t81 - t152 + (t169 + t346) * t264 + t322) * t264 + t321 * t263;
t83 = -t175 * t339 - t300;
t13 = (t83 - t153 + (t171 + t345) * t264 + t320) * t264 + t319 * t263;
t48 = t263 * t81 - t264 * (-t263 * (-t176 * t274 + t346) - t168 * t264);
t49 = t263 * t83 - t264 * (-t263 * (-t178 * t274 + t345) - t170 * t264);
t50 = t263 * t85 - t264 * t84;
t51 = t263 * t87 - t264 * t86;
t2 = (t51 / 0.2e1 - t13 / 0.2e1 + t50 / 0.2e1 - t12 / 0.2e1) * t264 + (t11 / 0.2e1 + t49 / 0.2e1 + t10 / 0.2e1 + t48 / 0.2e1) * t263;
t60 = 0.2e1 * t120;
t306 = qJD(4) * t2 + qJD(5) * t60;
t303 = rSges(6,2) * t272 - t359 - t368;
t287 = (-t200 * t264 + t201 * t263) * t246;
t278 = t285 + t402;
t277 = -t285 + (t272 * t426 + t274 * t427) * (t378 + t379);
t59 = t120 + t107;
t276 = t59 * qJD(5) + ((t12 + t13) * t264 / 0.2e1 + (t10 + t11 + t48 + t49) * t379 + (t418 * t274 + t416 * t272 + t409 * t264 + (t233 + t234) * t263) * t378 + (t409 * t263 - t272 * t415 + t274 * t417 - t342 - t343 + t50 + t51) * t377) * qJD(4);
t249 = -rSges(5,2) * t272 + t360;
t185 = t303 * t264;
t183 = t303 * t263;
t141 = -t200 * t263 - t201 * t264;
t110 = t422 * t308;
t58 = 0.2e1 * t107;
t55 = t59 * qJD(4);
t53 = t58 * qJD(4);
t35 = t385 / 0.2e1;
t34 = t386 / 0.2e1;
t23 = t285 - t414;
t20 = t285 + t382 + t391;
t17 = t374 + t384 + t393;
t16 = t35 - t386 / 0.2e1;
t15 = t35 + t34;
t14 = t34 - t385 / 0.2e1;
t5 = t278 - t366;
t4 = t278 + t366;
t3 = t277 + t366 - t402;
t1 = [qJD(3) * t17 + qJD(4) * t20 + t361 * t69, 0, t17 * qJD(1) + t4 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t57 * t396 + t62 * t397 + t374 / 0.2e1) * qJD(3), t20 * qJD(1) + t4 * qJD(3) + ((t134 * t185 + t135 * t183) * t396 + ((-t144 * t264 - t145 * t263) * t249 + t287) * t397) * t398 + t276, qJD(3) * t15 + t407 + t55; 0, 0, 0, (-t110 * t396 + t141 * t397) * t398, 0; t5 * qJD(4) + t16 * qJD(5) + (-t384 / 0.4e1 - t393 / 0.4e1 - t374 / 0.4e1) * t400, 0, qJD(4) * t23 + t361 * t79, t5 * qJD(1) + t23 * qJD(3) + ((t147 * t183 - t185 * t286) * t396 + ((-t148 * t264 - t149 * t263) * t249 + t287) * t397) * t398 + t276, qJD(1) * t16 + t410 + t55; t277 * qJD(1) + t3 * qJD(3) + (-t382 / 0.4e1 - t391 / 0.4e1) * t400 + t306, 0, t3 * qJD(1) + t306 + (t277 + t414) * qJD(3), (m(5) * (t246 * t249 * t308 + (t263 * (rSges(5,1) * t338 - t310) + t264 * (rSges(5,1) * t336 + t263 * rSges(5,3) - t232)) * t141) + m(6) * (-t110 * ((-pkin(3) * t263 + t260 + t286) * t263 + ((-pkin(3) + t302) * t264 + (rSges(6,3) - pkin(7)) * t263 + t309) * t264) - t430 * t183 - t429 * t185) + (t405 * t261 + (t403 * t264 + (t404 - t406) * t263) * t264) * t378 + (t406 * t262 + (t404 * t263 + (t403 - t405) * t264) * t263) * t377) * qJD(4) + t307 * t2, t307 * t60; t14 * qJD(3) - t407 + t53, 0, t14 * qJD(1) - t410 + t53, m(6) * (-t183 * t264 + t185 * t263) * qJD(4) + t307 * t58, 0;];
Cq = t1;
