% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:02
% EndTime: 2019-12-31 17:14:08
% DurationCPUTime: 3.92s
% Computational Cost: add. (13318->283), mult. (13885->362), div. (0->0), fcn. (12489->6), ass. (0->190)
t271 = qJ(1) + qJ(2);
t268 = cos(t271);
t272 = sin(qJ(3));
t337 = t268 * t272;
t267 = sin(t271);
t274 = cos(qJ(3));
t368 = rSges(5,1) + pkin(3);
t302 = t368 * t272;
t357 = rSges(5,3) + qJ(4);
t428 = -t357 * t274 + t302;
t433 = t428 * t267;
t434 = t337 * t433;
t269 = Icges(5,5) * t272;
t292 = Icges(5,3) * t274 - t269;
t354 = Icges(4,4) * t272;
t431 = Icges(4,2) * t274 + t292 + t354;
t270 = Icges(4,4) * t274;
t237 = Icges(4,1) * t272 + t270;
t353 = Icges(5,5) * t274;
t432 = Icges(5,1) * t272 + t237 - t353;
t234 = -Icges(4,2) * t272 + t270;
t430 = t234 + t432;
t238 = Icges(4,1) * t274 - t354;
t405 = Icges(5,1) * t274 + t269;
t429 = t238 + t405;
t265 = t267 ^ 2;
t266 = t268 ^ 2;
t305 = t265 + t266;
t406 = t357 * t272 + t368 * t274;
t264 = t268 * pkin(6);
t401 = pkin(2) + t406;
t140 = t268 * rSges(5,2) - t401 * t267 + t264;
t365 = sin(qJ(1)) * pkin(1);
t131 = t140 - t365;
t426 = t131 - t140;
t424 = (-Icges(4,6) + Icges(5,6)) * t274 + (-Icges(5,4) - Icges(4,5)) * t272;
t186 = Icges(5,4) * t267 + t268 * t405;
t188 = Icges(4,5) * t267 + t238 * t268;
t423 = -t431 * t268 + t186 + t188;
t185 = -Icges(5,4) * t268 + t267 * t405;
t339 = t267 * t272;
t251 = Icges(4,4) * t339;
t338 = t267 * t274;
t187 = Icges(4,1) * t338 - Icges(4,5) * t268 - t251;
t422 = -Icges(4,2) * t338 - t292 * t267 + t185 + t187 - t251;
t336 = t268 * t274;
t250 = Icges(5,5) * t336;
t178 = Icges(5,6) * t267 + Icges(5,3) * t337 + t250;
t184 = Icges(4,6) * t267 + t234 * t268;
t421 = -Icges(5,1) * t337 - t237 * t268 + t178 - t184 + t250;
t230 = Icges(5,3) * t272 + t353;
t177 = -Icges(5,6) * t268 + t230 * t267;
t183 = Icges(4,4) * t338 - Icges(4,2) * t339 - Icges(4,6) * t268;
t420 = t432 * t267 - t177 + t183;
t141 = (rSges(5,2) + pkin(6)) * t267 + t401 * t268;
t154 = -t268 * t302 + t357 * t336;
t67 = t140 * t433 + t141 * t154;
t358 = rSges(4,1) * t274;
t301 = pkin(2) + t358;
t306 = rSges(4,2) * t339 + t268 * rSges(4,3);
t151 = -t301 * t267 + t264 + t306;
t253 = rSges(4,2) * t337;
t152 = -t253 + t301 * t268 + (rSges(4,3) + pkin(6)) * t267;
t241 = rSges(4,1) * t272 + rSges(4,2) * t274;
t203 = t241 * t267;
t205 = t241 * t268;
t77 = t151 * t203 - t152 * t205;
t419 = -m(4) * t77 - m(5) * t67;
t396 = m(4) / 0.2e1;
t395 = m(5) / 0.2e1;
t364 = cos(qJ(1)) * pkin(1);
t367 = m(3) * (t364 * (-rSges(3,1) * t267 - rSges(3,2) * t268) + (t268 * rSges(3,1) - t267 * rSges(3,2)) * t365);
t232 = Icges(5,4) * t274 + Icges(5,6) * t272;
t343 = t232 * t267;
t181 = -Icges(5,2) * t268 + t343;
t168 = t267 * t181;
t92 = t177 * t337 + t185 * t336 + t168;
t415 = t268 * t92;
t303 = t141 * t337;
t73 = -t140 * t339 + t303;
t414 = t73 * m(5) * qJD(2);
t413 = (t429 - t431) * t274 + (t230 - t430) * t272;
t410 = (t177 * t272 + t185 * t274) * t267;
t146 = t151 - t365;
t147 = t152 + t364;
t68 = -t152 * t146 + t147 * t151;
t408 = t424 * t267;
t407 = t424 * t268;
t404 = -t423 * t272 + t421 * t274;
t403 = t422 * t272 + t420 * t274;
t132 = t141 + t364;
t172 = t428 * t268;
t363 = ((-t132 + t141) * t172 + t426 * t433) * t395 + ((-t147 + t152) * t268 + (t146 - t151) * t267) * t241 * t396;
t65 = t131 * t433 + t132 * t154;
t74 = t146 * t203 - t147 * t205;
t402 = (t67 + t65) * t395 + (t77 + t74) * t396;
t285 = (-t230 / 0.2e1 + t430 / 0.2e1) * t274 + (-t431 / 0.2e1 + t429 / 0.2e1) * t272;
t400 = 4 * qJD(1);
t398 = 2 * qJD(3);
t392 = m(4) * t68;
t390 = m(4) * t74;
t117 = t132 * t337;
t384 = m(5) * (-t426 * t339 + t117 - t303);
t383 = m(5) * (t117 + t303 + (-t131 - t140) * t339);
t42 = -t141 * t131 + t132 * t140;
t382 = m(5) * t42;
t321 = t154 * t339 + t434;
t323 = t131 * t336 + t132 * t338;
t381 = m(5) * (t321 + t323);
t322 = t140 * t336 + t141 * t338;
t379 = m(5) * (t321 + t322);
t298 = t172 * t339 - t434;
t378 = m(5) * (t298 + t323);
t377 = m(5) * (t298 + t322);
t376 = m(5) * t65;
t372 = -t267 / 0.2e1;
t371 = t267 / 0.2e1;
t370 = -t268 / 0.2e1;
t360 = m(5) * qJD(3);
t359 = m(5) * qJD(4);
t347 = t181 * t268;
t346 = t183 * t272;
t231 = Icges(4,5) * t274 - Icges(4,6) * t272;
t344 = t231 * t268;
t342 = t232 * t268;
t331 = t272 * t274;
t324 = (-t154 - t172) * t433;
t320 = -t172 * t336 - t338 * t433;
t179 = Icges(4,5) * t338 - Icges(4,6) * t339 - Icges(4,3) * t268;
t319 = -t267 * t179 - t187 * t336;
t180 = Icges(4,3) * t267 + t344;
t318 = t267 * t180 + t188 * t336;
t309 = t305 * t331;
t70 = -t131 * t339 + t117;
t304 = m(5) * t70 * qJD(1);
t182 = Icges(5,2) * t267 + t342;
t93 = t178 * t337 + t267 * t182 + t186 * t336;
t157 = t188 * t338;
t300 = t180 * t268 - t157;
t299 = t184 * t272 - t179;
t297 = -t178 * t339 + t182 * t268 - t186 * t338;
t287 = t93 + t347;
t286 = (-t203 * t268 + t205 * t267) * t241;
t278 = t285 + t402;
t277 = -t285 + ((t177 + t183) * t274 + (-t185 + t187) * t272) * (t371 + t372);
t19 = (-t287 + t93) * t268 + (t297 + t92 - t168) * t267;
t94 = -t183 * t337 - t319;
t95 = -t184 * t337 + t318;
t20 = (t299 * t268 - t318 + t95) * t268 + (t299 * t267 + t300 + t94) * t267;
t88 = -t347 + t410;
t21 = -t415 + (t287 + t88 - t410) * t267;
t91 = -t184 * t339 - t300;
t22 = (t91 - t157 + (t180 + t346) * t268 + t319) * t268 + t318 * t267;
t57 = -t267 * t297 - t268 * t88;
t58 = t267 * t91 - t268 * (-(-t187 * t274 + t346) * t267 - t179 * t268);
t59 = t267 * t93 - t415;
t60 = t267 * t95 - t268 * t94;
t276 = ((t21 + t22) * t268 / 0.2e1 + (t19 + t20 + t57 + t58) * t372 + (t231 * t267 + t413 * t268 + t421 * t272 + t423 * t274 + t343) * t371 + (t413 * t267 - t420 * t272 + t422 * t274 - t342 - t344 + t59 + t60) * t370) * qJD(3);
t245 = -rSges(4,2) * t272 + t358;
t208 = t305 * t272;
t174 = t309 - t331;
t173 = t406 * t268;
t171 = t406 * t267;
t98 = t154 * t268 - t265 * t428;
t85 = t305 * t406;
t62 = t208 * t85 + t320;
t56 = t377 / 0.2e1;
t49 = t378 / 0.2e1;
t47 = t379 / 0.2e1;
t43 = t381 / 0.2e1;
t37 = t383 / 0.2e1;
t36 = t384 / 0.2e1;
t27 = t285 - t419;
t26 = t285 + t376 + t390;
t23 = t367 + t382 + t392;
t14 = t56 - t379 / 0.2e1;
t13 = t56 + t47;
t12 = t47 - t377 / 0.2e1;
t11 = t49 - t381 / 0.2e1;
t10 = t49 + t43;
t9 = t43 - t378 / 0.2e1;
t8 = t37 - t384 / 0.2e1;
t7 = t37 + t36;
t6 = t36 - t383 / 0.2e1;
t5 = t278 + t363;
t4 = t278 - t363;
t3 = t277 + t363 - t402;
t2 = (t59 / 0.2e1 - t21 / 0.2e1 + t60 / 0.2e1 - t22 / 0.2e1) * t268 + (t19 / 0.2e1 + t57 / 0.2e1 + t20 / 0.2e1 + t58 / 0.2e1) * t267;
t1 = t2 * qJD(3);
t15 = [t23 * qJD(2) + t26 * qJD(3) + t70 * t359, t23 * qJD(1) + t5 * qJD(3) + t7 * qJD(4) + 0.2e1 * (t367 / 0.2e1 + t42 * t395 + t68 * t396) * qJD(2), t26 * qJD(1) + t5 * qJD(2) + t276 + t10 * qJD(4) + (((-t146 * t268 - t147 * t267) * t245 + t286) * t396 + (-t131 * t173 - t132 * t171 + t324) * t395) * t398, t7 * qJD(2) + t10 * qJD(3) + t304; t4 * qJD(3) + t8 * qJD(4) + (-t367 / 0.4e1 - t392 / 0.4e1 - t382 / 0.4e1) * t400, t27 * qJD(3) + t73 * t359, t4 * qJD(1) + t27 * qJD(2) + t276 + t13 * qJD(4) + (((-t151 * t268 - t152 * t267) * t245 + t286) * t396 + (-t140 * t173 - t141 * t171 + t324) * t395) * t398, t8 * qJD(1) + t13 * qJD(3) + t414; t277 * qJD(1) + t3 * qJD(2) + t1 + t11 * qJD(4) + (-t390 / 0.4e1 - t376 / 0.4e1) * t400, t3 * qJD(1) + t14 * qJD(4) + t1 + (t277 + t419) * qJD(2), (m(4) * ((t267 * (rSges(4,1) * t338 - t306) + t268 * (rSges(4,1) * t336 + t267 * rSges(4,3) - t253)) * (-t203 * t267 - t205 * t268) + t305 * t245 * t241) + m(5) * (t171 * t433 + t172 * t173 + t85 * t98) + (t407 * t265 + (t403 * t268 + (t404 - t408) * t267) * t268) * t371 + (t408 * t266 + (t404 * t267 + (t403 - t407) * t268) * t267) * t370) * qJD(3) + t62 * t359 + (qJD(1) + qJD(2)) * t2, t11 * qJD(1) + t14 * qJD(2) + t62 * t360 + (-t208 * t274 - t174 + t309) * t359; t6 * qJD(2) + t9 * qJD(3) - t304, t6 * qJD(1) + t12 * qJD(3) - t414, t9 * qJD(1) + t12 * qJD(2) + (-t274 * t98 + (-t171 * t267 - t173 * t268 + t85) * t272 - t62 + t320) * t360 + t174 * t359, t174 * t360;];
Cq = t15;
