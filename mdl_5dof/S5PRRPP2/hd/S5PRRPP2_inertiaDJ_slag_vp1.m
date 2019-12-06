% Calculate time derivative of joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:44
% EndTime: 2019-12-05 16:09:15
% DurationCPUTime: 14.88s
% Computational Cost: add. (12032->578), mult. (21422->834), div. (0->0), fcn. (20236->8), ass. (0->277)
t433 = Icges(5,4) - Icges(6,5);
t432 = Icges(5,1) + Icges(6,1);
t424 = Icges(6,4) + Icges(5,5);
t431 = Icges(5,2) + Icges(6,3);
t423 = Icges(5,6) - Icges(6,6);
t236 = qJ(3) + pkin(8);
t233 = cos(t236);
t430 = t433 * t233;
t232 = sin(t236);
t429 = t433 * t232;
t239 = cos(pkin(7));
t238 = sin(pkin(7));
t244 = cos(qJ(2));
t345 = t238 * t244;
t202 = t232 * t345 + t239 * t233;
t203 = -t239 * t232 + t233 * t345;
t242 = sin(qJ(2));
t346 = t238 * t242;
t125 = Icges(6,5) * t203 + Icges(6,6) * t346 + Icges(6,3) * t202;
t131 = Icges(5,4) * t203 - Icges(5,2) * t202 + Icges(5,6) * t346;
t422 = t125 - t131;
t342 = t239 * t244;
t204 = t232 * t342 - t238 * t233;
t205 = t238 * t232 + t233 * t342;
t343 = t239 * t242;
t126 = Icges(6,5) * t205 + Icges(6,6) * t343 + Icges(6,3) * t204;
t132 = Icges(5,4) * t205 - Icges(5,2) * t204 + Icges(5,6) * t343;
t421 = t126 - t132;
t133 = Icges(6,1) * t203 + Icges(6,4) * t346 + Icges(6,5) * t202;
t135 = Icges(5,1) * t203 - Icges(5,4) * t202 + Icges(5,5) * t346;
t420 = t133 + t135;
t134 = Icges(6,1) * t205 + Icges(6,4) * t343 + Icges(6,5) * t204;
t136 = Icges(5,1) * t205 - Icges(5,4) * t204 + Icges(5,5) * t343;
t419 = t134 + t136;
t428 = t431 * t232 - t430;
t427 = t432 * t233 - t429;
t426 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t128 = Icges(5,5) * t205 - Icges(5,6) * t204 + Icges(5,3) * t343;
t130 = Icges(6,4) * t205 + Icges(6,2) * t343 + Icges(6,6) * t204;
t243 = cos(qJ(3));
t241 = sin(qJ(3));
t341 = t241 * t244;
t221 = t238 * t243 - t239 * t341;
t340 = t243 * t244;
t347 = t238 * t241;
t222 = t239 * t340 + t347;
t160 = Icges(4,5) * t222 + Icges(4,6) * t221 + Icges(4,3) * t343;
t415 = t128 + t130 + t160;
t127 = Icges(5,5) * t203 - Icges(5,6) * t202 + Icges(5,3) * t346;
t129 = Icges(6,4) * t203 + Icges(6,2) * t346 + Icges(6,6) * t202;
t219 = -t238 * t341 - t239 * t243;
t344 = t239 * t241;
t220 = t238 * t340 - t344;
t159 = Icges(4,5) * t220 + Icges(4,6) * t219 + Icges(4,3) * t346;
t412 = -t129 - t127 - t159;
t425 = Icges(4,5) * t243 - Icges(4,6) * t241 - t423 * t232 + t424 * t233;
t414 = t242 * t425 + t244 * t426;
t402 = t414 * t244;
t407 = t242 * t428 + t244 * t423;
t418 = t427 * t242 - t244 * t424;
t162 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t343;
t164 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t343;
t417 = t162 * t219 + t164 * t220 + t202 * t421 + t203 * t419 + t346 * t415;
t161 = Icges(4,4) * t220 + Icges(4,2) * t219 + Icges(4,6) * t346;
t163 = Icges(4,1) * t220 + Icges(4,4) * t219 + Icges(4,5) * t346;
t416 = t161 * t221 + t163 * t222 + t204 * t422 + t205 * t420 - t343 * t412;
t411 = -t162 * t241 + t164 * t243 + t421 * t232 + t419 * t233;
t410 = -t161 * t241 + t163 * t243 + t422 * t232 + t420 * t233;
t317 = qJD(3) * t242;
t409 = (t431 * t233 + t429) * t317 + (-t423 * t242 + t244 * t428) * qJD(2);
t408 = (-t432 * t232 - t430) * t317 + (t424 * t242 + t244 * t427) * qJD(2);
t386 = t161 * t219 + t163 * t220 + t422 * t202 + t420 * t203 - t412 * t346;
t385 = t162 * t221 + t164 * t222 + t421 * t204 + t419 * t205 + t415 * t343;
t356 = Icges(4,4) * t243;
t276 = -Icges(4,2) * t241 + t356;
t207 = -Icges(4,6) * t244 + t242 * t276;
t357 = Icges(4,4) * t241;
t280 = Icges(4,1) * t243 - t357;
t208 = -Icges(4,5) * t244 + t242 * t280;
t405 = t407 * t202 + t418 * t203 + t207 * t219 + t208 * t220 + t414 * t346;
t404 = t407 * t204 + t418 * t205 + t207 * t221 + t208 * t222 + t414 * t343;
t403 = ((Icges(4,5) * t241 + Icges(4,6) * t243 + t424 * t232 + t423 * t233) * t317 + (t242 * t426 - t244 * t425) * qJD(2)) * t244;
t401 = t417 * t239;
t400 = t416 * t238;
t399 = t242 * t410 + t244 * t412;
t398 = t242 * t411 - t244 * t415;
t397 = t207 * t241 - t208 * t243 - t232 * t407 - t233 * t418;
t396 = rSges(6,1) + pkin(4);
t395 = rSges(6,3) + qJ(5);
t234 = t238 ^ 2;
t235 = t239 ^ 2;
t388 = t234 + t235;
t321 = qJD(2) * t244;
t171 = (-Icges(4,2) * t243 - t357) * t317 + (Icges(4,6) * t242 + t244 * t276) * qJD(2);
t172 = (-Icges(4,1) * t241 - t356) * t317 + (Icges(4,5) * t242 + t244 * t280) * qJD(2);
t319 = qJD(3) * t233;
t298 = t244 * t319;
t322 = qJD(2) * t242;
t306 = t238 * t322;
t176 = -t238 * t298 + (qJD(3) * t239 + t306) * t232;
t177 = -qJD(3) * t202 - t233 * t306;
t302 = t241 * t322;
t182 = -qJD(3) * t220 + t238 * t302;
t301 = t243 * t322;
t183 = qJD(3) * t219 - t238 * t301;
t305 = t238 * t321;
t91 = Icges(6,4) * t177 + Icges(6,2) * t305 - Icges(6,6) * t176;
t255 = t129 * t321 + t242 * t91;
t87 = Icges(6,5) * t177 + Icges(6,6) * t305 - Icges(6,3) * t176;
t95 = Icges(6,1) * t177 + Icges(6,4) * t305 - Icges(6,5) * t176;
t22 = -t176 * t125 + t177 * t133 + t202 * t87 + t203 * t95 + t238 * t255;
t304 = t239 * t322;
t320 = qJD(3) * t232;
t178 = t232 * t304 - t238 * t320 - t239 * t298;
t179 = -qJD(3) * t204 - t233 * t304;
t303 = t239 * t321;
t92 = Icges(6,4) * t179 + Icges(6,2) * t303 - Icges(6,6) * t178;
t254 = t130 * t321 + t242 * t92;
t88 = Icges(6,5) * t179 + Icges(6,6) * t303 - Icges(6,3) * t178;
t96 = Icges(6,1) * t179 + Icges(6,4) * t303 - Icges(6,5) * t178;
t23 = -t176 * t126 + t177 * t134 + t202 * t88 + t203 * t96 + t238 * t254;
t89 = Icges(5,5) * t177 + Icges(5,6) * t176 + Icges(5,3) * t305;
t257 = t127 * t321 + t242 * t89;
t93 = Icges(5,4) * t177 + Icges(5,2) * t176 + Icges(5,6) * t305;
t97 = Icges(5,1) * t177 + Icges(5,4) * t176 + Icges(5,5) * t305;
t24 = t176 * t131 + t177 * t135 - t202 * t93 + t203 * t97 + t238 * t257;
t90 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t303;
t256 = t128 * t321 + t242 * t90;
t94 = Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t303;
t98 = Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t303;
t25 = t176 * t132 + t177 * t136 - t202 * t94 + t203 * t98 + t238 * t256;
t112 = Icges(4,4) * t183 + Icges(4,2) * t182 + Icges(4,6) * t305;
t114 = Icges(4,1) * t183 + Icges(4,4) * t182 + Icges(4,5) * t305;
t110 = Icges(4,5) * t183 + Icges(4,6) * t182 + Icges(4,3) * t305;
t252 = t110 * t242 + t159 * t321;
t32 = t219 * t112 + t220 * t114 + t182 * t161 + t183 * t163 + t238 * t252;
t184 = -qJD(3) * t222 + t239 * t302;
t185 = qJD(3) * t221 - t239 * t301;
t113 = Icges(4,4) * t185 + Icges(4,2) * t184 + Icges(4,6) * t303;
t115 = Icges(4,1) * t185 + Icges(4,4) * t184 + Icges(4,5) * t303;
t111 = Icges(4,5) * t185 + Icges(4,6) * t184 + Icges(4,3) * t303;
t251 = t111 * t242 + t160 * t321;
t33 = t219 * t113 + t220 * t115 + t182 * t162 + t183 * t164 + t238 * t251;
t394 = (-t219 * t171 - t220 * t172 + t407 * t176 - t177 * t418 - t182 * t207 - t183 * t208 - t409 * t202 - t408 * t203) * t244 + ((t23 + t25 + t33) * t239 + (t22 + t24 + t32 + t403) * t238) * t242 + (((-t402 + t386) * t238 + t401) * t244 + t405 * t242) * qJD(2);
t26 = -t178 * t125 + t179 * t133 + t204 * t87 + t205 * t95 + t239 * t255;
t27 = -t178 * t126 + t179 * t134 + t204 * t88 + t205 * t96 + t239 * t254;
t28 = t178 * t131 + t179 * t135 - t204 * t93 + t205 * t97 + t239 * t257;
t29 = t178 * t132 + t179 * t136 - t204 * t94 + t205 * t98 + t239 * t256;
t34 = t221 * t112 + t222 * t114 + t184 * t161 + t185 * t163 + t239 * t252;
t35 = t221 * t113 + t222 * t115 + t184 * t162 + t185 * t164 + t239 * t251;
t393 = (-t221 * t171 - t222 * t172 + t407 * t178 - t179 * t418 - t184 * t207 - t185 * t208 - t409 * t204 - t408 * t205) * t244 + ((t27 + t29 + t35 + t403) * t239 + (t26 + t28 + t34) * t238) * t242 + (((-t402 + t385) * t239 + t400) * t244 + t404 * t242) * qJD(2);
t392 = (-qJD(2) * t410 + t110 + t89 + t91) * t244 + (t112 * t241 - t114 * t243 + (-t95 - t97) * t233 + (-t87 + t93) * t232 + (t161 * t243 + t163 * t241 + t420 * t232 - t422 * t233) * qJD(3) + t412 * qJD(2)) * t242;
t391 = (-qJD(2) * t411 + t111 + t90 + t92) * t244 + (t113 * t241 - t115 * t243 + (-t96 - t98) * t233 + (-t88 + t94) * t232 + (t162 * t243 + t164 * t241 + t419 * t232 - t421 * t233) * qJD(3) - t415 * qJD(2)) * t242;
t387 = qJD(2) * (t242 * rSges(3,1) + rSges(3,2) * t244);
t384 = t397 * t242 + t402;
t383 = t399 * t238 + t398 * t239;
t369 = pkin(3) * t243;
t187 = -qJ(4) * t244 + t242 * t369;
t380 = 2 * m(4);
t379 = 2 * m(5);
t378 = 2 * m(6);
t377 = 0.2e1 * qJD(2);
t376 = m(5) / 0.2e1;
t375 = m(6) / 0.2e1;
t367 = rSges(6,2) * t305 + qJD(5) * t202 - t395 * t176 + t396 * t177;
t360 = rSges(6,2) * t303 + qJD(5) * t204 - t395 * t178 + t396 * t179;
t102 = t179 * rSges(5,1) + t178 * rSges(5,2) + rSges(5,3) * t303;
t318 = qJD(3) * t241;
t315 = pkin(3) * t318;
t245 = -qJD(2) * t187 + qJD(4) * t242 - t244 * t315;
t316 = qJD(3) * t243;
t314 = pkin(3) * t316;
t119 = t238 * t314 + t239 * t245;
t336 = -t102 - t119;
t118 = t238 * t245 - t239 * t314;
t246 = qJ(4) * t242 + t244 * t369;
t165 = -pkin(3) * t344 + t238 * t246;
t335 = t118 * t343 + t165 * t303;
t285 = rSges(6,1) * t233 + rSges(6,3) * t232;
t334 = (pkin(4) * t321 + qJ(5) * t317) * t233 + (qJ(5) * t321 + (-pkin(4) * qJD(3) + qJD(5)) * t242) * t232 + (-rSges(6,1) * t232 + rSges(6,3) * t233) * t317 + (rSges(6,2) * t242 + t244 * t285) * qJD(2);
t333 = rSges(6,2) * t346 + t395 * t202 + t396 * t203;
t332 = rSges(6,2) * t343 + t395 * t204 + t396 * t205;
t140 = rSges(5,1) * t205 - rSges(5,2) * t204 + rSges(5,3) * t343;
t166 = pkin(3) * t347 + t239 * t246;
t331 = -t140 - t166;
t286 = rSges(5,1) * t233 - rSges(5,2) * t232;
t150 = (-rSges(5,1) * t232 - rSges(5,2) * t233) * t317 + (rSges(5,3) * t242 + t244 * t286) * qJD(2);
t169 = qJD(2) * t246 - qJD(4) * t244 - t242 * t315;
t330 = -t150 - t169;
t329 = t244 * t165 + t187 * t346;
t287 = rSges(4,1) * t243 - rSges(4,2) * t241;
t173 = (-rSges(4,1) * t241 - rSges(4,2) * t243) * t317 + (rSges(4,3) * t242 + t244 * t287) * qJD(2);
t289 = pkin(2) * t244 + pkin(6) * t242;
t226 = t289 * qJD(2);
t328 = -t173 - t226;
t195 = -rSges(5,3) * t244 + t242 * t286;
t327 = -t187 - t195;
t326 = -rSges(6,2) * t244 + (pkin(4) * t233 + qJ(5) * t232 + t285) * t242;
t229 = t242 * pkin(2) - pkin(6) * t244;
t325 = t388 * qJD(2) * t229;
t210 = -rSges(4,3) * t244 + t242 * t287;
t324 = -t210 - t229;
t323 = t388 * t289;
t313 = -t119 - t360;
t312 = t244 * t118 + t169 * t346 + t187 * t305;
t311 = -t169 - t334;
t310 = -t166 - t332;
t309 = -t226 + t330;
t308 = -t187 - t326;
t307 = -t229 + t327;
t300 = t242 * t321;
t299 = t233 * t317;
t297 = t238 * t326;
t296 = t333 * t244;
t295 = t331 * t244;
t294 = t238 * t118 + t239 * t119 - t325;
t293 = -t226 + t311;
t292 = t238 * t165 + t239 * t166 + t323;
t291 = -t229 + t308;
t290 = t310 * t244;
t167 = rSges(4,1) * t220 + rSges(4,2) * t219 + rSges(4,3) * t346;
t168 = rSges(4,1) * t222 + rSges(4,2) * t221 + rSges(4,3) * t343;
t264 = t167 * t239 - t168 * t238;
t247 = qJD(2) * (-Icges(3,5) * t242 - Icges(3,6) * t244);
t237 = t242 ^ 2;
t214 = t239 * t247;
t213 = t238 * t247;
t175 = t324 * t239;
t174 = t324 * t238;
t158 = t388 * t387;
t152 = t166 * t322;
t148 = t165 * t343;
t138 = rSges(5,1) * t203 - rSges(5,2) * t202 + rSges(5,3) * t346;
t123 = t328 * t239;
t122 = t328 * t238;
t121 = t307 * t239;
t120 = t307 * t238;
t117 = t185 * rSges(4,1) + t184 * rSges(4,2) + rSges(4,3) * t303;
t116 = t183 * rSges(4,1) + t182 * rSges(4,2) + rSges(4,3) * t305;
t108 = -t168 * t244 - t210 * t343;
t107 = t167 * t244 + t210 * t346;
t100 = t177 * rSges(5,1) + t176 * rSges(5,2) + rSges(5,3) * t305;
t86 = t291 * t239;
t85 = t291 * t238;
t78 = t264 * t242;
t77 = t309 * t239;
t76 = t309 * t238;
t71 = t167 * t238 + t168 * t239 + t323;
t68 = t327 * t343 + t295;
t67 = t138 * t244 + t195 * t346 + t329;
t62 = t293 * t239;
t61 = t293 * t238;
t56 = t116 * t238 + t117 * t239 - t325;
t55 = -t173 * t343 - t117 * t244 + (t168 * t242 - t210 * t342) * qJD(2);
t54 = t173 * t346 + t116 * t244 + (-t167 * t242 + t210 * t345) * qJD(2);
t45 = t148 + (t138 * t239 + t238 * t331) * t242;
t44 = t308 * t343 + t290;
t43 = t242 * t297 + t296 + t329;
t42 = t138 * t238 + t140 * t239 + t292;
t41 = (t116 * t239 - t117 * t238) * t242 + t264 * t321;
t40 = t148 + (t238 * t310 + t333 * t239) * t242;
t39 = t238 * t333 + t332 * t239 + t292;
t38 = t100 * t238 + t102 * t239 + t294;
t37 = t140 * t322 + t152 + t336 * t244 + (t242 * t330 + t321 * t327) * t239;
t36 = t150 * t346 + t100 * t244 + (t195 * t345 + (-t138 - t165) * t242) * qJD(2) + t312;
t21 = (t100 * t242 + t138 * t321) * t239 + (qJD(2) * t295 + t242 * t336) * t238 + t335;
t20 = t238 * t367 + t360 * t239 + t294;
t15 = t152 + t332 * t322 + t313 * t244 + (t242 * t311 + t308 * t321) * t239;
t14 = t367 * t244 + t334 * t346 + (t244 * t297 + (-t165 - t333) * t242) * qJD(2) + t312;
t13 = (qJD(2) * t296 + t242 * t367) * t239 + (qJD(2) * t290 + t242 * t313) * t238 + t335;
t12 = t238 * t35 - t239 * t34;
t11 = t238 * t33 - t239 * t32;
t10 = t238 * t29 - t239 * t28;
t9 = t238 * t27 - t239 * t26;
t8 = t238 * t25 - t239 * t24;
t7 = -t22 * t239 + t23 * t238;
t1 = [0; -m(3) * t158 + m(4) * t56 + m(5) * t38 + m(6) * t20; (t20 * t39 + t61 * t85 + t62 * t86) * t378 + (t120 * t76 + t121 * t77 + t38 * t42) * t379 + (t122 * t174 + t123 * t175 + t56 * t71) * t380 + 0.2e1 * m(3) * (-t158 + t387) * t388 * (rSges(3,1) * t244 - rSges(3,2) * t242) + (-t235 * t213 - t11 - t7 - t8) * t239 + (t234 * t214 + t10 + t12 + t9 + (-t238 * t213 + t239 * t214) * t239) * t238; m(4) * t41 + m(5) * t21 + m(6) * t13; m(4) * (t107 * t123 + t108 * t122 + t174 * t55 + t175 * t54 + t41 * t71 + t56 * t78) + m(5) * (t120 * t37 + t121 * t36 + t21 * t42 + t38 * t45 + t67 * t77 + t68 * t76) + m(6) * (t13 * t39 + t14 * t86 + t15 * t85 + t20 * t40 + t43 * t62 + t44 * t61) + ((t9 / 0.2e1 + t10 / 0.2e1 + t12 / 0.2e1) * t239 + (t7 / 0.2e1 + t8 / 0.2e1 + t11 / 0.2e1) * t238) * t242 + ((t385 * t238 - t416 * t239) * t342 / 0.2e1 + (t398 * t238 - t399 * t239) * t242 / 0.2e1) * qJD(2) + ((t417 * t238 - t386 * t239) * t321 + t393) * t238 / 0.2e1 - t394 * t239 / 0.2e1 - (-t238 * t391 + t239 * t392) * t244 / 0.2e1; (t21 * t45 + t36 * t67 + t37 * t68) * t379 + (t13 * t40 + t14 * t43 + t15 * t44) * t378 + (t107 * t54 + t108 * t55 + t41 * t78) * t380 + t394 * t346 + t393 * t343 + (t383 * t322 + (t238 * t386 + t401) * t305 + (t239 * t385 + t400) * t303) * t242 + (t403 * t244 + t384 * t322 - t405 * t305 - t404 * t303 + ((-t171 * t241 + t172 * t243 - t207 * t316 - t208 * t318 + t409 * t232 + t408 * t233 + t407 * t319 - t320 * t418) * t244 + t391 * t239 + t392 * t238) * t242 + ((-t397 * t244 - t383) * t244 + (t384 + t402) * t242) * qJD(2)) * t244; (m(5) + m(6)) * t322; m(6) * (-t20 * t244 + t343 * t62 + t346 * t61) + m(5) * (-t244 * t38 + t343 * t77 + t346 * t76) + ((t242 * t39 + t342 * t86 + t345 * t85) * t375 + (t120 * t345 + t121 * t342 + t242 * t42) * t376) * t377; m(5) * (-t21 * t244 + t343 * t36 + t346 * t37) + m(6) * (-t13 * t244 + t14 * t343 + t15 * t346) + ((t242 * t45 + t342 * t67 + t345 * t68) * t376 + (t242 * t40 + t342 * t43 + t345 * t44) * t375) * t377; 0.4e1 * (t376 + t375) * (-0.1e1 + t388) * t300; (t232 * t321 + t299) * m(6); m(6) * (t39 * t299 - t176 * t85 - t178 * t86 + t202 * t61 + t204 * t62 + (t20 * t242 + t321 * t39) * t232); m(6) * (t40 * t299 + t204 * t14 + t202 * t15 - t176 * t44 - t178 * t43 + (t13 * t242 + t321 * t40) * t232); m(6) * ((-t176 * t238 - t178 * t239 - t298) * t242 + (t232 * t237 + (t202 * t238 + t204 * t239 - t232 * t244) * t244) * qJD(2)); (-t202 * t176 - t204 * t178 + (t232 * t300 + t237 * t319) * t232) * t378;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
