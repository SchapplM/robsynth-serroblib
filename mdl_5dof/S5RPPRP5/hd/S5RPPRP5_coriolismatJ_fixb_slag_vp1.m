% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:22
% EndTime: 2019-12-31 17:53:29
% DurationCPUTime: 4.29s
% Computational Cost: add. (6486->307), mult. (16548->439), div. (0->0), fcn. (18228->6), ass. (0->191)
t246 = cos(pkin(7));
t319 = sin(pkin(7));
t338 = sin(qJ(4));
t339 = cos(qJ(4));
t232 = -t246 * t338 + t319 * t339;
t247 = sin(qJ(1));
t219 = t232 * t247;
t248 = cos(qJ(1));
t221 = t232 * t248;
t240 = t248 * qJ(2);
t276 = t319 * qJ(3) + pkin(1);
t375 = (pkin(2) + pkin(3)) * t246 + t276;
t254 = -t375 * t247 + t240;
t231 = t246 * t339 + t319 * t338;
t220 = t231 * t247;
t340 = rSges(6,1) + pkin(4);
t399 = rSges(6,3) + qJ(5);
t265 = t399 * t219 - t340 * t220;
t410 = (-rSges(6,2) - pkin(6)) * t248 + t254 + t265;
t252 = t375 * t248 + (-pkin(6) + qJ(2)) * t247;
t222 = t231 * t248;
t377 = t247 * rSges(6,2) + t399 * t221 - t340 * t222;
t74 = t252 - t377;
t33 = t410 * t219 - t74 * t221;
t414 = m(6) * t33 * qJD(1);
t413 = Icges(5,1) + Icges(6,1);
t412 = Icges(5,2) + Icges(6,3);
t196 = rSges(5,1) * t232 - rSges(5,2) * t231;
t384 = t247 ^ 2 + t248 ^ 2;
t107 = t384 * t196;
t368 = m(6) / 0.2e1;
t369 = m(5) / 0.2e1;
t292 = t399 * t231 + t340 * t232;
t106 = t292 * t247;
t109 = t292 * t248;
t389 = t106 * t247 + t109 * t248;
t320 = t107 * t369 + t389 * t368;
t175 = rSges(5,1) * t219 - rSges(5,2) * t220;
t177 = rSges(5,1) * t221 - rSges(5,2) * t222;
t267 = -t247 * t175 - t177 * t248;
t95 = -t340 * t219 - t220 * t399;
t96 = t340 * t221 + t399 * t222;
t274 = t247 * t95 - t248 * t96;
t327 = t267 * t369 + t274 * t368;
t14 = t327 - t320;
t411 = t14 * qJD(1);
t283 = t248 * t319;
t284 = t247 * t319;
t400 = m(6) * (t219 * t283 - t221 * t284);
t123 = -t400 / 0.2e1;
t124 = t400 / 0.2e1;
t409 = Icges(6,4) + Icges(5,5);
t408 = Icges(5,6) - Icges(6,6);
t316 = Icges(6,5) * t221;
t151 = Icges(6,1) * t222 - Icges(6,4) * t247 - t316;
t210 = Icges(5,4) * t221;
t154 = Icges(5,1) * t222 - Icges(5,5) * t247 + t210;
t407 = t412 * t222 - t151 - t154 - t210 + t316;
t205 = Icges(6,5) * t219;
t149 = Icges(6,1) * t220 + Icges(6,4) * t248 - t205;
t208 = Icges(5,4) * t219;
t153 = -Icges(5,1) * t220 - Icges(5,5) * t248 - t208;
t406 = t412 * t220 - t149 + t153 + t205 - t208;
t207 = Icges(6,5) * t222;
t139 = -Icges(6,6) * t247 - Icges(6,3) * t221 + t207;
t318 = Icges(5,4) * t222;
t148 = Icges(5,2) * t221 - Icges(5,6) * t247 + t318;
t405 = t413 * t221 + t139 - t148 + t207 - t318;
t206 = Icges(6,5) * t220;
t138 = -Icges(6,6) * t248 + Icges(6,3) * t219 - t206;
t209 = Icges(5,4) * t220;
t146 = Icges(5,2) * t219 + Icges(5,6) * t248 + t209;
t404 = t413 * t219 - t138 - t146 + t206 - t209;
t403 = t219 * t138 + t220 * t149;
t310 = t138 * t221 + t222 * t149;
t402 = t219 * t146 - t220 * t153;
t308 = t221 * t146 - t153 * t222;
t157 = t247 * t219 + t248 * t221;
t390 = m(6) * t157;
t135 = -t390 / 0.2e1;
t134 = t390 / 0.2e1;
t142 = Icges(5,5) * t222 + Icges(5,6) * t221 - Icges(5,3) * t247;
t307 = t221 * t148 + t222 * t154;
t277 = t247 * t142 - t307;
t140 = Icges(5,5) * t220 + Icges(5,6) * t219 + Icges(5,3) * t248;
t42 = t140 * t248 + t402;
t398 = t277 - t42;
t145 = Icges(6,4) * t222 - Icges(6,2) * t247 - Icges(6,6) * t221;
t309 = -t221 * t139 + t222 * t151;
t278 = t247 * t145 - t309;
t143 = Icges(6,4) * t220 + Icges(6,2) * t248 - Icges(6,6) * t219;
t40 = t143 * t248 + t403;
t397 = t278 - t40;
t385 = -t220 * rSges(5,1) - t219 * rSges(5,2);
t394 = (-rSges(5,3) - pkin(6)) * t248 + t254 + t385;
t160 = -t220 * t248 + t222 * t247;
t289 = -m(6) * t160 / 0.2e1;
t324 = m(6) * qJD(5);
t260 = t384 * t319;
t251 = (-t157 * t246 + t260 * t231) * t368;
t253 = m(6) * (t220 * t284 + t222 * t283 - t232 * t246);
t58 = t251 - t253 / 0.2e1;
t79 = -t219 * t220 - t221 * t222 + t231 * t232;
t392 = -qJD(2) * t289 - t58 * qJD(3) + t79 * t324;
t388 = -t219 * t139 + t220 * t151;
t387 = t219 * t148 + t220 * t154;
t386 = -t409 * t231 - t408 * t232;
t383 = -t407 * t247 + t406 * t248;
t382 = -t405 * t247 + t404 * t248;
t381 = (t409 * t219 - t408 * t220) * t248 + (-t409 * t221 + t408 * t222) * t247;
t229 = Icges(6,5) * t232;
t180 = Icges(6,3) * t231 + t229;
t317 = Icges(5,4) * t232;
t186 = -Icges(5,2) * t231 + t317;
t187 = -Icges(6,1) * t231 + t229;
t189 = -Icges(5,1) * t231 - t317;
t379 = -t186 + t189 + t180 + t187;
t315 = Icges(6,5) * t231;
t179 = Icges(6,3) * t232 - t315;
t230 = Icges(5,4) * t231;
t185 = -Icges(5,2) * t232 - t230;
t188 = Icges(6,1) * t232 + t315;
t190 = Icges(5,1) * t232 - t230;
t378 = -t185 - t190 + t179 - t188;
t328 = (t95 * t283 + t96 * t284) * t368 + (-t175 * t283 + t177 * t284) * t369;
t376 = t319 * rSges(4,3) + (rSges(4,1) + pkin(2)) * t246 + t276;
t374 = (t190 / 0.2e1 + t185 / 0.2e1 + t188 / 0.2e1 - t179 / 0.2e1) * t231 - (t189 / 0.2e1 - t186 / 0.2e1 + t187 / 0.2e1 + t180 / 0.2e1) * t232;
t373 = -0.2e1 * t260;
t371 = 0.4e1 * qJD(1);
t370 = 2 * qJD(4);
t161 = t248 * rSges(4,2) - t376 * t247 + t240;
t162 = (rSges(4,2) + qJ(2)) * t247 + t376 * t248;
t366 = m(4) * (-t161 * t284 + t162 * t283);
t365 = m(4) * (t161 * t248 + t162 * t247);
t264 = t222 * rSges(5,1) + t221 * rSges(5,2) - t247 * rSges(5,3);
t99 = t252 + t264;
t364 = m(5) * (-t175 * t394 + t177 * t99);
t363 = m(5) * (t99 * t283 - t284 * t394);
t53 = t99 * t247 + t394 * t248;
t362 = m(5) * t53;
t357 = m(6) * (-t219 * t96 + t220 * t74 - t221 * t95 + t222 * t410);
t35 = t74 * t247 + t410 * t248;
t356 = m(6) * (-t106 * t221 + t219 * t109 + t35 * t231);
t355 = m(6) * (t410 * t95 + t74 * t96);
t353 = m(6) * (t74 * t283 - t284 * t410);
t352 = m(6) * t35;
t349 = m(6) * (t106 * t283 - t109 * t284);
t345 = -t247 / 0.2e1;
t342 = t248 / 0.2e1;
t337 = m(3) * ((rSges(3,2) * t284 + t248 * rSges(3,3) + t240) * t248 + (-rSges(3,2) * t283 + (rSges(3,3) + qJ(2)) * t247) * t247);
t325 = m(6) * qJD(4);
t77 = (t368 + t369 + m(4) / 0.2e1) * t373;
t311 = t77 * qJD(1);
t293 = -t340 * t231 + t399 * t232;
t193 = -rSges(5,1) * t231 - rSges(5,2) * t232;
t126 = t160 * t325 / 0.2e1;
t108 = t293 * t248;
t105 = t293 * t247;
t76 = (m(6) / 0.4e1 + m(5) / 0.4e1 + m(4) / 0.4e1) * t373 + (m(6) + m(5) + m(4)) * t260 / 0.2e1;
t65 = 0.2e1 * t135;
t64 = t134 + t135;
t63 = 0.2e1 * t134;
t59 = t349 / 0.2e1;
t57 = t251 + t253 / 0.2e1;
t56 = 0.2e1 * t124;
t55 = t124 + t123;
t54 = 0.2e1 * t123;
t37 = (-t248 * rSges(6,2) + t265) * t247 + t377 * t248;
t30 = t247 * t277 + t248 * (-t140 * t247 + t308);
t29 = t247 * t278 + t248 * (-t143 * t247 + t310);
t28 = -(t248 * t142 + t387) * t247 + t248 * t42;
t27 = -(t248 * t145 + t388) * t247 + t248 * t40;
t26 = t37 * t157 + t389 * t231;
t23 = t356 / 0.2e1;
t21 = t357 / 0.2e1;
t20 = t353 + t363 + t366;
t17 = t337 + t352 + t362 + t365;
t15 = t320 + t327;
t12 = t59 - t328;
t11 = -t349 / 0.2e1 + t328;
t10 = t59 + t328;
t9 = t355 + t364 - t374;
t8 = t308 * t248 + (t398 + t402) * t247;
t7 = t310 * t248 + (t397 + t403) * t247;
t6 = (t307 + t398) * t248 + t387 * t247;
t5 = (t309 + t397) * t248 + t388 * t247;
t4 = t21 - t356 / 0.2e1;
t3 = t23 + t21;
t2 = t23 - t357 / 0.2e1;
t1 = (-t30 / 0.2e1 + t8 / 0.2e1 - t29 / 0.2e1 + t7 / 0.2e1) * t248 + (-t6 / 0.2e1 - t28 / 0.2e1 - t5 / 0.2e1 - t27 / 0.2e1) * t247;
t13 = [t17 * qJD(2) + t20 * qJD(3) + t9 * qJD(4) + t33 * t324, qJD(1) * t17 + qJD(3) * t76 + qJD(4) * t15 + qJD(5) * t64, qJD(1) * t20 + qJD(2) * t76 + qJD(4) * t10 + qJD(5) * t55, t9 * qJD(1) + t15 * qJD(2) + t10 * qJD(3) + t3 * qJD(5) + ((t105 * t74 + t106 * t96 + t108 * t410 + t109 * t95) * t368 + ((-t175 * t248 + t177 * t247) * t196 + t53 * t193) * t369) * t370 + ((-t378 * t221 + t379 * t222 + t407 * t231 + t405 * t232 - t386 * t247) * t345 + (t6 + t28 + t5 + t27) * t247 / 0.2e1 + (-t378 * t219 + t379 * t220 + t406 * t231 + t404 * t232 + t386 * t248 + t29 + t30) * t342 - (t8 + t7) * t248 / 0.2e1) * qJD(4), t64 * qJD(2) + t55 * qJD(3) + t3 * qJD(4) + t414; t77 * qJD(3) + t14 * qJD(4) + t63 * qJD(5) + (-t352 / 0.4e1 - t362 / 0.4e1 - t337 / 0.4e1 - t365 / 0.4e1) * t371, 0, t311, t411 + 0.2e1 * ((-t105 * t248 + t108 * t247) * qJD(4) / 0.2e1 + t160 * qJD(5) / 0.4e1) * m(6), qJD(1) * t63 + t126; -t77 * qJD(2) + t11 * qJD(4) + t56 * qJD(5) + (-t353 / 0.4e1 - t363 / 0.4e1 - t366 / 0.4e1) * t371, -t311, 0, t11 * qJD(1) + ((t193 * t260 - t246 * t267) * t369 + (t105 * t284 + t108 * t283 - t246 * t274) * t368) * t370 + t57 * qJD(5), qJD(1) * t56 + qJD(4) * t57; -t14 * qJD(2) + t12 * qJD(3) + t1 * qJD(4) + t2 * qJD(5) + (-t355 / 0.4e1 - t364 / 0.4e1) * t371 + t374 * qJD(1), qJD(5) * t289 - t411, qJD(1) * t12 + qJD(5) * t58, t1 * qJD(1) + (m(5) * (t193 * t107 + t267 * (-t247 * (t248 * rSges(5,3) - t385) - t248 * t264)) + m(6) * (t106 * t105 + t109 * t108 + t274 * t37) + (-t383 * t221 + t382 * t222 - t381 * t247) * t345 + (-t383 * t219 + t382 * t220 + t381 * t248) * t342) * qJD(4) + t26 * t324, t2 * qJD(1) + t26 * t325 - t392; t65 * qJD(2) + t54 * qJD(3) + t4 * qJD(4) - t414, qJD(1) * t65 + t126, qJD(1) * t54 - qJD(4) * t58, t4 * qJD(1) + (-t219 * t105 + t220 * t106 - t221 * t108 + t222 * t109 + t231 * t274 + t232 * t37 - t26) * t325 + t392, t79 * t325;];
Cq = t13;
