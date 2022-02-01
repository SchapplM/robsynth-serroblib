% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:09
% EndTime: 2022-01-23 09:14:15
% DurationCPUTime: 3.37s
% Computational Cost: add. (19489->287), mult. (11935->398), div. (0->0), fcn. (10892->9), ass. (0->200)
t232 = qJ(1) + pkin(8);
t225 = sin(t232);
t227 = cos(t232);
t231 = pkin(9) + qJ(4);
t228 = qJ(5) + t231;
t219 = sin(t228);
t220 = cos(t228);
t182 = rSges(6,1) * t219 + rSges(6,2) * t220;
t224 = sin(t231);
t337 = pkin(4) * t224;
t241 = t182 + t337;
t372 = t241 * t227;
t373 = t241 * t225;
t366 = t225 * t373 + t227 * t372;
t226 = cos(t231);
t189 = rSges(5,1) * t224 + rSges(5,2) * t226;
t222 = t225 ^ 2;
t223 = t227 ^ 2;
t273 = t222 + t223;
t371 = t273 * t189;
t319 = -m(5) * t371 / 0.2e1 - m(6) * t366 / 0.2e1;
t359 = m(6) / 0.2e1;
t360 = m(5) / 0.2e1;
t172 = t189 * t225;
t173 = t189 * t227;
t99 = t172 * t225 + t173 * t227;
t385 = t360 * t99;
t332 = t366 * t359 + t385;
t26 = t332 - t319;
t386 = t26 * qJD(1);
t384 = t182 * t273;
t208 = Icges(6,4) * t220;
t179 = -Icges(6,2) * t219 + t208;
t180 = Icges(6,1) * t219 + t208;
t383 = t180 + t179;
t221 = cos(pkin(9)) * pkin(3) + pkin(2);
t336 = pkin(4) * t226;
t196 = t221 + t336;
t234 = -pkin(6) - qJ(3);
t230 = pkin(7) - t234;
t206 = t230 * t227;
t305 = t219 * t225;
t275 = rSges(6,2) * t305 + t227 * rSges(6,3);
t326 = rSges(6,1) * t220;
t335 = sin(qJ(1)) * pkin(1);
t368 = t206 + (-t196 - t326) * t225 - t335 + t275;
t300 = t220 * t227;
t304 = t219 * t227;
t246 = rSges(6,1) * t300 - rSges(6,2) * t304 + t225 * rSges(6,3);
t338 = cos(qJ(1)) * pkin(1);
t378 = t227 * t196 + t230 * t225;
t95 = -t246 - t338 - t378;
t56 = -t225 * t95 + t227 * t368;
t342 = t225 / 0.2e1;
t340 = -t227 / 0.2e1;
t381 = t227 / 0.2e1;
t157 = t182 * t227;
t316 = Icges(5,4) * t224;
t185 = Icges(5,2) * t226 + t316;
t188 = Icges(5,1) * t226 - t316;
t379 = (t188 / 0.2e1 - t185 / 0.2e1) * t224;
t105 = (-t230 - t234) * t227 + (t196 - t221) * t225;
t207 = t225 * t234;
t369 = t207 + t378;
t106 = -t227 * t221 + t369;
t291 = t227 * t234;
t34 = (-t106 + t369) * t225 + (-t196 * t225 + t105 + t206 + t291) * t227;
t377 = m(6) * t34;
t376 = m(6) * qJD(5);
t183 = -rSges(6,2) * t219 + t326;
t250 = Icges(6,5) * t219 + Icges(6,6) * t220;
t150 = t250 * t225;
t151 = t227 * t250;
t315 = Icges(6,4) * t219;
t181 = Icges(6,1) * t220 - t315;
t139 = Icges(6,5) * t225 + t181 * t227;
t178 = Icges(6,2) * t220 + t315;
t281 = -t178 * t227 + t139;
t193 = Icges(6,4) * t305;
t301 = t220 * t225;
t138 = Icges(6,1) * t301 - Icges(6,5) * t227 - t193;
t282 = -Icges(6,2) * t301 + t138 - t193;
t137 = Icges(6,6) * t225 + t179 * t227;
t283 = -t180 * t227 - t137;
t136 = Icges(6,4) * t301 - Icges(6,2) * t305 - Icges(6,6) * t227;
t284 = t180 * t225 + t136;
t363 = (-t281 * t225 + t227 * t282) * t219 + (t283 * t225 + t227 * t284) * t220;
t333 = (-t222 * t151 + (t225 * t150 + t363) * t227) * t342 + (-t223 * t150 + (t227 * t151 + t363) * t225) * t340;
t156 = t182 * t225;
t370 = t225 * t156 + t227 * t157;
t82 = t225 * (rSges(6,1) * t301 - t275) + t227 * t246;
t48 = t105 * t225 + t106 * t227 + t82;
t6 = t333 + m(6) * (t366 * t183 - t370 * t48);
t375 = t6 * qJD(5);
t327 = rSges(5,1) * t226;
t267 = t221 + t327;
t298 = t224 * t225;
t274 = rSges(5,2) * t298 + t227 * rSges(5,3);
t101 = -t225 * t267 + t274 - t291 - t335;
t297 = t224 * t227;
t266 = -rSges(5,2) * t297 + t225 * rSges(5,3);
t102 = t227 * t267 - t207 + t266 + t338;
t374 = t101 * t227 + t102 * t225;
t215 = Icges(5,4) * t226;
t186 = -Icges(5,2) * t224 + t215;
t187 = Icges(5,1) * t224 + t215;
t257 = t383 * t220 / 0.2e1 + (t181 / 0.2e1 - t178 / 0.2e1) * t219;
t107 = t139 * t301;
t134 = Icges(6,5) * t301 - Icges(6,6) * t305 - Icges(6,3) * t227;
t177 = Icges(6,5) * t220 - Icges(6,6) * t219;
t309 = t177 * t227;
t135 = Icges(6,3) * t225 + t309;
t259 = t137 * t219 - t134;
t261 = t227 * t135 - t107;
t287 = t225 * t135 + t139 * t300;
t288 = -t225 * t134 - t138 * t300;
t310 = t136 * t219;
t58 = -t137 * t305 - t261;
t59 = -t136 * t304 - t288;
t60 = -t137 * t304 + t287;
t269 = ((t58 - t107 + (t135 + t310) * t227 + t288) * t227 + t287 * t225) * t340 + (t225 * t60 - t227 * t59) * t381 + ((t225 * t259 + t261 + t58 + t59) * t225 + (t225 * (-t138 * t220 + t310) - t287 + t60 + (t134 + t259) * t227) * t227) * t342;
t148 = Icges(5,5) * t225 + t188 * t227;
t277 = -t185 * t227 + t148;
t199 = Icges(5,4) * t298;
t296 = t226 * t225;
t147 = Icges(5,1) * t296 - Icges(5,5) * t227 - t199;
t278 = -Icges(5,2) * t296 + t147 - t199;
t146 = Icges(5,6) * t225 + t186 * t227;
t279 = -t187 * t227 - t146;
t145 = Icges(5,4) * t296 - Icges(5,2) * t298 - Icges(5,6) * t227;
t280 = t187 * t225 + t145;
t364 = (-t277 * t225 + t227 * t278) * t224 + (t279 * t225 + t227 * t280) * t226;
t362 = 0.4e1 * qJD(1);
t361 = 2 * qJD(4);
t320 = rSges(4,3) + qJ(3);
t358 = m(4) * ((t320 * t227 - t335) * t227 + (t225 * t320 + t338) * t225);
t357 = m(5) * (t101 * t172 - t102 * t173);
t356 = m(5) * t374;
t354 = t82 * t377;
t353 = t48 * t377;
t244 = t56 * t183;
t351 = m(6) * (-t156 * t372 + t157 * t373 - t244);
t350 = m(6) * (-t244 + (t225 * t372 - t227 * t373) * t182);
t349 = m(6) * (t368 * t373 + t372 * t95);
t348 = m(6) * (t156 * t368 + t157 * t95);
t347 = m(6) * t56;
t343 = -t225 / 0.2e1;
t328 = m(6) * qJD(4);
t299 = t224 * t145;
t295 = t226 * t227;
t25 = t34 * t359;
t290 = t25 * qJD(2);
t245 = t370 * t359;
t255 = m(6) * t384;
t67 = t245 + t255 / 0.2e1;
t289 = t67 * qJD(1);
t143 = Icges(5,5) * t296 - Icges(5,6) * t298 - Icges(5,3) * t227;
t286 = -t225 * t143 - t147 * t295;
t252 = Icges(5,5) * t226 - Icges(5,6) * t224;
t144 = Icges(5,3) * t225 + t227 * t252;
t285 = t225 * t144 + t148 * t295;
t268 = -t183 - t336;
t113 = t148 * t296;
t260 = t227 * t144 - t113;
t258 = t224 * t146 - t143;
t256 = t354 / 0.2e1 + t269;
t251 = -Icges(5,5) * t224 - Icges(5,6) * t226;
t243 = t25 * qJD(4);
t239 = (-t178 + t181) * t220 - t383 * t219;
t242 = -t269 + (t177 * t225 + t219 * t283 + t220 * t281 + t227 * t239) * t342 + (-t219 * t284 + t220 * t282 + t225 * t239 - t309) * t340;
t240 = -t257 + (t342 + t343) * (t136 * t220 + t138 * t219);
t190 = -rSges(5,2) * t224 + t327;
t167 = t251 * t227;
t166 = t251 * t225;
t133 = t268 * t227;
t131 = t268 * t225;
t87 = t370 * t376;
t72 = -t273 * t337 - t370;
t66 = t245 - t255 / 0.2e1;
t65 = -t146 * t297 + t285;
t64 = -t145 * t297 - t286;
t63 = -t146 * t298 - t260;
t44 = t225 * t65 - t227 * t64;
t43 = t225 * t63 - t227 * (-t225 * (-t226 * t147 + t299) - t227 * t143);
t36 = t257 + t348;
t32 = t350 / 0.2e1;
t30 = t351 / 0.2e1;
t27 = t319 + t332;
t24 = t347 + t356 + t358;
t23 = t25 * qJD(1);
t20 = (t187 / 0.2e1 + t186 / 0.2e1) * t226 + t379 + t357 + t349 + t257;
t18 = (t63 - t113 + (t144 + t299) * t227 + t286) * t227 + t285 * t225;
t17 = (t227 * t258 - t285 + t65) * t227 + (t225 * t258 + t260 + t64) * t225;
t8 = m(6) * (t183 * t384 - t370 * t82) + t333;
t7 = t8 * qJD(5);
t4 = t30 - t350 / 0.2e1 + t256;
t3 = t32 - t351 / 0.2e1 + t256;
t2 = t30 + t32 - t354 / 0.2e1 + t242;
t1 = t353 + (t44 / 0.2e1 - t18 / 0.2e1) * t227 + (t17 / 0.2e1 + t43 / 0.2e1) * t225 + t269;
t5 = [t24 * qJD(3) + t20 * qJD(4) + t36 * qJD(5), -t243, qJD(1) * t24 + qJD(4) * t27 + qJD(5) * t66, t20 * qJD(1) - t290 + t27 * qJD(3) + t2 * qJD(5) + ((-t131 * t95 + t133 * t368) * t359 + (-t374 * t190 + (-t172 * t227 + t173 * t225) * t189) * t360) * t361 + ((t224 * t279 + t226 * t277) * t342 + t18 * t381 - t353 + t242 + (t222 / 0.2e1 + t223 / 0.2e1) * t252 + (t17 + t43) * t343 + (-t224 * t280 + t226 * t278 + t44) * t340) * qJD(4), t36 * qJD(1) + t66 * qJD(3) + t2 * qJD(4) + t242 * qJD(5) + (-t244 + (-t156 * t227 + t157 * t225) * t182) * t376; t243, 0, 0, t23 + (t72 * t359 - t385) * t361 - t87, -t328 * t370 - t87; t26 * qJD(4) + t67 * qJD(5) + (-t347 / 0.4e1 - t356 / 0.4e1 - t358 / 0.4e1) * t362, 0, 0, t386 + (-t131 * t227 + t133 * t225) * t328, t289; (t240 - (t187 + t186) * t226 / 0.2e1 - t379) * qJD(1) + t290 - t26 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t349 / 0.4e1 - t357 / 0.4e1) * t362, t23, -t386, t1 * qJD(1) + (m(5) * (t190 * t371 - t99 * (t225 * (rSges(5,1) * t296 - t274) + t227 * (rSges(5,1) * t295 + t266))) + (t222 * t167 + (-t225 * t166 + t364) * t227) * t342 + (t223 * t166 + (-t227 * t167 + t364) * t225) * t340 + m(6) * (-t131 * t373 - t133 * t372 + t48 * t72) + t333) * qJD(4) + t375, t4 * qJD(1) + t6 * qJD(4) + t375; (t240 - t348) * qJD(1) - t67 * qJD(3) + t3 * qJD(4) + t269 * qJD(5), 0, -t289, t3 * qJD(1) + ((t72 * t82 + (-t131 * t225 - t133 * t227) * t182) * m(6) + t333) * qJD(4) + t7, qJD(1) * t269 + qJD(4) * t8 + t7;];
Cq = t5;
