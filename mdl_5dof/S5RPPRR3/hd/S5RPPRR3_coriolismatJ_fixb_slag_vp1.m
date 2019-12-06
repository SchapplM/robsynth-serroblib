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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:41:38
% EndTime: 2019-12-05 17:41:46
% DurationCPUTime: 3.26s
% Computational Cost: add. (18915->288), mult. (11629->395), div. (0->0), fcn. (10604->9), ass. (0->190)
t226 = pkin(9) + qJ(4);
t223 = qJ(5) + t226;
t214 = sin(t223);
t215 = cos(t223);
t179 = rSges(6,1) * t214 + rSges(6,2) * t215;
t219 = sin(t226);
t227 = qJ(1) + pkin(8);
t220 = sin(t227);
t282 = t219 * t220;
t203 = pkin(4) * t282;
t127 = t179 * t220 + t203;
t340 = m(6) / 0.2e1;
t221 = cos(t226);
t186 = rSges(5,1) * t219 + rSges(5,2) * t221;
t217 = t220 ^ 2;
t222 = cos(t227);
t218 = t222 ^ 2;
t256 = t217 + t218;
t356 = t256 * t186;
t318 = pkin(4) * t219;
t244 = (t179 + t318) * t222;
t359 = t222 * t244;
t304 = -m(5) * t356 / 0.2e1 + (-t127 * t220 - t359) * t340;
t285 = t215 * t220;
t289 = t214 * t220;
t153 = rSges(6,1) * t289 + rSges(6,2) * t285;
t122 = t203 + t153;
t341 = m(5) / 0.2e1;
t169 = t186 * t220;
t170 = t186 * t222;
t97 = t169 * t220 + t170 * t222;
t372 = t341 * t97;
t314 = (t122 * t220 + t359) * t340 + t372;
t25 = t314 - t304;
t373 = t25 * qJD(1);
t284 = t215 * t222;
t288 = t214 * t222;
t132 = Icges(6,5) * t284 - Icges(6,6) * t288 + Icges(6,3) * t220;
t208 = Icges(6,4) * t215;
t351 = Icges(6,2) * t214 - t208;
t133 = Icges(6,6) * t222 + t351 * t220;
t174 = Icges(6,5) * t215 - Icges(6,6) * t214;
t295 = t174 * t220;
t272 = t222 * (Icges(6,3) * t222 - t295) + t133 * t289;
t134 = Icges(6,4) * t284 - Icges(6,2) * t288 + Icges(6,6) * t220;
t191 = Icges(6,4) * t288;
t136 = Icges(6,1) * t284 + Icges(6,5) * t220 - t191;
t357 = (t134 * t214 - t136 * t215) * t222;
t371 = t220 * t132 + t272 - t357;
t370 = t133 * t288;
t369 = t179 * t256;
t368 = -t220 / 0.2e1;
t324 = t220 / 0.2e1;
t367 = -t221 / 0.2e1;
t321 = t222 / 0.2e1;
t366 = (Icges(5,1) * t221 / 0.2e1 - Icges(5,4) * t219 + Icges(5,2) * t367) * t219;
t265 = -Icges(6,2) * t284 + t136 - t191;
t190 = Icges(6,4) * t289;
t135 = -Icges(6,1) * t285 + Icges(6,5) * t222 + t190;
t266 = Icges(6,2) * t285 + t135 + t190;
t354 = Icges(6,1) * t214 + t208;
t267 = t222 * t354 + t134;
t268 = -t220 * t354 + t133;
t365 = -(t265 * t220 + t266 * t222) * t214 - (t267 * t220 + t268 * t222) * t215;
t281 = t219 * t222;
t200 = Icges(5,4) * t281;
t276 = t221 * t222;
t145 = Icges(5,1) * t276 + Icges(5,5) * t220 - t200;
t261 = -Icges(5,2) * t276 + t145 - t200;
t199 = Icges(5,4) * t282;
t279 = t220 * t221;
t144 = -Icges(5,1) * t279 + Icges(5,5) * t222 + t199;
t262 = Icges(5,2) * t279 + t144 + t199;
t143 = Icges(5,4) * t276 - Icges(5,2) * t281 + Icges(5,6) * t220;
t210 = Icges(5,4) * t221;
t353 = Icges(5,1) * t219 + t210;
t263 = t222 * t353 + t143;
t350 = Icges(5,2) * t219 - t210;
t142 = Icges(5,6) * t222 + t350 * t220;
t264 = -t220 * t353 + t142;
t364 = -(t261 * t220 + t262 * t222) * t219 - (t263 * t220 + t264 * t222) * t221;
t310 = rSges(6,1) * t215;
t180 = -rSges(6,2) * t214 + t310;
t292 = t180 * t222;
t293 = t180 * t220;
t239 = Icges(6,5) * t214 + Icges(6,6) * t215;
t147 = t220 * t239;
t148 = t239 * t222;
t315 = (t218 * t147 + (-t222 * t148 - t365) * t220) * t321 + (-t217 * t148 + (t220 * t147 + t365) * t222) * t324;
t154 = t179 * t222;
t355 = t220 * t153 + t222 * t154;
t245 = -rSges(6,1) * t284 + rSges(6,2) * t288;
t121 = t222 * (t220 * rSges(6,3) - t245);
t258 = -rSges(6,2) * t289 - t222 * rSges(6,3);
t137 = -rSges(6,1) * t285 - t258;
t216 = cos(pkin(9)) * pkin(3) + pkin(2);
t317 = pkin(4) * t221;
t196 = t216 + t317;
t173 = t222 * t196;
t229 = -pkin(6) - qJ(3);
t207 = t220 * t229;
t360 = t222 * t229;
t48 = -t222 * (t222 * t216 - t173 - t207) + t121 + (-t137 + (t196 - t216) * t220 - t360) * t220;
t6 = t315 + m(6) * (t127 * t293 + t244 * t292 - t355 * t48);
t363 = t6 * qJD(5);
t362 = -rSges(4,3) - qJ(3);
t311 = rSges(5,1) * t221;
t252 = t216 + t311;
t257 = -rSges(5,2) * t282 - t222 * rSges(5,3);
t319 = sin(qJ(1)) * pkin(1);
t100 = t252 * t220 + t257 + t319 + t360;
t251 = rSges(5,2) * t281 - t220 * rSges(5,3);
t316 = cos(qJ(1)) * pkin(1);
t101 = -t252 * t222 + t207 + t251 - t316;
t361 = t222 * t100 + t101 * t220;
t277 = t221 * t145;
t358 = (t219 * t143 - t277) * t222;
t300 = Icges(6,4) * t214;
t175 = Icges(6,2) * t215 + t300;
t178 = Icges(6,1) * t215 - t300;
t348 = (-t351 + t354) * t214 + (t175 - t178) * t215;
t248 = (t354 / 0.2e1 - t351 / 0.2e1) * t215 + (t178 / 0.2e1 - t175 / 0.2e1) * t214;
t250 = -t135 * t215 - t132;
t296 = t133 * t214;
t57 = -t135 * t285 + t272;
t58 = t222 * t132 + t134 * t289 - t136 * t285;
t5 = (t220 * t58 + t222 * t57) * t368 + ((t357 + t371) * t222 + ((t250 - t296) * t222 + t58 + t370) * t220) * t324 + ((t58 + (-t132 + t296) * t222 - t370) * t222 + (t250 * t220 + t371 - t57) * t220) * t321;
t343 = 0.4e1 * qJD(1);
t342 = 2 * qJD(4);
t339 = m(4) * (-t222 * (t362 * t222 + t319) + (-t362 * t220 + t316) * t220);
t338 = m(5) * (-t100 * t169 + t101 * t170);
t337 = m(5) * t361;
t225 = -pkin(7) + t229;
t92 = t319 + t222 * t225 + (t196 + t310) * t220 + t258;
t93 = -t316 - t173 + (-rSges(6,3) + t225) * t220 + t245;
t246 = t92 * t292 + t93 * t293;
t332 = m(6) * (t127 * t154 - t153 * t244 + t246);
t331 = m(6) * ((-t122 * t222 + t220 * t244) * t179 + t246);
t330 = m(6) * (-t122 * t92 + t244 * t93);
t329 = m(6) * (-t153 * t92 + t154 * t93);
t328 = m(6) * (-t220 * t93 - t222 * t92);
t322 = -t222 / 0.2e1;
t313 = m(6) * qJD(4);
t312 = m(6) * qJD(5);
t283 = t219 * t142;
t278 = t221 * t144;
t234 = m(6) * t355;
t235 = -m(6) * t369 / 0.2e1;
t66 = -t234 / 0.2e1 + t235;
t273 = t66 * qJD(1);
t241 = Icges(5,5) * t221 - Icges(5,6) * t219;
t140 = Icges(5,3) * t222 - t241 * t220;
t270 = t222 * t140 + t142 * t282;
t269 = t220 * t140 + t144 * t276;
t253 = t180 + t317;
t141 = Icges(5,5) * t276 - Icges(5,6) * t281 + Icges(5,3) * t220;
t249 = -t141 - t278;
t240 = Icges(5,5) * t219 + Icges(5,6) * t221;
t63 = t222 * t141 + t143 * t282 - t220 * t277;
t233 = -t5 + (-t267 * t214 + t265 * t215 - t348 * t222 + t295) * t324 + (t174 * t222 - t268 * t214 + t266 * t215 + t348 * t220) * t321;
t232 = -t248 + (t321 + t322) * (t134 * t215 + t136 * t214);
t187 = -rSges(5,2) * t219 + t311;
t164 = t240 * t222;
t163 = t220 * t240;
t130 = t253 * t222;
t128 = t253 * t220;
t88 = t355 * t312;
t84 = -t220 * t137 + t121;
t72 = -t256 * t318 - t355;
t67 = t234 / 0.2e1 + t235;
t65 = t220 * t141 - t358;
t64 = -t142 * t281 + t269;
t62 = -t220 * t278 + t270;
t42 = t220 * t65 + t222 * t64;
t41 = t220 * t63 + t222 * t62;
t32 = t248 + t329;
t29 = t331 / 0.2e1;
t28 = t332 / 0.2e1;
t26 = t304 + t314;
t23 = t328 - t337 + t339;
t19 = (t353 / 0.2e1 - t350 / 0.2e1) * t221 + t366 + t338 + t330 + t248;
t16 = (t63 + (-t141 + t283) * t222 - t269) * t222 + (t249 * t220 + t270 - t62) * t220;
t15 = (t270 + t65 + t358) * t222 + (-t64 + (t249 - t283) * t222 + t63 + t269) * t220;
t8 = m(6) * (t180 * t369 - t355 * t84) + t315;
t7 = t8 * qJD(5);
t4 = t28 - t331 / 0.2e1 + t5;
t3 = t29 - t332 / 0.2e1 + t5;
t2 = t28 + t29 + t233;
t1 = (t16 / 0.2e1 + t42 / 0.2e1) * t222 + (-t41 / 0.2e1 + t15 / 0.2e1) * t220 + t5;
t9 = [qJD(3) * t23 + qJD(4) * t19 + qJD(5) * t32, 0, qJD(1) * t23 + qJD(4) * t26 + qJD(5) * t67, t19 * qJD(1) + t26 * qJD(3) + t2 * qJD(5) + ((t361 * t187 + (-t169 * t222 + t170 * t220) * t186) * t341 + (t128 * t93 + t130 * t92 + (-t122 + t127) * t244) * t340) * t342 + (t15 * t368 + (-t264 * t219 + t262 * t221) * t321 + t233 + (t218 / 0.2e1 + t217 / 0.2e1) * t241 + (-t263 * t219 + t261 * t221 + t41) * t324 + (t16 + t42) * t322) * qJD(4), t32 * qJD(1) + t67 * qJD(3) + t2 * qJD(4) + t233 * qJD(5) + ((-t153 * t222 + t154 * t220) * t179 + t246) * t312; 0, 0, 0, (t72 * t340 - t372) * t342 - t88, -t313 * t355 - t88; t25 * qJD(4) - t66 * qJD(5) + (-t339 / 0.4e1 + t337 / 0.4e1 - t328 / 0.4e1) * t343, 0, 0, t373 + (t128 * t222 - t130 * t220) * t313, -t273; (t232 + (-t350 + t353) * t367 - t366) * qJD(1) - t25 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t330 / 0.4e1 - t338 / 0.4e1) * t343, 0, -t373, t1 * qJD(1) + (m(5) * (t187 * t356 - (t222 * (rSges(5,1) * t276 - t251) - t220 * (-rSges(5,1) * t279 - t257)) * t97) + (t218 * t163 + (-t222 * t164 - t364) * t220) * t321 + (-t217 * t164 + (t220 * t163 + t364) * t222) * t324 + m(6) * (t127 * t128 + t130 * t244 + t48 * t72) + t315) * qJD(4) + t363, t4 * qJD(1) + t6 * qJD(4) + t363; (t232 - t329) * qJD(1) + t66 * qJD(3) + t3 * qJD(4) + t5 * qJD(5), 0, t273, t3 * qJD(1) + ((t72 * t84 + (t128 * t220 + t130 * t222) * t179) * m(6) + t315) * qJD(4) + t7, qJD(1) * t5 + qJD(4) * t8 + t7;];
Cq = t9;
