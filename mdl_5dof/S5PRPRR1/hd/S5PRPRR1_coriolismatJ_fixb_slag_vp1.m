% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:39
% EndTime: 2019-12-05 15:42:47
% DurationCPUTime: 3.18s
% Computational Cost: add. (18813->267), mult. (11519->380), div. (0->0), fcn. (10500->7), ass. (0->185)
t224 = pkin(8) + qJ(2);
t218 = sin(t224);
t220 = cos(t224);
t223 = pkin(9) + qJ(4);
t221 = qJ(5) + t223;
t212 = sin(t221);
t213 = cos(t221);
t176 = rSges(6,1) * t212 + rSges(6,2) * t213;
t217 = sin(t223);
t323 = pkin(4) * t217;
t228 = t176 + t323;
t355 = t228 * t220;
t356 = t228 * t218;
t351 = t218 * t356 + t220 * t355;
t219 = cos(t223);
t183 = rSges(5,1) * t217 + rSges(5,2) * t219;
t215 = t218 ^ 2;
t216 = t220 ^ 2;
t262 = t215 + t216;
t354 = t262 * t183;
t306 = -m(5) * t354 / 0.2e1 - m(6) * t351 / 0.2e1;
t344 = m(6) / 0.2e1;
t345 = m(5) / 0.2e1;
t168 = t183 * t218;
t169 = t183 * t220;
t95 = t168 * t218 + t169 * t220;
t366 = t345 * t95;
t318 = t351 * t344 + t366;
t23 = t318 - t306;
t367 = t23 * qJD(2);
t365 = t176 * t262;
t202 = Icges(6,4) * t213;
t173 = -Icges(6,2) * t212 + t202;
t174 = Icges(6,1) * t212 + t202;
t364 = t173 + t174;
t327 = t218 / 0.2e1;
t325 = -t220 / 0.2e1;
t362 = t220 / 0.2e1;
t303 = Icges(5,4) * t217;
t179 = Icges(5,2) * t219 + t303;
t182 = Icges(5,1) * t219 - t303;
t361 = (t182 / 0.2e1 - t179 / 0.2e1) * t217;
t360 = m(6) * qJD(5);
t315 = rSges(6,1) * t213;
t177 = -rSges(6,2) * t212 + t315;
t236 = Icges(6,5) * t212 + Icges(6,6) * t213;
t146 = t236 * t218;
t147 = t220 * t236;
t302 = Icges(6,4) * t212;
t175 = Icges(6,1) * t213 - t302;
t135 = Icges(6,5) * t218 + t175 * t220;
t172 = Icges(6,2) * t213 + t302;
t269 = -t172 * t220 + t135;
t293 = t212 * t218;
t189 = Icges(6,4) * t293;
t289 = t213 * t218;
t134 = Icges(6,1) * t289 - Icges(6,5) * t220 - t189;
t270 = -Icges(6,2) * t289 + t134 - t189;
t133 = Icges(6,6) * t218 + t173 * t220;
t271 = -t174 * t220 - t133;
t132 = Icges(6,4) * t289 - Icges(6,2) * t293 - Icges(6,6) * t220;
t272 = t174 * t218 + t132;
t348 = (-t269 * t218 + t270 * t220) * t212 + (t271 * t218 + t272 * t220) * t213;
t319 = (-t215 * t147 + (t218 * t146 + t348) * t220) * t327 + (-t216 * t146 + (t220 * t147 + t348) * t218) * t325;
t152 = t176 * t218;
t153 = t176 * t220;
t353 = t218 * t152 + t220 * t153;
t214 = cos(pkin(9)) * pkin(3) + pkin(2);
t322 = pkin(4) * t219;
t192 = t214 + t322;
t320 = -pkin(6) - qJ(3);
t261 = -pkin(7) + t320;
t199 = t218 * t261;
t201 = t218 * t320;
t258 = t220 * t320;
t264 = -t218 * t192 - t220 * t261;
t136 = rSges(6,1) * t289 - rSges(6,2) * t293 - t220 * rSges(6,3);
t292 = t212 * t220;
t252 = -rSges(6,2) * t292 + rSges(6,3) * t218;
t288 = t213 * t220;
t82 = t218 * t136 + t220 * (rSges(6,1) * t288 + t252);
t47 = -t218 * (t218 * t214 + t258 + t264) + (-t199 + t201 + (t192 - t214) * t220) * t220 + t82;
t6 = t319 + m(6) * (t351 * t177 - t353 * t47);
t359 = t6 * qJD(5);
t93 = -t136 + t264;
t94 = -t199 + (t192 + t315) * t220 + t252;
t358 = t94 * t218 + t220 * t93;
t316 = rSges(5,1) * t219;
t254 = t214 + t316;
t286 = t217 * t218;
t263 = rSges(5,2) * t286 + t220 * rSges(5,3);
t98 = -t254 * t218 - t258 + t263;
t285 = t217 * t220;
t253 = -rSges(5,2) * t285 + rSges(5,3) * t218;
t99 = t254 * t220 - t201 + t253;
t357 = t99 * t218 + t220 * t98;
t209 = Icges(5,4) * t219;
t180 = -Icges(5,2) * t217 + t209;
t181 = Icges(5,1) * t217 + t209;
t243 = t364 * t213 / 0.2e1 + (t175 / 0.2e1 - t172 / 0.2e1) * t212;
t103 = t135 * t289;
t130 = Icges(6,5) * t289 - Icges(6,6) * t293 - Icges(6,3) * t220;
t171 = Icges(6,5) * t213 - Icges(6,6) * t212;
t297 = t171 * t220;
t131 = Icges(6,3) * t218 + t297;
t245 = t133 * t212 - t130;
t247 = t220 * t131 - t103;
t276 = t218 * t131 + t135 * t288;
t277 = -t218 * t130 - t134 * t288;
t298 = t132 * t212;
t58 = -t133 * t293 - t247;
t59 = -t132 * t292 - t277;
t60 = -t133 * t292 + t276;
t256 = ((t58 - t103 + (t131 + t298) * t220 + t277) * t220 + t276 * t218) * t325 + (t218 * t60 - t220 * t59) * t362 + ((t245 * t218 + t247 + t58 + t59) * t218 + ((-t134 * t213 + t298) * t218 - t276 + t60 + (t130 + t245) * t220) * t220) * t327;
t144 = Icges(5,5) * t218 + t182 * t220;
t265 = -t179 * t220 + t144;
t195 = Icges(5,4) * t286;
t284 = t218 * t219;
t143 = Icges(5,1) * t284 - Icges(5,5) * t220 - t195;
t266 = -Icges(5,2) * t284 + t143 - t195;
t142 = Icges(5,6) * t218 + t180 * t220;
t267 = -t181 * t220 - t142;
t141 = Icges(5,4) * t284 - Icges(5,2) * t286 - Icges(5,6) * t220;
t268 = t181 * t218 + t141;
t349 = (-t265 * t218 + t220 * t266) * t217 + (t267 * t218 + t220 * t268) * t219;
t347 = 0.4e1 * qJD(2);
t346 = 2 * qJD(4);
t343 = m(4) * t262 * (rSges(4,3) + qJ(3));
t342 = m(5) * (t168 * t98 - t169 * t99);
t341 = m(5) * t357;
t231 = t358 * t177;
t336 = m(6) * (-t152 * t355 + t153 * t356 - t231);
t335 = m(6) * (-t231 + (t218 * t355 - t220 * t356) * t176);
t334 = m(6) * (-t355 * t94 + t356 * t93);
t333 = m(6) * (t152 * t93 - t153 * t94);
t332 = m(6) * t358;
t328 = -t218 / 0.2e1;
t317 = m(6) * qJD(4);
t287 = t217 * t141;
t283 = t219 * t220;
t232 = t353 * t344;
t241 = m(6) * t365;
t66 = t232 + t241 / 0.2e1;
t278 = t66 * qJD(2);
t139 = Icges(5,5) * t284 - Icges(5,6) * t286 - Icges(5,3) * t220;
t275 = -t218 * t139 - t143 * t283;
t238 = Icges(5,5) * t219 - Icges(5,6) * t217;
t140 = Icges(5,3) * t218 + t238 * t220;
t274 = t218 * t140 + t144 * t283;
t255 = -t177 - t322;
t107 = t144 * t284;
t246 = t220 * t140 - t107;
t244 = t217 * t142 - t139;
t237 = -Icges(5,5) * t217 - Icges(5,6) * t219;
t226 = (-t172 + t175) * t213 - t364 * t212;
t229 = -t256 + (t171 * t218 + t271 * t212 + t269 * t213 + t226 * t220) * t327 + (-t272 * t212 + t270 * t213 + t226 * t218 - t297) * t325;
t227 = -t243 + (t327 + t328) * (t132 * t213 + t134 * t212);
t185 = -rSges(5,2) * t217 + t316;
t163 = t237 * t220;
t162 = t237 * t218;
t129 = t255 * t220;
t127 = t255 * t218;
t87 = t353 * t360;
t72 = -t262 * t323 - t353;
t65 = t232 - t241 / 0.2e1;
t64 = -t142 * t285 + t274;
t63 = -t141 * t285 - t275;
t62 = -t142 * t286 - t246;
t42 = t218 * t64 - t220 * t63;
t41 = t218 * t62 - t220 * (-(-t219 * t143 + t287) * t218 - t220 * t139);
t36 = t243 + t333;
t29 = t335 / 0.2e1;
t28 = t336 / 0.2e1;
t27 = t332 + t341 + t343;
t24 = t306 + t318;
t19 = (t181 / 0.2e1 + t180 / 0.2e1) * t219 + t361 + t342 + t334 + t243;
t16 = (t62 - t107 + (t140 + t287) * t220 + t275) * t220 + t274 * t218;
t15 = (t244 * t220 - t274 + t64) * t220 + (t244 * t218 + t246 + t63) * t218;
t8 = m(6) * (t177 * t365 - t353 * t82) + t319;
t7 = t8 * qJD(5);
t4 = t28 - t335 / 0.2e1 + t256;
t3 = t29 - t336 / 0.2e1 + t256;
t2 = t28 + t29 + t229;
t1 = (t42 / 0.2e1 - t16 / 0.2e1) * t220 + (t15 / 0.2e1 + t41 / 0.2e1) * t218 + t256;
t5 = [0, 0, 0, (t72 * t344 - t366) * t346 - t87, -t317 * t353 - t87; 0, qJD(3) * t27 + qJD(4) * t19 + qJD(5) * t36, qJD(2) * t27 + qJD(4) * t24 + qJD(5) * t65, t19 * qJD(2) + t24 * qJD(3) + t2 * qJD(5) + ((t127 * t94 + t129 * t93) * t344 + (-t357 * t185 + (-t168 * t220 + t169 * t218) * t183) * t345) * t346 + (t16 * t362 + (t217 * t267 + t219 * t265) * t327 + t229 + (t216 / 0.2e1 + t215 / 0.2e1) * t238 + (t15 + t41) * t328 + (-t217 * t268 + t219 * t266 + t42) * t325) * qJD(4), t36 * qJD(2) + t65 * qJD(3) + t2 * qJD(4) + t229 * qJD(5) + (-t231 + (-t152 * t220 + t153 * t218) * t176) * t360; 0, t23 * qJD(4) + t66 * qJD(5) + (-t332 / 0.4e1 - t341 / 0.4e1 - t343 / 0.4e1) * t347, 0, t367 + (-t127 * t220 + t129 * t218) * t317, t278; 0, (t227 - (t181 + t180) * t219 / 0.2e1 - t361) * qJD(2) - t23 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t334 / 0.4e1 - t342 / 0.4e1) * t347, -t367, t1 * qJD(2) + (m(5) * (t185 * t354 - (t218 * (rSges(5,1) * t284 - t263) + t220 * (rSges(5,1) * t283 + t253)) * t95) + (t215 * t163 + (-t218 * t162 + t349) * t220) * t327 + (t216 * t162 + (-t220 * t163 + t349) * t218) * t325 + m(6) * (-t127 * t356 - t129 * t355 + t47 * t72) + t319) * qJD(4) + t359, t4 * qJD(2) + t6 * qJD(4) + t359; 0, (t227 - t333) * qJD(2) - t66 * qJD(3) + t3 * qJD(4) + t256 * qJD(5), -t278, t3 * qJD(2) + ((t72 * t82 + (-t127 * t218 - t129 * t220) * t176) * m(6) + t319) * qJD(4) + t7, qJD(2) * t256 + qJD(4) * t8 + t7;];
Cq = t5;
