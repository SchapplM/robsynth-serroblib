% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:50
% DurationCPUTime: 2.79s
% Computational Cost: add. (27928->236), mult. (13707->305), div. (0->0), fcn. (11584->10), ass. (0->163)
t246 = qJ(1) + qJ(2);
t245 = qJ(3) + t246;
t240 = sin(t245);
t241 = cos(t245);
t205 = rSges(4,1) * t240 + rSges(4,2) * t241;
t242 = sin(t246);
t320 = pkin(2) * t242;
t201 = t205 + t320;
t322 = sin(qJ(1)) * pkin(1);
t191 = t201 + t322;
t206 = -t241 * rSges(4,1) + t240 * rSges(4,2);
t243 = cos(t246);
t319 = pkin(2) * t243;
t202 = t206 - t319;
t321 = cos(qJ(1)) * pkin(1);
t192 = t202 - t321;
t118 = -t206 * t191 + t192 * t205;
t127 = -t206 * t201 + t202 * t205;
t358 = m(6) / 0.2e1;
t359 = m(5) / 0.2e1;
t360 = m(4) / 0.2e1;
t239 = pkin(9) + t245;
t236 = sin(t239);
t237 = cos(t239);
t249 = cos(qJ(5));
t314 = rSges(6,1) * t249;
t266 = pkin(4) + t314;
t247 = sin(qJ(5));
t291 = t236 * t247;
t273 = -rSges(6,2) * t291 - t237 * rSges(6,3);
t318 = pkin(3) * t240;
t150 = -t237 * pkin(8) + t266 * t236 + t273 + t318;
t140 = t150 + t320;
t134 = t140 + t322;
t288 = t237 * t247;
t222 = rSges(6,2) * t288;
t317 = pkin(3) * t241;
t151 = -t317 + t222 - t266 * t237 + (-rSges(6,3) - pkin(8)) * t236;
t141 = t151 - t319;
t135 = t141 - t321;
t50 = -t151 * t134 + t135 * t150;
t54 = -t151 * t140 + t141 * t150;
t189 = rSges(5,1) * t236 + rSges(5,2) * t237 + t318;
t181 = t189 + t320;
t172 = t181 + t322;
t190 = -t237 * rSges(5,1) + t236 * rSges(5,2) - t317;
t182 = t190 - t319;
t173 = t182 - t321;
t91 = -t190 * t172 + t173 * t189;
t96 = -t190 * t181 + t182 * t189;
t269 = (-t54 + t50) * t358 + (-t96 + t91) * t359 + (-t127 + t118) * t360;
t270 = (t54 + t50) * t358 + (t96 + t91) * t359 + (t127 + t118) * t360;
t1 = t270 - t269;
t383 = t1 * qJD(1);
t228 = rSges(6,1) * t247 + rSges(6,2) * t249;
t200 = t228 * t237;
t199 = t228 * t236;
t308 = t140 * t199;
t64 = t141 * t200 - t308;
t382 = m(6) * t64;
t330 = m(3) * (-t321 * (rSges(3,1) * t242 + rSges(3,2) * t243) - (-t243 * rSges(3,1) + t242 * rSges(3,2)) * t322);
t312 = Icges(6,4) * t247;
t224 = Icges(6,2) * t249 + t312;
t227 = Icges(6,1) * t249 - t312;
t244 = Icges(6,4) * t249;
t373 = Icges(6,2) * t247 - t244;
t375 = Icges(6,1) * t247 + t244;
t264 = (-t373 / 0.2e1 + t375 / 0.2e1) * t249 + (-t224 / 0.2e1 + t227 / 0.2e1) * t247;
t374 = -t134 * t199 + t135 * t200;
t305 = t150 * t199;
t66 = t151 * t200 - t305;
t346 = m(6) * (t66 + t374);
t381 = t264 + t346 / 0.2e1;
t345 = m(6) * (t66 + t64);
t380 = t264 + t345 / 0.2e1;
t378 = m(6) * t66;
t230 = -rSges(6,2) * t247 + t314;
t377 = m(6) * t230;
t287 = t237 * t249;
t177 = Icges(6,4) * t287 - Icges(6,2) * t288 + t236 * Icges(6,6);
t220 = Icges(6,4) * t288;
t179 = Icges(6,1) * t287 + t236 * Icges(6,5) - t220;
t376 = (t177 * t247 - t179 * t249) * t237;
t77 = -t182 * t172 + t173 * t181;
t103 = -t202 * t191 + t192 * t201;
t371 = (-t373 + t375) * t247 + (t224 - t227) * t249;
t274 = -Icges(6,2) * t287 + t179 - t220;
t276 = t237 * t375 + t177;
t370 = t274 * t247 + t276 * t249;
t219 = Icges(6,4) * t291;
t290 = t236 * t249;
t178 = -Icges(6,1) * t290 + Icges(6,5) * t237 + t219;
t275 = Icges(6,2) * t290 + t178 + t219;
t176 = Icges(6,6) * t237 + t373 * t236;
t277 = -t236 * t375 + t176;
t369 = -t275 * t247 - t277 * t249;
t368 = t236 ^ 2;
t367 = t237 ^ 2;
t365 = 0.4e1 * qJD(1);
t362 = 2 * qJD(3);
t353 = m(5) * t77;
t351 = m(5) * t91;
t350 = m(5) * t96;
t347 = m(6) * (t64 + t374);
t344 = m(6) * (t374 - t64);
t343 = m(6) * (t374 - t66);
t342 = m(6) * (-t308 + t305 + (t141 - t151) * t200);
t48 = -t141 * t134 + t135 * t140;
t341 = m(6) * t48;
t339 = m(6) * t50;
t338 = m(6) * t54;
t336 = m(6) * t374;
t333 = t236 / 0.2e1;
t332 = -t237 / 0.2e1;
t331 = t237 / 0.2e1;
t328 = m(4) * t103;
t326 = m(4) * t118;
t325 = m(4) * t127;
t175 = Icges(6,5) * t287 - Icges(6,6) * t288 + Icges(6,3) * t236;
t265 = -t178 * t249 - t175;
t223 = Icges(6,5) * t249 - Icges(6,6) * t247;
t292 = t223 * t236;
t174 = Icges(6,3) * t237 - t292;
t279 = t236 * t174 + t178 * t287;
t280 = t237 * t174 + t176 * t291;
t302 = t176 * t247;
t84 = t237 * t175 + t177 * t291 - t179 * t290;
t85 = -t176 * t288 + t279;
t86 = t175 * t236 - t376;
t17 = (t280 + t86 + t376) * t237 + (-t85 + (t265 - t302) * t237 + t84 + t279) * t236;
t83 = -t178 * t290 + t280;
t18 = (t84 + (-t175 + t302) * t237 - t279) * t237 + (t265 * t236 + t280 - t83) * t236;
t46 = t236 * t84 + t237 * t83;
t47 = t236 * t86 + t237 * t85;
t5 = (t18 / 0.2e1 + t47 / 0.2e1) * t237 + (-t46 / 0.2e1 + t17 / 0.2e1) * t236;
t316 = t5 * qJD(5);
t261 = t347 / 0.2e1 + t264;
t257 = Icges(6,5) * t247 + Icges(6,6) * t249;
t193 = t236 * t257;
t254 = -t236 * t17 / 0.2e1 + (t18 + t47) * t332 + (t223 * t237 + t371 * t236 - t277 * t247 + t275 * t249) * t331 + (-t371 * t237 - t276 * t247 + t274 * t249 + t292 + t46) * t333;
t253 = -t264 + (t331 + t332) * (t249 * t177 + t247 * t179);
t194 = t257 * t237;
t138 = -t199 * t236 - t200 * t237;
t62 = t264 + t378;
t59 = t264 + t382;
t53 = t264 + t336;
t42 = t342 / 0.2e1;
t40 = t343 / 0.2e1;
t38 = t344 / 0.2e1;
t29 = t325 + t338 + t350;
t28 = t326 + t339 + t351;
t19 = t328 + t330 + t341 + t353;
t14 = -t342 / 0.2e1 + t380;
t13 = t42 + t380;
t12 = -t343 / 0.2e1 + t381;
t11 = t40 + t381;
t10 = -t344 / 0.2e1 + t261;
t9 = t38 + t261;
t8 = t42 - t345 / 0.2e1 + t253;
t7 = t40 - t346 / 0.2e1 + t253;
t6 = t38 - t347 / 0.2e1 + t253;
t2 = t269 + t270;
t3 = [qJD(2) * t19 + qJD(3) * t28 + qJD(5) * t53, t19 * qJD(1) + t2 * qJD(3) + t9 * qJD(5) + 0.2e1 * (t330 / 0.2e1 + t103 * t360 + t48 * t358 + t77 * t359) * qJD(2), t28 * qJD(1) + t2 * qJD(2) + t11 * qJD(5) + (t118 * t360 + t50 * t358 + t91 * t359) * t362, 0, t53 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + ((t134 * t237 + t135 * t236) * t377 + t254) * qJD(5); t1 * qJD(3) + t10 * qJD(5) + (-t330 / 0.4e1 - t328 / 0.4e1 - t353 / 0.4e1 - t341 / 0.4e1) * t365, qJD(3) * t29 + qJD(5) * t59, t383 + t29 * qJD(2) + t13 * qJD(5) + (t127 * t360 + t54 * t358 + t96 * t359) * t362, 0, t10 * qJD(1) + t59 * qJD(2) + t13 * qJD(3) + ((t140 * t237 + t141 * t236) * t377 + t254) * qJD(5); -t1 * qJD(2) + t12 * qJD(5) + (-t326 / 0.4e1 - t351 / 0.4e1 - t339 / 0.4e1) * t365, -t383 + t14 * qJD(5) + 0.4e1 * (-t338 / 0.4e1 - t350 / 0.4e1 - t325 / 0.4e1) * qJD(2), qJD(5) * t62, 0, t12 * qJD(1) + t14 * qJD(2) + t62 * qJD(3) + ((t150 * t237 + t151 * t236) * t377 + t254) * qJD(5); 0, 0, 0, 0, m(6) * t138 * qJD(5); (t253 - t336) * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + t316, t6 * qJD(1) + (t253 - t382) * qJD(2) + t8 * qJD(3) + t316, t7 * qJD(1) + t8 * qJD(2) + (t253 - t378) * qJD(3) + t316, 0, (m(6) * ((t237 * (rSges(6,1) * t287 + t236 * rSges(6,3) - t222) - t236 * (-rSges(6,1) * t290 - t273)) * t138 + (t367 + t368) * t230 * t228) + (t367 * t193 + (t370 * t236 + (-t194 - t369) * t237) * t236) * t331 + (-t368 * t194 + (t369 * t237 + (t193 - t370) * t236) * t237) * t333) * qJD(5) + (qJD(1) + qJD(2) + qJD(3)) * t5;];
Cq = t3;
