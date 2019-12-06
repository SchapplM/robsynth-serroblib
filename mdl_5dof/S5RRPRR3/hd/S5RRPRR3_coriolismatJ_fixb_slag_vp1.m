% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:19
% EndTime: 2019-12-05 18:30:23
% DurationCPUTime: 2.44s
% Computational Cost: add. (25408->222), mult. (12475->292), div. (0->0), fcn. (10650->10), ass. (0->151)
t329 = m(6) / 0.2e1;
t330 = m(5) / 0.2e1;
t226 = qJ(1) + qJ(2);
t222 = pkin(9) + t226;
t221 = qJ(4) + t222;
t216 = sin(t221);
t217 = cos(t221);
t229 = cos(qJ(5));
t290 = rSges(6,1) * t229;
t248 = pkin(4) + t290;
t227 = sin(qJ(5));
t270 = t216 * t227;
t253 = -rSges(6,2) * t270 - t217 * rSges(6,3);
t141 = -t217 * pkin(8) + t248 * t216 + t253;
t219 = sin(t222);
t223 = sin(t226);
t296 = pkin(2) * t223;
t244 = pkin(3) * t219 + t296;
t130 = t141 + t244;
t298 = sin(qJ(1)) * pkin(1);
t122 = t130 + t298;
t267 = t217 * t227;
t202 = rSges(6,2) * t267;
t142 = t202 - t248 * t217 + (-rSges(6,3) - pkin(8)) * t216;
t220 = cos(t222);
t224 = cos(t226);
t295 = pkin(2) * t224;
t243 = -pkin(3) * t220 - t295;
t131 = t142 + t243;
t297 = cos(qJ(1)) * pkin(1);
t123 = t131 - t297;
t51 = -t142 * t122 + t123 * t141;
t53 = -t142 * t130 + t131 * t141;
t183 = rSges(5,1) * t216 + rSges(5,2) * t217;
t165 = t183 + t244;
t156 = t165 + t298;
t184 = -t217 * rSges(5,1) + t216 * rSges(5,2);
t166 = t184 + t243;
t157 = t166 - t297;
t89 = -t184 * t156 + t157 * t183;
t92 = -t184 * t165 + t166 * t183;
t292 = (-t53 + t51) * t329 + (-t92 + t89) * t330;
t293 = (t53 + t51) * t329 + (t92 + t89) * t330;
t5 = t292 - t293;
t353 = t5 * qJD(1);
t328 = m(4) * (-t297 * (rSges(4,1) * t219 + rSges(4,2) * t220 + t296) - (-t220 * rSges(4,1) + t219 * rSges(4,2) - t295) * t298);
t208 = rSges(6,1) * t227 + rSges(6,2) * t229;
t180 = t208 * t217;
t179 = t208 * t216;
t284 = t130 * t179;
t64 = t131 * t180 - t284;
t352 = m(6) * t64;
t301 = m(3) * (-t297 * (rSges(3,1) * t223 + rSges(3,2) * t224) - (-t224 * rSges(3,1) + t223 * rSges(3,2)) * t298);
t288 = Icges(6,4) * t227;
t204 = Icges(6,2) * t229 + t288;
t207 = Icges(6,1) * t229 - t288;
t225 = Icges(6,4) * t229;
t343 = Icges(6,2) * t227 - t225;
t345 = Icges(6,1) * t227 + t225;
t245 = (-t343 / 0.2e1 + t345 / 0.2e1) * t229 + (-t204 / 0.2e1 + t207 / 0.2e1) * t227;
t344 = -t122 * t179 + t123 * t180;
t281 = t141 * t179;
t67 = t142 * t180 - t281;
t317 = m(6) * (t67 + t344);
t351 = t245 + t317 / 0.2e1;
t315 = m(6) * (t67 + t64);
t350 = t245 + t315 / 0.2e1;
t348 = m(6) * t67;
t210 = -rSges(6,2) * t227 + t290;
t347 = m(6) * t210;
t266 = t217 * t229;
t161 = Icges(6,4) * t266 - Icges(6,2) * t267 + Icges(6,6) * t216;
t200 = Icges(6,4) * t267;
t163 = Icges(6,1) * t266 + Icges(6,5) * t216 - t200;
t346 = (t161 * t227 - t163 * t229) * t217;
t70 = -t166 * t156 + t157 * t165;
t341 = (t204 - t207) * t229 + (-t343 + t345) * t227;
t254 = -Icges(6,2) * t266 + t163 - t200;
t256 = t217 * t345 + t161;
t340 = t227 * t254 + t229 * t256;
t199 = Icges(6,4) * t270;
t269 = t216 * t229;
t162 = -Icges(6,1) * t269 + Icges(6,5) * t217 + t199;
t255 = Icges(6,2) * t269 + t162 + t199;
t160 = Icges(6,6) * t217 + t343 * t216;
t257 = -t216 * t345 + t160;
t339 = -t227 * t255 - t229 * t257;
t338 = t216 ^ 2;
t337 = t217 ^ 2;
t336 = 0.4e1 * qJD(1);
t333 = 2 * qJD(4);
t324 = m(5) * t70;
t322 = m(5) * t89;
t321 = m(5) * t92;
t318 = m(6) * (t64 + t344);
t316 = m(6) * (t344 - t64);
t314 = m(6) * (t344 - t67);
t313 = m(6) * (-t284 + t281 + (t131 - t142) * t180);
t48 = -t131 * t122 + t123 * t130;
t312 = m(6) * t48;
t310 = m(6) * t51;
t309 = m(6) * t53;
t307 = m(6) * t344;
t304 = t216 / 0.2e1;
t303 = -t217 / 0.2e1;
t302 = t217 / 0.2e1;
t159 = Icges(6,5) * t266 - Icges(6,6) * t267 + Icges(6,3) * t216;
t246 = -t162 * t229 - t159;
t203 = Icges(6,5) * t229 - Icges(6,6) * t227;
t271 = t203 * t216;
t158 = Icges(6,3) * t217 - t271;
t258 = t216 * t158 + t162 * t266;
t259 = t217 * t158 + t160 * t270;
t278 = t160 * t227;
t78 = t217 * t159 + t161 * t270 - t163 * t269;
t79 = -t160 * t267 + t258;
t80 = t159 * t216 - t346;
t17 = (t259 + t80 + t346) * t217 + (-t79 + (t246 - t278) * t217 + t78 + t258) * t216;
t77 = -t162 * t269 + t259;
t18 = (t78 + (-t159 + t278) * t217 - t258) * t217 + (t246 * t216 + t259 - t77) * t216;
t46 = t216 * t78 + t217 * t77;
t47 = t216 * t80 + t217 * t79;
t2 = (t18 / 0.2e1 + t47 / 0.2e1) * t217 + (-t46 / 0.2e1 + t17 / 0.2e1) * t216;
t294 = t2 * qJD(5);
t241 = t318 / 0.2e1 + t245;
t237 = Icges(6,5) * t227 + Icges(6,6) * t229;
t173 = t216 * t237;
t233 = -t216 * t17 / 0.2e1 + (t18 + t47) * t303 + (t203 * t217 + t341 * t216 - t257 * t227 + t255 * t229) * t302 + (-t341 * t217 - t256 * t227 + t254 * t229 + t271 + t46) * t304;
t232 = -t245 + (t302 + t303) * (t229 * t161 + t227 * t163);
t174 = t237 * t217;
t129 = -t179 * t216 - t180 * t217;
t59 = t245 + t348;
t52 = t245 + t352;
t50 = t245 + t307;
t40 = t313 / 0.2e1;
t38 = t314 / 0.2e1;
t33 = t316 / 0.2e1;
t31 = t309 + t321;
t30 = t310 + t322;
t19 = t301 + t312 + t324 + t328;
t14 = -t313 / 0.2e1 + t350;
t13 = t40 + t350;
t12 = -t314 / 0.2e1 + t351;
t11 = t38 + t351;
t10 = -t316 / 0.2e1 + t241;
t9 = t33 + t241;
t8 = t40 - t315 / 0.2e1 + t232;
t7 = t38 - t317 / 0.2e1 + t232;
t6 = t33 - t318 / 0.2e1 + t232;
t4 = t292 + t293;
t1 = [qJD(2) * t19 + qJD(4) * t30 + qJD(5) * t50, t19 * qJD(1) + t4 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t328 / 0.2e1 + t301 / 0.2e1 + t48 * t329 + t70 * t330) * qJD(2), 0, t30 * qJD(1) + t4 * qJD(2) + t11 * qJD(5) + (t51 * t329 + t89 * t330) * t333, t50 * qJD(1) + t9 * qJD(2) + t11 * qJD(4) + ((t122 * t217 + t123 * t216) * t347 + t233) * qJD(5); -t5 * qJD(4) + t10 * qJD(5) + (-t301 / 0.4e1 - t328 / 0.4e1 - t324 / 0.4e1 - t312 / 0.4e1) * t336, qJD(4) * t31 + qJD(5) * t52, 0, -t353 + t31 * qJD(2) + t13 * qJD(5) + (t53 * t329 + t92 * t330) * t333, t10 * qJD(1) + t52 * qJD(2) + t13 * qJD(4) + ((t130 * t217 + t131 * t216) * t347 + t233) * qJD(5); 0, 0, 0, 0, m(6) * t129 * qJD(5); t5 * qJD(2) + t12 * qJD(5) + (-t322 / 0.4e1 - t310 / 0.4e1) * t336, t353 + t14 * qJD(5) + 0.4e1 * (-t309 / 0.4e1 - t321 / 0.4e1) * qJD(2), 0, qJD(5) * t59, t12 * qJD(1) + t14 * qJD(2) + t59 * qJD(4) + ((t141 * t217 + t142 * t216) * t347 + t233) * qJD(5); (t232 - t307) * qJD(1) + t6 * qJD(2) + t7 * qJD(4) + t294, t6 * qJD(1) + (t232 - t352) * qJD(2) + t8 * qJD(4) + t294, 0, t7 * qJD(1) + t8 * qJD(2) + (t232 - t348) * qJD(4) + t294, (m(6) * ((t217 * (rSges(6,1) * t266 + t216 * rSges(6,3) - t202) - t216 * (-rSges(6,1) * t269 - t253)) * t129 + (t337 + t338) * t210 * t208) + (t337 * t173 + (t340 * t216 + (-t174 - t339) * t217) * t216) * t302 + (-t338 * t174 + (t339 * t217 + (t173 - t340) * t216) * t217) * t304) * qJD(5) + (qJD(1) + qJD(2) + qJD(4)) * t2;];
Cq = t1;
