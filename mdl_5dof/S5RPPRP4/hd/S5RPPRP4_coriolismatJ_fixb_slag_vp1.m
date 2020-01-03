% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:11
% DurationCPUTime: 3.90s
% Computational Cost: add. (5457->219), mult. (11458->321), div. (0->0), fcn. (13464->6), ass. (0->160)
t356 = Icges(5,6) + Icges(6,6);
t195 = sin(qJ(4));
t196 = cos(qJ(4));
t252 = Icges(6,4) * t196;
t214 = -Icges(6,2) * t195 + t252;
t254 = Icges(5,4) * t196;
t216 = -Icges(5,2) * t195 + t254;
t355 = t214 + t216;
t349 = Icges(5,5) + Icges(6,5);
t353 = Icges(5,3) + Icges(6,3);
t210 = Icges(6,5) * t196 - Icges(6,6) * t195;
t212 = Icges(5,5) * t196 - Icges(5,6) * t195;
t352 = t210 + t212;
t253 = Icges(6,4) * t195;
t218 = Icges(6,1) * t196 - t253;
t255 = Icges(5,4) * t195;
t220 = Icges(5,1) * t196 - t255;
t351 = t218 + t220;
t256 = sin(pkin(7));
t257 = cos(pkin(7));
t281 = sin(qJ(1));
t282 = cos(qJ(1));
t168 = -t281 * t256 - t282 * t257;
t169 = t282 * t256 - t281 * t257;
t332 = t356 * t168 + t355 * t169;
t350 = t332 * t195;
t348 = -t355 * t168 + t356 * t169;
t347 = -t352 * t168 + t353 * t169;
t328 = -t351 * t168 + t349 * t169;
t341 = t353 * t168 + t352 * t169;
t259 = rSges(6,3) + qJ(5) + pkin(6);
t317 = t259 * t168;
t329 = t349 * t168 + t351 * t169;
t344 = t348 * t195;
t262 = rSges(6,2) * t196;
t311 = (rSges(6,1) + pkin(4)) * t195 + t262;
t201 = t311 * t168;
t120 = t201 * t282;
t222 = -t195 * rSges(6,1) - t262;
t242 = t169 * t195;
t126 = pkin(4) * t242 - t169 * t222;
t223 = t195 * rSges(5,1) + rSges(5,2) * t196;
t304 = m(6) / 0.2e1;
t305 = m(5) / 0.2e1;
t272 = (t126 * t281 + t120) * t304 + (t168 * t282 + t169 * t281) * t223 * t305;
t145 = t223 * t169;
t146 = t223 * t168;
t93 = t311 * t169;
t274 = (-t281 * t93 - t120) * t304 + (-t281 * t145 - t282 * t146) * t305;
t324 = t272 - t274;
t342 = t324 * qJD(1);
t313 = t281 * t168 - t282 * t169;
t337 = m(6) * t313;
t132 = -t337 / 0.2e1;
t130 = t337 / 0.2e1;
t340 = t349 * t195 + t356 * t196;
t243 = t168 * t196;
t325 = t347 * t169 - t328 * t243;
t241 = t169 * t196;
t339 = t347 * t168 + t328 * t241;
t165 = t168 * pkin(6);
t276 = pkin(4) * t196;
t192 = pkin(3) + t276;
t225 = rSges(6,1) * t241 + t169 * t192 + t317;
t260 = t195 * rSges(6,2);
t263 = rSges(6,1) * t196;
t186 = t260 - t263;
t327 = -t186 + t192;
t258 = t165 - t317 + (pkin(3) - t327) * t169;
t338 = m(6) * ((pkin(3) + t260) * t169 + t165 - t225 - t258);
t336 = t344 * t169 - t339;
t334 = t344 * t168 + t325;
t286 = -t168 / 0.2e1;
t323 = t168 / 0.2e1;
t284 = t169 / 0.2e1;
t164 = t168 * rSges(5,3);
t261 = t195 * rSges(5,2);
t264 = rSges(5,1) * t196;
t187 = t261 - t264;
t116 = t187 * t169 - t164;
t233 = rSges(5,1) * t241 + t164;
t321 = m(5) * (rSges(5,2) * t242 - t116 - t233);
t224 = -t281 * pkin(1) + t282 * qJ(2);
t200 = -t281 * pkin(2) + t224;
t316 = t165 + t200;
t67 = t164 + (pkin(3) - t187) * t169 + t316;
t69 = (pkin(3) - t261) * t169 + t233 + t316;
t320 = m(5) * (t67 - t69);
t91 = t201 * t168;
t319 = t340 * t169;
t318 = t340 * t168;
t213 = Icges(6,2) * t196 + t253;
t215 = Icges(5,2) * t196 + t255;
t315 = -t213 - t215;
t217 = Icges(6,1) * t195 + t252;
t219 = Icges(5,1) * t195 + t254;
t314 = -t217 - t219;
t310 = (-t220 / 0.2e1 + t215 / 0.2e1 - t218 / 0.2e1 + t213 / 0.2e1) * t195 + (-t219 / 0.2e1 - t216 / 0.2e1 - t217 / 0.2e1 - t214 / 0.2e1) * t196;
t167 = t168 ^ 2;
t309 = t169 ^ 2;
t308 = 0.4e1 * qJD(1);
t307 = 2 * qJD(4);
t306 = 4 * qJD(4);
t303 = m(6) / 0.4e1;
t166 = t169 * pkin(6);
t244 = t168 * t195;
t221 = rSges(6,2) * t244 + t259 * t169;
t227 = -t192 - t263;
t24 = t258 * t169 + (-t166 + (pkin(3) + t227) * t168 + t221) * t168;
t302 = t24 * t168 * t338;
t202 = t282 * pkin(1) + t281 * qJ(2);
t198 = t282 * pkin(2) + t202;
t301 = m(4) * ((-t168 * rSges(4,1) - t169 * rSges(4,2) + t198) * t281 + (t169 * rSges(4,1) - t168 * rSges(4,2) + t200) * t282);
t226 = rSges(5,2) * t244 + t169 * rSges(5,3);
t45 = t169 * t116 + t168 * (-rSges(5,1) * t243 + t226);
t300 = t45 * t168 * t321;
t68 = t166 + (-pkin(3) - t264) * t168 + t198 + t226;
t299 = m(5) * (-t145 * t67 + t146 * t68);
t63 = t67 * t282;
t298 = m(5) * (t68 * t281 + t63);
t60 = t327 * t169 + t200 + t317;
t61 = t227 * t168 + t198 + t221;
t294 = m(6) * (t201 * t61 - t60 * t93);
t57 = t60 * t282;
t292 = m(6) * (t61 * t281 + t57);
t290 = m(6) * (-t169 * t93 - t91);
t288 = m(6) * (t126 * t169 + t91);
t280 = m(3) * ((t282 * rSges(3,3) + t224) * t282 + (t281 * rSges(3,3) + t202) * t281);
t62 = -rSges(6,2) * t242 + t200 + t225;
t273 = -t60 + t62;
t266 = m(6) * qJD(1);
t265 = m(6) * qJD(4);
t245 = t168 * t169;
t20 = (t321 + t338) * t323;
t236 = t20 * qJD(3);
t235 = t20 * qJD(4);
t230 = -t212 / 0.2e1 - t210 / 0.2e1;
t228 = -t186 + t276;
t208 = -t325 * t169 / 0.2e1 + t334 * t284 + (t336 + t339) * t286;
t207 = ((-t344 + t341) * t168 - t325 + t334) * t323 + (t341 * t168 + (t329 * t196 - t350) * t169) * t286 + (-t168 * t350 + t329 * t243 + t336) * t284 + (t347 * t323 + (t328 * t196 - t344) * t284) * t169;
t129 = t228 * t168;
t127 = t228 * t169;
t59 = t145 * t169 + t146 * t168;
t51 = t60 * t168;
t50 = 0.2e1 * t132;
t49 = t130 + t132;
t48 = 0.2e1 * t130;
t47 = t288 / 0.2e1;
t42 = t290 / 0.2e1;
t41 = -t311 * t309 + (-pkin(4) * t195 + t222) * t167;
t26 = t61 * t169 + t51;
t15 = t20 * qJD(1);
t14 = t42 - t288 / 0.2e1;
t13 = t47 + t42;
t12 = t47 - t290 / 0.2e1;
t11 = t280 + t292 + t298 + t301;
t10 = t272 + t274;
t7 = t294 + t299 - t310;
t1 = t208 * t168 + t207 * t169 - t300 - t302;
t2 = [(-t68 * t320 / 0.4e1 + t273 * t61 * t303) * t308 + t11 * qJD(2) + t7 * qJD(4) + m(6) * t26 * qJD(5), qJD(1) * t11 + qJD(4) * t10 + qJD(5) * t49, -t235, t7 * qJD(1) + t10 * qJD(2) - t236 + (t127 * t61 + t129 * t60 + (t126 - t93) * t201) * t265 + t13 * qJD(5) + (t300 / 0.4e1 + t302 / 0.4e1) * t306 + ((m(5) * (t146 * t223 - t187 * t68) + t230 * t169 - t207) * t169 + (m(5) * (-t145 * t223 - t187 * t67) + t230 * t168 - t208) * t168 + ((t169 * t314 - t332) * t286 + (t168 * t314 + t348) * t284) * t195 + ((t169 * t315 + t329) * t286 + (t168 * t315 - t328) * t284) * t196) * qJD(4), t49 * qJD(2) + t13 * qJD(4) + t26 * t266; -t324 * qJD(4) + t48 * qJD(5) + (-t280 / 0.4e1 - t301 / 0.4e1 - t298 / 0.4e1 - t292 / 0.4e1) * t308 + 0.2e1 * ((-t282 * t69 + t63) * t305 + (-t282 * t62 + t57) * t304) * qJD(1), 0, 0, -t342 + ((-t127 * t282 + t129 * t281) * t304 - t313 * t187 * t305) * t307, t48 * qJD(1); t235, 0, 0, t15 + (-m(5) * t59 / 0.2e1 + t41 * t304) * t307, 0; t273 * t126 * t266 + t324 * qJD(2) + t236 + t1 * qJD(4) + t12 * qJD(5) + (-t294 / 0.4e1 - t299 / 0.4e1) * t308 + (-t320 * t145 + t310) * qJD(1), t342, t15, t1 * qJD(1) + ((t167 * t319 - t245 * t318) * t286 + (-t245 * t319 + t309 * t318) * t284) * qJD(4) + (m(5) * (t45 * t59 - (t167 + t309) * t187 * t223) / 0.4e1 + (t126 * t127 + t129 * t201 - t24 * t41) * t303) * t306, t12 * qJD(1); (-t168 * t62 - t26 + t51) * t266 + t50 * qJD(2) + t14 * qJD(4), t50 * qJD(1), 0, t14 * qJD(1) + (-t127 * t168 + t129 * t169) * t265, 0;];
Cq = t2;
