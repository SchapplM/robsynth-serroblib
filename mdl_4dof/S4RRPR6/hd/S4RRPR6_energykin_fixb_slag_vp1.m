% Calculate kinetic energy for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:22
% EndTime: 2019-12-31 17:04:22
% DurationCPUTime: 0.93s
% Computational Cost: add. (606->152), mult. (735->244), div. (0->0), fcn. (622->8), ass. (0->94)
t317 = Icges(3,3) + Icges(4,3);
t258 = qJ(2) + pkin(7);
t252 = sin(t258);
t253 = cos(t258);
t260 = sin(qJ(2));
t262 = cos(qJ(2));
t316 = Icges(3,5) * t262 + Icges(4,5) * t253 - Icges(3,6) * t260 - Icges(4,6) * t252;
t261 = sin(qJ(1));
t263 = cos(qJ(1));
t315 = t316 * t261 - t317 * t263;
t314 = t317 * t261 + t316 * t263;
t313 = Icges(3,5) * t260 + Icges(4,5) * t252 + Icges(3,6) * t262 + Icges(4,6) * t253;
t300 = Icges(4,4) * t252;
t235 = Icges(4,2) * t253 + t300;
t299 = Icges(4,4) * t253;
t236 = Icges(4,1) * t252 + t299;
t302 = Icges(3,4) * t260;
t243 = Icges(3,2) * t262 + t302;
t301 = Icges(3,4) * t262;
t244 = Icges(3,1) * t260 + t301;
t312 = -t235 * t252 + t236 * t253 - t243 * t260 + t244 * t262;
t278 = -Icges(4,2) * t252 + t299;
t217 = Icges(4,6) * t261 + t278 * t263;
t281 = Icges(4,1) * t253 - t300;
t219 = Icges(4,5) * t261 + t281 * t263;
t279 = -Icges(3,2) * t260 + t301;
t225 = Icges(3,6) * t261 + t279 * t263;
t282 = Icges(3,1) * t262 - t302;
t227 = Icges(3,5) * t261 + t282 * t263;
t311 = -t217 * t252 + t219 * t253 - t225 * t260 + t227 * t262;
t216 = -Icges(4,6) * t263 + t278 * t261;
t218 = -Icges(4,5) * t263 + t281 * t261;
t224 = -Icges(3,6) * t263 + t279 * t261;
t226 = -Icges(3,5) * t263 + t282 * t261;
t310 = t216 * t252 - t218 * t253 + t224 * t260 - t226 * t262;
t306 = pkin(2) * t260;
t304 = t262 * pkin(2);
t254 = qJ(4) + t258;
t249 = sin(t254);
t298 = Icges(5,4) * t249;
t250 = cos(t254);
t297 = Icges(5,4) * t250;
t212 = -qJ(3) * t263 + t304 * t261;
t213 = qJ(3) * t261 + t304 * t263;
t291 = qJD(2) * t263;
t292 = qJD(2) * t261;
t296 = t212 * t292 + t213 * t291;
t248 = pkin(1) * t261 - pkin(5) * t263;
t295 = -t212 - t248;
t294 = pkin(3) * t253;
t290 = qJD(2) + qJD(4);
t239 = qJD(1) * (pkin(1) * t263 + pkin(5) * t261);
t287 = qJD(1) * t213 - qJD(3) * t263 + t239;
t286 = rSges(3,1) * t262 - rSges(3,2) * t260;
t285 = rSges(4,1) * t253 - rSges(4,2) * t252;
t284 = rSges(5,1) * t250 - rSges(5,2) * t249;
t283 = qJD(2) * (-rSges(4,1) * t252 - rSges(4,2) * t253 - t306);
t280 = Icges(5,1) * t250 - t298;
t277 = -Icges(5,2) * t249 + t297;
t274 = Icges(5,5) * t250 - Icges(5,6) * t249;
t267 = qJD(2) * (-pkin(3) * t252 - t306);
t240 = t290 * t261;
t241 = t290 * t263;
t266 = (Icges(5,5) * t249 + Icges(5,6) * t250) * qJD(1) - (-Icges(5,3) * t263 + t274 * t261) * t241 + (Icges(5,3) * t261 + t274 * t263) * t240;
t205 = -Icges(5,6) * t263 + t277 * t261;
t206 = Icges(5,6) * t261 + t277 * t263;
t207 = -Icges(5,5) * t263 + t280 * t261;
t208 = Icges(5,5) * t261 + t280 * t263;
t231 = Icges(5,2) * t250 + t298;
t232 = Icges(5,1) * t249 + t297;
t265 = (-t206 * t249 + t208 * t250) * t240 - (-t205 * t249 + t207 * t250) * t241 + (-t231 * t249 + t232 * t250) * qJD(1);
t255 = qJD(3) * t261;
t247 = rSges(2,1) * t263 - rSges(2,2) * t261;
t246 = rSges(2,1) * t261 + rSges(2,2) * t263;
t245 = rSges(3,1) * t260 + rSges(3,2) * t262;
t233 = rSges(5,1) * t249 + rSges(5,2) * t250;
t229 = rSges(3,3) * t261 + t286 * t263;
t228 = -rSges(3,3) * t263 + t286 * t261;
t221 = rSges(4,3) * t261 + t285 * t263;
t220 = -rSges(4,3) * t263 + t285 * t261;
t211 = rSges(5,3) * t261 + t284 * t263;
t210 = -rSges(5,3) * t263 + t284 * t261;
t200 = pkin(6) * t261 + t294 * t263;
t199 = -pkin(6) * t263 + t294 * t261;
t198 = qJD(1) * t229 - t245 * t292 + t239;
t197 = -t245 * t291 + (-t228 - t248) * qJD(1);
t196 = (t228 * t261 + t229 * t263) * qJD(2);
t195 = qJD(1) * t221 + t261 * t283 + t287;
t194 = t255 + t263 * t283 + (-t220 + t295) * qJD(1);
t193 = (t220 * t261 + t221 * t263) * qJD(2) + t296;
t192 = -t233 * t240 + t261 * t267 + (t200 + t211) * qJD(1) + t287;
t191 = -t233 * t241 + t255 + t263 * t267 + (-t199 - t210 + t295) * qJD(1);
t190 = t210 * t240 + t211 * t241 + (t199 * t261 + t200 * t263) * qJD(2) + t296;
t1 = m(3) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(4) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(5) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + t240 * (t266 * t261 + t265 * t263) / 0.2e1 - t241 * (t265 * t261 - t266 * t263) / 0.2e1 + (m(2) * (t246 ^ 2 + t247 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t314 * t261 ^ 2 + (t310 * t263 + (t311 - t315) * t261) * t263) * qJD(2) + (t313 * t261 + t312 * t263) * qJD(1)) * t292 / 0.2e1 - ((t315 * t263 ^ 2 + (t311 * t261 + (t310 - t314) * t263) * t261) * qJD(2) + (t312 * t261 - t313 * t263) * qJD(1)) * t291 / 0.2e1 + ((t206 * t250 + t208 * t249) * t240 - (t205 * t250 + t207 * t249) * t241 + ((-t216 * t253 - t218 * t252 - t224 * t262 - t226 * t260) * t263 + (t217 * t253 + t219 * t252 + t225 * t262 + t227 * t260) * t261) * qJD(2) + (t231 * t250 + t232 * t249 + t235 * t253 + t236 * t252 + t243 * t262 + t244 * t260) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
