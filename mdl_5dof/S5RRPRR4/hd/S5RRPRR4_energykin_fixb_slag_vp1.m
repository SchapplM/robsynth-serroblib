% Calculate kinetic energy for
% S5RRPRR4
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:01
% EndTime: 2022-01-20 10:48:02
% DurationCPUTime: 0.62s
% Computational Cost: add. (803->128), mult. (514->211), div. (0->0), fcn. (408->10), ass. (0->81)
t263 = qJ(1) + qJ(2);
t257 = sin(t263);
t298 = pkin(2) * t257;
t266 = cos(qJ(4));
t297 = pkin(4) * t266;
t295 = pkin(1) * qJD(1);
t264 = sin(qJ(4));
t294 = Icges(5,4) * t264;
t293 = Icges(5,4) * t266;
t262 = qJ(4) + qJ(5);
t256 = sin(t262);
t292 = Icges(6,4) * t256;
t258 = cos(t262);
t291 = Icges(6,4) * t258;
t267 = cos(qJ(1));
t254 = t267 * t295;
t259 = cos(t263);
t261 = qJD(1) + qJD(2);
t290 = t261 * pkin(2) * t259 + t254;
t255 = pkin(9) + t263;
t251 = sin(t255);
t289 = qJD(4) * t251;
t252 = cos(t255);
t288 = qJD(4) * t252;
t287 = qJD(4) + qJD(5);
t286 = pkin(4) * qJD(4) * t264;
t265 = sin(qJ(1));
t285 = t265 * t295;
t284 = t261 * (pkin(3) * t252 + pkin(7) * t251) + t290;
t283 = -pkin(3) * t251 + pkin(7) * t252 - t298;
t282 = rSges(5,1) * t266 - rSges(5,2) * t264;
t281 = rSges(6,1) * t258 - rSges(6,2) * t256;
t280 = Icges(5,1) * t266 - t294;
t279 = Icges(6,1) * t258 - t292;
t278 = -Icges(5,2) * t264 + t293;
t277 = -Icges(6,2) * t256 + t291;
t276 = Icges(5,5) * t266 - Icges(5,6) * t264;
t275 = Icges(6,5) * t258 - Icges(6,6) * t256;
t228 = -Icges(5,6) * t252 + t278 * t251;
t230 = -Icges(5,5) * t252 + t280 * t251;
t274 = t228 * t264 - t230 * t266;
t229 = Icges(5,6) * t251 + t278 * t252;
t231 = Icges(5,5) * t251 + t280 * t252;
t273 = -t229 * t264 + t231 * t266;
t245 = Icges(5,2) * t266 + t294;
t246 = Icges(5,1) * t264 + t293;
t272 = -t245 * t264 + t246 * t266;
t237 = t287 * t251;
t238 = t287 * t252;
t271 = -(-Icges(6,3) * t252 + t275 * t251) * t238 + (Icges(6,3) * t251 + t275 * t252) * t237 + (Icges(6,5) * t256 + Icges(6,6) * t258) * t261;
t220 = -Icges(6,6) * t252 + t277 * t251;
t221 = Icges(6,6) * t251 + t277 * t252;
t222 = -Icges(6,5) * t252 + t279 * t251;
t223 = Icges(6,5) * t251 + t279 * t252;
t241 = Icges(6,2) * t258 + t292;
t242 = Icges(6,1) * t256 + t291;
t270 = (-t221 * t256 + t223 * t258) * t237 - (-t220 * t256 + t222 * t258) * t238 + (-t241 * t256 + t242 * t258) * t261;
t250 = rSges(2,1) * t267 - rSges(2,2) * t265;
t249 = rSges(2,1) * t265 + rSges(2,2) * t267;
t248 = rSges(5,1) * t264 + rSges(5,2) * t266;
t244 = Icges(5,5) * t264 + Icges(5,6) * t266;
t243 = rSges(6,1) * t256 + rSges(6,2) * t258;
t235 = t254 + t261 * (rSges(3,1) * t259 - rSges(3,2) * t257);
t234 = -t285 - t261 * (rSges(3,1) * t257 + rSges(3,2) * t259);
t233 = rSges(5,3) * t251 + t282 * t252;
t232 = -rSges(5,3) * t252 + t282 * t251;
t227 = Icges(5,3) * t251 + t276 * t252;
t226 = -Icges(5,3) * t252 + t276 * t251;
t225 = rSges(6,3) * t251 + t281 * t252;
t224 = -rSges(6,3) * t252 + t281 * t251;
t217 = t261 * (rSges(4,1) * t252 - rSges(4,2) * t251) + t290;
t216 = -t285 + (-rSges(4,1) * t251 - rSges(4,2) * t252 - t298) * t261;
t215 = pkin(8) * t251 + t297 * t252;
t214 = -pkin(8) * t252 + t297 * t251;
t213 = qJD(3) + (t232 * t251 + t233 * t252) * qJD(4);
t212 = t233 * t261 - t248 * t289 + t284;
t211 = -t285 - t248 * t288 + (-t232 + t283) * t261;
t210 = -t251 * t286 - t237 * t243 + (t215 + t225) * t261 + t284;
t209 = -t252 * t286 - t285 - t238 * t243 + (-t214 - t224 + t283) * t261;
t208 = t224 * t237 + t225 * t238 + qJD(3) + (t214 * t251 + t215 * t252) * qJD(4);
t1 = m(3) * (t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + m(5) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + ((t251 * t244 + t272 * t252) * t261 + (t251 ^ 2 * t227 + (t274 * t252 + (-t226 + t273) * t251) * t252) * qJD(4)) * t289 / 0.2e1 - ((-t252 * t244 + t272 * t251) * t261 + (t252 ^ 2 * t226 + (t273 * t251 + (-t227 + t274) * t252) * t251) * qJD(4)) * t288 / 0.2e1 + m(6) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + t237 * (t271 * t251 + t270 * t252) / 0.2e1 - t238 * (t270 * t251 - t271 * t252) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t261 ^ 2 / 0.2e1 + (((t229 * t266 + t231 * t264) * t251 - (t228 * t266 + t230 * t264) * t252) * qJD(4) + (t221 * t258 + t223 * t256) * t237 - (t220 * t258 + t222 * t256) * t238 + (t258 * t241 + t256 * t242 + t266 * t245 + t264 * t246) * t261) * t261 / 0.2e1 + (m(2) * (t249 ^ 2 + t250 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
