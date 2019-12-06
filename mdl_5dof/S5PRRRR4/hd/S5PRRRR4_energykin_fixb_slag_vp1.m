% Calculate kinetic energy for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:39
% EndTime: 2019-12-05 17:07:39
% DurationCPUTime: 0.54s
% Computational Cost: add. (783->120), mult. (489->202), div. (0->0), fcn. (396->8), ass. (0->77)
t261 = cos(qJ(4));
t290 = pkin(4) * t261;
t288 = pkin(2) * qJD(2);
t260 = sin(qJ(4));
t287 = Icges(5,4) * t260;
t286 = Icges(5,4) * t261;
t259 = qJ(4) + qJ(5);
t255 = sin(t259);
t285 = Icges(6,4) * t255;
t256 = cos(t259);
t284 = Icges(6,4) * t256;
t257 = pkin(9) + qJ(2);
t253 = cos(t257);
t248 = t253 * t288;
t254 = qJ(3) + t257;
t249 = sin(t254);
t250 = cos(t254);
t258 = qJD(2) + qJD(3);
t283 = t258 * (pkin(3) * t250 + pkin(7) * t249) + t248;
t282 = qJD(4) * t249;
t281 = qJD(4) * t250;
t280 = qJD(4) + qJD(5);
t279 = pkin(4) * qJD(4) * t260;
t252 = sin(t257);
t278 = t252 * t288;
t277 = rSges(5,1) * t261 - rSges(5,2) * t260;
t276 = rSges(6,1) * t256 - rSges(6,2) * t255;
t275 = Icges(5,1) * t261 - t287;
t274 = Icges(6,1) * t256 - t285;
t273 = -Icges(5,2) * t260 + t286;
t272 = -Icges(6,2) * t255 + t284;
t271 = Icges(5,5) * t261 - Icges(5,6) * t260;
t270 = Icges(6,5) * t256 - Icges(6,6) * t255;
t226 = -Icges(5,6) * t250 + t273 * t249;
t228 = -Icges(5,5) * t250 + t275 * t249;
t269 = t226 * t260 - t228 * t261;
t227 = Icges(5,6) * t249 + t273 * t250;
t229 = Icges(5,5) * t249 + t275 * t250;
t268 = -t227 * t260 + t229 * t261;
t245 = Icges(5,2) * t261 + t287;
t246 = Icges(5,1) * t260 + t286;
t267 = -t245 * t260 + t246 * t261;
t235 = t280 * t249;
t236 = t280 * t250;
t266 = -(-Icges(6,3) * t250 + t270 * t249) * t236 + (Icges(6,3) * t249 + t270 * t250) * t235 + (Icges(6,5) * t255 + Icges(6,6) * t256) * t258;
t218 = -Icges(6,6) * t250 + t272 * t249;
t219 = Icges(6,6) * t249 + t272 * t250;
t220 = -Icges(6,5) * t250 + t274 * t249;
t221 = Icges(6,5) * t249 + t274 * t250;
t241 = Icges(6,2) * t256 + t285;
t242 = Icges(6,1) * t255 + t284;
t265 = (-t219 * t255 + t221 * t256) * t235 - (-t218 * t255 + t220 * t256) * t236 + (-t241 * t255 + t242 * t256) * t258;
t264 = qJD(1) ^ 2;
t263 = qJD(2) ^ 2;
t247 = rSges(5,1) * t260 + rSges(5,2) * t261;
t244 = Icges(5,5) * t260 + Icges(5,6) * t261;
t243 = rSges(6,1) * t255 + rSges(6,2) * t256;
t239 = rSges(3,1) * t253 - rSges(3,2) * t252;
t238 = rSges(3,1) * t252 + rSges(3,2) * t253;
t237 = pkin(3) * t249 - pkin(7) * t250;
t233 = t248 + t258 * (rSges(4,1) * t250 - rSges(4,2) * t249);
t232 = -t278 - t258 * (rSges(4,1) * t249 + rSges(4,2) * t250);
t231 = rSges(5,3) * t249 + t277 * t250;
t230 = -rSges(5,3) * t250 + t277 * t249;
t225 = Icges(5,3) * t249 + t271 * t250;
t224 = -Icges(5,3) * t250 + t271 * t249;
t223 = rSges(6,3) * t249 + t276 * t250;
t222 = -rSges(6,3) * t250 + t276 * t249;
t215 = pkin(8) * t249 + t290 * t250;
t214 = -pkin(8) * t250 + t290 * t249;
t213 = qJD(1) + (t230 * t249 + t231 * t250) * qJD(4);
t212 = t231 * t258 - t247 * t282 + t283;
t211 = -t278 - t247 * t281 + (-t230 - t237) * t258;
t210 = -t249 * t279 - t235 * t243 + (t215 + t223) * t258 + t283;
t209 = -t250 * t279 - t278 - t236 * t243 + (-t214 - t222 - t237) * t258;
t208 = t222 * t235 + t223 * t236 + qJD(1) + (t214 * t249 + t215 * t250) * qJD(4);
t1 = m(2) * t264 / 0.2e1 + m(3) * (t264 + (t238 ^ 2 + t239 ^ 2) * t263) / 0.2e1 + t263 * Icges(3,3) / 0.2e1 + m(4) * (t232 ^ 2 + t233 ^ 2 + t264) / 0.2e1 + t258 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + ((t244 * t249 + t267 * t250) * t258 + (t225 * t249 ^ 2 + (t269 * t250 + (-t224 + t268) * t249) * t250) * qJD(4)) * t282 / 0.2e1 - ((-t244 * t250 + t267 * t249) * t258 + (t224 * t250 ^ 2 + (t268 * t249 + (-t225 + t269) * t250) * t249) * qJD(4)) * t281 / 0.2e1 + m(6) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + t235 * (t266 * t249 + t265 * t250) / 0.2e1 - t236 * (t265 * t249 - t266 * t250) / 0.2e1 + (((t227 * t261 + t229 * t260) * t249 - (t226 * t261 + t228 * t260) * t250) * qJD(4) + (t219 * t256 + t221 * t255) * t235 - (t218 * t256 + t220 * t255) * t236 + (t241 * t256 + t242 * t255 + t245 * t261 + t246 * t260) * t258) * t258 / 0.2e1;
T = t1;
