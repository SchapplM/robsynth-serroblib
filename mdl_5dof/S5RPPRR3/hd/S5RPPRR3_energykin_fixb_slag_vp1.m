% Calculate kinetic energy for
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:06
% EndTime: 2022-01-23 09:14:07
% DurationCPUTime: 0.58s
% Computational Cost: add. (783->136), mult. (559->221), div. (0->0), fcn. (452->10), ass. (0->83)
t291 = sin(qJ(1));
t326 = t291 * pkin(1);
t289 = cos(pkin(9));
t324 = t289 * pkin(3);
t286 = pkin(9) + qJ(4);
t279 = sin(t286);
t323 = Icges(5,4) * t279;
t281 = cos(t286);
t322 = Icges(5,4) * t281;
t283 = qJ(5) + t286;
t275 = sin(t283);
t321 = Icges(6,4) * t275;
t276 = cos(t283);
t320 = Icges(6,4) * t276;
t292 = cos(qJ(1));
t278 = qJD(1) * t292 * pkin(1);
t287 = qJ(1) + pkin(8);
t280 = sin(t287);
t282 = cos(t287);
t318 = qJD(1) * (t282 * pkin(2) + t280 * qJ(3)) + t278;
t317 = pkin(4) * t281;
t315 = qJD(4) * t280;
t314 = qJD(4) * t282;
t313 = qJD(4) + qJD(5);
t312 = pkin(4) * qJD(4) * t279;
t311 = -t280 * pkin(2) + t282 * qJ(3) - t326;
t310 = pkin(6) * t282 - t324 * t280 + t311;
t288 = sin(pkin(9));
t309 = rSges(4,1) * t289 - rSges(4,2) * t288;
t308 = rSges(5,1) * t281 - rSges(5,2) * t279;
t307 = rSges(6,1) * t276 - rSges(6,2) * t275;
t306 = Icges(5,1) * t281 - t323;
t305 = Icges(6,1) * t276 - t321;
t304 = -Icges(5,2) * t279 + t322;
t303 = -Icges(6,2) * t275 + t320;
t302 = Icges(5,5) * t281 - Icges(5,6) * t279;
t301 = Icges(6,5) * t276 - Icges(6,6) * t275;
t251 = -Icges(5,6) * t282 + t304 * t280;
t253 = -Icges(5,5) * t282 + t306 * t280;
t300 = t251 * t279 - t253 * t281;
t252 = Icges(5,6) * t280 + t304 * t282;
t254 = Icges(5,5) * t280 + t306 * t282;
t299 = -t252 * t279 + t254 * t281;
t267 = Icges(5,2) * t281 + t323;
t268 = Icges(5,1) * t279 + t322;
t298 = -t267 * t279 + t268 * t281;
t297 = -qJD(3) * t282 + qJD(1) * (pkin(6) * t280 + t324 * t282) + t318;
t264 = t313 * t280;
t265 = t313 * t282;
t296 = (Icges(6,5) * t275 + Icges(6,6) * t276) * qJD(1) - (-Icges(6,3) * t282 + t301 * t280) * t265 + (Icges(6,3) * t280 + t301 * t282) * t264;
t243 = -Icges(6,6) * t282 + t303 * t280;
t244 = Icges(6,6) * t280 + t303 * t282;
t245 = -Icges(6,5) * t282 + t305 * t280;
t246 = Icges(6,5) * t280 + t305 * t282;
t260 = Icges(6,2) * t276 + t321;
t261 = Icges(6,1) * t275 + t320;
t295 = (-t244 * t275 + t246 * t276) * t264 - (-t243 * t275 + t245 * t276) * t265 + (-t260 * t275 + t261 * t276) * qJD(1);
t293 = qJD(2) ^ 2;
t274 = qJD(3) * t280;
t273 = t292 * rSges(2,1) - t291 * rSges(2,2);
t272 = t291 * rSges(2,1) + t292 * rSges(2,2);
t269 = t279 * rSges(5,1) + t281 * rSges(5,2);
t266 = Icges(5,5) * t279 + Icges(5,6) * t281;
t262 = t275 * rSges(6,1) + t276 * rSges(6,2);
t258 = t278 + qJD(1) * (t282 * rSges(3,1) - t280 * rSges(3,2));
t257 = (-t280 * rSges(3,1) - t282 * rSges(3,2) - t326) * qJD(1);
t256 = t280 * rSges(5,3) + t308 * t282;
t255 = -t282 * rSges(5,3) + t308 * t280;
t250 = Icges(5,3) * t280 + t302 * t282;
t249 = -Icges(5,3) * t282 + t302 * t280;
t248 = t280 * rSges(6,3) + t307 * t282;
t247 = -t282 * rSges(6,3) + t307 * t280;
t238 = pkin(7) * t280 + t317 * t282;
t237 = -pkin(7) * t282 + t317 * t280;
t236 = qJD(1) * t280 * rSges(4,3) + (qJD(1) * t309 - qJD(3)) * t282 + t318;
t235 = t274 + (t282 * rSges(4,3) - t309 * t280 + t311) * qJD(1);
t234 = qJD(2) + (t255 * t280 + t256 * t282) * qJD(4);
t233 = qJD(1) * t256 - t269 * t315 + t297;
t232 = -t269 * t314 + t274 + (-t255 + t310) * qJD(1);
t231 = -t280 * t312 - t264 * t262 + (t238 + t248) * qJD(1) + t297;
t230 = -t282 * t312 - t265 * t262 + t274 + (-t237 - t247 + t310) * qJD(1);
t229 = t264 * t247 + t265 * t248 + qJD(2) + (t237 * t280 + t238 * t282) * qJD(4);
t1 = m(3) * (t257 ^ 2 + t258 ^ 2 + t293) / 0.2e1 + m(4) * (t235 ^ 2 + t236 ^ 2 + t293) / 0.2e1 + m(5) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + ((t280 * t266 + t298 * t282) * qJD(1) + (t280 ^ 2 * t250 + (t300 * t282 + (-t249 + t299) * t280) * t282) * qJD(4)) * t315 / 0.2e1 - ((-t282 * t266 + t298 * t280) * qJD(1) + (t282 ^ 2 * t249 + (t299 * t280 + (-t250 + t300) * t282) * t280) * qJD(4)) * t314 / 0.2e1 + m(6) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t264 * (t296 * t280 + t295 * t282) / 0.2e1 - t265 * (t295 * t280 - t296 * t282) / 0.2e1 + (((t281 * t252 + t279 * t254) * t280 - (t281 * t251 + t279 * t253) * t282) * qJD(4) + (t276 * t244 + t275 * t246) * t264 - (t276 * t243 + t275 * t245) * t265 + (t276 * t260 + t275 * t261 + t281 * t267 + t279 * t268) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t272 ^ 2 + t273 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t289 ^ 2 + (Icges(4,1) * t288 + 0.2e1 * Icges(4,4) * t289) * t288) * qJD(1) ^ 2 / 0.2e1;
T = t1;
