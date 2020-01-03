% Calculate kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:04
% EndTime: 2020-01-03 11:20:04
% DurationCPUTime: 0.59s
% Computational Cost: add. (609->129), mult. (605->212), div. (0->0), fcn. (578->10), ass. (0->64)
t313 = pkin(4) * qJD(1);
t281 = cos(pkin(9));
t282 = cos(pkin(8));
t300 = t281 * t282;
t280 = sin(pkin(8));
t277 = pkin(9) + qJ(5);
t273 = sin(t277);
t275 = cos(t277);
t278 = qJ(1) + pkin(7);
t276 = cos(t278);
t274 = sin(t278);
t304 = t274 * t282;
t255 = -t273 * t304 - t276 * t275;
t256 = -t276 * t273 + t275 * t304;
t302 = t276 * t282;
t257 = t273 * t302 - t274 * t275;
t258 = -t274 * t273 - t275 * t302;
t303 = t276 * t280;
t305 = t274 * t280;
t288 = (Icges(6,5) * t256 + Icges(6,6) * t255 + Icges(6,3) * t305) * t274 - (Icges(6,5) * t258 + Icges(6,6) * t257 - Icges(6,3) * t303) * t276;
t312 = t288 * t280;
t311 = -qJD(3) + t300 * t313 + (qJD(1) * pkin(6) - qJD(5) * (-t282 * rSges(6,3) + (rSges(6,1) * t275 - rSges(6,2) * t273) * t280)) * t280;
t307 = pkin(1) * qJD(1);
t279 = sin(pkin(9));
t306 = t274 * t279;
t301 = t279 * t282;
t264 = -t276 * pkin(2) - t274 * qJ(3);
t289 = pkin(3) * t282 + qJ(4) * t280;
t298 = t289 * t276 - t264;
t284 = sin(qJ(1));
t271 = t284 * t307;
t297 = qJD(1) * (t274 * pkin(2) - t276 * qJ(3)) + t271;
t295 = qJD(4) * t280;
t293 = qJD(5) * t280;
t292 = qJD(1) * t289 * t274 + t297;
t285 = cos(qJ(1));
t272 = t285 * t307;
t291 = -qJD(3) * t276 + t272;
t268 = -qJD(4) * t282 + qJD(2);
t290 = rSges(4,1) * t282 - rSges(4,2) * t280;
t286 = qJD(2) ^ 2;
t269 = -qJD(5) * t282 + qJD(1);
t267 = -t285 * rSges(2,1) + t284 * rSges(2,2);
t266 = t284 * rSges(2,1) + t285 * rSges(2,2);
t265 = t274 * t295;
t260 = t272 - qJD(1) * (-t276 * rSges(3,1) + t274 * rSges(3,2));
t259 = t271 + qJD(1) * (t274 * rSges(3,1) + t276 * rSges(3,2));
t253 = -Icges(6,5) * t282 + (Icges(6,1) * t275 - Icges(6,4) * t273) * t280;
t252 = -Icges(6,6) * t282 + (Icges(6,4) * t275 - Icges(6,2) * t273) * t280;
t251 = -Icges(6,3) * t282 + (Icges(6,5) * t275 - Icges(6,6) * t273) * t280;
t250 = (t274 * rSges(4,3) + t290 * t276 - t264) * qJD(1) + t291;
t249 = -qJD(1) * t276 * rSges(4,3) + (qJD(1) * t290 - qJD(3)) * t274 + t297;
t248 = t258 * rSges(6,1) + t257 * rSges(6,2) - rSges(6,3) * t303;
t247 = t256 * rSges(6,1) + t255 * rSges(6,2) + rSges(6,3) * t305;
t246 = Icges(6,1) * t258 + Icges(6,4) * t257 - Icges(6,5) * t303;
t245 = Icges(6,1) * t256 + Icges(6,4) * t255 + Icges(6,5) * t305;
t244 = Icges(6,4) * t258 + Icges(6,2) * t257 - Icges(6,6) * t303;
t243 = Icges(6,4) * t256 + Icges(6,2) * t255 + Icges(6,6) * t305;
t240 = t265 + (-(-t276 * t300 - t306) * rSges(5,1) - (-t274 * t281 + t276 * t301) * rSges(5,2) + rSges(5,3) * t303 + t298) * qJD(1) + t291;
t239 = -qJD(3) * t274 - t276 * t295 + qJD(1) * ((t274 * t300 - t276 * t279) * rSges(5,1) + (-t274 * t301 - t276 * t281) * rSges(5,2) + rSges(5,3) * t305) + t292;
t238 = (t247 * t276 + t248 * t274) * t293 + t268;
t237 = -t269 * t248 + t265 + t272 + (pkin(4) * t306 + t298) * qJD(1) + t311 * t276;
t236 = t269 * t247 + (-t279 * t313 - t295) * t276 + t311 * t274 + t292;
t1 = m(3) * (t259 ^ 2 + t260 ^ 2 + t286) / 0.2e1 + m(4) * (t249 ^ 2 + t250 ^ 2 + t286) / 0.2e1 + m(5) * (t239 ^ 2 + t240 ^ 2 + t268 ^ 2) / 0.2e1 + m(6) * (t236 ^ 2 + t237 ^ 2 + t238 ^ 2) / 0.2e1 + t269 * ((-t282 * t251 + (-t252 * t273 + t253 * t275) * t280) * t269 + (((-t243 * t273 + t245 * t275) * t274 - (-t244 * t273 + t246 * t275) * t276) * t280 - t288 * t282) * t293) / 0.2e1 + t274 * ((t251 * t305 + t255 * t252 + t256 * t253) * t269 + (-(t255 * t244 + t256 * t246) * t276 + (t255 * t243 + t256 * t245 + t312) * t274) * t293) * t293 / 0.2e1 - t276 * ((-t251 * t303 + t257 * t252 + t258 * t253) * t269 + ((t257 * t243 + t258 * t245) * t274 + (-t257 * t244 - t258 * t246 - t312) * t276) * t293) * t293 / 0.2e1 + (m(2) * (t266 ^ 2 + t267 ^ 2) + Icges(2,3) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t282 ^ 2 + ((Icges(4,1) + Icges(5,1) * t281 ^ 2 + (-0.2e1 * Icges(5,4) * t281 + Icges(5,2) * t279) * t279) * t280 + 0.2e1 * (-Icges(5,5) * t281 + Icges(5,6) * t279 + Icges(4,4)) * t282) * t280) * qJD(1) ^ 2 / 0.2e1;
T = t1;
