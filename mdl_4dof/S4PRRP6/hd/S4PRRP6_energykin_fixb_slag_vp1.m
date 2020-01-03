% Calculate kinetic energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:11
% EndTime: 2019-12-31 16:30:12
% DurationCPUTime: 1.00s
% Computational Cost: add. (428->108), mult. (1074->180), div. (0->0), fcn. (1106->6), ass. (0->65)
t309 = Icges(4,1) + Icges(5,1);
t308 = Icges(4,4) - Icges(5,5);
t307 = Icges(5,4) + Icges(4,5);
t306 = Icges(4,2) + Icges(5,3);
t305 = -Icges(4,6) + Icges(5,6);
t304 = -Icges(4,3) - Icges(5,2);
t257 = cos(pkin(6));
t288 = t257 ^ 2;
t256 = sin(pkin(6));
t289 = t256 ^ 2;
t290 = t288 + t289;
t303 = qJD(2) * t290;
t302 = rSges(5,1) + pkin(3);
t301 = rSges(5,3) + qJ(4);
t300 = t256 * t257;
t260 = cos(qJ(3));
t258 = sin(qJ(3));
t261 = cos(qJ(2));
t282 = t258 * t261;
t245 = t256 * t282 + t257 * t260;
t281 = t260 * t261;
t246 = t256 * t281 - t257 * t258;
t259 = sin(qJ(2));
t284 = t256 * t259;
t299 = t306 * t245 - t308 * t246 + t305 * t284;
t247 = -t256 * t260 + t257 * t282;
t248 = t256 * t258 + t257 * t281;
t283 = t257 * t259;
t298 = t306 * t247 - t308 * t248 + t305 * t283;
t297 = t305 * t245 + t307 * t246 - t304 * t284;
t296 = t305 * t247 + t307 * t248 - t304 * t283;
t295 = -t308 * t245 + t309 * t246 + t307 * t284;
t294 = -t308 * t247 + t309 * t248 + t307 * t283;
t293 = t305 * t261 + (-t306 * t258 + t308 * t260) * t259;
t292 = t304 * t261 + (t305 * t258 + t307 * t260) * t259;
t291 = t307 * t261 + (t308 * t258 - t309 * t260) * t259;
t287 = qJD(2) ^ 2;
t280 = rSges(5,2) * t284 + t301 * t245 + t302 * t246;
t279 = -rSges(5,2) * t283 - t301 * t247 - t302 * t248;
t278 = -rSges(5,2) * t261 + (t301 * t258 + t302 * t260) * t259;
t277 = qJD(2) * t256;
t276 = qJD(2) * t257;
t275 = qJD(3) * t259;
t274 = qJD(3) * t261;
t254 = pkin(2) * t259 - pkin(5) * t261;
t273 = t254 * t277;
t272 = t254 * t276;
t271 = qJD(1) + (pkin(2) * t261 + pkin(5) * t259) * t303;
t267 = Icges(3,5) * t261 - Icges(3,6) * t259;
t253 = rSges(3,1) * t259 + rSges(3,2) * t261;
t251 = t256 * t275 - t276;
t250 = t257 * t275 + t277;
t244 = -rSges(4,3) * t261 + (rSges(4,1) * t260 - rSges(4,2) * t258) * t259;
t230 = Icges(3,3) * t256 + t267 * t257;
t229 = -Icges(3,3) * t257 + t267 * t256;
t226 = rSges(4,1) * t248 - rSges(4,2) * t247 + rSges(4,3) * t283;
t224 = rSges(4,1) * t246 - rSges(4,2) * t245 + rSges(4,3) * t284;
t210 = qJD(1) + (rSges(3,1) * t261 - rSges(3,2) * t259) * t303;
t209 = t224 * t274 + t244 * t251 - t272;
t208 = -t226 * t274 - t244 * t250 - t273;
t207 = t224 * t250 - t226 * t251 + t271;
t206 = qJD(4) * t247 + t278 * t251 + t280 * t274 - t272;
t205 = qJD(4) * t245 - t278 * t250 + t279 * t274 - t273;
t204 = qJD(4) * t258 * t259 + t280 * t250 + t279 * t251 + t271;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t290 * t287 * t253 ^ 2 + t210 ^ 2) / 0.2e1 + t287 * t256 * (-t229 * t300 + t289 * t230) / 0.2e1 - t287 * t257 * (t288 * t229 - t230 * t300) / 0.2e1 + m(4) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(5) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + ((t293 * t247 + t291 * t248 - t292 * t283) * t274 + (t299 * t247 + t295 * t248 + t297 * t283) * t251 + (t298 * t247 + t294 * t248 + t296 * t283) * t250) * t250 / 0.2e1 + ((t293 * t245 + t291 * t246 - t292 * t284) * t274 + (t299 * t245 + t295 * t246 + t297 * t284) * t251 + (t298 * t245 + t294 * t246 + t296 * t284) * t250) * t251 / 0.2e1 - ((-t296 * t250 - t297 * t251 + t292 * t274) * t261 + ((t293 * t258 + t291 * t260) * t274 + (t299 * t258 + t295 * t260) * t251 + (t298 * t258 + t294 * t260) * t250) * t259) * t274 / 0.2e1;
T = t1;
