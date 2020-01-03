% Calculate kinetic energy for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:38
% EndTime: 2019-12-31 21:50:40
% DurationCPUTime: 1.27s
% Computational Cost: add. (936->140), mult. (771->220), div. (0->0), fcn. (641->8), ass. (0->93)
t388 = Icges(5,4) - Icges(6,5);
t387 = Icges(5,1) + Icges(6,1);
t386 = Icges(5,2) + Icges(6,3);
t313 = qJ(3) + qJ(4);
t308 = sin(t313);
t385 = t388 * t308;
t310 = cos(t313);
t384 = t388 * t310;
t383 = Icges(6,4) + Icges(5,5);
t382 = Icges(5,6) - Icges(6,6);
t381 = t386 * t308 - t384;
t380 = t387 * t310 - t385;
t379 = rSges(6,1) + pkin(4);
t378 = rSges(6,3) + qJ(5);
t377 = Icges(6,2) + Icges(5,3);
t314 = qJ(1) + qJ(2);
t309 = sin(t314);
t311 = cos(t314);
t376 = -t381 * t309 - t382 * t311;
t375 = -t382 * t309 + t381 * t311;
t374 = t380 * t309 - t383 * t311;
t373 = t383 * t309 + t380 * t311;
t372 = -t386 * t310 - t385;
t371 = t387 * t308 + t384;
t370 = -t382 * t308 + t383 * t310;
t369 = t378 * t308 + t379 * t310;
t345 = qJD(3) + qJD(4);
t288 = t345 * t309;
t289 = t345 * t311;
t312 = qJD(1) + qJD(2);
t368 = (t372 * t308 + t371 * t310) * t312 + (t376 * t308 - t374 * t310) * t289 + (t375 * t308 + t373 * t310) * t288;
t367 = (t383 * t308 + t382 * t310) * t312 + (-t370 * t309 + t377 * t311) * t289 + (t377 * t309 + t370 * t311) * t288;
t317 = cos(qJ(3));
t362 = pkin(3) * t317;
t360 = pkin(1) * qJD(1);
t315 = sin(qJ(3));
t359 = Icges(4,4) * t315;
t358 = Icges(4,4) * t317;
t257 = -pkin(8) * t311 + t362 * t309;
t258 = pkin(8) * t309 + t362 * t311;
t346 = qJD(3) * t311;
t347 = qJD(3) * t309;
t353 = t257 * t347 + t258 * t346;
t299 = pkin(2) * t309 - pkin(7) * t311;
t352 = -t257 - t299;
t351 = -rSges(6,2) * t311 + t369 * t309;
t350 = rSges(6,2) * t309 + t369 * t311;
t318 = cos(qJ(1));
t307 = t318 * t360;
t349 = t312 * (pkin(2) * t311 + pkin(7) * t309) + t307;
t348 = t379 * t308 - t378 * t310;
t344 = pkin(3) * qJD(3) * t315;
t316 = sin(qJ(1));
t343 = t316 * t360;
t342 = t312 * t258 + t349;
t341 = rSges(4,1) * t317 - rSges(4,2) * t315;
t340 = rSges(5,1) * t310 - rSges(5,2) * t308;
t337 = Icges(4,1) * t317 - t359;
t334 = -Icges(4,2) * t315 + t358;
t331 = Icges(4,5) * t317 - Icges(4,6) * t315;
t277 = -Icges(4,6) * t311 + t334 * t309;
t279 = -Icges(4,5) * t311 + t337 * t309;
t328 = t277 * t315 - t279 * t317;
t278 = Icges(4,6) * t309 + t334 * t311;
t280 = Icges(4,5) * t309 + t337 * t311;
t327 = -t278 * t315 + t280 * t317;
t301 = Icges(4,2) * t317 + t359;
t302 = Icges(4,1) * t315 + t358;
t326 = -t301 * t315 + t302 * t317;
t325 = qJD(5) * t308 - t344;
t305 = rSges(2,1) * t318 - rSges(2,2) * t316;
t304 = rSges(2,1) * t316 + rSges(2,2) * t318;
t303 = rSges(4,1) * t315 + rSges(4,2) * t317;
t300 = Icges(4,5) * t315 + Icges(4,6) * t317;
t298 = rSges(5,1) * t308 + rSges(5,2) * t310;
t284 = t307 + t312 * (rSges(3,1) * t311 - rSges(3,2) * t309);
t283 = -t343 - t312 * (rSges(3,1) * t309 + rSges(3,2) * t311);
t282 = rSges(4,3) * t309 + t341 * t311;
t281 = -rSges(4,3) * t311 + t341 * t309;
t276 = Icges(4,3) * t309 + t331 * t311;
t275 = -Icges(4,3) * t311 + t331 * t309;
t274 = rSges(5,3) * t309 + t340 * t311;
t272 = -rSges(5,3) * t311 + t340 * t309;
t253 = (t281 * t309 + t282 * t311) * qJD(3);
t252 = t282 * t312 - t303 * t347 + t349;
t251 = -t343 - t303 * t346 + (-t281 - t299) * t312;
t250 = t274 * t312 - t288 * t298 - t309 * t344 + t342;
t249 = -t311 * t344 - t343 - t289 * t298 + (-t272 + t352) * t312;
t248 = t272 * t288 + t274 * t289 + t353;
t247 = -t348 * t288 + t325 * t309 + t350 * t312 + t342;
t246 = -t343 + t325 * t311 - t348 * t289 + (-t351 + t352) * t312;
t245 = -qJD(5) * t310 + t351 * t288 + t350 * t289 + t353;
t1 = m(3) * (t283 ^ 2 + t284 ^ 2) / 0.2e1 + t312 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + ((t309 * t300 + t326 * t311) * t312 + (t309 ^ 2 * t276 + (t328 * t311 + (-t275 + t327) * t309) * t311) * qJD(3)) * t347 / 0.2e1 - ((-t311 * t300 + t326 * t309) * t312 + (t311 ^ 2 * t275 + (t327 * t309 + (-t276 + t328) * t311) * t309) * qJD(3)) * t346 / 0.2e1 + m(5) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(6) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + (t367 * t309 + t368 * t311) * t288 / 0.2e1 - (t368 * t309 - t367 * t311) * t289 / 0.2e1 + (m(2) * (t304 ^ 2 + t305 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t278 * t317 + t280 * t315) * t309 - (t277 * t317 + t315 * t279) * t311) * qJD(3) - (t374 * t308 + t376 * t310) * t289 + (t373 * t308 - t375 * t310) * t288 + (t317 * t301 + t315 * t302 + t371 * t308 - t372 * t310) * t312) * t312 / 0.2e1;
T = t1;
