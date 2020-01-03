% Calculate kinetic energy for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:22
% EndTime: 2019-12-31 20:49:23
% DurationCPUTime: 1.45s
% Computational Cost: add. (888->134), mult. (747->200), div. (0->0), fcn. (617->8), ass. (0->87)
t393 = Icges(5,4) - Icges(6,5);
t392 = Icges(5,1) + Icges(6,1);
t391 = Icges(5,2) + Icges(6,3);
t310 = qJ(3) + pkin(8);
t306 = cos(t310);
t390 = t393 * t306;
t305 = sin(t310);
t389 = t393 * t305;
t388 = Icges(6,4) + Icges(5,5);
t387 = Icges(5,6) - Icges(6,6);
t386 = t391 * t305 - t390;
t385 = t392 * t306 - t389;
t384 = rSges(6,1) + pkin(4);
t383 = rSges(6,3) + qJ(5);
t311 = qJ(1) + qJ(2);
t307 = sin(t311);
t308 = cos(t311);
t382 = t386 * t307 + t387 * t308;
t381 = -t387 * t307 + t386 * t308;
t380 = -t385 * t307 + t388 * t308;
t379 = t388 * t307 + t385 * t308;
t378 = -t391 * t306 - t389;
t377 = t392 * t305 + t390;
t376 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t313 = sin(qJ(3));
t315 = cos(qJ(3));
t375 = Icges(4,5) * t315 - Icges(4,6) * t313 - t387 * t305 + t388 * t306;
t374 = t383 * t305 + t384 * t306;
t373 = t375 * t307 - t376 * t308;
t372 = t376 * t307 + t375 * t308;
t371 = Icges(4,5) * t313 + Icges(4,6) * t315 + t388 * t305 + t387 * t306;
t359 = Icges(4,4) * t313;
t297 = Icges(4,2) * t315 + t359;
t358 = Icges(4,4) * t315;
t298 = Icges(4,1) * t313 + t358;
t370 = -t297 * t313 + t298 * t315 + t378 * t305 + t377 * t306;
t334 = -Icges(4,2) * t313 + t358;
t276 = Icges(4,6) * t307 + t334 * t308;
t337 = Icges(4,1) * t315 - t359;
t278 = Icges(4,5) * t307 + t337 * t308;
t369 = -t276 * t313 + t278 * t315 + t381 * t305 + t379 * t306;
t275 = -Icges(4,6) * t308 + t334 * t307;
t277 = -Icges(4,5) * t308 + t337 * t307;
t368 = t275 * t313 - t277 * t315 - t382 * t305 + t380 * t306;
t363 = pkin(3) * t313;
t362 = pkin(3) * t315;
t360 = pkin(1) * qJD(1);
t255 = -qJ(4) * t308 + t362 * t307;
t256 = qJ(4) * t307 + t362 * t308;
t347 = qJD(3) * t308;
t348 = qJD(3) * t307;
t353 = t255 * t348 + t256 * t347;
t295 = pkin(2) * t307 - pkin(7) * t308;
t352 = -t255 - t295;
t351 = -rSges(6,2) * t308 + t374 * t307;
t350 = rSges(6,2) * t307 + t374 * t308;
t316 = cos(qJ(1));
t304 = t316 * t360;
t309 = qJD(1) + qJD(2);
t349 = t309 * (pkin(2) * t308 + pkin(7) * t307) + t304;
t314 = sin(qJ(1));
t346 = t314 * t360;
t343 = qJD(4) * t307 - t346;
t342 = rSges(4,1) * t315 - rSges(4,2) * t313;
t341 = rSges(5,1) * t306 - rSges(5,2) * t305;
t338 = qJD(3) * (-rSges(5,1) * t305 - rSges(5,2) * t306 - t363);
t319 = -qJD(4) * t308 + t309 * t256 + t349;
t318 = (t383 * t306 - t363) * qJD(3) + (-t384 * qJD(3) + qJD(5)) * t305;
t301 = rSges(2,1) * t316 - rSges(2,2) * t314;
t300 = rSges(2,1) * t314 + rSges(2,2) * t316;
t299 = rSges(4,1) * t313 + rSges(4,2) * t315;
t282 = t304 + t309 * (rSges(3,1) * t308 - rSges(3,2) * t307);
t281 = -t346 - t309 * (rSges(3,1) * t307 + rSges(3,2) * t308);
t280 = rSges(4,3) * t307 + t342 * t308;
t279 = -rSges(4,3) * t308 + t342 * t307;
t272 = rSges(5,3) * t307 + t341 * t308;
t270 = -rSges(5,3) * t308 + t341 * t307;
t251 = (t279 * t307 + t280 * t308) * qJD(3);
t250 = t280 * t309 - t299 * t348 + t349;
t249 = -t346 - t299 * t347 + (-t279 - t295) * t309;
t248 = t272 * t309 + t307 * t338 + t319;
t247 = t308 * t338 + (-t270 + t352) * t309 + t343;
t246 = (t270 * t307 + t272 * t308) * qJD(3) + t353;
t245 = t318 * t307 + t350 * t309 + t319;
t244 = (-t351 + t352) * t309 + t318 * t308 + t343;
t243 = -qJD(5) * t306 + (t351 * t307 + t350 * t308) * qJD(3) + t353;
t1 = m(3) * (t281 ^ 2 + t282 ^ 2) / 0.2e1 + t309 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(5) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(6) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + (m(2) * (t300 ^ 2 + t301 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t315 * t297 + t313 * t298 + t377 * t305 - t378 * t306) * t309 + ((-t275 * t315 - t277 * t313 + t380 * t305 + t382 * t306) * t308 + (t276 * t315 + t278 * t313 + t379 * t305 - t381 * t306) * t307) * qJD(3)) * t309 / 0.2e1 + ((t371 * t307 + t370 * t308) * t309 + (t372 * t307 ^ 2 + (t368 * t308 + (t369 - t373) * t307) * t308) * qJD(3)) * t348 / 0.2e1 - ((t370 * t307 - t371 * t308) * t309 + (t373 * t308 ^ 2 + (t369 * t307 + (t368 - t372) * t308) * t307) * qJD(3)) * t347 / 0.2e1;
T = t1;
