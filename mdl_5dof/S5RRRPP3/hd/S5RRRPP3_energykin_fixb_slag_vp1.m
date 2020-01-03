% Calculate kinetic energy for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:06
% EndTime: 2019-12-31 20:53:08
% DurationCPUTime: 1.27s
% Computational Cost: add. (673->110), mult. (743->172), div. (0->0), fcn. (613->6), ass. (0->75)
t387 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t386 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t385 = -Icges(4,2) - Icges(6,2) - Icges(5,3);
t310 = cos(qJ(3));
t384 = t387 * t310;
t308 = sin(qJ(3));
t383 = t387 * t308;
t382 = -Icges(5,4) + Icges(4,5) + Icges(6,5);
t381 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t380 = t386 * t310 - t383;
t379 = t385 * t308 + t384;
t378 = rSges(6,3) + qJ(5);
t377 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t307 = qJ(1) + qJ(2);
t304 = sin(t307);
t305 = cos(t307);
t376 = t379 * t304 + t381 * t305;
t375 = -t381 * t304 + t379 * t305;
t374 = -t380 * t304 + t382 * t305;
t373 = t382 * t304 + t380 * t305;
t372 = t386 * t308 + t384;
t371 = t385 * t310 - t383;
t370 = t381 * t308 + t382 * t310;
t369 = rSges(6,1) + pkin(4);
t368 = rSges(6,2) * t308 + t378 * t310;
t367 = t370 * t304 - t377 * t305;
t366 = t377 * t304 + t370 * t305;
t365 = t382 * t308 - t381 * t310;
t364 = t371 * t308 + t372 * t310;
t363 = -t375 * t308 + t373 * t310;
t362 = t376 * t308 + t374 * t310;
t357 = pkin(1) * qJD(1);
t349 = t369 * t304 + t368 * t305;
t348 = t368 * t304 - t369 * t305;
t332 = pkin(3) * t310 + qJ(4) * t308;
t280 = t332 * t304;
t285 = t304 * pkin(2) - t305 * pkin(7);
t347 = -t280 - t285;
t311 = cos(qJ(1));
t303 = t311 * t357;
t306 = qJD(1) + qJD(2);
t346 = t306 * (t305 * pkin(2) + t304 * pkin(7)) + t303;
t345 = qJD(3) * t304;
t344 = qJD(3) * t305;
t343 = qJD(4) * t308;
t309 = sin(qJ(1));
t342 = t309 * t357;
t295 = t308 * pkin(3) - t310 * qJ(4);
t339 = qJD(3) * (t308 * rSges(5,2) + t310 * rSges(5,3) - t295);
t281 = t332 * t305;
t338 = t306 * t281 + t304 * t343 + t346;
t337 = t305 * t343 - t342;
t336 = -qJD(4) * t310 + t280 * t345 + t281 * t344;
t335 = rSges(4,1) * t310 - rSges(4,2) * t308;
t334 = -rSges(5,2) * t310 + rSges(5,3) * t308;
t313 = qJD(5) * t310 + (t310 * rSges(6,2) - t378 * t308 - t295) * qJD(3);
t300 = t311 * rSges(2,1) - t309 * rSges(2,2);
t298 = t309 * rSges(2,1) + t311 * rSges(2,2);
t297 = t308 * rSges(4,1) + t310 * rSges(4,2);
t278 = t303 + t306 * (t305 * rSges(3,1) - t304 * rSges(3,2));
t277 = -t342 - t306 * (t304 * rSges(3,1) + t305 * rSges(3,2));
t276 = -t305 * rSges(5,1) + t304 * t334;
t274 = t304 * rSges(5,1) + t305 * t334;
t272 = t304 * rSges(4,3) + t305 * t335;
t271 = -t305 * rSges(4,3) + t304 * t335;
t250 = (t271 * t304 + t272 * t305) * qJD(3);
t249 = t306 * t272 - t297 * t345 + t346;
t248 = -t342 - t297 * t344 + (-t271 - t285) * t306;
t247 = t306 * t274 + t304 * t339 + t338;
t246 = t305 * t339 + (-t276 + t347) * t306 + t337;
t245 = (t274 * t305 + t276 * t304) * qJD(3) + t336;
t244 = t304 * t313 + t306 * t349 + t338;
t243 = (t347 - t348) * t306 + t313 * t305 + t337;
t242 = qJD(5) * t308 + (t304 * t348 + t305 * t349) * qJD(3) + t336;
t1 = m(3) * (t277 ^ 2 + t278 ^ 2) / 0.2e1 + t306 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(5) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + m(6) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + (m(2) * (t298 ^ 2 + t300 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t372 * t308 - t371 * t310) * t306 + ((t374 * t308 - t376 * t310) * t305 + (t373 * t308 + t375 * t310) * t304) * qJD(3)) * t306 / 0.2e1 + ((t365 * t304 + t364 * t305) * t306 + (t366 * t304 ^ 2 + (t362 * t305 + (t363 - t367) * t304) * t305) * qJD(3)) * t345 / 0.2e1 - ((t364 * t304 - t365 * t305) * t306 + (t367 * t305 ^ 2 + (t363 * t304 + (t362 - t366) * t305) * t304) * qJD(3)) * t344 / 0.2e1;
T = t1;
