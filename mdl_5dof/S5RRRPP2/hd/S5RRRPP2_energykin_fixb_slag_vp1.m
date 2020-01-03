% Calculate kinetic energy for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:28
% EndTime: 2019-12-31 20:51:29
% DurationCPUTime: 1.25s
% Computational Cost: add. (672->110), mult. (740->172), div. (0->0), fcn. (610->6), ass. (0->75)
t386 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t385 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t384 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t309 = cos(qJ(3));
t383 = t386 * t309;
t307 = sin(qJ(3));
t382 = t386 * t307;
t381 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t380 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t379 = t384 * t307 - t383;
t378 = t385 * t309 - t382;
t377 = rSges(6,1) + pkin(4);
t376 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t306 = qJ(1) + qJ(2);
t303 = sin(t306);
t304 = cos(t306);
t375 = t379 * t303 + t380 * t304;
t374 = -t380 * t303 + t379 * t304;
t373 = -t378 * t303 + t381 * t304;
t372 = t381 * t303 + t378 * t304;
t371 = -t384 * t309 - t382;
t370 = t385 * t307 + t383;
t369 = -t380 * t307 + t381 * t309;
t368 = rSges(6,3) + qJ(5);
t367 = rSges(6,2) * t307 + t377 * t309;
t366 = t369 * t303 - t376 * t304;
t365 = t376 * t303 + t369 * t304;
t364 = -t381 * t307 - t380 * t309;
t363 = t371 * t307 + t370 * t309;
t362 = t374 * t307 + t372 * t309;
t361 = -t375 * t307 + t373 * t309;
t355 = pkin(1) * qJD(1);
t348 = t367 * t303 + t368 * t304;
t347 = -t368 * t303 + t367 * t304;
t331 = pkin(3) * t309 + qJ(4) * t307;
t279 = t331 * t303;
t284 = t303 * pkin(2) - t304 * pkin(7);
t346 = -t279 - t284;
t310 = cos(qJ(1));
t302 = t310 * t355;
t305 = qJD(1) + qJD(2);
t345 = t305 * (t304 * pkin(2) + t303 * pkin(7)) + t302;
t344 = qJD(3) * t303;
t343 = qJD(3) * t304;
t342 = qJD(4) * t307;
t308 = sin(qJ(1));
t341 = t308 * t355;
t294 = t307 * pkin(3) - t309 * qJ(4);
t338 = qJD(3) * (-t307 * rSges(5,1) + t309 * rSges(5,3) - t294);
t280 = t331 * t304;
t337 = t305 * t280 + t303 * t342 + t345;
t336 = t304 * t342 - t341;
t335 = -qJD(4) * t309 + t279 * t344 + t280 * t343;
t334 = rSges(4,1) * t309 - rSges(4,2) * t307;
t333 = rSges(5,1) * t309 + rSges(5,3) * t307;
t312 = qJD(3) * (t309 * rSges(6,2) - t377 * t307 - t294);
t299 = t310 * rSges(2,1) - t308 * rSges(2,2);
t298 = t308 * rSges(2,1) + t310 * rSges(2,2);
t297 = t307 * rSges(4,1) + t309 * rSges(4,2);
t277 = t302 + t305 * (t304 * rSges(3,1) - t303 * rSges(3,2));
t276 = -t341 - t305 * (t303 * rSges(3,1) + t304 * rSges(3,2));
t275 = t303 * rSges(4,3) + t304 * t334;
t274 = t303 * rSges(5,2) + t304 * t333;
t272 = -t304 * rSges(4,3) + t303 * t334;
t271 = -t304 * rSges(5,2) + t303 * t333;
t249 = (t272 * t303 + t275 * t304) * qJD(3);
t248 = t305 * t275 - t297 * t344 + t345;
t247 = -t341 - t297 * t343 + (-t272 - t284) * t305;
t246 = t305 * t274 + t303 * t338 + t337;
t245 = t304 * t338 + (-t271 + t346) * t305 + t336;
t244 = (t271 * t303 + t274 * t304) * qJD(3) + t335;
t243 = qJD(5) * t304 + t303 * t312 + t305 * t347 + t337;
t242 = -qJD(5) * t303 + t304 * t312 + (t346 - t348) * t305 + t336;
t241 = (t303 * t348 + t304 * t347) * qJD(3) + t335;
t1 = m(3) * (t276 ^ 2 + t277 ^ 2) / 0.2e1 + t305 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + m(5) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + m(6) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + (m(2) * (t298 ^ 2 + t299 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t370 * t307 - t371 * t309) * t305 + ((t373 * t307 + t375 * t309) * t304 + (t372 * t307 - t374 * t309) * t303) * qJD(3)) * t305 / 0.2e1 + ((-t364 * t303 + t363 * t304) * t305 + (t365 * t303 ^ 2 + (t361 * t304 + (t362 - t366) * t303) * t304) * qJD(3)) * t344 / 0.2e1 - ((t363 * t303 + t364 * t304) * t305 + (t366 * t304 ^ 2 + (t362 * t303 + (t361 - t365) * t304) * t303) * qJD(3)) * t343 / 0.2e1;
T = t1;
