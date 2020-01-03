% Calculate kinetic energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:06
% EndTime: 2019-12-31 18:16:07
% DurationCPUTime: 1.12s
% Computational Cost: add. (321->117), mult. (750->173), div. (0->0), fcn. (620->4), ass. (0->71)
t388 = -Icges(4,4) + Icges(6,4) + Icges(5,5);
t387 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t386 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t313 = cos(qJ(3));
t385 = t388 * t313;
t311 = sin(qJ(3));
t384 = t388 * t311;
t383 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t382 = -Icges(4,6) + Icges(5,6) - Icges(6,6);
t381 = -t386 * t313 + t384;
t380 = t387 * t311 - t385;
t379 = rSges(6,1) + pkin(4);
t378 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t312 = sin(qJ(1));
t314 = cos(qJ(1));
t377 = t381 * t312 + t382 * t314;
t376 = t382 * t312 - t381 * t314;
t375 = t380 * t312 + t383 * t314;
t374 = t383 * t312 - t380 * t314;
t373 = t386 * t311 + t385;
t372 = t387 * t313 + t384;
t371 = -t383 * t311 + t382 * t313;
t370 = -rSges(6,3) - qJ(5);
t369 = -rSges(6,2) * t313 + t379 * t311;
t368 = -t371 * t312 + t378 * t314;
t367 = t378 * t312 + t371 * t314;
t366 = t382 * t311 + t383 * t313;
t365 = t372 * t311 - t373 * t313;
t364 = t374 * t311 - t376 * t313;
t363 = -t375 * t311 + t377 * t313;
t334 = pkin(3) * t311 - qJ(4) * t313;
t286 = t334 * t314;
t345 = qJD(3) * t314;
t351 = qJD(4) * t311 - t286 * t345;
t350 = t369 * t312 + t370 * t314;
t349 = t370 * t312 - t369 * t314;
t302 = t313 * pkin(3) + t311 * qJ(4);
t310 = qJD(2) * t312;
t346 = qJD(3) * t312;
t348 = t302 * t346 + t310;
t290 = qJD(1) * (t314 * pkin(1) + t312 * qJ(2));
t347 = qJD(1) * t314 * pkin(6) + t290;
t344 = qJD(4) * t313;
t341 = t311 * rSges(6,2) + t379 * t313;
t300 = t312 * pkin(1) - t314 * qJ(2);
t340 = -pkin(6) * t312 - t300;
t285 = t334 * t312;
t339 = qJD(1) * t285 + t314 * t344 + t347;
t338 = t286 + t340;
t337 = rSges(4,1) * t311 + rSges(4,2) * t313;
t336 = rSges(5,1) * t311 - rSges(5,3) * t313;
t306 = t314 * rSges(2,1) - t312 * rSges(2,2);
t305 = t313 * rSges(4,1) - t311 * rSges(4,2);
t304 = t313 * rSges(5,1) + t311 * rSges(5,3);
t301 = t312 * rSges(2,1) + t314 * rSges(2,2);
t283 = t312 * rSges(4,3) - t337 * t314;
t282 = t312 * rSges(5,2) - t336 * t314;
t280 = t314 * rSges(4,3) + t337 * t312;
t279 = t314 * rSges(5,2) + t336 * t312;
t258 = t290 - qJD(2) * t314 + qJD(1) * (-t314 * rSges(3,2) + t312 * rSges(3,3));
t257 = t310 + (t312 * rSges(3,2) + t314 * rSges(3,3) - t300) * qJD(1);
t256 = (-t280 * t312 + t283 * t314) * qJD(3);
t255 = qJD(1) * t280 + (-qJD(3) * t305 - qJD(2)) * t314 + t347;
t254 = t305 * t346 + t310 + (-t283 + t340) * qJD(1);
t253 = (t282 * t314 + (-t279 - t285) * t312) * qJD(3) + t351;
t252 = qJD(1) * t279 + (-qJD(2) + (-t302 - t304) * qJD(3)) * t314 + t339;
t251 = (qJD(3) * t304 - t344) * t312 + (-t282 + t338) * qJD(1) + t348;
t250 = -qJD(5) * t312 + t350 * qJD(1) + (-qJD(2) + (-t302 - t341) * qJD(3)) * t314 + t339;
t249 = -qJD(5) * t314 + (t341 * qJD(3) - t344) * t312 + (t338 - t349) * qJD(1) + t348;
t248 = (t349 * t314 + (-t285 - t350) * t312) * qJD(3) + t351;
t1 = m(3) * (t257 ^ 2 + t258 ^ 2) / 0.2e1 + m(4) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + m(5) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + m(6) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + (m(2) * (t301 ^ 2 + t306 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1 + (((t377 * t311 + t375 * t313) * t314 + (t376 * t311 + t374 * t313) * t312) * qJD(3) + (t373 * t311 + t372 * t313) * qJD(1)) * qJD(1) / 0.2e1 + ((t367 * t312 ^ 2 + (t363 * t314 + (-t364 + t368) * t312) * t314) * qJD(3) + (t366 * t312 - t365 * t314) * qJD(1)) * t346 / 0.2e1 + ((t368 * t314 ^ 2 + (t364 * t312 + (-t363 + t367) * t314) * t312) * qJD(3) + (t365 * t312 + t366 * t314) * qJD(1)) * t345 / 0.2e1;
T = t1;
