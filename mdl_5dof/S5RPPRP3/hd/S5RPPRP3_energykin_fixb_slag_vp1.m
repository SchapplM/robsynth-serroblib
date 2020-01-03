% Calculate kinetic energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:43
% EndTime: 2019-12-31 17:50:43
% DurationCPUTime: 0.85s
% Computational Cost: add. (456->102), mult. (509->151), div. (0->0), fcn. (402->6), ass. (0->62)
t346 = Icges(5,4) + Icges(6,4);
t345 = Icges(5,1) + Icges(6,1);
t344 = -Icges(5,2) - Icges(6,2);
t282 = cos(qJ(4));
t343 = t346 * t282;
t280 = sin(qJ(4));
t342 = t346 * t280;
t341 = Icges(5,5) + Icges(6,5);
t340 = Icges(5,6) + Icges(6,6);
t339 = t344 * t282 - t342;
t338 = t345 * t280 + t343;
t337 = rSges(6,1) + pkin(4);
t336 = Icges(5,3) + Icges(6,3);
t278 = qJ(1) + pkin(7);
t276 = sin(t278);
t277 = cos(t278);
t335 = t339 * t276 - t340 * t277;
t334 = t340 * t276 + t339 * t277;
t333 = t338 * t276 + t341 * t277;
t332 = t341 * t276 - t338 * t277;
t331 = t344 * t280 + t343;
t330 = t345 * t282 - t342;
t329 = t341 * t280 + t340 * t282;
t328 = rSges(6,3) + qJ(5);
t327 = rSges(6,2) * t282 + t337 * t280;
t326 = t329 * t276 + t336 * t277;
t325 = t336 * t276 - t329 * t277;
t324 = -t340 * t280 + t341 * t282;
t323 = t330 * t280 + t331 * t282;
t322 = t332 * t280 + t334 * t282;
t321 = -t333 * t280 + t335 * t282;
t281 = sin(qJ(1));
t316 = t281 * pkin(1);
t309 = t327 * t276 + t328 * t277;
t308 = t328 * t276 - t327 * t277;
t283 = cos(qJ(1));
t275 = qJD(1) * t283 * pkin(1);
t307 = qJD(1) * (t277 * pkin(2) + t276 * qJ(3)) + t275;
t306 = qJD(4) * t276;
t305 = qJD(1) * t277 * pkin(6) + t307;
t302 = -t276 * pkin(2) + t277 * qJ(3) - t316;
t301 = -t280 * rSges(6,2) + t337 * t282;
t300 = rSges(5,1) * t280 + rSges(5,2) * t282;
t286 = -t276 * pkin(6) + t302;
t284 = qJD(2) ^ 2;
t274 = qJD(3) * t276;
t272 = t283 * rSges(2,1) - t281 * rSges(2,2);
t271 = t282 * rSges(5,1) - t280 * rSges(5,2);
t269 = t281 * rSges(2,1) + t283 * rSges(2,2);
t260 = t275 + qJD(1) * (t277 * rSges(3,1) - t276 * rSges(3,2));
t259 = (-t276 * rSges(3,1) - t277 * rSges(3,2) - t316) * qJD(1);
t256 = t276 * rSges(5,3) - t300 * t277;
t254 = t277 * rSges(5,3) + t300 * t276;
t240 = -qJD(3) * t277 + qJD(1) * (-t277 * rSges(4,2) + t276 * rSges(4,3)) + t307;
t239 = t274 + (t276 * rSges(4,2) + t277 * rSges(4,3) + t302) * qJD(1);
t238 = qJD(2) + (-t254 * t276 + t256 * t277) * qJD(4);
t237 = qJD(1) * t254 + (-qJD(4) * t271 - qJD(3)) * t277 + t305;
t236 = t271 * t306 + t274 + (-t256 + t286) * qJD(1);
t235 = qJD(5) * t276 + t309 * qJD(1) + (-t301 * qJD(4) - qJD(3)) * t277 + t305;
t234 = qJD(5) * t277 + t274 + t301 * t306 + (t286 - t308) * qJD(1);
t233 = qJD(2) + (-t309 * t276 + t308 * t277) * qJD(4);
t1 = m(3) * (t259 ^ 2 + t260 ^ 2 + t284) / 0.2e1 + m(4) * (t239 ^ 2 + t240 ^ 2 + t284) / 0.2e1 + m(5) * (t236 ^ 2 + t237 ^ 2 + t238 ^ 2) / 0.2e1 + m(6) * (t233 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + (((t335 * t280 + t333 * t282) * t277 + (-t334 * t280 + t332 * t282) * t276) * qJD(4) + (-t331 * t280 + t330 * t282) * qJD(1)) * qJD(1) / 0.2e1 + ((t325 * t276 ^ 2 + (t321 * t277 + (-t322 + t326) * t276) * t277) * qJD(4) + (t324 * t276 - t323 * t277) * qJD(1)) * t306 / 0.2e1 + ((t326 * t277 ^ 2 + (t322 * t276 + (-t321 + t325) * t277) * t276) * qJD(4) + (t323 * t276 + t324 * t277) * qJD(1)) * qJD(4) * t277 / 0.2e1 + (m(2) * (t269 ^ 2 + t272 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
