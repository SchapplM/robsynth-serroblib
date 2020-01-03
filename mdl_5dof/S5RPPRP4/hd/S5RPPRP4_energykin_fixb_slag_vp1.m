% Calculate kinetic energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:02
% EndTime: 2019-12-31 17:52:03
% DurationCPUTime: 0.89s
% Computational Cost: add. (457->85), mult. (948->132), div. (0->0), fcn. (1056->6), ass. (0->49)
t348 = -Icges(5,5) - Icges(6,5);
t347 = Icges(5,6) + Icges(6,6);
t320 = sin(pkin(7));
t321 = cos(pkin(7));
t325 = sin(qJ(1));
t326 = cos(qJ(1));
t267 = -t325 * t320 - t326 * t321;
t346 = t267 ^ 2;
t268 = t326 * t320 - t325 * t321;
t345 = t268 ^ 2;
t344 = rSges(6,1) + pkin(4);
t288 = sin(qJ(4));
t289 = cos(qJ(4));
t343 = t348 * t288 - t347 * t289;
t341 = Icges(5,3) + Icges(6,3);
t340 = t347 * t288 + t348 * t289;
t339 = rSges(6,3) + qJ(5);
t338 = t268 * t267;
t337 = rSges(6,2) * t288 - t344 * t289;
t336 = t343 * qJD(1);
t335 = -t267 * t341 + t340 * t268;
t334 = t340 * t267 + t268 * t341;
t315 = -t339 * t267 + t337 * t268;
t314 = t337 * t267 + t339 * t268;
t313 = qJD(4) * (-t288 * rSges(5,1) - t289 * rSges(5,2));
t278 = t325 * pkin(1) - t326 * qJ(2);
t310 = -t325 * pkin(2) - t278;
t309 = -qJD(2) * t326 + qJD(1) * (t326 * pkin(1) + t325 * qJ(2));
t308 = -rSges(5,1) * t289 + rSges(5,2) * t288;
t306 = qJD(4) * (t289 * rSges(6,2) + t344 * t288);
t293 = t268 * pkin(3) + t267 * pkin(6) + t310;
t292 = qJD(1) * t326 * pkin(2) + t309;
t291 = qJD(1) * (-t267 * pkin(3) + t268 * pkin(6)) + t292;
t286 = qJD(2) * t325;
t280 = t326 * rSges(2,1) - t325 * rSges(2,2);
t279 = t325 * rSges(2,1) + t326 * rSges(2,2);
t266 = qJD(1) * (t326 * rSges(3,1) + t325 * rSges(3,3)) + t309;
t265 = t286 + (-t325 * rSges(3,1) + t326 * rSges(3,3) - t278) * qJD(1);
t262 = t268 * rSges(5,3) + t308 * t267;
t260 = -t267 * rSges(5,3) + t308 * t268;
t246 = qJD(1) * (-t267 * rSges(4,1) - t268 * rSges(4,2)) + t292;
t245 = t286 + (t268 * rSges(4,1) - t267 * rSges(4,2) + t310) * qJD(1);
t242 = qJD(1) * t262 - t268 * t313 + t291;
t241 = -t267 * t313 + t286 + (-t260 + t293) * qJD(1);
t240 = -qJD(3) + (t260 * t268 + t262 * t267) * qJD(4);
t239 = t314 * qJD(1) - qJD(5) * t267 + t268 * t306 + t291;
t238 = qJD(5) * t268 + t286 + t267 * t306 + (t293 - t315) * qJD(1);
t237 = -qJD(3) + (t314 * t267 + t315 * t268) * qJD(4);
t1 = m(3) * (t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + m(5) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(6) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + (((Icges(5,2) + Icges(6,2)) * t289 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t288 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t289) * t288) * qJD(1) + (t345 + t346) * t343 * qJD(4)) * qJD(1) / 0.2e1 - ((-t334 * t338 + t335 * t346) * qJD(4) - t267 * t336) * qJD(4) * t267 / 0.2e1 + ((t334 * t345 - t335 * t338) * qJD(4) + t268 * t336) * qJD(4) * t268 / 0.2e1 + (m(2) * (t279 ^ 2 + t280 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
