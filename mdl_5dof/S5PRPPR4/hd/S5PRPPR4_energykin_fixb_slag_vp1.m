% Calculate kinetic energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:45
% EndTime: 2019-12-31 17:36:46
% DurationCPUTime: 0.33s
% Computational Cost: add. (415->102), mult. (535->164), div. (0->0), fcn. (540->6), ass. (0->52)
t288 = cos(pkin(8));
t304 = pkin(4) * t288;
t286 = pkin(7) + qJ(2);
t284 = sin(t286);
t285 = cos(t286);
t275 = t284 * pkin(2) - t285 * qJ(3);
t287 = sin(pkin(8));
t294 = pkin(3) * t288 + qJ(4) * t287;
t303 = -t294 * t284 - t275;
t283 = qJD(3) * t284;
t301 = qJD(4) * t287;
t302 = t285 * t301 + t283;
t300 = qJD(5) * t284;
t299 = qJD(5) * t285;
t298 = t284 * qJD(2);
t274 = qJD(2) * (t285 * pkin(2) + t284 * qJ(3));
t297 = qJD(2) * t294 * t285 + t284 * t301 + t274;
t282 = -qJD(4) * t288 + qJD(1);
t296 = rSges(4,1) * t288 - rSges(4,2) * t287;
t295 = rSges(5,1) * t288 + rSges(5,3) * t287;
t289 = sin(qJ(5));
t290 = cos(qJ(5));
t279 = t287 * t290 - t288 * t289;
t293 = t287 * t289 + t288 * t290;
t292 = qJD(1) ^ 2;
t291 = qJD(2) ^ 2;
t277 = t285 * rSges(3,1) - t284 * rSges(3,2);
t276 = t284 * rSges(3,1) + t285 * rSges(3,2);
t271 = t293 * t285;
t270 = t279 * t285;
t269 = t293 * t284;
t268 = t279 * t284;
t267 = t279 * rSges(6,1) - rSges(6,2) * t293;
t266 = Icges(6,1) * t279 - Icges(6,4) * t293;
t265 = Icges(6,4) * t279 - Icges(6,2) * t293;
t264 = Icges(6,5) * t279 - Icges(6,6) * t293;
t263 = rSges(4,3) * t298 + t274 + (qJD(2) * t296 - qJD(3)) * t285;
t262 = t283 + (t285 * rSges(4,3) - t284 * t296 - t275) * qJD(2);
t261 = t271 * rSges(6,1) + t270 * rSges(6,2) - t284 * rSges(6,3);
t260 = t269 * rSges(6,1) + t268 * rSges(6,2) + t285 * rSges(6,3);
t259 = Icges(6,1) * t271 + Icges(6,4) * t270 - Icges(6,5) * t284;
t258 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t285;
t257 = Icges(6,4) * t271 + Icges(6,2) * t270 - Icges(6,6) * t284;
t256 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t285;
t255 = Icges(6,5) * t271 + Icges(6,6) * t270 - Icges(6,3) * t284;
t254 = Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t285;
t253 = rSges(5,2) * t298 + (qJD(2) * t295 - qJD(3)) * t285 + t297;
t252 = (t285 * rSges(5,2) - t284 * t295 + t303) * qJD(2) + t302;
t251 = (-t260 * t284 - t261 * t285) * qJD(5) + t282;
t250 = t267 * t300 - qJD(3) * t285 + (-t284 * pkin(6) + t285 * t304 + t261) * qJD(2) + t297;
t249 = t267 * t299 + (-t285 * pkin(6) - t284 * t304 - t260 + t303) * qJD(2) + t302;
t1 = m(2) * t292 / 0.2e1 + m(3) * (t292 + (t276 ^ 2 + t277 ^ 2) * t291) / 0.2e1 + m(4) * (t262 ^ 2 + t263 ^ 2 + t292) / 0.2e1 + m(5) * (t252 ^ 2 + t253 ^ 2 + t282 ^ 2) / 0.2e1 + m(6) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 - ((-t284 * t264 + t270 * t265 + t271 * t266) * qJD(2) + (-(-t284 * t255 + t270 * t257 + t271 * t259) * t284 + (-t284 * t254 + t270 * t256 + t271 * t258) * t285) * qJD(5)) * t300 / 0.2e1 + ((t285 * t264 + t268 * t265 + t269 * t266) * qJD(2) + (-(t285 * t255 + t268 * t257 + t269 * t259) * t284 + (t285 * t254 + t268 * t256 + t269 * t258) * t285) * qJD(5)) * t299 / 0.2e1 + qJD(2) * ((-t265 * t293 + t279 * t266) * qJD(2) + (-(-t257 * t293 + t279 * t259) * t284 + (-t256 * t293 + t279 * t258) * t285) * qJD(5)) / 0.2e1 + (Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t288 ^ 2 + ((Icges(4,1) + Icges(5,1)) * t287 + 0.2e1 * (Icges(4,4) - Icges(5,5)) * t288) * t287) * t291 / 0.2e1;
T = t1;
