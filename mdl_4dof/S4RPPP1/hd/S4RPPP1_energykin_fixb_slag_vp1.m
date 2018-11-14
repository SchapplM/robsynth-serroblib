% Calculate kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:24
% EndTime: 2018-11-14 13:45:24
% DurationCPUTime: 0.19s
% Computational Cost: add. (284->77), mult. (340->103), div. (0->0), fcn. (259->10), ass. (0->45)
t272 = cos(pkin(4));
t297 = t272 ^ 2;
t295 = -0.2e1 * t272;
t294 = 0.2e1 * t272;
t292 = rSges(5,3) + qJ(4);
t270 = sin(pkin(4));
t273 = sin(qJ(1));
t291 = t270 * t273;
t274 = cos(qJ(1));
t290 = t270 * t274;
t269 = sin(pkin(6));
t286 = pkin(4) - pkin(6);
t280 = cos(t286) / 0.2e1;
t285 = pkin(4) + pkin(6);
t283 = cos(t285);
t278 = t280 + t283 / 0.2e1;
t249 = t269 * t273 - t274 * t278;
t271 = cos(pkin(6));
t279 = sin(t285) / 0.2e1;
t282 = sin(t286);
t277 = t279 - t282 / 0.2e1;
t250 = t271 * t273 + t274 * t277;
t260 = pkin(1) * t273 - qJ(2) * t290;
t289 = -pkin(2) * t250 - qJ(3) * t249 - t260;
t251 = t269 * t274 + t273 * t278;
t287 = qJD(2) * t270;
t263 = t273 * t287;
t288 = qJD(3) * t251 + t263;
t258 = t279 + t282 / 0.2e1;
t253 = qJD(2) * t272 - qJD(3) * t258;
t284 = t270 * (rSges(5,1) + pkin(3));
t281 = -t274 * t287 + qJD(1) * (pkin(1) * t274 + qJ(2) * t291);
t252 = t271 * t274 - t273 * t277;
t276 = qJD(1) * (pkin(2) * t252 + qJ(3) * t251) + qJD(3) * t249 + t281;
t262 = rSges(2,1) * t274 - rSges(2,2) * t273;
t261 = rSges(2,1) * t273 + rSges(2,2) * t274;
t259 = t280 - t283 / 0.2e1;
t246 = qJD(4) * t259 + t253;
t243 = qJD(1) * (rSges(3,1) * t252 - rSges(3,2) * t251 + rSges(3,3) * t291) + t281;
t242 = t263 + (-rSges(3,1) * t250 + rSges(3,2) * t249 + rSges(3,3) * t290 - t260) * qJD(1);
t241 = qJD(1) * (rSges(4,1) * t291 - rSges(4,2) * t252 + rSges(4,3) * t251) + t276;
t240 = (rSges(4,1) * t290 + rSges(4,2) * t250 - rSges(4,3) * t249 + t289) * qJD(1) + t288;
t239 = qJD(4) * t250 + (t251 * rSges(5,2) + t252 * t292 + t273 * t284) * qJD(1) + t276;
t238 = qJD(4) * t252 + (-t249 * rSges(5,2) - t250 * t292 + t274 * t284 + t289) * qJD(1) + t288;
t1 = m(3) * (qJD(2) ^ 2 * t297 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + m(4) * (t240 ^ 2 + t241 ^ 2 + t253 ^ 2) / 0.2e1 + m(5) * (t238 ^ 2 + t239 ^ 2 + t246 ^ 2) / 0.2e1 + (m(2) * (t261 ^ 2 + t262 ^ 2) + Icges(2,3) + (Icges(4,1) + Icges(5,1) + Icges(3,3)) * t297 + (Icges(4,4) * t295 + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t259 + (Icges(3,5) + Icges(5,5)) * t294) * t259 + (Icges(3,6) * t294 - 0.2e1 * Icges(5,6) * t259 + (Icges(3,2) + Icges(5,2) + Icges(4,3)) * t258 + 0.2e1 * (Icges(3,4) + Icges(4,6)) * t259 + (Icges(5,4) + Icges(4,5)) * t295) * t258) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
