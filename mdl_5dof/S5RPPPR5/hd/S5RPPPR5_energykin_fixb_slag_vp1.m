% Calculate kinetic energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:20
% EndTime: 2019-12-31 17:46:20
% DurationCPUTime: 0.33s
% Computational Cost: add. (359->98), mult. (574->159), div. (0->0), fcn. (608->8), ass. (0->57)
t289 = cos(qJ(1));
t288 = sin(qJ(1));
t262 = cos(pkin(8));
t287 = pkin(4) * t262;
t286 = cos(pkin(7));
t285 = sin(pkin(7));
t260 = pkin(8) + qJ(5);
t257 = sin(t260);
t284 = Icges(6,4) * t257;
t258 = cos(t260);
t283 = Icges(6,4) * t258;
t247 = t289 * t285 - t288 * t286;
t259 = qJD(2) * t288;
t281 = qJD(4) * t247 + t259;
t246 = -t288 * t285 - t289 * t286;
t280 = qJD(5) * t246;
t279 = qJD(5) * t247;
t249 = t288 * pkin(1) - t289 * qJ(2);
t278 = -t288 * pkin(2) - t249;
t277 = -qJD(2) * t289 + qJD(1) * (t289 * pkin(1) + t288 * qJ(2));
t261 = sin(pkin(8));
t276 = rSges(5,1) * t262 - rSges(5,2) * t261;
t275 = -rSges(6,1) * t258 + rSges(6,2) * t257;
t274 = -Icges(6,1) * t258 + t284;
t273 = Icges(6,2) * t257 - t283;
t272 = -Icges(6,5) * t258 + Icges(6,6) * t257;
t228 = -Icges(6,6) * t246 + t273 * t247;
t230 = -Icges(6,5) * t246 + t274 * t247;
t271 = -t228 * t257 + t230 * t258;
t229 = Icges(6,6) * t247 + t273 * t246;
t231 = Icges(6,5) * t247 + t274 * t246;
t270 = t229 * t257 - t231 * t258;
t243 = -Icges(6,2) * t258 - t284;
t244 = -Icges(6,1) * t257 - t283;
t269 = t243 * t257 - t244 * t258;
t268 = t247 * pkin(3) + t246 * qJ(4) + t278;
t267 = qJD(1) * t289 * pkin(2) + t277;
t266 = qJD(1) * (-t246 * pkin(3) + t247 * qJ(4)) - qJD(4) * t246 + t267;
t264 = qJD(3) ^ 2;
t251 = t289 * rSges(2,1) - t288 * rSges(2,2);
t250 = t288 * rSges(2,1) + t289 * rSges(2,2);
t245 = -rSges(6,1) * t257 - rSges(6,2) * t258;
t242 = -Icges(6,5) * t257 - Icges(6,6) * t258;
t239 = qJD(1) * (t289 * rSges(3,1) + t288 * rSges(3,3)) + t277;
t238 = t259 + (-t288 * rSges(3,1) + t289 * rSges(3,3) - t249) * qJD(1);
t235 = qJD(1) * (-rSges(4,1) * t246 - rSges(4,2) * t247) + t267;
t234 = t259 + (t247 * rSges(4,1) - t246 * rSges(4,2) + t278) * qJD(1);
t233 = rSges(6,3) * t247 + t275 * t246;
t232 = -rSges(6,3) * t246 + t275 * t247;
t227 = Icges(6,3) * t247 + t272 * t246;
t226 = -Icges(6,3) * t246 + t272 * t247;
t225 = qJD(1) * (rSges(5,3) * t247 - t276 * t246) + t266;
t224 = (t246 * rSges(5,3) + t276 * t247 + t268) * qJD(1) + t281;
t223 = -qJD(3) + (t232 * t247 + t233 * t246) * qJD(5);
t222 = -t245 * t279 + (pkin(6) * t247 - t287 * t246 + t233) * qJD(1) + t266;
t221 = -t245 * t280 + (pkin(6) * t246 + t287 * t247 - t232 + t268) * qJD(1) + t281;
t1 = m(3) * (t238 ^ 2 + t239 ^ 2) / 0.2e1 + m(4) * (t234 ^ 2 + t235 ^ 2 + t264) / 0.2e1 + m(5) * (t224 ^ 2 + t225 ^ 2 + t264) / 0.2e1 + m(6) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + ((t247 * t242 + t269 * t246) * qJD(1) + (t247 ^ 2 * t227 + (t271 * t246 + (-t226 + t270) * t247) * t246) * qJD(5)) * t279 / 0.2e1 - ((-t246 * t242 + t269 * t247) * qJD(1) + (t246 ^ 2 * t226 + (t270 * t247 + (-t227 + t271) * t246) * t247) * qJD(5)) * t280 / 0.2e1 + qJD(1) * ((-t258 * t243 - t257 * t244) * qJD(1) + ((-t229 * t258 - t231 * t257) * t247 - (-t228 * t258 - t257 * t230) * t246) * qJD(5)) / 0.2e1 + (m(2) * (t250 ^ 2 + t251 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3) + t262 ^ 2 * Icges(5,2) + (Icges(5,1) * t261 + 0.2e1 * Icges(5,4) * t262) * t261) * qJD(1) ^ 2 / 0.2e1;
T = t1;
