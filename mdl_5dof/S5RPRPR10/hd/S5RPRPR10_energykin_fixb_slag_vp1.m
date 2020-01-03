% Calculate kinetic energy for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:54
% EndTime: 2019-12-31 18:25:54
% DurationCPUTime: 0.50s
% Computational Cost: add. (467->79), mult. (527->131), div. (0->0), fcn. (544->8), ass. (0->50)
t280 = qJ(3) + pkin(8);
t274 = sin(t280);
t275 = cos(t280);
t284 = sin(qJ(1));
t285 = cos(qJ(1));
t238 = -t284 * t274 - t285 * t275;
t293 = t238 ^ 2;
t239 = t285 * t274 - t284 * t275;
t292 = t239 ^ 2;
t258 = sin(qJ(5));
t260 = cos(qJ(5));
t291 = -Icges(6,5) * t258 - Icges(6,6) * t260;
t289 = t239 * t238;
t257 = qJD(1) - qJD(3);
t288 = t291 * t257;
t281 = qJD(5) * (-t258 * rSges(6,1) - t260 * rSges(6,2));
t279 = t285 * pkin(2);
t278 = t284 * pkin(2);
t259 = sin(qJ(3));
t277 = t285 * t259;
t276 = t284 * t259;
t273 = -qJD(2) * t285 + qJD(1) * (t285 * pkin(1) + t284 * qJ(2));
t272 = -rSges(6,1) * t260 + rSges(6,2) * t258;
t269 = -Icges(6,5) * t260 + Icges(6,6) * t258;
t265 = qJD(1) * t279 + t273;
t261 = cos(qJ(3));
t253 = t261 * pkin(3) + pkin(2);
t264 = t257 * (pkin(3) * t276 + t285 * t253 - t279) + t265;
t247 = t284 * pkin(1) - t285 * qJ(2);
t255 = qJD(2) * t284;
t263 = t255 + (-t278 - t247) * qJD(1);
t249 = t285 * rSges(2,1) - t284 * rSges(2,2);
t248 = t284 * rSges(2,1) + t285 * rSges(2,2);
t241 = t284 * t261 - t277;
t240 = -t285 * t261 - t276;
t237 = -pkin(3) * t277 + t284 * t253 - t278;
t235 = qJD(1) * (t285 * rSges(3,1) + t284 * rSges(3,3)) + t273;
t234 = t255 + (-t284 * rSges(3,1) + t285 * rSges(3,3) - t247) * qJD(1);
t233 = t257 * (-t240 * rSges(4,1) + t241 * rSges(4,2)) + t265;
t232 = -t257 * (t241 * rSges(4,1) + t240 * rSges(4,2)) + t263;
t231 = t239 * rSges(6,3) + t272 * t238;
t230 = -t238 * rSges(6,3) + t272 * t239;
t225 = Icges(6,3) * t239 + t269 * t238;
t224 = -Icges(6,3) * t238 + t269 * t239;
t223 = t257 * (-t238 * rSges(5,1) - t239 * rSges(5,2)) + t264;
t222 = (t239 * rSges(5,1) - t238 * rSges(5,2) - t237) * t257 + t263;
t221 = -qJD(4) + (t230 * t239 + t231 * t238) * qJD(5);
t220 = -t239 * t281 + (-t238 * pkin(4) + t239 * pkin(7) + t231) * t257 + t264;
t219 = -t238 * t281 + (t239 * pkin(4) + t238 * pkin(7) - t230 - t237) * t257 + t263;
t1 = m(3) * (t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(4) * (t232 ^ 2 + t233 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(6) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + qJD(5) * t239 * (t239 * t288 + (-t224 * t289 + t292 * t225) * qJD(5)) / 0.2e1 - qJD(5) * t238 * (-t238 * t288 + (t293 * t224 - t225 * t289) * qJD(5)) / 0.2e1 + t257 * ((t260 ^ 2 * Icges(6,2) + (Icges(6,1) * t258 + 0.2e1 * Icges(6,4) * t260) * t258) * t257 + (t292 + t293) * t291 * qJD(5)) / 0.2e1 + (Icges(4,3) + Icges(5,3)) * t257 ^ 2 / 0.2e1 + (m(2) * (t248 ^ 2 + t249 ^ 2) + Icges(2,3) + Icges(3,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
