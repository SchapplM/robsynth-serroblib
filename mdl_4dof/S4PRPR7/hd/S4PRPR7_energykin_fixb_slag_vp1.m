% Calculate kinetic energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:36
% EndTime: 2019-12-31 16:25:36
% DurationCPUTime: 0.65s
% Computational Cost: add. (290->100), mult. (733->176), div. (0->0), fcn. (710->6), ass. (0->56)
t297 = Icges(4,1) + Icges(3,3);
t257 = sin(qJ(2));
t259 = cos(qJ(2));
t296 = (-Icges(4,4) + Icges(3,5)) * t259 + (Icges(4,5) - Icges(3,6)) * t257;
t255 = cos(pkin(6));
t288 = t255 ^ 2;
t254 = sin(pkin(6));
t289 = t254 ^ 2;
t278 = t288 + t289;
t290 = qJD(2) * t278;
t295 = t254 * t255;
t294 = t296 * t254 - t297 * t255;
t293 = t297 * t254 + t296 * t255;
t287 = qJD(2) ^ 2;
t286 = t254 * t259;
t285 = t255 * t259;
t256 = sin(qJ(4));
t284 = t256 * t257;
t258 = cos(qJ(4));
t283 = t257 * t258;
t282 = qJD(3) * t257;
t281 = qJD(4) * t257;
t280 = qJD(4) * t259;
t279 = qJD(1) + (pkin(2) * t259 + qJ(3) * t257) * t290;
t249 = t257 * pkin(2) - t259 * qJ(3);
t275 = qJD(2) * (t257 * rSges(4,2) + t259 * rSges(4,3) - t249);
t274 = qJD(2) * (-pkin(5) * t257 - t249);
t253 = t255 * t282;
t252 = t254 * t282;
t251 = t257 * rSges(3,1) + t259 * rSges(3,2);
t247 = -qJD(2) * t255 + t254 * t280;
t246 = qJD(2) * t254 + t255 * t280;
t245 = t254 * t284 - t255 * t258;
t244 = t254 * t283 + t255 * t256;
t243 = t254 * t258 + t255 * t284;
t242 = -t254 * t256 + t255 * t283;
t241 = t257 * rSges(5,3) + (-rSges(5,1) * t256 - rSges(5,2) * t258) * t259;
t240 = Icges(5,5) * t257 + (-Icges(5,1) * t256 - Icges(5,4) * t258) * t259;
t239 = Icges(5,6) * t257 + (-Icges(5,4) * t256 - Icges(5,2) * t258) * t259;
t238 = Icges(5,3) * t257 + (-Icges(5,5) * t256 - Icges(5,6) * t258) * t259;
t223 = t255 * t275 + t253;
t222 = t254 * t275 + t252;
t221 = t245 * rSges(5,1) + t244 * rSges(5,2) + rSges(5,3) * t286;
t220 = t243 * rSges(5,1) + t242 * rSges(5,2) + rSges(5,3) * t285;
t219 = Icges(5,1) * t245 + Icges(5,4) * t244 + Icges(5,5) * t286;
t218 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t285;
t217 = Icges(5,4) * t245 + Icges(5,2) * t244 + Icges(5,6) * t286;
t216 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t285;
t215 = Icges(5,5) * t245 + Icges(5,6) * t244 + Icges(5,3) * t286;
t214 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t285;
t213 = qJD(1) + (rSges(3,1) * t259 - rSges(3,2) * t257) * t290;
t212 = -qJD(3) * t259 + t279 + (-rSges(4,2) * t259 + rSges(4,3) * t257) * t290;
t211 = -t221 * t281 + t247 * t241 + t255 * t274 + t253;
t210 = t220 * t281 - t246 * t241 + t254 * t274 + t252;
t209 = -t247 * t220 + t246 * t221 + (pkin(5) * t290 - qJD(3)) * t259 + t279;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t278 * t287 * t251 ^ 2 + t213 ^ 2) / 0.2e1 + m(4) * (t212 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(5) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + t246 * ((t214 * t285 + t242 * t216 + t243 * t218) * t246 + (t215 * t285 + t242 * t217 + t243 * t219) * t247 + (t238 * t285 + t242 * t239 + t243 * t240) * t281) / 0.2e1 + t247 * ((t214 * t286 + t244 * t216 + t245 * t218) * t246 + (t215 * t286 + t244 * t217 + t245 * t219) * t247 + (t238 * t286 + t244 * t239 + t245 * t240) * t281) / 0.2e1 + ((t214 * t246 + t215 * t247 + t238 * t281) * t257 + ((-t216 * t258 - t218 * t256) * t246 + (-t217 * t258 - t219 * t256) * t247 + (-t239 * t258 - t240 * t256) * t281) * t259) * t281 / 0.2e1 + (t293 * t289 - t294 * t295) * t254 * t287 / 0.2e1 - (t294 * t288 - t293 * t295) * t255 * t287 / 0.2e1;
T = t1;
