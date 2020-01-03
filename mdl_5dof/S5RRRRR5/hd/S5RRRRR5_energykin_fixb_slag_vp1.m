% Calculate kinetic energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:02
% EndTime: 2020-01-03 12:13:03
% DurationCPUTime: 0.58s
% Computational Cost: add. (821->125), mult. (513->213), div. (0->0), fcn. (408->10), ass. (0->83)
t265 = qJ(1) + qJ(2);
t262 = qJ(3) + t265;
t252 = sin(t262);
t253 = cos(t262);
t264 = qJ(4) + qJ(5);
t258 = sin(t264);
t260 = cos(t264);
t294 = Icges(6,4) * t260;
t281 = -Icges(6,2) * t258 + t294;
t220 = -Icges(6,6) * t253 + t281 * t252;
t221 = -Icges(6,6) * t252 - t281 * t253;
t295 = Icges(6,4) * t258;
t283 = Icges(6,1) * t260 - t295;
t222 = -Icges(6,5) * t253 + t283 * t252;
t223 = -Icges(6,5) * t252 - t283 * t253;
t289 = -qJD(4) - qJD(5);
t237 = t289 * t252;
t238 = t289 * t253;
t241 = Icges(6,2) * t260 + t295;
t242 = Icges(6,1) * t258 + t294;
t263 = qJD(1) + qJD(2);
t257 = qJD(3) + t263;
t304 = (t241 * t258 - t242 * t260) * t257 + (t220 * t258 - t222 * t260) * t238 + (t221 * t258 - t223 * t260) * t237;
t301 = pkin(2) * t263;
t268 = cos(qJ(4));
t300 = pkin(4) * t268;
t298 = pkin(1) * qJD(1);
t266 = sin(qJ(4));
t297 = Icges(5,4) * t266;
t296 = Icges(5,4) * t268;
t267 = sin(qJ(1));
t255 = t267 * t298;
t259 = sin(t265);
t293 = t259 * t301 + t255;
t269 = cos(qJ(1));
t256 = t269 * t298;
t261 = cos(t265);
t292 = t261 * t301 + t256;
t291 = qJD(4) * t252;
t290 = qJD(4) * t253;
t288 = pkin(4) * qJD(4) * t266;
t287 = t257 * (pkin(3) * t252 - pkin(8) * t253) + t293;
t286 = rSges(5,1) * t268 - rSges(5,2) * t266;
t285 = rSges(6,1) * t260 - rSges(6,2) * t258;
t284 = Icges(5,1) * t268 - t297;
t282 = -Icges(5,2) * t266 + t296;
t280 = Icges(5,5) * t268 - Icges(5,6) * t266;
t279 = Icges(6,5) * t260 - Icges(6,6) * t258;
t228 = -Icges(5,6) * t253 + t282 * t252;
t230 = -Icges(5,5) * t253 + t284 * t252;
t276 = -t228 * t266 + t230 * t268;
t229 = -Icges(5,6) * t252 - t282 * t253;
t231 = -Icges(5,5) * t252 - t284 * t253;
t275 = t229 * t266 - t231 * t268;
t245 = Icges(5,2) * t268 + t297;
t246 = Icges(5,1) * t266 + t296;
t273 = t245 * t266 - t246 * t268;
t272 = -(-Icges(6,3) * t253 + t279 * t252) * t238 - (-Icges(6,3) * t252 - t279 * t253) * t237 - (Icges(6,5) * t258 + Icges(6,6) * t260) * t257;
t251 = -rSges(2,1) * t269 + rSges(2,2) * t267;
t250 = rSges(2,1) * t267 + rSges(2,2) * t269;
t249 = rSges(5,1) * t266 + rSges(5,2) * t268;
t244 = Icges(5,5) * t266 + Icges(5,6) * t268;
t243 = rSges(6,1) * t258 + rSges(6,2) * t260;
t239 = -pkin(3) * t253 - pkin(8) * t252;
t235 = t256 - t263 * (-rSges(3,1) * t261 + rSges(3,2) * t259);
t234 = t255 + t263 * (rSges(3,1) * t259 + rSges(3,2) * t261);
t233 = -rSges(5,3) * t252 - t286 * t253;
t232 = -rSges(5,3) * t253 + t286 * t252;
t227 = -Icges(5,3) * t252 - t280 * t253;
t226 = -Icges(5,3) * t253 + t280 * t252;
t225 = -rSges(6,3) * t252 - t285 * t253;
t224 = -rSges(6,3) * t253 + t285 * t252;
t217 = -pkin(9) * t252 - t300 * t253;
t216 = -pkin(9) * t253 + t300 * t252;
t215 = -t257 * (-rSges(4,1) * t253 + rSges(4,2) * t252) + t292;
t214 = t257 * (rSges(4,1) * t252 + rSges(4,2) * t253) + t293;
t213 = (t232 * t252 - t233 * t253) * qJD(4);
t212 = -t249 * t291 + (-t233 - t239) * t257 + t292;
t211 = t232 * t257 + t249 * t290 + t287;
t210 = -t252 * t288 + t237 * t243 + (-t217 - t225 - t239) * t257 + t292;
t209 = t253 * t288 - t238 * t243 + (t216 + t224) * t257 + t287;
t208 = -t224 * t237 + t225 * t238 + (t216 * t252 - t217 * t253) * qJD(4);
t1 = m(3) * (t234 ^ 2 + t235 ^ 2) / 0.2e1 + t263 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t214 ^ 2 + t215 ^ 2) / 0.2e1 + t257 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 - ((-t253 * t244 - t273 * t252) * t257 + (t253 ^ 2 * t226 + (t275 * t252 + (t227 - t276) * t253) * t252) * qJD(4)) * t290 / 0.2e1 - ((-t252 * t244 + t273 * t253) * t257 + (t252 ^ 2 * t227 + (t276 * t253 + (t226 - t275) * t252) * t253) * qJD(4)) * t291 / 0.2e1 + m(6) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + t238 * (-t304 * t252 + t272 * t253) / 0.2e1 + t237 * (t272 * t252 + t304 * t253) / 0.2e1 + ((-(t228 * t268 + t230 * t266) * t253 - (t229 * t268 + t231 * t266) * t252) * qJD(4) + (t220 * t260 + t222 * t258) * t238 + (t221 * t260 + t223 * t258) * t237 + (t241 * t260 + t242 * t258 + t245 * t268 + t246 * t266) * t257) * t257 / 0.2e1 + (m(2) * (t250 ^ 2 + t251 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
