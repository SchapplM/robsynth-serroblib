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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:58:08
% EndTime: 2019-12-05 18:58:08
% DurationCPUTime: 0.53s
% Computational Cost: add. (821->125), mult. (513->213), div. (0->0), fcn. (408->10), ass. (0->83)
t262 = qJ(1) + qJ(2);
t259 = qJ(3) + t262;
t251 = sin(t259);
t252 = cos(t259);
t261 = qJ(4) + qJ(5);
t255 = sin(t261);
t257 = cos(t261);
t293 = Icges(6,4) * t257;
t281 = -Icges(6,2) * t255 + t293;
t221 = Icges(6,6) * t252 - t281 * t251;
t222 = Icges(6,6) * t251 + t281 * t252;
t294 = Icges(6,4) * t255;
t283 = Icges(6,1) * t257 - t294;
t223 = Icges(6,5) * t252 - t283 * t251;
t224 = Icges(6,5) * t251 + t283 * t252;
t290 = qJD(4) + qJD(5);
t238 = t290 * t251;
t239 = t290 * t252;
t242 = Icges(6,2) * t257 + t294;
t243 = Icges(6,1) * t255 + t293;
t260 = qJD(1) + qJD(2);
t254 = qJD(3) + t260;
t303 = (t242 * t255 - t243 * t257) * t254 + (t221 * t255 - t223 * t257) * t239 + (t222 * t255 - t224 * t257) * t238;
t300 = pkin(2) * t260;
t265 = cos(qJ(4));
t299 = pkin(4) * t265;
t297 = pkin(1) * qJD(1);
t263 = sin(qJ(4));
t296 = Icges(5,4) * t263;
t295 = Icges(5,4) * t265;
t292 = qJD(4) * t251;
t291 = qJD(4) * t252;
t289 = pkin(4) * qJD(4) * t263;
t264 = sin(qJ(1));
t288 = t264 * t297;
t266 = cos(qJ(1));
t287 = t266 * t297;
t286 = rSges(5,1) * t265 - rSges(5,2) * t263;
t285 = rSges(6,1) * t257 - rSges(6,2) * t255;
t284 = Icges(5,1) * t265 - t296;
t282 = -Icges(5,2) * t263 + t295;
t280 = Icges(5,5) * t265 - Icges(5,6) * t263;
t279 = Icges(6,5) * t257 - Icges(6,6) * t255;
t229 = Icges(5,6) * t252 - t282 * t251;
t231 = Icges(5,5) * t252 - t284 * t251;
t276 = -t229 * t263 + t231 * t265;
t230 = Icges(5,6) * t251 + t282 * t252;
t232 = Icges(5,5) * t251 + t284 * t252;
t275 = t230 * t263 - t232 * t265;
t246 = Icges(5,2) * t265 + t296;
t247 = Icges(5,1) * t263 + t295;
t273 = t246 * t263 - t247 * t265;
t256 = sin(t262);
t272 = -t256 * t300 - t288;
t258 = cos(t262);
t271 = -t258 * t300 - t287;
t270 = (Icges(6,3) * t252 - t279 * t251) * t239 + (Icges(6,3) * t251 + t279 * t252) * t238 + (Icges(6,5) * t255 + Icges(6,6) * t257) * t254;
t269 = t254 * (-pkin(3) * t251 + pkin(8) * t252) + t272;
t250 = rSges(2,1) * t266 - rSges(2,2) * t264;
t249 = -rSges(2,1) * t264 - rSges(2,2) * t266;
t248 = rSges(5,1) * t263 + rSges(5,2) * t265;
t245 = Icges(5,5) * t263 + Icges(5,6) * t265;
t244 = rSges(6,1) * t255 + rSges(6,2) * t257;
t240 = pkin(3) * t252 + pkin(8) * t251;
t236 = -t287 - t260 * (rSges(3,1) * t258 - rSges(3,2) * t256);
t235 = -t288 + t260 * (-rSges(3,1) * t256 - rSges(3,2) * t258);
t234 = t251 * rSges(5,3) + t286 * t252;
t233 = t252 * rSges(5,3) - t286 * t251;
t228 = Icges(5,3) * t251 + t280 * t252;
t227 = Icges(5,3) * t252 - t280 * t251;
t226 = t251 * rSges(6,3) + t285 * t252;
t225 = t252 * rSges(6,3) - t285 * t251;
t218 = pkin(9) * t251 + t299 * t252;
t217 = pkin(9) * t252 - t299 * t251;
t216 = -t254 * (rSges(4,1) * t252 - rSges(4,2) * t251) + t271;
t215 = t254 * (-rSges(4,1) * t251 - rSges(4,2) * t252) + t272;
t214 = (-t233 * t251 + t234 * t252) * qJD(4);
t213 = t248 * t292 + (-t234 - t240) * t254 + t271;
t212 = t233 * t254 - t248 * t291 + t269;
t211 = t251 * t289 + t238 * t244 + (-t218 - t226 - t240) * t254 + t271;
t210 = -t252 * t289 - t239 * t244 + (t217 + t225) * t254 + t269;
t209 = -t225 * t238 + t226 * t239 + (-t217 * t251 + t218 * t252) * qJD(4);
t1 = m(3) * (t235 ^ 2 + t236 ^ 2) / 0.2e1 + t260 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t215 ^ 2 + t216 ^ 2) / 0.2e1 + t254 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t212 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + ((t252 * t245 + t273 * t251) * t254 + (t252 ^ 2 * t227 + (t275 * t251 + (t228 - t276) * t252) * t251) * qJD(4)) * t291 / 0.2e1 + ((t251 * t245 - t273 * t252) * t254 + (t251 ^ 2 * t228 + (t276 * t252 + (t227 - t275) * t251) * t252) * qJD(4)) * t292 / 0.2e1 + m(6) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + t239 * (t303 * t251 + t270 * t252) / 0.2e1 + t238 * (t270 * t251 - t303 * t252) / 0.2e1 + (((t229 * t265 + t231 * t263) * t252 + (t230 * t265 + t232 * t263) * t251) * qJD(4) + (t221 * t257 + t223 * t255) * t239 + (t222 * t257 + t224 * t255) * t238 + (t257 * t242 + t255 * t243 + t265 * t246 + t263 * t247) * t254) * t254 / 0.2e1 + (m(2) * (t249 ^ 2 + t250 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
