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
% m [6x1]
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:01:38
% EndTime: 2022-01-20 12:01:39
% DurationCPUTime: 0.66s
% Computational Cost: add. (821->125), mult. (513->213), div. (0->0), fcn. (408->10), ass. (0->83)
t259 = qJD(1) + qJD(2);
t296 = pkin(2) * t259;
t264 = cos(qJ(4));
t295 = pkin(4) * t264;
t293 = pkin(1) * qJD(1);
t262 = sin(qJ(4));
t292 = Icges(5,4) * t262;
t291 = Icges(5,4) * t264;
t260 = qJ(4) + qJ(5);
t254 = sin(t260);
t290 = Icges(6,4) * t254;
t256 = cos(t260);
t289 = Icges(6,4) * t256;
t261 = qJ(1) + qJ(2);
t265 = cos(qJ(1));
t252 = t265 * t293;
t257 = cos(t261);
t288 = t257 * t296 + t252;
t258 = qJ(3) + t261;
t249 = sin(t258);
t287 = qJD(4) * t249;
t250 = cos(t258);
t286 = qJD(4) * t250;
t285 = qJD(4) + qJD(5);
t284 = pkin(4) * qJD(4) * t262;
t263 = sin(qJ(1));
t283 = t263 * t293;
t253 = qJD(3) + t259;
t282 = t253 * (pkin(3) * t250 + pkin(8) * t249) + t288;
t281 = rSges(5,1) * t264 - rSges(5,2) * t262;
t280 = rSges(6,1) * t256 - rSges(6,2) * t254;
t279 = Icges(5,1) * t264 - t292;
t278 = Icges(6,1) * t256 - t290;
t277 = -Icges(5,2) * t262 + t291;
t276 = -Icges(6,2) * t254 + t289;
t275 = Icges(5,5) * t264 - Icges(5,6) * t262;
t274 = Icges(6,5) * t256 - Icges(6,6) * t254;
t226 = -Icges(5,6) * t250 + t277 * t249;
t228 = -Icges(5,5) * t250 + t279 * t249;
t273 = t226 * t262 - t228 * t264;
t227 = Icges(5,6) * t249 + t277 * t250;
t229 = Icges(5,5) * t249 + t279 * t250;
t272 = -t227 * t262 + t229 * t264;
t243 = Icges(5,2) * t264 + t292;
t244 = Icges(5,1) * t262 + t291;
t271 = -t243 * t262 + t244 * t264;
t255 = sin(t261);
t270 = -t255 * t296 - t283;
t235 = t285 * t249;
t236 = t285 * t250;
t269 = -(-Icges(6,3) * t250 + t274 * t249) * t236 + (Icges(6,3) * t249 + t274 * t250) * t235 + (Icges(6,5) * t254 + Icges(6,6) * t256) * t253;
t218 = -Icges(6,6) * t250 + t276 * t249;
t219 = Icges(6,6) * t249 + t276 * t250;
t220 = -Icges(6,5) * t250 + t278 * t249;
t221 = Icges(6,5) * t249 + t278 * t250;
t239 = Icges(6,2) * t256 + t290;
t240 = Icges(6,1) * t254 + t289;
t268 = (-t219 * t254 + t221 * t256) * t235 - (-t218 * t254 + t220 * t256) * t236 + (-t239 * t254 + t240 * t256) * t253;
t248 = rSges(2,1) * t265 - rSges(2,2) * t263;
t247 = rSges(2,1) * t263 + rSges(2,2) * t265;
t246 = rSges(5,1) * t262 + rSges(5,2) * t264;
t242 = Icges(5,5) * t262 + Icges(5,6) * t264;
t241 = rSges(6,1) * t254 + rSges(6,2) * t256;
t237 = pkin(3) * t249 - pkin(8) * t250;
t233 = t252 + t259 * (rSges(3,1) * t257 - rSges(3,2) * t255);
t232 = -t283 - t259 * (rSges(3,1) * t255 + rSges(3,2) * t257);
t231 = rSges(5,3) * t249 + t281 * t250;
t230 = -rSges(5,3) * t250 + t281 * t249;
t225 = Icges(5,3) * t249 + t275 * t250;
t224 = -Icges(5,3) * t250 + t275 * t249;
t223 = rSges(6,3) * t249 + t280 * t250;
t222 = -rSges(6,3) * t250 + t280 * t249;
t215 = pkin(9) * t249 + t295 * t250;
t214 = -pkin(9) * t250 + t295 * t249;
t213 = t253 * (rSges(4,1) * t250 - rSges(4,2) * t249) + t288;
t212 = -t253 * (rSges(4,1) * t249 + rSges(4,2) * t250) + t270;
t211 = (t230 * t249 + t231 * t250) * qJD(4);
t210 = t231 * t253 - t246 * t287 + t282;
t209 = -t246 * t286 + (-t230 - t237) * t253 + t270;
t208 = -t249 * t284 - t235 * t241 + (t215 + t223) * t253 + t282;
t207 = -t250 * t284 - t236 * t241 + (-t214 - t222 - t237) * t253 + t270;
t206 = t222 * t235 + t223 * t236 + (t214 * t249 + t215 * t250) * qJD(4);
t1 = m(3) * (t232 ^ 2 + t233 ^ 2) / 0.2e1 + t259 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t212 ^ 2 + t213 ^ 2) / 0.2e1 + t253 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + ((t249 * t242 + t271 * t250) * t253 + (t249 ^ 2 * t225 + (t273 * t250 + (-t224 + t272) * t249) * t250) * qJD(4)) * t287 / 0.2e1 - ((-t250 * t242 + t271 * t249) * t253 + (t250 ^ 2 * t224 + (t272 * t249 + (-t225 + t273) * t250) * t249) * qJD(4)) * t286 / 0.2e1 + m(6) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + t235 * (t269 * t249 + t268 * t250) / 0.2e1 - t236 * (t268 * t249 - t269 * t250) / 0.2e1 + (((t227 * t264 + t229 * t262) * t249 - (t226 * t264 + t228 * t262) * t250) * qJD(4) + (t219 * t256 + t221 * t254) * t235 - (t218 * t256 + t220 * t254) * t236 + (t256 * t239 + t254 * t240 + t264 * t243 + t262 * t244) * t253) * t253 / 0.2e1 + (m(2) * (t247 ^ 2 + t248 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
