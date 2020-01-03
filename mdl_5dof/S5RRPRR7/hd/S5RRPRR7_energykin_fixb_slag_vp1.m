% Calculate kinetic energy for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:27
% EndTime: 2019-12-31 20:15:27
% DurationCPUTime: 0.56s
% Computational Cost: add. (580->128), mult. (517->213), div. (0->0), fcn. (412->8), ass. (0->78)
t269 = qJ(1) + qJ(2);
t263 = sin(t269);
t265 = cos(t269);
t268 = qJ(4) + qJ(5);
t264 = cos(t268);
t262 = sin(t268);
t300 = Icges(6,4) * t262;
t285 = Icges(6,2) * t264 + t300;
t227 = Icges(6,6) * t265 + t285 * t263;
t228 = Icges(6,6) * t263 - t285 * t265;
t299 = Icges(6,4) * t264;
t287 = Icges(6,1) * t262 + t299;
t229 = Icges(6,5) * t265 + t287 * t263;
t230 = Icges(6,5) * t263 - t287 * t265;
t296 = qJD(4) + qJD(5);
t246 = t296 * t263;
t247 = t296 * t265;
t249 = -Icges(6,2) * t262 + t299;
t250 = Icges(6,1) * t264 - t300;
t267 = qJD(1) + qJD(2);
t309 = (t227 * t264 + t229 * t262) * t247 + (t228 * t264 + t230 * t262) * t246 + (t249 * t264 + t250 * t262) * t267;
t270 = sin(qJ(4));
t305 = pkin(4) * t270;
t303 = pkin(1) * qJD(1);
t302 = Icges(5,4) * t270;
t272 = cos(qJ(4));
t301 = Icges(5,4) * t272;
t273 = cos(qJ(1));
t261 = t273 * t303;
t298 = t267 * (pkin(2) * t265 + qJ(3) * t263) + t261;
t297 = qJD(4) * t263;
t295 = pkin(4) * qJD(4) * t272;
t271 = sin(qJ(1));
t294 = t271 * t303;
t293 = t267 * t265 * pkin(7) + t298;
t251 = pkin(2) * t263 - qJ(3) * t265;
t292 = -pkin(7) * t263 - t251;
t291 = qJD(3) * t263 - t294;
t290 = rSges(5,1) * t270 + rSges(5,2) * t272;
t289 = rSges(6,1) * t262 + rSges(6,2) * t264;
t288 = Icges(5,1) * t270 + t301;
t286 = Icges(5,2) * t272 + t302;
t284 = Icges(5,5) * t270 + Icges(5,6) * t272;
t283 = Icges(6,5) * t262 + Icges(6,6) * t264;
t235 = Icges(5,6) * t265 + t286 * t263;
t237 = Icges(5,5) * t265 + t288 * t263;
t280 = -t235 * t272 - t237 * t270;
t236 = Icges(5,6) * t263 - t286 * t265;
t238 = Icges(5,5) * t263 - t288 * t265;
t279 = t236 * t272 + t238 * t270;
t254 = -Icges(5,2) * t270 + t301;
t255 = Icges(5,1) * t272 - t302;
t277 = t254 * t272 + t255 * t270;
t276 = (Icges(6,3) * t265 + t283 * t263) * t247 + (Icges(6,3) * t263 - t283 * t265) * t246 + (Icges(6,5) * t264 - Icges(6,6) * t262) * t267;
t259 = rSges(2,1) * t273 - rSges(2,2) * t271;
t258 = rSges(5,1) * t272 - rSges(5,2) * t270;
t257 = rSges(2,1) * t271 + rSges(2,2) * t273;
t253 = Icges(5,5) * t272 - Icges(5,6) * t270;
t252 = rSges(6,1) * t264 - rSges(6,2) * t262;
t244 = pkin(8) * t265 + t263 * t305;
t243 = pkin(8) * t263 - t265 * t305;
t242 = t261 + t267 * (rSges(3,1) * t265 - rSges(3,2) * t263);
t241 = -t294 - t267 * (rSges(3,1) * t263 + rSges(3,2) * t265);
t240 = rSges(5,3) * t263 - t290 * t265;
t239 = rSges(5,3) * t265 + t290 * t263;
t234 = Icges(5,3) * t263 - t284 * t265;
t233 = Icges(5,3) * t265 + t284 * t263;
t232 = rSges(6,3) * t263 - t289 * t265;
t231 = rSges(6,3) * t265 + t289 * t263;
t224 = -qJD(3) * t265 + t267 * (-rSges(4,2) * t265 + rSges(4,3) * t263) + t298;
t223 = (rSges(4,2) * t263 + rSges(4,3) * t265 - t251) * t267 + t291;
t222 = (-t239 * t263 + t240 * t265) * qJD(4);
t221 = t239 * t267 + (-qJD(4) * t258 - qJD(3)) * t265 + t293;
t220 = t258 * t297 + (-t240 + t292) * t267 + t291;
t219 = -t247 * t252 + (t231 + t244) * t267 + (-qJD(3) - t295) * t265 + t293;
t218 = t263 * t295 + t246 * t252 + (-t232 - t243 + t292) * t267 + t291;
t217 = -t231 * t246 + t232 * t247 + (t243 * t265 - t244 * t263) * qJD(4);
t1 = m(3) * (t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(4) * (t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(5) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + qJD(4) * t265 * ((t265 * t253 + t277 * t263) * t267 + (t265 ^ 2 * t233 + (t279 * t263 + (t234 - t280) * t265) * t263) * qJD(4)) / 0.2e1 + ((t263 * t253 - t277 * t265) * t267 + (t263 ^ 2 * t234 + (t280 * t265 + (t233 - t279) * t263) * t265) * qJD(4)) * t297 / 0.2e1 + m(6) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + t247 * (t309 * t263 + t276 * t265) / 0.2e1 + t246 * (t276 * t263 - t309 * t265) / 0.2e1 + (Icges(3,3) + Icges(4,1)) * t267 ^ 2 / 0.2e1 + (((-t235 * t270 + t237 * t272) * t265 + (-t236 * t270 + t238 * t272) * t263) * qJD(4) + (-t227 * t262 + t229 * t264) * t247 + (-t228 * t262 + t230 * t264) * t246 + (-t249 * t262 + t250 * t264 - t254 * t270 + t255 * t272) * t267) * t267 / 0.2e1 + (m(2) * (t257 ^ 2 + t259 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
