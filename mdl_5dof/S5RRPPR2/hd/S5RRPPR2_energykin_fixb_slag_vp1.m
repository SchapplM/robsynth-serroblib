% Calculate kinetic energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:25
% EndTime: 2022-01-20 10:05:25
% DurationCPUTime: 0.39s
% Computational Cost: add. (611->109), mult. (528->181), div. (0->0), fcn. (498->10), ass. (0->57)
t254 = qJ(1) + qJ(2);
t250 = sin(t254);
t279 = pkin(2) * t250;
t278 = pkin(1) * qJD(1);
t249 = pkin(8) + t254;
t246 = sin(t249);
t255 = sin(pkin(9));
t277 = t246 * t255;
t247 = cos(t249);
t276 = t247 * t255;
t256 = cos(pkin(9));
t257 = sin(qJ(5));
t275 = t256 * t257;
t259 = cos(qJ(5));
t274 = t256 * t259;
t260 = cos(qJ(1));
t248 = t260 * t278;
t251 = cos(t254);
t253 = qJD(1) + qJD(2);
t273 = t253 * pkin(2) * t251 + t248;
t272 = qJD(5) * t255;
t258 = sin(qJ(1));
t271 = t258 * t278;
t270 = t253 * (t247 * pkin(3) + t246 * qJ(4)) + t273;
t268 = -t246 * pkin(3) + t247 * qJ(4) - t279;
t267 = qJD(4) * t246 - t271;
t266 = rSges(5,1) * t256 - rSges(5,2) * t255;
t231 = -t246 * t275 - t247 * t259;
t232 = t246 * t274 - t247 * t257;
t233 = t246 * t259 - t247 * t275;
t234 = t246 * t257 + t247 * t274;
t265 = (Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t277) * t246 + (Icges(6,5) * t234 + Icges(6,6) * t233 + Icges(6,3) * t276) * t247;
t264 = t265 * t255;
t263 = -t253 * (pkin(4) * t256 + pkin(7) * t255) + (-t256 * rSges(6,3) + (rSges(6,1) * t259 - rSges(6,2) * t257) * t255) * t272;
t261 = qJD(3) ^ 2;
t244 = -qJD(5) * t256 + t253;
t243 = t260 * rSges(2,1) - t258 * rSges(2,2);
t242 = t258 * rSges(2,1) + t260 * rSges(2,2);
t237 = -Icges(6,5) * t256 + (Icges(6,1) * t259 - Icges(6,4) * t257) * t255;
t236 = -Icges(6,6) * t256 + (Icges(6,4) * t259 - Icges(6,2) * t257) * t255;
t235 = -Icges(6,3) * t256 + (Icges(6,5) * t259 - Icges(6,6) * t257) * t255;
t230 = t248 + t253 * (t251 * rSges(3,1) - t250 * rSges(3,2));
t229 = -t271 - t253 * (t250 * rSges(3,1) + t251 * rSges(3,2));
t228 = t253 * (t247 * rSges(4,1) - t246 * rSges(4,2)) + t273;
t227 = -t271 + (-t246 * rSges(4,1) - t247 * rSges(4,2) - t279) * t253;
t226 = t234 * rSges(6,1) + t233 * rSges(6,2) + rSges(6,3) * t276;
t225 = t232 * rSges(6,1) + t231 * rSges(6,2) + rSges(6,3) * t277;
t224 = Icges(6,1) * t234 + Icges(6,4) * t233 + Icges(6,5) * t276;
t223 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t277;
t222 = Icges(6,4) * t234 + Icges(6,2) * t233 + Icges(6,6) * t276;
t221 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t277;
t218 = t253 * t246 * rSges(5,3) + (t253 * t266 - qJD(4)) * t247 + t270;
t217 = (t247 * rSges(5,3) - t266 * t246 + t268) * t253 + t267;
t216 = qJD(3) + (t225 * t247 - t226 * t246) * t272;
t215 = t244 * t226 + (-qJD(4) - t263) * t247 + t270;
t214 = -t244 * t225 + t263 * t246 + t268 * t253 + t267;
t1 = m(3) * (t229 ^ 2 + t230 ^ 2) / 0.2e1 + m(4) * (t227 ^ 2 + t228 ^ 2 + t261) / 0.2e1 + m(5) * (t217 ^ 2 + t218 ^ 2 + t261) / 0.2e1 + m(6) * (t214 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 + t244 * ((-t256 * t235 + (-t236 * t257 + t237 * t259) * t255) * t244 + (((-t222 * t257 + t224 * t259) * t247 + (-t221 * t257 + t223 * t259) * t246) * t255 - t265 * t256) * t272) / 0.2e1 + (m(2) * (t242 ^ 2 + t243 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (t247 * ((t233 * t236 + t234 * t237 + t235 * t276) * t244 + ((t233 * t221 + t234 * t223) * t246 + (t233 * t222 + t234 * t224 + t264) * t247) * t272) + t246 * ((t231 * t236 + t232 * t237 + t235 * t277) * t244 + ((t231 * t222 + t232 * t224) * t247 + (t231 * t221 + t232 * t223 + t264) * t246) * t272)) * t272 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,2) * t256 ^ 2 + (Icges(5,1) * t255 + 0.2e1 * Icges(5,4) * t256) * t255) * t253 ^ 2 / 0.2e1;
T = t1;
