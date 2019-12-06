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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:20:02
% EndTime: 2019-12-05 18:20:03
% DurationCPUTime: 0.44s
% Computational Cost: add. (611->110), mult. (528->183), div. (0->0), fcn. (498->10), ass. (0->57)
t253 = sin(pkin(9));
t252 = qJ(1) + qJ(2);
t247 = pkin(8) + t252;
t245 = sin(t247);
t246 = cos(t247);
t257 = cos(qJ(5));
t254 = cos(pkin(9));
t255 = sin(qJ(5));
t272 = t254 * t255;
t230 = t245 * t272 + t246 * t257;
t271 = t254 * t257;
t231 = -t245 * t271 + t246 * t255;
t232 = t245 * t257 - t246 * t272;
t233 = t245 * t255 + t246 * t271;
t273 = t246 * t253;
t274 = t245 * t253;
t262 = (Icges(6,5) * t231 + Icges(6,6) * t230 - Icges(6,3) * t274) * t245 - (Icges(6,5) * t233 + Icges(6,6) * t232 + Icges(6,3) * t273) * t246;
t281 = t262 * t253;
t251 = qJD(1) + qJD(2);
t270 = qJD(5) * t253;
t280 = -t251 * (pkin(4) * t254 + pkin(7) * t253) + (-rSges(6,3) * t254 + (rSges(6,1) * t257 - rSges(6,2) * t255) * t253) * t270;
t248 = sin(t252);
t277 = pkin(2) * t248;
t249 = cos(t252);
t276 = pkin(2) * t249;
t275 = pkin(1) * qJD(1);
t256 = sin(qJ(1));
t269 = t256 * t275;
t258 = cos(qJ(1));
t268 = t258 * t275;
t266 = -pkin(3) * t246 - qJ(4) * t245 - t276;
t265 = qJD(4) * t246 - t268;
t263 = -rSges(5,1) * t254 + rSges(5,2) * t253;
t261 = t251 * (-pkin(3) * t245 + qJ(4) * t246) + qJD(4) * t245 - t269;
t259 = qJD(3) ^ 2;
t242 = -t254 * qJD(5) + t251;
t241 = rSges(2,1) * t258 - rSges(2,2) * t256;
t240 = -rSges(2,1) * t256 - rSges(2,2) * t258;
t236 = -Icges(6,5) * t254 + (Icges(6,1) * t257 - Icges(6,4) * t255) * t253;
t235 = -Icges(6,6) * t254 + (Icges(6,4) * t257 - Icges(6,2) * t255) * t253;
t234 = -Icges(6,3) * t254 + (Icges(6,5) * t257 - Icges(6,6) * t255) * t253;
t229 = -t268 - t251 * (rSges(3,1) * t249 - rSges(3,2) * t248);
t228 = -t269 + t251 * (-rSges(3,1) * t248 - rSges(3,2) * t249);
t227 = -t268 + (-rSges(4,1) * t246 + rSges(4,2) * t245 - t276) * t251;
t226 = -t269 + (-rSges(4,1) * t245 - rSges(4,2) * t246 - t277) * t251;
t225 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t273;
t224 = rSges(6,1) * t231 + rSges(6,2) * t230 - rSges(6,3) * t274;
t223 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t273;
t222 = Icges(6,1) * t231 + Icges(6,4) * t230 - Icges(6,5) * t274;
t221 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t273;
t220 = Icges(6,4) * t231 + Icges(6,2) * t230 - Icges(6,6) * t274;
t217 = (-t245 * rSges(5,3) + t263 * t246 + t266) * t251 + t265;
t216 = (t246 * rSges(5,3) + t263 * t245 - t277) * t251 + t261;
t215 = qJD(3) + (-t224 * t246 - t225 * t245) * t270;
t214 = -t242 * t225 + t280 * t246 + t266 * t251 + t265;
t213 = t242 * t224 + t280 * t245 - t251 * t277 + t261;
t1 = m(3) * (t228 ^ 2 + t229 ^ 2) / 0.2e1 + m(4) * (t226 ^ 2 + t227 ^ 2 + t259) / 0.2e1 + m(5) * (t216 ^ 2 + t217 ^ 2 + t259) / 0.2e1 + m(6) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + t242 * ((-t254 * t234 + (-t235 * t255 + t236 * t257) * t253) * t242 + ((-(-t220 * t255 + t222 * t257) * t245 + (-t221 * t255 + t223 * t257) * t246) * t253 + t262 * t254) * t270) / 0.2e1 - t245 * ((t230 * t235 + t231 * t236 - t234 * t274) * t242 + ((t221 * t230 + t223 * t231) * t246 + (-t230 * t220 - t231 * t222 + t281) * t245) * t270) * t270 / 0.2e1 + t246 * ((t232 * t235 + t233 * t236 + t234 * t273) * t242 + (-(t220 * t232 + t222 * t233) * t245 + (t232 * t221 + t233 * t223 - t281) * t246) * t270) * t270 / 0.2e1 + (m(2) * (t240 ^ 2 + t241 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,2) * t254 ^ 2 + (Icges(5,1) * t253 + 0.2e1 * Icges(5,4) * t254) * t253) * t251 ^ 2 / 0.2e1;
T = t1;
