% Calculate kinetic energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:18
% EndTime: 2019-12-31 16:24:19
% DurationCPUTime: 0.82s
% Computational Cost: add. (488->132), mult. (853->228), div. (0->0), fcn. (866->8), ass. (0->68)
t256 = cos(pkin(6));
t284 = t256 ^ 2;
t254 = sin(pkin(6));
t285 = t254 ^ 2;
t287 = t284 + t285;
t286 = qJD(2) * t287;
t283 = qJD(2) ^ 2;
t255 = cos(pkin(7));
t282 = t255 * pkin(3);
t258 = sin(qJ(2));
t281 = t254 * t258;
t259 = cos(qJ(2));
t280 = t254 * t259;
t279 = t256 * t258;
t278 = t256 * t259;
t276 = qJD(3) * t258;
t275 = qJD(4) * t258;
t274 = qJD(4) * t259;
t245 = t258 * pkin(2) - t259 * qJ(3);
t271 = qJD(2) * (pkin(5) * t259 - t282 * t258 - t245);
t253 = sin(pkin(7));
t270 = qJD(2) * (t259 * rSges(4,3) - (rSges(4,1) * t255 - rSges(4,2) * t253) * t258 - t245);
t267 = Icges(3,5) * t259 - Icges(3,6) * t258;
t264 = -qJD(3) * t259 + qJD(1) + (pkin(2) * t259 + qJ(3) * t258) * t286;
t252 = pkin(7) + qJ(4);
t251 = cos(t252);
t250 = sin(t252);
t248 = t256 * t276;
t247 = t254 * t276;
t246 = t258 * rSges(3,1) + t259 * rSges(3,2);
t243 = -qJD(2) * t256 + t254 * t275;
t242 = qJD(2) * t254 + t256 * t275;
t241 = t254 * t253 + t255 * t278;
t240 = -t253 * t278 + t254 * t255;
t239 = -t256 * t253 + t255 * t280;
t238 = -t253 * t280 - t256 * t255;
t234 = t254 * t250 + t251 * t278;
t233 = -t250 * t278 + t254 * t251;
t232 = -t256 * t250 + t251 * t280;
t231 = -t250 * t280 - t256 * t251;
t226 = Icges(3,3) * t254 + t267 * t256;
t225 = -Icges(3,3) * t256 + t267 * t254;
t224 = -t259 * rSges(5,3) + (rSges(5,1) * t251 - rSges(5,2) * t250) * t258;
t223 = -Icges(5,5) * t259 + (Icges(5,1) * t251 - Icges(5,4) * t250) * t258;
t222 = -Icges(5,6) * t259 + (Icges(5,4) * t251 - Icges(5,2) * t250) * t258;
t221 = -Icges(5,3) * t259 + (Icges(5,5) * t251 - Icges(5,6) * t250) * t258;
t219 = Icges(4,1) * t241 + Icges(4,4) * t240 + Icges(4,5) * t279;
t218 = Icges(4,1) * t239 + Icges(4,4) * t238 + Icges(4,5) * t281;
t217 = Icges(4,4) * t241 + Icges(4,2) * t240 + Icges(4,6) * t279;
t216 = Icges(4,4) * t239 + Icges(4,2) * t238 + Icges(4,6) * t281;
t215 = Icges(4,5) * t241 + Icges(4,6) * t240 + Icges(4,3) * t279;
t214 = Icges(4,5) * t239 + Icges(4,6) * t238 + Icges(4,3) * t281;
t213 = t256 * t270 + t248;
t212 = t254 * t270 + t247;
t211 = qJD(1) + (rSges(3,1) * t259 - rSges(3,2) * t258) * t286;
t210 = t234 * rSges(5,1) + t233 * rSges(5,2) + rSges(5,3) * t279;
t209 = t232 * rSges(5,1) + t231 * rSges(5,2) + rSges(5,3) * t281;
t208 = Icges(5,1) * t234 + Icges(5,4) * t233 + Icges(5,5) * t279;
t207 = Icges(5,1) * t232 + Icges(5,4) * t231 + Icges(5,5) * t281;
t206 = Icges(5,4) * t234 + Icges(5,2) * t233 + Icges(5,6) * t279;
t205 = Icges(5,4) * t232 + Icges(5,2) * t231 + Icges(5,6) * t281;
t204 = Icges(5,5) * t234 + Icges(5,6) * t233 + Icges(5,3) * t279;
t203 = Icges(5,5) * t232 + Icges(5,6) * t231 + Icges(5,3) * t281;
t202 = (t254 * (t239 * rSges(4,1) + t238 * rSges(4,2) + rSges(4,3) * t281) + t256 * (t241 * rSges(4,1) + t240 * rSges(4,2) + rSges(4,3) * t279)) * qJD(2) + t264;
t201 = t209 * t274 + t243 * t224 + t256 * t271 + t248;
t200 = -t210 * t274 - t242 * t224 + t254 * t271 + t247;
t199 = t242 * t209 - t243 * t210 + t264 + (pkin(5) * t258 + t282 * t259) * t286;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t287 * t283 * t246 ^ 2 + t211 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(5) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + t242 * ((t204 * t279 + t233 * t206 + t234 * t208) * t242 + (t203 * t279 + t233 * t205 + t234 * t207) * t243 - (t221 * t279 + t233 * t222 + t234 * t223) * t274) / 0.2e1 + t243 * ((t204 * t281 + t231 * t206 + t232 * t208) * t242 + (t203 * t281 + t231 * t205 + t232 * t207) * t243 - (t221 * t281 + t231 * t222 + t232 * t223) * t274) / 0.2e1 - ((-t203 * t243 - t204 * t242 + t221 * t274) * t259 + ((-t206 * t250 + t208 * t251) * t242 + (-t205 * t250 + t207 * t251) * t243 - (-t222 * t250 + t223 * t251) * t274) * t258) * t274 / 0.2e1 + (t285 * t226 + (t215 * t279 + t240 * t217 + t241 * t219) * t254 + (-t214 * t279 - t240 * t216 - t241 * t218 - t254 * t225) * t256) * t254 * t283 / 0.2e1 - (t284 * t225 - (t214 * t281 + t238 * t216 + t239 * t218) * t256 + (t215 * t281 + t238 * t217 + t239 * t219 - t256 * t226) * t254) * t256 * t283 / 0.2e1;
T = t1;
