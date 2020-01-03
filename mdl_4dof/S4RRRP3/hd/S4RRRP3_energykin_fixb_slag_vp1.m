% Calculate kinetic energy for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:00
% EndTime: 2019-12-31 17:14:00
% DurationCPUTime: 0.70s
% Computational Cost: add. (425->85), mult. (475->138), div. (0->0), fcn. (387->6), ass. (0->59)
t291 = Icges(4,4) - Icges(5,5);
t290 = Icges(4,1) + Icges(5,1);
t289 = Icges(4,2) + Icges(5,3);
t230 = cos(qJ(3));
t288 = t291 * t230;
t228 = sin(qJ(3));
t287 = t291 * t228;
t286 = Icges(5,4) + Icges(4,5);
t285 = Icges(4,6) - Icges(5,6);
t284 = t289 * t228 - t288;
t283 = t290 * t230 - t287;
t282 = rSges(5,1) + pkin(3);
t281 = rSges(5,3) + qJ(4);
t280 = Icges(5,2) + Icges(4,3);
t227 = qJ(1) + qJ(2);
t224 = sin(t227);
t225 = cos(t227);
t279 = t284 * t224 + t285 * t225;
t278 = -t285 * t224 + t284 * t225;
t277 = -t283 * t224 + t286 * t225;
t276 = t286 * t224 + t283 * t225;
t275 = -t289 * t230 - t287;
t274 = t290 * t228 + t288;
t273 = -t285 * t228 + t286 * t230;
t272 = t281 * t228 + t282 * t230;
t271 = t273 * t224 - t280 * t225;
t270 = t280 * t224 + t273 * t225;
t269 = t286 * t228 + t285 * t230;
t268 = t275 * t228 + t274 * t230;
t267 = t278 * t228 + t276 * t230;
t266 = -t279 * t228 + t277 * t230;
t261 = pkin(1) * qJD(1);
t256 = -t225 * rSges(5,2) + t272 * t224;
t255 = t224 * rSges(5,2) + t272 * t225;
t231 = cos(qJ(1));
t223 = t231 * t261;
t226 = qJD(1) + qJD(2);
t254 = t226 * (t225 * pkin(2) + t224 * pkin(6)) + t223;
t253 = qJD(3) * t224;
t252 = qJD(3) * t225;
t229 = sin(qJ(1));
t251 = t229 * t261;
t248 = rSges(4,1) * t230 - rSges(4,2) * t228;
t233 = t281 * qJD(3) * t230 + (-t282 * qJD(3) + qJD(4)) * t228;
t222 = t231 * rSges(2,1) - t229 * rSges(2,2);
t221 = t229 * rSges(2,1) + t231 * rSges(2,2);
t220 = t228 * rSges(4,1) + t230 * rSges(4,2);
t211 = t224 * pkin(2) - t225 * pkin(6);
t207 = t223 + t226 * (t225 * rSges(3,1) - t224 * rSges(3,2));
t206 = -t251 - t226 * (t224 * rSges(3,1) + t225 * rSges(3,2));
t205 = t224 * rSges(4,3) + t225 * t248;
t203 = -t225 * rSges(4,3) + t224 * t248;
t189 = (t203 * t224 + t205 * t225) * qJD(3);
t188 = t226 * t205 - t220 * t253 + t254;
t187 = -t251 - t220 * t252 + (-t203 - t211) * t226;
t186 = t224 * t233 + t226 * t255 + t254;
t185 = -t251 + (-t211 - t256) * t226 + t233 * t225;
t184 = -qJD(4) * t230 + (t224 * t256 + t225 * t255) * qJD(3);
t1 = m(3) * (t206 ^ 2 + t207 ^ 2) / 0.2e1 + t226 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(5) * (t184 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + ((t274 * t228 - t275 * t230) * t226 + ((t277 * t228 + t279 * t230) * t225 + (t276 * t228 - t278 * t230) * t224) * qJD(3)) * t226 / 0.2e1 + (m(2) * (t221 ^ 2 + t222 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t269 * t224 + t268 * t225) * t226 + (t270 * t224 ^ 2 + (t266 * t225 + (t267 - t271) * t224) * t225) * qJD(3)) * t253 / 0.2e1 - ((t268 * t224 - t269 * t225) * t226 + (t271 * t225 ^ 2 + (t267 * t224 + (t266 - t270) * t225) * t224) * qJD(3)) * t252 / 0.2e1;
T = t1;
