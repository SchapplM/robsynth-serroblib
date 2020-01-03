% Calculate kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:43
% EndTime: 2019-12-31 16:26:43
% DurationCPUTime: 0.69s
% Computational Cost: add. (418->78), mult. (453->127), div. (0->0), fcn. (376->4), ass. (0->53)
t282 = Icges(4,4) + Icges(5,4);
t281 = Icges(4,1) + Icges(5,1);
t280 = Icges(4,2) + Icges(5,2);
t224 = cos(qJ(3));
t279 = t282 * t224;
t223 = sin(qJ(3));
t278 = t282 * t223;
t277 = Icges(4,5) + Icges(5,5);
t276 = Icges(4,6) + Icges(5,6);
t275 = -t280 * t223 + t279;
t274 = t281 * t224 - t278;
t273 = rSges(5,1) + pkin(3);
t272 = Icges(4,3) + Icges(5,3);
t221 = pkin(6) + qJ(2);
t219 = sin(t221);
t220 = cos(t221);
t271 = t275 * t219 - t276 * t220;
t270 = t276 * t219 + t275 * t220;
t269 = -t274 * t219 + t277 * t220;
t268 = t277 * t219 + t274 * t220;
t267 = t280 * t224 + t278;
t266 = t281 * t223 + t279;
t265 = -t276 * t223 + t277 * t224;
t264 = rSges(5,3) + qJ(4);
t263 = -rSges(5,2) * t223 + t273 * t224;
t262 = t265 * t219 - t272 * t220;
t261 = t272 * t219 + t265 * t220;
t260 = t277 * t223 + t276 * t224;
t259 = -t267 * t223 + t266 * t224;
t258 = -t270 * t223 + t268 * t224;
t257 = t271 * t223 + t269 * t224;
t247 = t263 * t219 - t264 * t220;
t246 = t264 * t219 + t263 * t220;
t245 = qJD(3) * t219;
t244 = qJD(3) * t220;
t241 = rSges(4,1) * t224 - rSges(4,2) * t223;
t239 = qJD(3) * (-rSges(5,2) * t224 - t273 * t223);
t226 = qJD(1) ^ 2;
t225 = qJD(2) ^ 2;
t217 = rSges(4,1) * t223 + rSges(4,2) * t224;
t209 = pkin(2) * t219 - pkin(5) * t220;
t208 = rSges(3,1) * t220 - rSges(3,2) * t219;
t207 = rSges(3,1) * t219 + rSges(3,2) * t220;
t206 = qJD(2) * (pkin(2) * t220 + pkin(5) * t219);
t205 = rSges(4,3) * t219 + t241 * t220;
t203 = -rSges(4,3) * t220 + t241 * t219;
t187 = qJD(2) * t205 - t217 * t245 + t206;
t186 = -t217 * t244 + (-t203 - t209) * qJD(2);
t185 = qJD(1) + (t203 * t219 + t205 * t220) * qJD(3);
t184 = t246 * qJD(2) - qJD(4) * t220 + t219 * t239 + t206;
t183 = qJD(4) * t219 + t220 * t239 + (-t209 - t247) * qJD(2);
t182 = qJD(1) + (t247 * t219 + t246 * t220) * qJD(3);
t1 = m(2) * t226 / 0.2e1 + m(3) * (t226 + (t207 ^ 2 + t208 ^ 2) * t225) / 0.2e1 + t225 * Icges(3,3) / 0.2e1 + m(4) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(5) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + (((t269 * t223 - t271 * t224) * t220 + (t268 * t223 + t270 * t224) * t219) * qJD(3) + (t266 * t223 + t267 * t224) * qJD(2)) * qJD(2) / 0.2e1 + ((t261 * t219 ^ 2 + (t257 * t220 + (t258 - t262) * t219) * t220) * qJD(3) + (t260 * t219 + t259 * t220) * qJD(2)) * t245 / 0.2e1 - ((t262 * t220 ^ 2 + (t258 * t219 + (t257 - t261) * t220) * t219) * qJD(3) + (t259 * t219 - t260 * t220) * qJD(2)) * t244 / 0.2e1;
T = t1;
