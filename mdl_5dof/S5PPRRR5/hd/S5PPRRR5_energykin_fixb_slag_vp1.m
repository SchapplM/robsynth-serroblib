% Calculate kinetic energy for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:37
% EndTime: 2019-12-31 17:35:38
% DurationCPUTime: 0.41s
% Computational Cost: add. (431->61), mult. (474->114), div. (0->0), fcn. (514->8), ass. (0->42)
t250 = sin(pkin(8));
t251 = cos(pkin(8));
t271 = qJ(3) + qJ(4);
t266 = sin(t271);
t267 = cos(t271);
t236 = -t250 * t266 - t251 * t267;
t282 = t236 ^ 2;
t237 = -t250 * t267 + t251 * t266;
t281 = t237 ^ 2;
t252 = sin(qJ(5));
t254 = cos(qJ(5));
t280 = -Icges(6,5) * t252 - Icges(6,6) * t254;
t278 = t236 * t237;
t249 = -qJD(3) - qJD(4);
t277 = t280 * t249;
t255 = cos(qJ(3));
t276 = t255 * pkin(3);
t253 = sin(qJ(3));
t273 = t250 * t253;
t272 = t251 * t253;
t248 = qJD(2) * t250;
t270 = qJD(3) * (-pkin(3) * t272 + t276 * t250) + t248;
t269 = qJD(2) * t251;
t268 = qJD(5) * (-t252 * rSges(6,1) - t254 * rSges(6,2));
t265 = -rSges(6,1) * t254 + rSges(6,2) * t252;
t262 = -Icges(6,5) * t254 + Icges(6,6) * t252;
t258 = -qJD(3) * (pkin(3) * t273 + t276 * t251) - t269;
t257 = qJD(1) ^ 2;
t239 = t250 * t255 - t272;
t238 = -t251 * t255 - t273;
t233 = -t269 - qJD(3) * (-t238 * rSges(4,1) + t239 * rSges(4,2));
t232 = t248 + qJD(3) * (t239 * rSges(4,1) + t238 * rSges(4,2));
t231 = t237 * rSges(6,3) + t265 * t236;
t230 = -t236 * rSges(6,3) + t265 * t237;
t225 = Icges(6,3) * t237 + t262 * t236;
t224 = -Icges(6,3) * t236 + t262 * t237;
t223 = t249 * (-t236 * rSges(5,1) - t237 * rSges(5,2)) + t258;
t222 = -t249 * (-t237 * rSges(5,1) + t236 * rSges(5,2)) + t270;
t221 = qJD(1) + (t230 * t237 + t231 * t236) * qJD(5);
t220 = -t237 * t268 + (-t236 * pkin(4) + t237 * pkin(7) + t231) * t249 + t258;
t219 = -t236 * t268 + (t237 * pkin(4) + t236 * pkin(7) - t230) * t249 + t270;
t1 = m(2) * t257 / 0.2e1 + m(3) * (t257 + (t250 ^ 2 + t251 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t232 ^ 2 + t233 ^ 2 + t257) / 0.2e1 + qJD(3) ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t222 ^ 2 + t223 ^ 2 + t257) / 0.2e1 + t249 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + qJD(5) * t237 * (t237 * t277 + (-t224 * t278 + t281 * t225) * qJD(5)) / 0.2e1 - qJD(5) * t236 * (-t236 * t277 + (t282 * t224 - t225 * t278) * qJD(5)) / 0.2e1 + t249 * ((t254 ^ 2 * Icges(6,2) + (Icges(6,1) * t252 + 0.2e1 * Icges(6,4) * t254) * t252) * t249 + (t281 + t282) * t280 * qJD(5)) / 0.2e1;
T = t1;
