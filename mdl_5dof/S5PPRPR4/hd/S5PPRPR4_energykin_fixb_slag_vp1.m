% Calculate kinetic energy for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:16
% EndTime: 2019-12-31 17:32:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (329->82), mult. (520->142), div. (0->0), fcn. (578->8), ass. (0->50)
t281 = cos(qJ(3));
t280 = sin(qJ(3));
t257 = cos(pkin(8));
t279 = t257 * pkin(4);
t254 = pkin(8) + qJ(5);
t252 = sin(t254);
t278 = Icges(6,4) * t252;
t253 = cos(t254);
t277 = Icges(6,4) * t253;
t258 = cos(pkin(7));
t275 = qJD(2) * t258;
t256 = sin(pkin(7));
t245 = -t256 * t280 - t258 * t281;
t274 = qJD(5) * t245;
t246 = -t256 * t281 + t258 * t280;
t273 = qJD(5) * t246;
t251 = qJD(2) * t256;
t272 = qJD(3) * (-t246 * pkin(3) - t245 * qJ(4)) + qJD(4) * t246 + t251;
t271 = -qJD(4) * t245 - t275;
t255 = sin(pkin(8));
t270 = rSges(5,1) * t257 - rSges(5,2) * t255;
t269 = -rSges(6,1) * t253 + rSges(6,2) * t252;
t268 = -Icges(6,1) * t253 + t278;
t267 = Icges(6,2) * t252 - t277;
t266 = -Icges(6,5) * t253 + Icges(6,6) * t252;
t229 = -Icges(6,6) * t245 + t267 * t246;
t231 = -Icges(6,5) * t245 + t268 * t246;
t265 = -t229 * t252 + t231 * t253;
t230 = Icges(6,6) * t246 + t267 * t245;
t232 = Icges(6,5) * t246 + t268 * t245;
t264 = t230 * t252 - t232 * t253;
t242 = -Icges(6,2) * t253 - t278;
t243 = -Icges(6,1) * t252 - t277;
t263 = t242 * t252 - t243 * t253;
t262 = qJD(1) ^ 2;
t244 = -t252 * rSges(6,1) - t253 * rSges(6,2);
t241 = -Icges(6,5) * t252 - Icges(6,6) * t253;
t238 = -t245 * pkin(3) + t246 * qJ(4);
t236 = -t275 - qJD(3) * (-t245 * rSges(4,1) - t246 * rSges(4,2));
t235 = t251 + qJD(3) * (-t246 * rSges(4,1) + t245 * rSges(4,2));
t234 = t246 * rSges(6,3) + t269 * t245;
t233 = -t245 * rSges(6,3) + t269 * t246;
t228 = Icges(6,3) * t246 + t266 * t245;
t227 = -Icges(6,3) * t245 + t266 * t246;
t226 = (-t246 * rSges(5,3) + t270 * t245 - t238) * qJD(3) + t271;
t225 = qJD(3) * (-t245 * rSges(5,3) - t270 * t246) + t272;
t224 = qJD(1) + (t233 * t246 + t234 * t245) * qJD(5);
t223 = -t244 * t273 + (-pkin(6) * t246 + t279 * t245 - t234 - t238) * qJD(3) + t271;
t222 = -t244 * t274 + (-pkin(6) * t245 - t279 * t246 + t233) * qJD(3) + t272;
t1 = m(2) * t262 / 0.2e1 + m(3) * (t262 + (t256 ^ 2 + t258 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t235 ^ 2 + t236 ^ 2 + t262) / 0.2e1 + m(5) * (t225 ^ 2 + t226 ^ 2 + t262) / 0.2e1 + m(6) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + (-(t246 * t241 + t263 * t245) * qJD(3) + (t246 ^ 2 * t228 + (t265 * t245 + (-t227 + t264) * t246) * t245) * qJD(5)) * t273 / 0.2e1 - (-(-t245 * t241 + t263 * t246) * qJD(3) + (t245 ^ 2 * t227 + (t264 * t246 + (-t228 + t265) * t245) * t246) * qJD(5)) * t274 / 0.2e1 - qJD(3) * (-(-t253 * t242 - t252 * t243) * qJD(3) + ((-t253 * t230 - t252 * t232) * t246 - (-t253 * t229 - t252 * t231) * t245) * qJD(5)) / 0.2e1 + (Icges(4,3) + t257 ^ 2 * Icges(5,2) + (Icges(5,1) * t255 + 0.2e1 * Icges(5,4) * t257) * t255) * qJD(3) ^ 2 / 0.2e1;
T = t1;
