% Calculate kinetic energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:45
% EndTime: 2019-12-31 16:44:46
% DurationCPUTime: 0.80s
% Computational Cost: add. (402->93), mult. (511->149), div. (0->0), fcn. (423->6), ass. (0->60)
t315 = Icges(4,4) - Icges(5,5);
t314 = Icges(4,1) + Icges(5,1);
t313 = Icges(4,2) + Icges(5,3);
t248 = pkin(6) + qJ(3);
t246 = cos(t248);
t312 = t315 * t246;
t245 = sin(t248);
t311 = t315 * t245;
t310 = Icges(5,4) + Icges(4,5);
t309 = Icges(4,6) - Icges(5,6);
t308 = t313 * t245 - t312;
t307 = t314 * t246 - t311;
t306 = rSges(5,1) + pkin(3);
t305 = rSges(5,3) + qJ(4);
t304 = Icges(5,2) + Icges(4,3);
t252 = sin(qJ(1));
t253 = cos(qJ(1));
t303 = t308 * t252 + t309 * t253;
t302 = -t309 * t252 + t308 * t253;
t301 = -t307 * t252 + t310 * t253;
t300 = t310 * t252 + t307 * t253;
t299 = -t313 * t246 - t311;
t298 = t314 * t245 + t312;
t297 = -t309 * t245 + t310 * t246;
t296 = t305 * t245 + t306 * t246;
t295 = t297 * t252 - t304 * t253;
t294 = t304 * t252 + t297 * t253;
t293 = t310 * t245 + t309 * t246;
t292 = t299 * t245 + t298 * t246;
t291 = t302 * t245 + t300 * t246;
t290 = -t303 * t245 + t301 * t246;
t250 = cos(pkin(6));
t285 = pkin(2) * t250;
t241 = pkin(1) * t252 - qJ(2) * t253;
t279 = pkin(5) * t253 - t285 * t252 - t241;
t278 = -rSges(5,2) * t253 + t296 * t252;
t277 = rSges(5,2) * t252 + t296 * t253;
t276 = qJD(3) * t252;
t275 = qJD(3) * t253;
t240 = qJD(1) * (pkin(1) * t253 + qJ(2) * t252);
t272 = -qJD(2) * t253 + qJD(1) * (pkin(5) * t252 + t285 * t253) + t240;
t249 = sin(pkin(6));
t271 = rSges(3,1) * t250 - rSges(3,2) * t249;
t270 = rSges(4,1) * t246 - rSges(4,2) * t245;
t255 = t305 * qJD(3) * t246 + (-t306 * qJD(3) + qJD(4)) * t245;
t247 = qJD(2) * t252;
t243 = rSges(2,1) * t253 - rSges(2,2) * t252;
t242 = rSges(2,1) * t252 + rSges(2,2) * t253;
t239 = rSges(4,1) * t245 + rSges(4,2) * t246;
t228 = rSges(4,3) * t252 + t270 * t253;
t226 = -rSges(4,3) * t253 + t270 * t252;
t210 = qJD(1) * t252 * rSges(3,3) + t240 + (qJD(1) * t271 - qJD(2)) * t253;
t209 = t247 + (t253 * rSges(3,3) - t271 * t252 - t241) * qJD(1);
t208 = (t226 * t252 + t228 * t253) * qJD(3);
t207 = qJD(1) * t228 - t239 * t276 + t272;
t206 = -t239 * t275 + t247 + (-t226 + t279) * qJD(1);
t205 = -qJD(4) * t246 + (t278 * t252 + t277 * t253) * qJD(3);
t204 = t277 * qJD(1) + t255 * t252 + t272;
t203 = t247 + t255 * t253 + (-t278 + t279) * qJD(1);
t1 = m(3) * (t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(4) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(5) * (t203 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + (((t301 * t245 + t303 * t246) * t253 + (t300 * t245 - t302 * t246) * t252) * qJD(3) + (t298 * t245 - t299 * t246) * qJD(1)) * qJD(1) / 0.2e1 + ((t294 * t252 ^ 2 + (t290 * t253 + (t291 - t295) * t252) * t253) * qJD(3) + (t293 * t252 + t292 * t253) * qJD(1)) * t276 / 0.2e1 - ((t295 * t253 ^ 2 + (t291 * t252 + (t290 - t294) * t253) * t252) * qJD(3) + (t292 * t252 - t293 * t253) * qJD(1)) * t275 / 0.2e1 + (m(2) * (t242 ^ 2 + t243 ^ 2) + Icges(2,3) + Icges(3,2) * t250 ^ 2 + (Icges(3,1) * t249 + 0.2e1 * Icges(3,4) * t250) * t249) * qJD(1) ^ 2 / 0.2e1;
T = t1;
