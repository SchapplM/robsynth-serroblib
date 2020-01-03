% Calculate kinetic energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:02
% EndTime: 2019-12-31 16:47:02
% DurationCPUTime: 0.75s
% Computational Cost: add. (206->88), mult. (483->138), div. (0->0), fcn. (395->4), ass. (0->54)
t302 = -Icges(4,4) + Icges(5,5);
t301 = Icges(4,1) + Icges(5,1);
t300 = Icges(4,2) + Icges(5,3);
t241 = cos(qJ(3));
t299 = t302 * t241;
t239 = sin(qJ(3));
t298 = t302 * t239;
t297 = Icges(5,4) + Icges(4,5);
t296 = Icges(4,6) - Icges(5,6);
t295 = -t300 * t241 + t298;
t294 = t301 * t239 - t299;
t293 = rSges(5,1) + pkin(3);
t292 = rSges(5,3) + qJ(4);
t291 = Icges(5,2) + Icges(4,3);
t240 = sin(qJ(1));
t242 = cos(qJ(1));
t290 = t295 * t240 - t296 * t242;
t289 = -t296 * t240 - t295 * t242;
t288 = t294 * t240 + t297 * t242;
t287 = t297 * t240 - t294 * t242;
t286 = t300 * t239 + t299;
t285 = t301 * t241 + t298;
t284 = t297 * t239 + t296 * t241;
t283 = t293 * t239 - t292 * t241;
t282 = t284 * t240 + t291 * t242;
t281 = t291 * t240 - t284 * t242;
t280 = -t296 * t239 + t297 * t241;
t279 = t285 * t239 - t286 * t241;
t278 = t287 * t239 - t289 * t241;
t277 = -t288 * t239 + t290 * t241;
t276 = (t292 * t239 + t293 * t241) * qJD(3) - qJD(4) * t241;
t267 = t242 * rSges(5,2) + t283 * t240;
t266 = t240 * rSges(5,2) - t283 * t242;
t224 = qJD(1) * (t242 * pkin(1) + t240 * qJ(2));
t265 = qJD(1) * t242 * pkin(5) + t224;
t263 = qJD(3) * t240;
t231 = t240 * pkin(1) - t242 * qJ(2);
t259 = -pkin(5) * t240 - t231;
t258 = rSges(4,1) * t239 + rSges(4,2) * t241;
t238 = qJD(2) * t240;
t236 = t242 * rSges(2,1) - t240 * rSges(2,2);
t235 = t241 * rSges(4,1) - t239 * rSges(4,2);
t232 = t240 * rSges(2,1) + t242 * rSges(2,2);
t221 = t240 * rSges(4,3) - t258 * t242;
t219 = t242 * rSges(4,3) + t258 * t240;
t205 = t224 - qJD(2) * t242 + qJD(1) * (-t242 * rSges(3,2) + t240 * rSges(3,3));
t204 = t238 + (t240 * rSges(3,2) + t242 * rSges(3,3) - t231) * qJD(1);
t203 = (-t219 * t240 + t221 * t242) * qJD(3);
t202 = qJD(1) * t219 + (-qJD(3) * t235 - qJD(2)) * t242 + t265;
t201 = t235 * t263 + t238 + (-t221 + t259) * qJD(1);
t200 = qJD(4) * t239 + (-t267 * t240 + t266 * t242) * qJD(3);
t199 = t267 * qJD(1) + (-qJD(2) - t276) * t242 + t265;
t198 = t238 + t276 * t240 + (t259 - t266) * qJD(1);
t1 = m(3) * (t204 ^ 2 + t205 ^ 2) / 0.2e1 + m(4) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(5) * (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) / 0.2e1 + (((t290 * t239 + t288 * t241) * t242 + (t289 * t239 + t287 * t241) * t240) * qJD(3) + (t286 * t239 + t285 * t241) * qJD(1)) * qJD(1) / 0.2e1 + ((t281 * t240 ^ 2 + (t277 * t242 + (-t278 + t282) * t240) * t242) * qJD(3) + (t240 * t280 - t242 * t279) * qJD(1)) * t263 / 0.2e1 + ((t282 * t242 ^ 2 + (t278 * t240 + (-t277 + t281) * t242) * t240) * qJD(3) + (t240 * t279 + t242 * t280) * qJD(1)) * qJD(3) * t242 / 0.2e1 + (m(2) * (t232 ^ 2 + t236 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
