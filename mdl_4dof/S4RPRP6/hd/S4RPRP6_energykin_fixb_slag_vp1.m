% Calculate kinetic energy for
% S4RPRP6
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:58
% EndTime: 2019-12-31 16:45:59
% DurationCPUTime: 0.77s
% Computational Cost: add. (211->88), mult. (478->138), div. (0->0), fcn. (390->4), ass. (0->54)
t294 = Icges(4,4) + Icges(5,4);
t293 = Icges(4,1) + Icges(5,1);
t292 = -Icges(4,2) - Icges(5,2);
t234 = cos(qJ(3));
t291 = t294 * t234;
t232 = sin(qJ(3));
t290 = t294 * t232;
t289 = Icges(4,5) + Icges(5,5);
t288 = Icges(4,6) + Icges(5,6);
t287 = t292 * t234 - t290;
t286 = t293 * t232 + t291;
t285 = rSges(5,1) + pkin(3);
t284 = Icges(4,3) + Icges(5,3);
t233 = sin(qJ(1));
t235 = cos(qJ(1));
t283 = t287 * t233 - t288 * t235;
t282 = t288 * t233 + t287 * t235;
t281 = t286 * t233 + t289 * t235;
t280 = t289 * t233 - t286 * t235;
t279 = t292 * t232 + t291;
t278 = t293 * t234 - t290;
t277 = t289 * t232 + t288 * t234;
t276 = rSges(5,3) + qJ(4);
t275 = rSges(5,2) * t234 + t285 * t232;
t274 = t277 * t233 + t284 * t235;
t273 = t284 * t233 - t277 * t235;
t272 = -t288 * t232 + t289 * t234;
t271 = t278 * t232 + t279 * t234;
t270 = t280 * t232 + t282 * t234;
t269 = -t281 * t232 + t283 * t234;
t258 = t275 * t233 + t276 * t235;
t257 = t276 * t233 - t275 * t235;
t217 = qJD(1) * (pkin(1) * t235 + qJ(2) * t233);
t256 = qJD(1) * t235 * pkin(5) + t217;
t255 = qJD(3) * t233;
t252 = -rSges(5,2) * t232 + t285 * t234;
t224 = pkin(1) * t233 - qJ(2) * t235;
t251 = -pkin(5) * t233 - t224;
t250 = rSges(4,1) * t232 + rSges(4,2) * t234;
t230 = qJD(2) * t233;
t228 = rSges(2,1) * t235 - rSges(2,2) * t233;
t227 = rSges(4,1) * t234 - rSges(4,2) * t232;
t225 = rSges(2,1) * t233 + rSges(2,2) * t235;
t214 = rSges(4,3) * t233 - t250 * t235;
t212 = rSges(4,3) * t235 + t250 * t233;
t198 = t217 - qJD(2) * t235 + qJD(1) * (-rSges(3,2) * t235 + rSges(3,3) * t233);
t197 = t230 + (rSges(3,2) * t233 + rSges(3,3) * t235 - t224) * qJD(1);
t196 = (-t212 * t233 + t214 * t235) * qJD(3);
t195 = qJD(1) * t212 + (-qJD(3) * t227 - qJD(2)) * t235 + t256;
t194 = t227 * t255 + t230 + (-t214 + t251) * qJD(1);
t193 = qJD(4) * t233 + t258 * qJD(1) + (-t252 * qJD(3) - qJD(2)) * t235 + t256;
t192 = qJD(4) * t235 + t230 + t252 * t255 + (t251 - t257) * qJD(1);
t191 = (-t258 * t233 + t257 * t235) * qJD(3);
t1 = m(3) * (t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(4) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + (((t283 * t232 + t281 * t234) * t235 + (-t282 * t232 + t280 * t234) * t233) * qJD(3) + (-t279 * t232 + t278 * t234) * qJD(1)) * qJD(1) / 0.2e1 + ((t273 * t233 ^ 2 + (t269 * t235 + (-t270 + t274) * t233) * t235) * qJD(3) + (t272 * t233 - t271 * t235) * qJD(1)) * t255 / 0.2e1 + ((t274 * t235 ^ 2 + (t270 * t233 + (-t269 + t273) * t235) * t233) * qJD(3) + (t271 * t233 + t272 * t235) * qJD(1)) * qJD(3) * t235 / 0.2e1 + (m(2) * (t225 ^ 2 + t228 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
