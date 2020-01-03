% Calculate kinetic energy for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:51
% EndTime: 2019-12-31 16:58:52
% DurationCPUTime: 0.96s
% Computational Cost: add. (298->98), mult. (711->156), div. (0->0), fcn. (598->4), ass. (0->66)
t328 = Icges(3,4) - Icges(5,4) - Icges(4,5);
t327 = Icges(3,1) + Icges(4,1) + Icges(5,1);
t326 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t255 = cos(qJ(2));
t325 = t328 * t255;
t253 = sin(qJ(2));
t324 = t328 * t253;
t323 = Icges(4,4) + Icges(3,5) - Icges(5,5);
t322 = Icges(3,6) - Icges(4,6) + Icges(5,6);
t321 = t326 * t253 - t325;
t320 = t327 * t255 - t324;
t319 = rSges(5,1) + pkin(3);
t318 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t254 = sin(qJ(1));
t256 = cos(qJ(1));
t317 = t321 * t254 + t322 * t256;
t316 = -t322 * t254 + t321 * t256;
t315 = -t320 * t254 + t323 * t256;
t314 = t323 * t254 + t320 * t256;
t313 = -t326 * t255 - t324;
t312 = t327 * t253 + t325;
t311 = -t322 * t253 + t323 * t255;
t310 = rSges(5,3) + qJ(4);
t309 = rSges(5,2) * t253 + t319 * t255;
t308 = t311 * t254 - t318 * t256;
t307 = t318 * t254 + t311 * t256;
t306 = -t323 * t253 - t322 * t255;
t305 = t313 * t253 + t312 * t255;
t304 = t316 * t253 + t314 * t255;
t303 = -t317 * t253 + t315 * t255;
t291 = t309 * t254 + t310 * t256;
t290 = -t310 * t254 + t309 * t256;
t277 = pkin(2) * t255 + qJ(3) * t253;
t230 = t277 * t254;
t250 = pkin(1) * t254 - pkin(5) * t256;
t289 = -t230 - t250;
t288 = qJD(2) * t254;
t287 = qJD(2) * t256;
t286 = qJD(3) * t253;
t231 = t277 * t256;
t234 = qJD(1) * (pkin(1) * t256 + pkin(5) * t254);
t285 = qJD(1) * t231 + t254 * t286 + t234;
t244 = pkin(2) * t253 - qJ(3) * t255;
t282 = qJD(2) * (-rSges(4,1) * t253 + rSges(4,3) * t255 - t244);
t281 = -qJD(3) * t255 + t230 * t288 + t231 * t287;
t280 = rSges(3,1) * t255 - rSges(3,2) * t253;
t279 = rSges(4,1) * t255 + rSges(4,3) * t253;
t258 = qJD(2) * (rSges(5,2) * t255 - t319 * t253 - t244);
t252 = t256 * t286;
t249 = rSges(2,1) * t256 - rSges(2,2) * t254;
t248 = rSges(2,1) * t254 + rSges(2,2) * t256;
t247 = rSges(3,1) * t253 + rSges(3,2) * t255;
t228 = rSges(3,3) * t254 + t280 * t256;
t227 = rSges(4,2) * t254 + t279 * t256;
t225 = -rSges(3,3) * t256 + t280 * t254;
t224 = -rSges(4,2) * t256 + t279 * t254;
t202 = qJD(1) * t228 - t247 * t288 + t234;
t201 = -t247 * t287 + (-t225 - t250) * qJD(1);
t200 = (t225 * t254 + t228 * t256) * qJD(2);
t199 = qJD(1) * t227 + t254 * t282 + t285;
t198 = t252 + t256 * t282 + (-t224 + t289) * qJD(1);
t197 = (t224 * t254 + t227 * t256) * qJD(2) + t281;
t196 = t290 * qJD(1) + qJD(4) * t256 + t254 * t258 + t285;
t195 = -qJD(4) * t254 + t252 + t256 * t258 + (t289 - t291) * qJD(1);
t194 = (t291 * t254 + t290 * t256) * qJD(2) + t281;
t1 = m(3) * (t200 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(4) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + m(5) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + (m(2) * (t248 ^ 2 + t249 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t315 * t253 + t317 * t255) * t256 + (t314 * t253 - t316 * t255) * t254) * qJD(2) + (t312 * t253 - t313 * t255) * qJD(1)) * qJD(1) / 0.2e1 + ((t307 * t254 ^ 2 + (t303 * t256 + (t304 - t308) * t254) * t256) * qJD(2) + (-t306 * t254 + t305 * t256) * qJD(1)) * t288 / 0.2e1 - ((t308 * t256 ^ 2 + (t304 * t254 + (t303 - t307) * t256) * t254) * qJD(2) + (t305 * t254 + t306 * t256) * qJD(1)) * t287 / 0.2e1;
T = t1;
