% Calculate kinetic energy for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:09
% EndTime: 2019-12-31 17:00:10
% DurationCPUTime: 1.01s
% Computational Cost: add. (299->98), mult. (714->156), div. (0->0), fcn. (601->4), ass. (0->66)
t329 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t328 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t327 = -Icges(3,2) - Icges(5,2) - Icges(4,3);
t256 = cos(qJ(2));
t326 = t329 * t256;
t254 = sin(qJ(2));
t325 = t329 * t254;
t324 = -Icges(4,4) + Icges(3,5) + Icges(5,5);
t323 = Icges(5,4) + Icges(4,5) - Icges(3,6);
t322 = t328 * t256 - t325;
t321 = t327 * t254 + t326;
t320 = rSges(5,3) + qJ(4);
t319 = Icges(4,1) + Icges(5,1) + Icges(3,3);
t255 = sin(qJ(1));
t257 = cos(qJ(1));
t318 = t321 * t255 + t323 * t257;
t317 = -t323 * t255 + t321 * t257;
t316 = -t322 * t255 + t324 * t257;
t315 = t324 * t255 + t322 * t257;
t314 = t328 * t254 + t326;
t313 = t327 * t256 - t325;
t312 = t323 * t254 + t324 * t256;
t311 = rSges(5,1) + pkin(3);
t310 = rSges(5,2) * t254 + t320 * t256;
t309 = t312 * t255 - t319 * t257;
t308 = t319 * t255 + t312 * t257;
t307 = t324 * t254 - t323 * t256;
t306 = t313 * t254 + t314 * t256;
t305 = -t317 * t254 + t315 * t256;
t304 = t318 * t254 + t316 * t256;
t292 = t311 * t255 + t310 * t257;
t291 = t310 * t255 - t311 * t257;
t278 = pkin(2) * t256 + qJ(3) * t254;
t231 = t278 * t255;
t251 = t255 * pkin(1) - t257 * pkin(5);
t290 = -t231 - t251;
t289 = qJD(2) * t255;
t288 = qJD(2) * t257;
t287 = qJD(3) * t254;
t232 = t278 * t257;
t235 = qJD(1) * (t257 * pkin(1) + t255 * pkin(5));
t286 = qJD(1) * t232 + t255 * t287 + t235;
t245 = t254 * pkin(2) - t256 * qJ(3);
t283 = qJD(2) * (t254 * rSges(4,2) + t256 * rSges(4,3) - t245);
t282 = -qJD(3) * t256 + t231 * t289 + t232 * t288;
t281 = rSges(3,1) * t256 - rSges(3,2) * t254;
t280 = -rSges(4,2) * t256 + rSges(4,3) * t254;
t259 = qJD(4) * t256 + (t256 * rSges(5,2) - t320 * t254 - t245) * qJD(2);
t253 = t257 * t287;
t250 = t257 * rSges(2,1) - t255 * rSges(2,2);
t248 = t255 * rSges(2,1) + t257 * rSges(2,2);
t247 = t254 * rSges(3,1) + t256 * rSges(3,2);
t229 = -t257 * rSges(4,1) + t280 * t255;
t227 = t255 * rSges(4,1) + t280 * t257;
t225 = t255 * rSges(3,3) + t281 * t257;
t224 = -t257 * rSges(3,3) + t281 * t255;
t203 = qJD(1) * t225 - t247 * t289 + t235;
t202 = -t247 * t288 + (-t224 - t251) * qJD(1);
t201 = (t224 * t255 + t225 * t257) * qJD(2);
t200 = qJD(1) * t227 + t255 * t283 + t286;
t199 = t253 + t257 * t283 + (-t229 + t290) * qJD(1);
t198 = (t227 * t257 + t229 * t255) * qJD(2) + t282;
t197 = t292 * qJD(1) + t259 * t255 + t286;
t196 = t253 + t259 * t257 + (t290 - t291) * qJD(1);
t195 = qJD(4) * t254 + (t291 * t255 + t292 * t257) * qJD(2) + t282;
t1 = m(3) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(4) * (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) / 0.2e1 + m(5) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + (m(2) * (t248 ^ 2 + t250 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t316 * t254 - t318 * t256) * t257 + (t315 * t254 + t317 * t256) * t255) * qJD(2) + (t314 * t254 - t313 * t256) * qJD(1)) * qJD(1) / 0.2e1 + ((t308 * t255 ^ 2 + (t304 * t257 + (t305 - t309) * t255) * t257) * qJD(2) + (t307 * t255 + t306 * t257) * qJD(1)) * t289 / 0.2e1 - ((t309 * t257 ^ 2 + (t305 * t255 + (t304 - t308) * t257) * t255) * qJD(2) + (t306 * t255 - t307 * t257) * qJD(1)) * t288 / 0.2e1;
T = t1;
