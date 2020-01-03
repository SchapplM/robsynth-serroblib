% Calculate kinetic energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:05
% EndTime: 2020-01-03 11:45:06
% DurationCPUTime: 1.05s
% Computational Cost: add. (687->99), mult. (505->151), div. (0->0), fcn. (398->8), ass. (0->68)
t347 = Icges(5,4) + Icges(6,4);
t346 = Icges(5,1) + Icges(6,1);
t345 = Icges(5,2) + Icges(6,2);
t282 = cos(qJ(4));
t344 = t347 * t282;
t280 = sin(qJ(4));
t343 = t347 * t280;
t342 = Icges(5,5) + Icges(6,5);
t341 = Icges(5,6) + Icges(6,6);
t340 = t345 * t280 - t344;
t339 = t346 * t282 - t343;
t338 = rSges(6,1) + pkin(4);
t337 = -Icges(5,3) - Icges(6,3);
t278 = qJ(1) + pkin(8);
t276 = qJ(3) + t278;
t269 = sin(t276);
t270 = cos(t276);
t336 = t340 * t269 + t341 * t270;
t335 = -t341 * t269 + t340 * t270;
t334 = t339 * t269 - t342 * t270;
t333 = t342 * t269 + t339 * t270;
t332 = t345 * t282 + t343;
t331 = t346 * t280 + t344;
t330 = -t341 * t280 + t342 * t282;
t329 = rSges(6,3) + qJ(5);
t328 = -rSges(6,2) * t280 + t338 * t282;
t327 = t330 * t269 + t337 * t270;
t326 = t337 * t269 - t330 * t270;
t325 = -t342 * t280 - t341 * t282;
t324 = t332 * t280 - t331 * t282;
t323 = t335 * t280 + t333 * t282;
t322 = t336 * t280 + t334 * t282;
t315 = pkin(1) * qJD(1);
t314 = pkin(2) * qJD(1);
t309 = t328 * t269 - t329 * t270;
t308 = t329 * t269 + t328 * t270;
t281 = sin(qJ(1));
t272 = t281 * t315;
t274 = sin(t278);
t307 = t274 * t314 + t272;
t283 = cos(qJ(1));
t273 = t283 * t315;
t275 = cos(t278);
t306 = t275 * t314 + t273;
t305 = qJD(4) * t269;
t304 = qJD(4) * t270;
t277 = qJD(1) + qJD(3);
t303 = t277 * (pkin(3) * t269 - pkin(7) * t270) + t307;
t300 = rSges(6,2) * t282 + t338 * t280;
t299 = rSges(5,1) * t282 - rSges(5,2) * t280;
t284 = qJD(2) ^ 2;
t266 = -rSges(2,1) * t283 + rSges(2,2) * t281;
t265 = rSges(2,1) * t281 + rSges(2,2) * t283;
t264 = rSges(5,1) * t280 + rSges(5,2) * t282;
t256 = -pkin(3) * t270 - pkin(7) * t269;
t254 = t273 - qJD(1) * (-rSges(3,1) * t275 + rSges(3,2) * t274);
t253 = t272 + qJD(1) * (rSges(3,1) * t274 + rSges(3,2) * t275);
t252 = -rSges(5,3) * t269 - t299 * t270;
t250 = -rSges(5,3) * t270 + t299 * t269;
t236 = -t277 * (-rSges(4,1) * t270 + rSges(4,2) * t269) + t306;
t235 = t277 * (rSges(4,1) * t269 + rSges(4,2) * t270) + t307;
t232 = qJD(2) + (t250 * t269 - t252 * t270) * qJD(4);
t231 = -t264 * t305 + (-t252 - t256) * t277 + t306;
t230 = t250 * t277 + t264 * t304 + t303;
t229 = -qJD(5) * t270 - t300 * t305 + (-t256 + t308) * t277 + t306;
t228 = -qJD(5) * t269 + t309 * t277 + t300 * t304 + t303;
t227 = qJD(2) + (t309 * t269 + t308 * t270) * qJD(4);
t1 = m(3) * (t253 ^ 2 + t254 ^ 2 + t284) / 0.2e1 + m(4) * (t235 ^ 2 + t236 ^ 2 + t284) / 0.2e1 + t277 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(6) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + ((t331 * t280 + t332 * t282) * t277 + ((-t334 * t280 + t336 * t282) * t270 + (t333 * t280 - t335 * t282) * t269) * qJD(4)) * t277 / 0.2e1 - ((t325 * t269 + t324 * t270) * t277 + (t326 * t269 ^ 2 + (t322 * t270 + (-t323 + t327) * t269) * t270) * qJD(4)) * t305 / 0.2e1 - ((-t324 * t269 + t325 * t270) * t277 + (t327 * t270 ^ 2 + (t323 * t269 + (-t322 + t326) * t270) * t269) * qJD(4)) * t304 / 0.2e1 + (m(2) * (t265 ^ 2 + t266 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
