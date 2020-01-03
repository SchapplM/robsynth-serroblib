% Calculate kinetic energy for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:25
% EndTime: 2019-12-31 19:49:26
% DurationCPUTime: 0.85s
% Computational Cost: add. (666->99), mult. (505->152), div. (0->0), fcn. (399->8), ass. (0->66)
t347 = Icges(5,4) - Icges(6,5);
t346 = Icges(5,1) + Icges(6,1);
t345 = Icges(5,2) + Icges(6,3);
t282 = cos(qJ(4));
t344 = t347 * t282;
t280 = sin(qJ(4));
t343 = t347 * t280;
t342 = Icges(6,4) + Icges(5,5);
t341 = Icges(5,6) - Icges(6,6);
t340 = t345 * t280 - t344;
t339 = t346 * t282 - t343;
t338 = rSges(6,1) + pkin(4);
t337 = rSges(6,3) + qJ(5);
t336 = Icges(6,2) + Icges(5,3);
t279 = qJ(1) + qJ(2);
t274 = pkin(8) + t279;
t271 = sin(t274);
t272 = cos(t274);
t335 = t340 * t271 + t341 * t272;
t334 = -t341 * t271 + t340 * t272;
t333 = -t339 * t271 + t342 * t272;
t332 = t342 * t271 + t339 * t272;
t331 = -t345 * t282 - t343;
t330 = t346 * t280 + t344;
t329 = -t341 * t280 + t342 * t282;
t328 = t337 * t280 + t338 * t282;
t327 = t329 * t271 - t336 * t272;
t326 = t336 * t271 + t329 * t272;
t325 = t342 * t280 + t341 * t282;
t324 = t331 * t280 + t330 * t282;
t323 = t334 * t280 + t332 * t282;
t322 = -t335 * t280 + t333 * t282;
t275 = sin(t279);
t316 = pkin(2) * t275;
t315 = pkin(1) * qJD(1);
t310 = -t272 * rSges(6,2) + t328 * t271;
t309 = t271 * rSges(6,2) + t328 * t272;
t283 = cos(qJ(1));
t273 = t283 * t315;
t276 = cos(t279);
t278 = qJD(1) + qJD(2);
t308 = t278 * pkin(2) * t276 + t273;
t307 = qJD(4) * t271;
t306 = qJD(4) * t272;
t281 = sin(qJ(1));
t305 = t281 * t315;
t304 = t278 * (t272 * pkin(3) + t271 * pkin(7)) + t308;
t301 = -t271 * pkin(3) + t272 * pkin(7) - t316;
t300 = rSges(5,1) * t282 - rSges(5,2) * t280;
t285 = t337 * qJD(4) * t282 + (-t338 * qJD(4) + qJD(5)) * t280;
t270 = t283 * rSges(2,1) - t281 * rSges(2,2);
t269 = t281 * rSges(2,1) + t283 * rSges(2,2);
t268 = t280 * rSges(5,1) + t282 * rSges(5,2);
t254 = t273 + t278 * (t276 * rSges(3,1) - t275 * rSges(3,2));
t253 = -t305 - t278 * (t275 * rSges(3,1) + t276 * rSges(3,2));
t252 = t271 * rSges(5,3) + t272 * t300;
t250 = -t272 * rSges(5,3) + t271 * t300;
t236 = t278 * (t272 * rSges(4,1) - t271 * rSges(4,2)) + t308;
t235 = -t305 + (-t271 * rSges(4,1) - t272 * rSges(4,2) - t316) * t278;
t234 = qJD(3) + (t250 * t271 + t252 * t272) * qJD(4);
t233 = t278 * t252 - t268 * t307 + t304;
t232 = -t305 - t268 * t306 + (-t250 + t301) * t278;
t231 = t271 * t285 + t278 * t309 + t304;
t230 = -t305 + t285 * t272 + (t301 - t310) * t278;
t229 = -qJD(5) * t282 + qJD(3) + (t271 * t310 + t272 * t309) * qJD(4);
t1 = m(3) * (t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + m(5) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + m(6) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t278 ^ 2 / 0.2e1 + ((t330 * t280 - t331 * t282) * t278 + ((t333 * t280 + t335 * t282) * t272 + (t332 * t280 - t334 * t282) * t271) * qJD(4)) * t278 / 0.2e1 + (m(2) * (t269 ^ 2 + t270 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t325 * t271 + t324 * t272) * t278 + (t326 * t271 ^ 2 + (t322 * t272 + (t323 - t327) * t271) * t272) * qJD(4)) * t307 / 0.2e1 - ((t324 * t271 - t325 * t272) * t278 + (t327 * t272 ^ 2 + (t323 * t271 + (t322 - t326) * t272) * t271) * qJD(4)) * t306 / 0.2e1;
T = t1;
