% Calculate kinetic energy for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:37
% EndTime: 2019-12-31 16:43:37
% DurationCPUTime: 0.75s
% Computational Cost: add. (407->86), mult. (476->135), div. (0->0), fcn. (387->6), ass. (0->57)
t292 = Icges(4,4) - Icges(5,5);
t291 = Icges(4,1) + Icges(5,1);
t290 = Icges(4,2) + Icges(5,3);
t229 = sin(qJ(3));
t289 = t292 * t229;
t231 = cos(qJ(3));
t288 = t292 * t231;
t287 = Icges(5,4) + Icges(4,5);
t286 = Icges(4,6) - Icges(5,6);
t285 = t290 * t229 - t288;
t284 = t291 * t231 - t289;
t283 = rSges(5,1) + pkin(3);
t282 = rSges(5,3) + qJ(4);
t281 = Icges(5,2) + Icges(4,3);
t228 = qJ(1) + pkin(6);
t226 = sin(t228);
t227 = cos(t228);
t280 = t285 * t226 + t286 * t227;
t279 = -t286 * t226 + t285 * t227;
t278 = -t284 * t226 + t287 * t227;
t277 = t287 * t226 + t284 * t227;
t276 = -t290 * t231 - t289;
t275 = t291 * t229 + t288;
t274 = -t286 * t229 + t287 * t231;
t273 = t282 * t229 + t283 * t231;
t272 = t274 * t226 - t281 * t227;
t271 = t281 * t226 + t274 * t227;
t270 = t287 * t229 + t286 * t231;
t269 = t276 * t229 + t275 * t231;
t268 = t279 * t229 + t277 * t231;
t267 = -t280 * t229 + t278 * t231;
t230 = sin(qJ(1));
t263 = t230 * pkin(1);
t257 = -t227 * rSges(5,2) + t273 * t226;
t256 = t226 * rSges(5,2) + t273 * t227;
t232 = cos(qJ(1));
t225 = qJD(1) * t232 * pkin(1);
t255 = qJD(1) * (t227 * pkin(2) + t226 * pkin(5)) + t225;
t254 = qJD(3) * t226;
t253 = qJD(3) * t227;
t250 = -t226 * pkin(2) + t227 * pkin(5) - t263;
t249 = rSges(4,1) * t231 - rSges(4,2) * t229;
t234 = t282 * qJD(3) * t231 + (-t283 * qJD(3) + qJD(4)) * t229;
t224 = t232 * rSges(2,1) - t230 * rSges(2,2);
t223 = t230 * rSges(2,1) + t232 * rSges(2,2);
t222 = t229 * rSges(4,1) + t231 * rSges(4,2);
t209 = t225 + qJD(1) * (t227 * rSges(3,1) - t226 * rSges(3,2));
t208 = (-t226 * rSges(3,1) - t227 * rSges(3,2) - t263) * qJD(1);
t207 = t226 * rSges(4,3) + t249 * t227;
t205 = -t227 * rSges(4,3) + t249 * t226;
t191 = qJD(1) * t207 - t222 * t254 + t255;
t190 = -t222 * t253 + (-t205 + t250) * qJD(1);
t189 = qJD(2) + (t205 * t226 + t207 * t227) * qJD(3);
t188 = t256 * qJD(1) + t234 * t226 + t255;
t187 = t234 * t227 + (t250 - t257) * qJD(1);
t186 = -qJD(4) * t231 + qJD(2) + (t257 * t226 + t256 * t227) * qJD(3);
t1 = m(3) * (qJD(2) ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(4) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + m(5) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + (((t278 * t229 + t280 * t231) * t227 + (t277 * t229 - t279 * t231) * t226) * qJD(3) + (t275 * t229 - t276 * t231) * qJD(1)) * qJD(1) / 0.2e1 + ((t271 * t226 ^ 2 + (t267 * t227 + (t268 - t272) * t226) * t227) * qJD(3) + (t270 * t226 + t269 * t227) * qJD(1)) * t254 / 0.2e1 - ((t272 * t227 ^ 2 + (t268 * t226 + (t267 - t271) * t227) * t226) * qJD(3) + (t269 * t226 - t270 * t227) * qJD(1)) * t253 / 0.2e1 + (m(2) * (t223 ^ 2 + t224 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
