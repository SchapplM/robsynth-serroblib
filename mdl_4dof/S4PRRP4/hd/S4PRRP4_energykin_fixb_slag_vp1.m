% Calculate kinetic energy for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:44
% EndTime: 2019-12-31 16:27:45
% DurationCPUTime: 0.69s
% Computational Cost: add. (397->78), mult. (454->127), div. (0->0), fcn. (377->4), ass. (0->53)
t285 = Icges(4,4) - Icges(5,5);
t284 = Icges(4,1) + Icges(5,1);
t283 = Icges(4,2) + Icges(5,3);
t228 = cos(qJ(3));
t282 = t285 * t228;
t227 = sin(qJ(3));
t281 = t285 * t227;
t280 = Icges(5,4) + Icges(4,5);
t279 = Icges(4,6) - Icges(5,6);
t278 = t283 * t227 - t282;
t277 = t284 * t228 - t281;
t276 = rSges(5,1) + pkin(3);
t275 = rSges(5,3) + qJ(4);
t274 = Icges(5,2) + Icges(4,3);
t226 = pkin(6) + qJ(2);
t224 = sin(t226);
t225 = cos(t226);
t273 = t278 * t224 + t279 * t225;
t272 = -t279 * t224 + t278 * t225;
t271 = -t277 * t224 + t280 * t225;
t270 = t280 * t224 + t277 * t225;
t269 = -t283 * t228 - t281;
t268 = t284 * t227 + t282;
t267 = -t279 * t227 + t280 * t228;
t266 = t275 * t227 + t276 * t228;
t265 = t267 * t224 - t274 * t225;
t264 = t274 * t224 + t267 * t225;
t263 = t280 * t227 + t279 * t228;
t262 = t269 * t227 + t268 * t228;
t261 = t272 * t227 + t270 * t228;
t260 = -t273 * t227 + t271 * t228;
t252 = -t225 * rSges(5,2) + t266 * t224;
t251 = t224 * rSges(5,2) + t266 * t225;
t250 = qJD(3) * t224;
t249 = qJD(3) * t225;
t246 = rSges(4,1) * t228 - rSges(4,2) * t227;
t231 = t275 * qJD(3) * t228 + (-t276 * qJD(3) + qJD(4)) * t227;
t230 = qJD(1) ^ 2;
t229 = qJD(2) ^ 2;
t223 = t227 * rSges(4,1) + t228 * rSges(4,2);
t214 = t224 * pkin(2) - t225 * pkin(5);
t213 = t225 * rSges(3,1) - t224 * rSges(3,2);
t212 = t224 * rSges(3,1) + t225 * rSges(3,2);
t211 = qJD(2) * (t225 * pkin(2) + t224 * pkin(5));
t208 = t224 * rSges(4,3) + t246 * t225;
t206 = -t225 * rSges(4,3) + t246 * t224;
t192 = qJD(2) * t208 - t223 * t250 + t211;
t191 = -t223 * t249 + (-t206 - t214) * qJD(2);
t190 = qJD(1) + (t206 * t224 + t208 * t225) * qJD(3);
t189 = t251 * qJD(2) + t231 * t224 + t211;
t188 = t231 * t225 + (-t214 - t252) * qJD(2);
t187 = -qJD(4) * t228 + qJD(1) + (t252 * t224 + t251 * t225) * qJD(3);
t1 = m(2) * t230 / 0.2e1 + m(3) * (t230 + (t212 ^ 2 + t213 ^ 2) * t229) / 0.2e1 + t229 * Icges(3,3) / 0.2e1 + m(4) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(5) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + (((t271 * t227 + t273 * t228) * t225 + (t270 * t227 - t272 * t228) * t224) * qJD(3) + (t268 * t227 - t269 * t228) * qJD(2)) * qJD(2) / 0.2e1 + ((t264 * t224 ^ 2 + (t260 * t225 + (t261 - t265) * t224) * t225) * qJD(3) + (t263 * t224 + t262 * t225) * qJD(2)) * t250 / 0.2e1 - ((t265 * t225 ^ 2 + (t261 * t224 + (t260 - t264) * t225) * t224) * qJD(3) + (t262 * t224 - t263 * t225) * qJD(2)) * t249 / 0.2e1;
T = t1;
