% Calculate kinetic energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:17
% EndTime: 2022-01-20 09:51:18
% DurationCPUTime: 0.35s
% Computational Cost: add. (500->95), mult. (318->156), div. (0->0), fcn. (232->10), ass. (0->59)
t247 = qJ(1) + qJ(2);
t242 = sin(t247);
t274 = pkin(2) * t242;
t249 = cos(pkin(9));
t273 = t249 * pkin(4);
t272 = pkin(1) * qJD(1);
t245 = pkin(9) + qJ(5);
t239 = sin(t245);
t271 = Icges(6,4) * t239;
t240 = cos(t245);
t270 = Icges(6,4) * t240;
t252 = cos(qJ(1));
t238 = t252 * t272;
t243 = cos(t247);
t246 = qJD(1) + qJD(2);
t268 = t246 * pkin(2) * t243 + t238;
t241 = pkin(8) + t247;
t236 = cos(t241);
t267 = qJD(5) * t236;
t251 = sin(qJ(1));
t266 = t251 * t272;
t235 = sin(t241);
t265 = t246 * (t236 * pkin(3) + t235 * qJ(4)) + t268;
t264 = -t235 * pkin(3) + t236 * qJ(4) - t274;
t263 = qJD(4) * t235 - t266;
t248 = sin(pkin(9));
t262 = rSges(5,1) * t249 - rSges(5,2) * t248;
t261 = rSges(6,1) * t240 - rSges(6,2) * t239;
t260 = Icges(6,1) * t240 - t271;
t259 = -Icges(6,2) * t239 + t270;
t258 = Icges(6,5) * t240 - Icges(6,6) * t239;
t217 = -Icges(6,6) * t236 + t259 * t235;
t219 = -Icges(6,5) * t236 + t260 * t235;
t257 = t217 * t239 - t219 * t240;
t218 = Icges(6,6) * t235 + t259 * t236;
t220 = Icges(6,5) * t235 + t260 * t236;
t256 = -t218 * t239 + t220 * t240;
t228 = Icges(6,2) * t240 + t271;
t229 = Icges(6,1) * t239 + t270;
t255 = -t228 * t239 + t229 * t240;
t253 = qJD(3) ^ 2;
t233 = t252 * rSges(2,1) - t251 * rSges(2,2);
t232 = t251 * rSges(2,1) + t252 * rSges(2,2);
t230 = t239 * rSges(6,1) + t240 * rSges(6,2);
t227 = Icges(6,5) * t239 + Icges(6,6) * t240;
t224 = t238 + t246 * (t243 * rSges(3,1) - t242 * rSges(3,2));
t223 = -t266 - t246 * (t242 * rSges(3,1) + t243 * rSges(3,2));
t222 = t235 * rSges(6,3) + t261 * t236;
t221 = -t236 * rSges(6,3) + t261 * t235;
t216 = Icges(6,3) * t235 + t258 * t236;
t215 = -Icges(6,3) * t236 + t258 * t235;
t214 = t246 * (t236 * rSges(4,1) - t235 * rSges(4,2)) + t268;
t213 = -t266 + (-t235 * rSges(4,1) - t236 * rSges(4,2) - t274) * t246;
t212 = t246 * t235 * rSges(5,3) + (t246 * t262 - qJD(4)) * t236 + t265;
t211 = (t236 * rSges(5,3) - t262 * t235 + t264) * t246 + t263;
t210 = qJD(3) + (t221 * t235 + t222 * t236) * qJD(5);
t209 = t246 * t222 + (t246 * t273 - qJD(4)) * t236 + (pkin(7) * t246 - qJD(5) * t230) * t235 + t265;
t208 = -t230 * t267 + (pkin(7) * t236 - t273 * t235 - t221 + t264) * t246 + t263;
t1 = m(3) * (t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(4) * (t213 ^ 2 + t214 ^ 2 + t253) / 0.2e1 + m(5) * (t211 ^ 2 + t212 ^ 2 + t253) / 0.2e1 + m(6) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + qJD(5) * t235 * ((t235 * t227 + t255 * t236) * t246 + (t235 ^ 2 * t216 + (t257 * t236 + (-t215 + t256) * t235) * t236) * qJD(5)) / 0.2e1 - ((-t236 * t227 + t255 * t235) * t246 + (t236 ^ 2 * t215 + (t256 * t235 + (-t216 + t257) * t236) * t235) * qJD(5)) * t267 / 0.2e1 + t246 * ((t240 * t228 + t239 * t229) * t246 + ((t240 * t218 + t239 * t220) * t235 - (t240 * t217 + t239 * t219) * t236) * qJD(5)) / 0.2e1 + (m(2) * (t232 ^ 2 + t233 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,2) * t249 ^ 2 + (Icges(5,1) * t248 + 0.2e1 * Icges(5,4) * t249) * t248) * t246 ^ 2 / 0.2e1;
T = t1;
