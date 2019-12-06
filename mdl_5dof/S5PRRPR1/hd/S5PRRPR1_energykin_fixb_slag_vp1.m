% Calculate kinetic energy for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:33
% EndTime: 2019-12-05 16:15:33
% DurationCPUTime: 0.29s
% Computational Cost: add. (480->87), mult. (293->147), div. (0->0), fcn. (220->8), ass. (0->54)
t246 = cos(pkin(9));
t266 = t246 * pkin(4);
t265 = pkin(2) * qJD(2);
t242 = pkin(9) + qJ(5);
t236 = sin(t242);
t264 = Icges(6,4) * t236;
t238 = cos(t242);
t263 = Icges(6,4) * t238;
t243 = pkin(8) + qJ(2);
t239 = cos(t243);
t232 = t239 * t265;
t240 = qJ(3) + t243;
t233 = sin(t240);
t234 = cos(t240);
t244 = qJD(2) + qJD(3);
t261 = t244 * (t234 * pkin(3) + t233 * qJ(4)) + t232;
t260 = qJD(5) * t234;
t237 = sin(t243);
t259 = t237 * t265;
t258 = qJD(4) * t233 - t259;
t245 = sin(pkin(9));
t257 = rSges(5,1) * t246 - rSges(5,2) * t245;
t256 = rSges(6,1) * t238 - rSges(6,2) * t236;
t255 = Icges(6,1) * t238 - t264;
t254 = -Icges(6,2) * t236 + t263;
t253 = Icges(6,5) * t238 - Icges(6,6) * t236;
t215 = -Icges(6,6) * t234 + t254 * t233;
t217 = -Icges(6,5) * t234 + t255 * t233;
t252 = t215 * t236 - t217 * t238;
t216 = Icges(6,6) * t233 + t254 * t234;
t218 = Icges(6,5) * t233 + t255 * t234;
t251 = -t216 * t236 + t218 * t238;
t226 = Icges(6,2) * t238 + t264;
t227 = Icges(6,1) * t236 + t263;
t250 = -t226 * t236 + t227 * t238;
t249 = qJD(1) ^ 2;
t248 = qJD(2) ^ 2;
t230 = t239 * rSges(3,1) - t237 * rSges(3,2);
t229 = t237 * rSges(3,1) + t239 * rSges(3,2);
t228 = t236 * rSges(6,1) + t238 * rSges(6,2);
t225 = Icges(6,5) * t236 + Icges(6,6) * t238;
t224 = t233 * pkin(3) - t234 * qJ(4);
t222 = t232 + t244 * (t234 * rSges(4,1) - t233 * rSges(4,2));
t221 = -t259 - t244 * (t233 * rSges(4,1) + t234 * rSges(4,2));
t220 = t233 * rSges(6,3) + t256 * t234;
t219 = -t234 * rSges(6,3) + t256 * t233;
t214 = Icges(6,3) * t233 + t253 * t234;
t213 = -Icges(6,3) * t234 + t253 * t233;
t212 = t244 * t233 * rSges(5,3) + (t244 * t257 - qJD(4)) * t234 + t261;
t211 = (t234 * rSges(5,3) - t257 * t233 - t224) * t244 + t258;
t210 = qJD(1) + (t219 * t233 + t220 * t234) * qJD(5);
t209 = t244 * t220 + (t244 * t266 - qJD(4)) * t234 + (pkin(7) * t244 - qJD(5) * t228) * t233 + t261;
t208 = -t228 * t260 + (pkin(7) * t234 - t266 * t233 - t219 - t224) * t244 + t258;
t1 = m(2) * t249 / 0.2e1 + m(3) * (t249 + (t229 ^ 2 + t230 ^ 2) * t248) / 0.2e1 + t248 * Icges(3,3) / 0.2e1 + m(4) * (t221 ^ 2 + t222 ^ 2 + t249) / 0.2e1 + m(5) * (t211 ^ 2 + t212 ^ 2 + t249) / 0.2e1 + m(6) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + qJD(5) * t233 * ((t233 * t225 + t250 * t234) * t244 + (t233 ^ 2 * t214 + (t252 * t234 + (-t213 + t251) * t233) * t234) * qJD(5)) / 0.2e1 - ((-t234 * t225 + t250 * t233) * t244 + (t234 ^ 2 * t213 + (t251 * t233 + (-t214 + t252) * t234) * t233) * qJD(5)) * t260 / 0.2e1 + t244 * ((t238 * t226 + t236 * t227) * t244 + ((t238 * t216 + t236 * t218) * t233 - (t238 * t215 + t236 * t217) * t234) * qJD(5)) / 0.2e1 + (Icges(4,3) + Icges(5,2) * t246 ^ 2 + (Icges(5,1) * t245 + 0.2e1 * Icges(5,4) * t246) * t245) * t244 ^ 2 / 0.2e1;
T = t1;
