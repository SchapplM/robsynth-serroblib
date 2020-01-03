% Calculate kinetic energy for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:27
% EndTime: 2019-12-31 19:26:27
% DurationCPUTime: 0.27s
% Computational Cost: add. (391->89), mult. (298->144), div. (0->0), fcn. (212->8), ass. (0->54)
t232 = qJ(1) + qJ(2);
t228 = sin(t232);
t255 = pkin(2) * t228;
t254 = pkin(1) * qJD(1);
t233 = sin(qJ(5));
t253 = Icges(6,4) * t233;
t235 = cos(qJ(5));
t252 = Icges(6,4) * t235;
t236 = cos(qJ(1));
t226 = t236 * t254;
t229 = cos(t232);
t231 = qJD(1) + qJD(2);
t251 = t231 * pkin(2) * t229 + t226;
t227 = pkin(8) + t232;
t224 = sin(t227);
t250 = qJD(5) * t224;
t234 = sin(qJ(1));
t249 = t234 * t254;
t225 = cos(t227);
t248 = t231 * (t225 * pkin(3) + t224 * qJ(4)) + t251;
t247 = -t224 * pkin(3) + t225 * qJ(4) - t255;
t246 = qJD(4) * t224 - t249;
t245 = rSges(6,1) * t233 + rSges(6,2) * t235;
t244 = Icges(6,1) * t233 + t252;
t243 = Icges(6,2) * t235 + t253;
t242 = Icges(6,5) * t233 + Icges(6,6) * t235;
t206 = Icges(6,6) * t225 + t243 * t224;
t208 = Icges(6,5) * t225 + t244 * t224;
t241 = -t206 * t235 - t208 * t233;
t207 = Icges(6,6) * t224 - t243 * t225;
t209 = Icges(6,5) * t224 - t244 * t225;
t240 = t207 * t235 + t209 * t233;
t217 = -Icges(6,2) * t233 + t252;
t218 = Icges(6,1) * t235 - t253;
t239 = t217 * t235 + t218 * t233;
t237 = qJD(3) ^ 2;
t222 = t236 * rSges(2,1) - t234 * rSges(2,2);
t221 = t235 * rSges(6,1) - t233 * rSges(6,2);
t220 = t234 * rSges(2,1) + t236 * rSges(2,2);
t216 = Icges(6,5) * t235 - Icges(6,6) * t233;
t213 = t226 + t231 * (t229 * rSges(3,1) - t228 * rSges(3,2));
t212 = -t249 - t231 * (t228 * rSges(3,1) + t229 * rSges(3,2));
t211 = t224 * rSges(6,3) - t245 * t225;
t210 = t225 * rSges(6,3) + t245 * t224;
t205 = Icges(6,3) * t224 - t242 * t225;
t204 = Icges(6,3) * t225 + t242 * t224;
t203 = t231 * (t225 * rSges(4,1) - t224 * rSges(4,2)) + t251;
t202 = -t249 + (-t224 * rSges(4,1) - t225 * rSges(4,2) - t255) * t231;
t201 = -qJD(4) * t225 + t231 * (-t225 * rSges(5,2) + t224 * rSges(5,3)) + t248;
t200 = (t224 * rSges(5,2) + t225 * rSges(5,3) + t247) * t231 + t246;
t199 = qJD(3) + (-t210 * t224 + t211 * t225) * qJD(5);
t198 = t231 * t210 + (pkin(7) * t231 - qJD(5) * t221 - qJD(4)) * t225 + t248;
t197 = t221 * t250 + (-pkin(7) * t224 - t211 + t247) * t231 + t246;
t1 = m(3) * (t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t203 ^ 2 + t237) / 0.2e1 + m(5) * (t200 ^ 2 + t201 ^ 2 + t237) / 0.2e1 + m(6) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + qJD(5) * t225 * ((t225 * t216 + t239 * t224) * t231 + (t225 ^ 2 * t204 + (t240 * t224 + (t205 - t241) * t225) * t224) * qJD(5)) / 0.2e1 + ((t224 * t216 - t239 * t225) * t231 + (t224 ^ 2 * t205 + (t241 * t225 + (t204 - t240) * t224) * t225) * qJD(5)) * t250 / 0.2e1 + t231 * ((-t233 * t217 + t235 * t218) * t231 + ((-t233 * t206 + t235 * t208) * t225 + (-t233 * t207 + t235 * t209) * t224) * qJD(5)) / 0.2e1 + (m(2) * (t220 ^ 2 + t222 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,1)) * t231 ^ 2 / 0.2e1;
T = t1;
