% Calculate kinetic energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:24
% EndTime: 2019-12-31 16:18:25
% DurationCPUTime: 0.51s
% Computational Cost: add. (422->90), mult. (575->159), div. (0->0), fcn. (574->6), ass. (0->55)
t230 = sin(pkin(6));
t231 = cos(pkin(6));
t259 = t230 * t231;
t255 = t231 ^ 2;
t256 = t230 ^ 2;
t258 = t255 + t256;
t257 = qJD(3) * t258;
t229 = pkin(7) + qJ(3);
t227 = sin(t229);
t253 = t227 * t230;
t252 = t227 * t231;
t232 = sin(qJ(4));
t251 = t230 * t232;
t233 = cos(qJ(4));
t250 = t230 * t233;
t249 = t231 * t232;
t248 = t231 * t233;
t247 = qJD(2) * t231;
t246 = qJD(3) * t230;
t245 = qJD(3) * t231;
t244 = qJD(4) * t227;
t228 = cos(t229);
t243 = qJD(4) * t228;
t239 = Icges(4,5) * t228 - Icges(4,6) * t227;
t235 = qJD(1) ^ 2;
t226 = qJD(2) * t230;
t225 = t227 * pkin(3) - t228 * pkin(5);
t224 = t227 * rSges(4,1) + t228 * rSges(4,2);
t223 = t230 * t244 - t245;
t222 = t231 * t244 + t246;
t221 = t228 * t248 + t251;
t220 = -t228 * t249 + t250;
t219 = t228 * t250 - t249;
t218 = -t228 * t251 - t248;
t217 = -t224 * t246 - t247;
t216 = -t224 * t245 + t226;
t211 = Icges(4,3) * t230 + t231 * t239;
t210 = -Icges(4,3) * t231 + t230 * t239;
t209 = -t228 * rSges(5,3) + (rSges(5,1) * t233 - rSges(5,2) * t232) * t227;
t208 = -Icges(5,5) * t228 + (Icges(5,1) * t233 - Icges(5,4) * t232) * t227;
t207 = -Icges(5,6) * t228 + (Icges(5,4) * t233 - Icges(5,2) * t232) * t227;
t206 = -Icges(5,3) * t228 + (Icges(5,5) * t233 - Icges(5,6) * t232) * t227;
t205 = t221 * rSges(5,1) + t220 * rSges(5,2) + rSges(5,3) * t252;
t204 = t219 * rSges(5,1) + t218 * rSges(5,2) + rSges(5,3) * t253;
t203 = Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t252;
t202 = Icges(5,1) * t219 + Icges(5,4) * t218 + Icges(5,5) * t253;
t201 = Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t252;
t200 = Icges(5,4) * t219 + Icges(5,2) * t218 + Icges(5,6) * t253;
t199 = Icges(5,5) * t221 + Icges(5,6) * t220 + Icges(5,3) * t252;
t198 = Icges(5,5) * t219 + Icges(5,6) * t218 + Icges(5,3) * t253;
t197 = qJD(1) + (rSges(4,1) * t228 - rSges(4,2) * t227) * t257;
t196 = -t205 * t243 - t222 * t209 - t225 * t246 - t247;
t195 = t204 * t243 + t223 * t209 - t225 * t245 + t226;
t194 = t222 * t204 - t223 * t205 + qJD(1) + (pkin(3) * t228 + pkin(5) * t227) * t257;
t1 = m(2) * t235 / 0.2e1 + m(3) * (t258 * qJD(2) ^ 2 + t235) / 0.2e1 + m(4) * (t197 ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + m(5) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + t222 * ((t199 * t252 + t220 * t201 + t221 * t203) * t222 + (t198 * t252 + t220 * t200 + t221 * t202) * t223 - (t206 * t252 + t220 * t207 + t221 * t208) * t243) / 0.2e1 + t223 * ((t199 * t253 + t218 * t201 + t219 * t203) * t222 + (t198 * t253 + t218 * t200 + t219 * t202) * t223 - (t206 * t253 + t218 * t207 + t219 * t208) * t243) / 0.2e1 - ((-t198 * t223 - t199 * t222 + t206 * t243) * t228 + ((-t201 * t232 + t203 * t233) * t222 + (-t200 * t232 + t202 * t233) * t223 - (-t207 * t232 + t208 * t233) * t243) * t227) * t243 / 0.2e1 + (t230 * (-t210 * t259 + t256 * t211) / 0.2e1 - t231 * (t255 * t210 - t211 * t259) / 0.2e1) * qJD(3) ^ 2;
T = t1;
