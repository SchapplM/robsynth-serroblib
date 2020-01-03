% Calculate kinetic energy for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:04
% EndTime: 2019-12-31 16:55:05
% DurationCPUTime: 0.52s
% Computational Cost: add. (299->117), mult. (488->197), div. (0->0), fcn. (400->6), ass. (0->69)
t221 = sin(qJ(1));
t223 = cos(qJ(1));
t219 = qJ(3) + qJ(4);
t218 = cos(t219);
t217 = sin(t219);
t247 = Icges(5,4) * t217;
t235 = Icges(5,2) * t218 + t247;
t185 = Icges(5,6) * t223 + t235 * t221;
t186 = Icges(5,6) * t221 - t235 * t223;
t246 = Icges(5,4) * t218;
t237 = Icges(5,1) * t217 + t246;
t187 = Icges(5,5) * t223 + t237 * t221;
t188 = Icges(5,5) * t221 - t237 * t223;
t202 = -Icges(5,2) * t217 + t246;
t203 = Icges(5,1) * t218 - t247;
t243 = qJD(3) + qJD(4);
t206 = t243 * t221;
t207 = t243 * t223;
t254 = (t185 * t218 + t187 * t217) * t207 + (t186 * t218 + t188 * t217) * t206 + (t202 * t218 + t203 * t217) * qJD(1);
t220 = sin(qJ(3));
t252 = pkin(3) * t220;
t249 = Icges(4,4) * t220;
t222 = cos(qJ(3));
t248 = Icges(4,4) * t222;
t205 = qJD(1) * (pkin(1) * t223 + qJ(2) * t221);
t245 = qJD(1) * t223 * pkin(5) + t205;
t244 = qJD(3) * t221;
t242 = pkin(3) * qJD(3) * t222;
t211 = pkin(1) * t221 - qJ(2) * t223;
t241 = -pkin(5) * t221 - t211;
t240 = rSges(4,1) * t220 + rSges(4,2) * t222;
t239 = rSges(5,1) * t217 + rSges(5,2) * t218;
t238 = Icges(4,1) * t220 + t248;
t236 = Icges(4,2) * t222 + t249;
t234 = Icges(4,5) * t220 + Icges(4,6) * t222;
t233 = Icges(5,5) * t217 + Icges(5,6) * t218;
t193 = Icges(4,6) * t223 + t236 * t221;
t195 = Icges(4,5) * t223 + t238 * t221;
t230 = -t193 * t222 - t195 * t220;
t194 = Icges(4,6) * t221 - t236 * t223;
t196 = Icges(4,5) * t221 - t238 * t223;
t229 = t194 * t222 + t196 * t220;
t209 = -Icges(4,2) * t220 + t248;
t210 = Icges(4,1) * t222 - t249;
t227 = t209 * t222 + t210 * t220;
t226 = qJD(1) * (Icges(5,5) * t218 - Icges(5,6) * t217) + (Icges(5,3) * t223 + t233 * t221) * t207 + (Icges(5,3) * t221 - t233 * t223) * t206;
t216 = qJD(2) * t221;
t214 = rSges(2,1) * t223 - rSges(2,2) * t221;
t213 = rSges(4,1) * t222 - rSges(4,2) * t220;
t212 = rSges(2,1) * t221 + rSges(2,2) * t223;
t208 = Icges(4,5) * t222 - Icges(4,6) * t220;
t204 = rSges(5,1) * t218 - rSges(5,2) * t217;
t200 = pkin(6) * t223 + t221 * t252;
t199 = pkin(6) * t221 - t223 * t252;
t198 = t221 * rSges(4,3) - t240 * t223;
t197 = t223 * rSges(4,3) + t240 * t221;
t192 = Icges(4,3) * t221 - t234 * t223;
t191 = Icges(4,3) * t223 + t234 * t221;
t190 = t221 * rSges(5,3) - t239 * t223;
t189 = t223 * rSges(5,3) + t239 * t221;
t182 = t205 - qJD(2) * t223 + qJD(1) * (-rSges(3,2) * t223 + rSges(3,3) * t221);
t181 = t216 + (rSges(3,2) * t221 + rSges(3,3) * t223 - t211) * qJD(1);
t180 = (-t197 * t221 + t198 * t223) * qJD(3);
t179 = qJD(1) * t197 + (-qJD(3) * t213 - qJD(2)) * t223 + t245;
t178 = t213 * t244 + t216 + (-t198 + t241) * qJD(1);
t177 = -t204 * t207 + (-qJD(2) - t242) * t223 + (t189 + t200) * qJD(1) + t245;
t176 = t221 * t242 + t206 * t204 + t216 + (-t190 - t199 + t241) * qJD(1);
t175 = -t206 * t189 + t207 * t190 + (t199 * t223 - t200 * t221) * qJD(3);
t1 = m(3) * (t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(4) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + qJD(3) * t223 * ((t223 * t208 + t227 * t221) * qJD(1) + (t223 ^ 2 * t191 + (t229 * t221 + (t192 - t230) * t223) * t221) * qJD(3)) / 0.2e1 + ((t221 * t208 - t227 * t223) * qJD(1) + (t221 ^ 2 * t192 + (t230 * t223 + (t191 - t229) * t221) * t223) * qJD(3)) * t244 / 0.2e1 + m(5) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t207 * (t254 * t221 + t226 * t223) / 0.2e1 + t206 * (t226 * t221 - t254 * t223) / 0.2e1 + (((-t193 * t220 + t195 * t222) * t223 + (-t220 * t194 + t222 * t196) * t221) * qJD(3) + (-t185 * t217 + t187 * t218) * t207 + (-t186 * t217 + t188 * t218) * t206 + (-t217 * t202 + t218 * t203 - t220 * t209 + t222 * t210) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t212 ^ 2 + t214 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
