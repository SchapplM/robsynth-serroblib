% Calculate kinetic energy for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:04
% EndTime: 2019-12-31 16:49:05
% DurationCPUTime: 0.50s
% Computational Cost: add. (526->115), mult. (485->194), div. (0->0), fcn. (396->8), ass. (0->72)
t214 = sin(qJ(1));
t245 = pkin(1) * t214;
t215 = cos(qJ(3));
t243 = pkin(3) * t215;
t213 = sin(qJ(3));
t241 = Icges(4,4) * t213;
t240 = Icges(4,4) * t215;
t212 = qJ(3) + qJ(4);
t209 = sin(t212);
t239 = Icges(5,4) * t209;
t210 = cos(t212);
t238 = Icges(5,4) * t210;
t216 = cos(qJ(1));
t206 = qJD(1) * t216 * pkin(1);
t211 = qJ(1) + pkin(7);
t207 = sin(t211);
t208 = cos(t211);
t237 = qJD(1) * (pkin(2) * t208 + pkin(5) * t207) + t206;
t236 = qJD(3) * t207;
t235 = qJD(3) * t208;
t234 = qJD(3) + qJD(4);
t233 = pkin(3) * qJD(3) * t213;
t232 = -pkin(2) * t207 + pkin(5) * t208 - t245;
t231 = rSges(4,1) * t215 - rSges(4,2) * t213;
t230 = rSges(5,1) * t210 - rSges(5,2) * t209;
t229 = Icges(4,1) * t215 - t241;
t228 = Icges(5,1) * t210 - t239;
t227 = -Icges(4,2) * t213 + t240;
t226 = -Icges(5,2) * t209 + t238;
t225 = Icges(4,5) * t215 - Icges(4,6) * t213;
t224 = Icges(5,5) * t210 - Icges(5,6) * t209;
t183 = -Icges(4,6) * t208 + t227 * t207;
t185 = -Icges(4,5) * t208 + t229 * t207;
t223 = t183 * t213 - t185 * t215;
t184 = Icges(4,6) * t207 + t227 * t208;
t186 = Icges(4,5) * t207 + t229 * t208;
t222 = -t184 * t213 + t186 * t215;
t200 = Icges(4,2) * t215 + t241;
t201 = Icges(4,1) * t213 + t240;
t221 = -t200 * t213 + t201 * t215;
t192 = t234 * t207;
t193 = t234 * t208;
t220 = qJD(1) * (Icges(5,5) * t209 + Icges(5,6) * t210) - (-Icges(5,3) * t208 + t224 * t207) * t193 + (Icges(5,3) * t207 + t224 * t208) * t192;
t175 = -Icges(5,6) * t208 + t226 * t207;
t176 = Icges(5,6) * t207 + t226 * t208;
t177 = -Icges(5,5) * t208 + t228 * t207;
t178 = Icges(5,5) * t207 + t228 * t208;
t196 = Icges(5,2) * t210 + t239;
t197 = Icges(5,1) * t209 + t238;
t219 = (-t176 * t209 + t178 * t210) * t192 - (-t175 * t209 + t177 * t210) * t193 + (-t196 * t209 + t197 * t210) * qJD(1);
t204 = rSges(2,1) * t216 - rSges(2,2) * t214;
t203 = rSges(2,1) * t214 + rSges(2,2) * t216;
t202 = rSges(4,1) * t213 + rSges(4,2) * t215;
t199 = Icges(4,5) * t213 + Icges(4,6) * t215;
t198 = rSges(5,1) * t209 + rSges(5,2) * t210;
t190 = t206 + qJD(1) * (rSges(3,1) * t208 - rSges(3,2) * t207);
t189 = (-rSges(3,1) * t207 - rSges(3,2) * t208 - t245) * qJD(1);
t188 = rSges(4,3) * t207 + t231 * t208;
t187 = -rSges(4,3) * t208 + t231 * t207;
t182 = Icges(4,3) * t207 + t225 * t208;
t181 = -Icges(4,3) * t208 + t225 * t207;
t180 = rSges(5,3) * t207 + t230 * t208;
t179 = -rSges(5,3) * t208 + t230 * t207;
t172 = pkin(6) * t207 + t243 * t208;
t171 = -pkin(6) * t208 + t243 * t207;
t170 = qJD(1) * t188 - t202 * t236 + t237;
t169 = -t202 * t235 + (-t187 + t232) * qJD(1);
t168 = qJD(2) + (t187 * t207 + t188 * t208) * qJD(3);
t167 = -t207 * t233 - t192 * t198 + (t172 + t180) * qJD(1) + t237;
t166 = -t208 * t233 - t193 * t198 + (-t171 - t179 + t232) * qJD(1);
t165 = t179 * t192 + t180 * t193 + qJD(2) + (t171 * t207 + t172 * t208) * qJD(3);
t1 = m(3) * (qJD(2) ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(4) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + ((t207 * t199 + t221 * t208) * qJD(1) + (t207 ^ 2 * t182 + (t223 * t208 + (-t181 + t222) * t207) * t208) * qJD(3)) * t236 / 0.2e1 - ((-t208 * t199 + t221 * t207) * qJD(1) + (t208 ^ 2 * t181 + (t222 * t207 + (-t182 + t223) * t208) * t207) * qJD(3)) * t235 / 0.2e1 + m(5) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t192 * (t220 * t207 + t219 * t208) / 0.2e1 - t193 * (t219 * t207 - t220 * t208) / 0.2e1 + (((t184 * t215 + t186 * t213) * t207 - (t183 * t215 + t185 * t213) * t208) * qJD(3) + (t176 * t210 + t178 * t209) * t192 - (t175 * t210 + t177 * t209) * t193 + (t210 * t196 + t209 * t197 + t215 * t200 + t213 * t201) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t203 ^ 2 + t204 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
