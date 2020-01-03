% Calculate kinetic energy for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:07
% EndTime: 2019-12-31 17:23:07
% DurationCPUTime: 0.49s
% Computational Cost: add. (544->114), mult. (484->197), div. (0->0), fcn. (396->8), ass. (0->74)
t214 = cos(qJ(3));
t243 = pkin(3) * t214;
t241 = pkin(1) * qJD(1);
t212 = sin(qJ(3));
t240 = Icges(4,4) * t212;
t239 = Icges(4,4) * t214;
t210 = qJ(3) + qJ(4);
t205 = sin(t210);
t238 = Icges(5,4) * t205;
t207 = cos(t210);
t237 = Icges(5,4) * t207;
t215 = cos(qJ(1));
t204 = t215 * t241;
t211 = qJ(1) + qJ(2);
t206 = sin(t211);
t208 = cos(t211);
t209 = qJD(1) + qJD(2);
t236 = t209 * (pkin(2) * t208 + pkin(6) * t206) + t204;
t235 = qJD(3) * t206;
t234 = qJD(3) * t208;
t233 = qJD(3) + qJD(4);
t232 = pkin(3) * qJD(3) * t212;
t213 = sin(qJ(1));
t231 = t213 * t241;
t230 = rSges(4,1) * t214 - rSges(4,2) * t212;
t229 = rSges(5,1) * t207 - rSges(5,2) * t205;
t228 = Icges(4,1) * t214 - t240;
t227 = Icges(5,1) * t207 - t238;
t226 = -Icges(4,2) * t212 + t239;
t225 = -Icges(5,2) * t205 + t237;
t224 = Icges(4,5) * t214 - Icges(4,6) * t212;
t223 = Icges(5,5) * t207 - Icges(5,6) * t205;
t181 = -Icges(4,6) * t208 + t226 * t206;
t183 = -Icges(4,5) * t208 + t228 * t206;
t222 = t181 * t212 - t183 * t214;
t182 = Icges(4,6) * t206 + t226 * t208;
t184 = Icges(4,5) * t206 + t228 * t208;
t221 = -t182 * t212 + t184 * t214;
t198 = Icges(4,2) * t214 + t240;
t199 = Icges(4,1) * t212 + t239;
t220 = -t198 * t212 + t199 * t214;
t190 = t233 * t206;
t191 = t233 * t208;
t219 = -(-Icges(5,3) * t208 + t223 * t206) * t191 + (Icges(5,3) * t206 + t223 * t208) * t190 + (Icges(5,5) * t205 + Icges(5,6) * t207) * t209;
t173 = -Icges(5,6) * t208 + t225 * t206;
t174 = Icges(5,6) * t206 + t225 * t208;
t175 = -Icges(5,5) * t208 + t227 * t206;
t176 = Icges(5,5) * t206 + t227 * t208;
t193 = Icges(5,2) * t207 + t238;
t194 = Icges(5,1) * t205 + t237;
t218 = (-t174 * t205 + t176 * t207) * t190 - (-t173 * t205 + t175 * t207) * t191 + (-t193 * t205 + t194 * t207) * t209;
t202 = rSges(2,1) * t215 - rSges(2,2) * t213;
t201 = rSges(2,1) * t213 + rSges(2,2) * t215;
t200 = rSges(4,1) * t212 + rSges(4,2) * t214;
t197 = Icges(4,5) * t212 + Icges(4,6) * t214;
t196 = pkin(2) * t206 - pkin(6) * t208;
t195 = rSges(5,1) * t205 + rSges(5,2) * t207;
t188 = t204 + t209 * (rSges(3,1) * t208 - rSges(3,2) * t206);
t187 = -t231 - t209 * (rSges(3,1) * t206 + rSges(3,2) * t208);
t186 = rSges(4,3) * t206 + t230 * t208;
t185 = -rSges(4,3) * t208 + t230 * t206;
t180 = Icges(4,3) * t206 + t224 * t208;
t179 = -Icges(4,3) * t208 + t224 * t206;
t178 = rSges(5,3) * t206 + t229 * t208;
t177 = -rSges(5,3) * t208 + t229 * t206;
t170 = pkin(7) * t206 + t243 * t208;
t169 = -pkin(7) * t208 + t243 * t206;
t168 = (t185 * t206 + t186 * t208) * qJD(3);
t167 = t186 * t209 - t200 * t235 + t236;
t166 = -t231 - t200 * t234 + (-t185 - t196) * t209;
t165 = -t206 * t232 - t190 * t195 + (t170 + t178) * t209 + t236;
t164 = -t208 * t232 - t231 - t191 * t195 + (-t169 - t177 - t196) * t209;
t163 = t177 * t190 + t178 * t191 + (t169 * t206 + t170 * t208) * qJD(3);
t1 = m(3) * (t187 ^ 2 + t188 ^ 2) / 0.2e1 + t209 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + ((t206 * t197 + t220 * t208) * t209 + (t206 ^ 2 * t180 + (t222 * t208 + (-t179 + t221) * t206) * t208) * qJD(3)) * t235 / 0.2e1 - ((-t208 * t197 + t220 * t206) * t209 + (t208 ^ 2 * t179 + (t221 * t206 + (-t180 + t222) * t208) * t206) * qJD(3)) * t234 / 0.2e1 + m(5) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + t190 * (t219 * t206 + t218 * t208) / 0.2e1 - t191 * (t218 * t206 - t219 * t208) / 0.2e1 + (((t182 * t214 + t184 * t212) * t206 - (t181 * t214 + t183 * t212) * t208) * qJD(3) + (t174 * t207 + t176 * t205) * t190 - (t173 * t207 + t175 * t205) * t191 + (t207 * t193 + t205 * t194 + t214 * t198 + t212 * t199) * t209) * t209 / 0.2e1 + (m(2) * (t201 ^ 2 + t202 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
