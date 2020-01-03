% Calculate kinetic energy for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:29
% EndTime: 2019-12-31 16:32:30
% DurationCPUTime: 0.44s
% Computational Cost: add. (516->107), mult. (463->186), div. (0->0), fcn. (386->6), ass. (0->68)
t212 = cos(qJ(3));
t238 = pkin(3) * t212;
t211 = sin(qJ(3));
t236 = Icges(4,4) * t211;
t235 = Icges(4,4) * t212;
t210 = qJ(3) + qJ(4);
t207 = sin(t210);
t234 = Icges(5,4) * t207;
t208 = cos(t210);
t233 = Icges(5,4) * t208;
t209 = pkin(7) + qJ(2);
t205 = sin(t209);
t232 = qJD(3) * t205;
t206 = cos(t209);
t231 = qJD(3) * t206;
t230 = qJD(3) + qJD(4);
t229 = pkin(3) * qJD(3) * t211;
t228 = rSges(4,1) * t212 - rSges(4,2) * t211;
t227 = rSges(5,1) * t208 - rSges(5,2) * t207;
t226 = Icges(4,1) * t212 - t236;
t225 = Icges(5,1) * t208 - t234;
t224 = -Icges(4,2) * t211 + t235;
t223 = -Icges(5,2) * t207 + t233;
t222 = Icges(4,5) * t212 - Icges(4,6) * t211;
t221 = Icges(5,5) * t208 - Icges(5,6) * t207;
t184 = -Icges(4,6) * t206 + t205 * t224;
t186 = -Icges(4,5) * t206 + t205 * t226;
t220 = t184 * t211 - t186 * t212;
t185 = Icges(4,6) * t205 + t206 * t224;
t187 = Icges(4,5) * t205 + t206 * t226;
t219 = -t185 * t211 + t187 * t212;
t201 = Icges(4,2) * t212 + t236;
t202 = Icges(4,1) * t211 + t235;
t218 = -t201 * t211 + t202 * t212;
t191 = t230 * t205;
t192 = t230 * t206;
t217 = (Icges(5,5) * t207 + Icges(5,6) * t208) * qJD(2) - (-Icges(5,3) * t206 + t205 * t221) * t192 + (Icges(5,3) * t205 + t206 * t221) * t191;
t176 = -Icges(5,6) * t206 + t205 * t223;
t177 = Icges(5,6) * t205 + t206 * t223;
t178 = -Icges(5,5) * t206 + t205 * t225;
t179 = Icges(5,5) * t205 + t206 * t225;
t197 = Icges(5,2) * t208 + t234;
t198 = Icges(5,1) * t207 + t233;
t216 = (-t177 * t207 + t179 * t208) * t191 - (-t176 * t207 + t178 * t208) * t192 + (-t197 * t207 + t198 * t208) * qJD(2);
t215 = qJD(1) ^ 2;
t214 = qJD(2) ^ 2;
t203 = rSges(4,1) * t211 + rSges(4,2) * t212;
t200 = Icges(4,5) * t211 + Icges(4,6) * t212;
t199 = rSges(5,1) * t207 + rSges(5,2) * t208;
t195 = pkin(2) * t205 - pkin(5) * t206;
t194 = rSges(3,1) * t206 - rSges(3,2) * t205;
t193 = rSges(3,1) * t205 + rSges(3,2) * t206;
t190 = qJD(2) * (pkin(2) * t206 + pkin(5) * t205);
t189 = rSges(4,3) * t205 + t206 * t228;
t188 = -rSges(4,3) * t206 + t205 * t228;
t183 = Icges(4,3) * t205 + t206 * t222;
t182 = -Icges(4,3) * t206 + t205 * t222;
t181 = rSges(5,3) * t205 + t206 * t227;
t180 = -rSges(5,3) * t206 + t205 * t227;
t173 = pkin(6) * t205 + t206 * t238;
t172 = -pkin(6) * t206 + t205 * t238;
t171 = qJD(2) * t189 - t203 * t232 + t190;
t170 = -t203 * t231 + (-t188 - t195) * qJD(2);
t169 = qJD(1) + (t188 * t205 + t189 * t206) * qJD(3);
t168 = -t205 * t229 - t191 * t199 + t190 + (t173 + t181) * qJD(2);
t167 = -t206 * t229 - t192 * t199 + (-t172 - t180 - t195) * qJD(2);
t166 = t180 * t191 + t181 * t192 + qJD(1) + (t172 * t205 + t173 * t206) * qJD(3);
t1 = m(2) * t215 / 0.2e1 + m(3) * (t215 + (t193 ^ 2 + t194 ^ 2) * t214) / 0.2e1 + t214 * Icges(3,3) / 0.2e1 + m(4) * (t169 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + ((t205 * t200 + t206 * t218) * qJD(2) + (t205 ^ 2 * t183 + (t220 * t206 + (-t182 + t219) * t205) * t206) * qJD(3)) * t232 / 0.2e1 - ((-t206 * t200 + t205 * t218) * qJD(2) + (t206 ^ 2 * t182 + (t219 * t205 + (-t183 + t220) * t206) * t205) * qJD(3)) * t231 / 0.2e1 + m(5) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + t191 * (t217 * t205 + t216 * t206) / 0.2e1 - t192 * (t216 * t205 - t217 * t206) / 0.2e1 + (((t185 * t212 + t187 * t211) * t205 - (t184 * t212 + t186 * t211) * t206) * qJD(3) + (t177 * t208 + t179 * t207) * t191 - (t176 * t208 + t178 * t207) * t192 + (t208 * t197 + t207 * t198 + t212 * t201 + t211 * t202) * qJD(2)) * qJD(2) / 0.2e1;
T = t1;
