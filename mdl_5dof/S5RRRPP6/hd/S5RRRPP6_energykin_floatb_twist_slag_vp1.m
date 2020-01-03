% Calculate kinetic energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:55
% EndTime: 2019-12-31 20:59:57
% DurationCPUTime: 2.33s
% Computational Cost: add. (1300->271), mult. (1916->381), div. (0->0), fcn. (1876->8), ass. (0->127)
t268 = Icges(5,1) + Icges(6,1);
t267 = -Icges(5,4) + Icges(6,5);
t266 = Icges(6,4) + Icges(5,5);
t265 = Icges(5,2) + Icges(6,3);
t264 = -Icges(6,6) + Icges(5,6);
t263 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t262 = rSges(6,1) + pkin(4);
t261 = rSges(6,3) + qJ(5);
t205 = qJ(3) + pkin(8);
t199 = sin(t205);
t200 = cos(t205);
t212 = cos(qJ(1));
t209 = sin(qJ(1));
t211 = cos(qJ(2));
t238 = t209 * t211;
t145 = t199 * t238 + t200 * t212;
t146 = -t199 * t212 + t200 * t238;
t208 = sin(qJ(2));
t241 = t208 * t209;
t260 = t265 * t145 + t267 * t146 - t264 * t241;
t237 = t211 * t212;
t147 = t199 * t237 - t209 * t200;
t148 = t209 * t199 + t200 * t237;
t240 = t208 * t212;
t259 = t265 * t147 + t267 * t148 - t264 * t240;
t258 = t267 * t145 + t268 * t146 + t266 * t241;
t257 = t267 * t147 + t268 * t148 + t266 * t240;
t256 = t264 * t211 + (t265 * t199 + t267 * t200) * t208;
t255 = -t266 * t211 + (t267 * t199 + t268 * t200) * t208;
t207 = sin(qJ(3));
t210 = cos(qJ(3));
t166 = -t207 * t238 - t210 * t212;
t242 = t207 * t212;
t167 = t210 * t238 - t242;
t254 = Icges(4,5) * t167 + Icges(4,6) * t166 - t264 * t145 + t266 * t146 - t263 * t241;
t168 = -t207 * t237 + t209 * t210;
t239 = t209 * t207;
t169 = t210 * t237 + t239;
t253 = Icges(4,5) * t169 + Icges(4,6) * t168 - t264 * t147 + t266 * t148 - t263 * t240;
t252 = t263 * t211 + (Icges(4,5) * t210 - Icges(4,6) * t207 - t264 * t199 + t266 * t200) * t208;
t247 = pkin(3) * t210;
t245 = Icges(2,4) * t209;
t244 = Icges(3,4) * t208;
t243 = Icges(3,4) * t211;
t236 = rSges(6,2) * t241 + t261 * t145 + t262 * t146;
t235 = rSges(6,2) * t240 + t261 * t147 + t262 * t148;
t234 = -rSges(6,2) * t211 + (t261 * t199 + t262 * t200) * t208;
t233 = qJD(3) * t208;
t232 = qJD(4) * t208;
t231 = V_base(5) * pkin(5) + V_base(1);
t195 = qJD(2) * t209 + V_base(4);
t201 = V_base(6) + qJD(1);
t228 = pkin(2) * t211 + pkin(7) * t208;
t194 = -qJD(2) * t212 + V_base(5);
t227 = rSges(3,1) * t211 - rSges(3,2) * t208;
t226 = Icges(3,1) * t211 - t244;
t225 = -Icges(3,2) * t208 + t243;
t224 = Icges(3,5) * t211 - Icges(3,6) * t208;
t192 = pkin(1) * t212 + t209 * pkin(6);
t223 = -V_base(4) * pkin(5) + t201 * t192 + V_base(2);
t191 = t209 * pkin(1) - pkin(6) * t212;
t222 = V_base(4) * t191 - t192 * V_base(5) + V_base(3);
t221 = qJ(4) * t208 + t211 * t247;
t171 = t228 * t209;
t190 = pkin(2) * t208 - pkin(7) * t211;
t220 = t194 * t190 + (-t171 - t191) * t201 + t231;
t219 = (-Icges(3,3) * t212 + t209 * t224) * t194 + (Icges(3,3) * t209 + t212 * t224) * t195 + (Icges(3,5) * t208 + Icges(3,6) * t211) * t201;
t172 = t228 * t212;
t218 = t201 * t172 - t190 * t195 + t223;
t135 = -qJ(4) * t211 + t208 * t247;
t164 = t209 * t233 + t194;
t217 = t164 * t135 + t212 * t232 + t220;
t216 = t195 * t171 - t172 * t194 + t222;
t128 = pkin(3) * t239 + t212 * t221;
t186 = -qJD(3) * t211 + t201;
t215 = t186 * t128 + t209 * t232 + t218;
t127 = -pkin(3) * t242 + t209 * t221;
t165 = t212 * t233 + t195;
t214 = -qJD(4) * t211 + t165 * t127 + t216;
t154 = -Icges(3,6) * t212 + t209 * t225;
t155 = Icges(3,6) * t209 + t212 * t225;
t157 = -Icges(3,5) * t212 + t209 * t226;
t158 = Icges(3,5) * t209 + t212 * t226;
t180 = Icges(3,2) * t211 + t244;
t183 = Icges(3,1) * t208 + t243;
t213 = (-t155 * t208 + t158 * t211) * t195 + (-t154 * t208 + t157 * t211) * t194 + (-t180 * t208 + t183 * t211) * t201;
t203 = Icges(2,4) * t212;
t189 = rSges(2,1) * t212 - t209 * rSges(2,2);
t188 = t209 * rSges(2,1) + rSges(2,2) * t212;
t187 = rSges(3,1) * t208 + rSges(3,2) * t211;
t185 = Icges(2,1) * t212 - t245;
t184 = Icges(2,1) * t209 + t203;
t182 = -Icges(2,2) * t209 + t203;
t181 = Icges(2,2) * t212 + t245;
t176 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t175 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t174 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t162 = t209 * rSges(3,3) + t212 * t227;
t161 = -rSges(3,3) * t212 + t209 * t227;
t160 = -rSges(4,3) * t211 + (rSges(4,1) * t210 - rSges(4,2) * t207) * t208;
t156 = -Icges(4,5) * t211 + (Icges(4,1) * t210 - Icges(4,4) * t207) * t208;
t153 = -Icges(4,6) * t211 + (Icges(4,4) * t210 - Icges(4,2) * t207) * t208;
t143 = -rSges(5,3) * t211 + (rSges(5,1) * t200 - rSges(5,2) * t199) * t208;
t134 = V_base(5) * rSges(2,3) - t188 * t201 + t231;
t133 = t189 * t201 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t132 = t188 * V_base(4) - t189 * V_base(5) + V_base(3);
t130 = t169 * rSges(4,1) + t168 * rSges(4,2) + rSges(4,3) * t240;
t129 = rSges(4,1) * t167 + rSges(4,2) * t166 + rSges(4,3) * t241;
t126 = Icges(4,1) * t169 + Icges(4,4) * t168 + Icges(4,5) * t240;
t125 = Icges(4,1) * t167 + Icges(4,4) * t166 + Icges(4,5) * t241;
t124 = Icges(4,4) * t169 + Icges(4,2) * t168 + Icges(4,6) * t240;
t123 = Icges(4,4) * t167 + Icges(4,2) * t166 + Icges(4,6) * t241;
t118 = t148 * rSges(5,1) - t147 * rSges(5,2) + rSges(5,3) * t240;
t116 = rSges(5,1) * t146 - rSges(5,2) * t145 + rSges(5,3) * t241;
t100 = t187 * t194 + (-t161 - t191) * t201 + t231;
t99 = t162 * t201 - t187 * t195 + t223;
t98 = t161 * t195 - t162 * t194 + t222;
t97 = -t129 * t186 + t160 * t164 + t220;
t96 = t130 * t186 - t160 * t165 + t218;
t95 = t129 * t165 - t130 * t164 + t216;
t94 = t143 * t164 + (-t116 - t127) * t186 + t217;
t93 = t118 * t186 + (-t135 - t143) * t165 + t215;
t92 = t116 * t165 + (-t118 - t128) * t164 + t214;
t91 = qJD(5) * t147 + t234 * t164 + (-t127 - t236) * t186 + t217;
t90 = qJD(5) * t145 + t235 * t186 + (-t135 - t234) * t165 + t215;
t89 = qJD(5) * t199 * t208 + t236 * t165 + (-t128 - t235) * t164 + t214;
t1 = m(1) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(2) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t195 * (t219 * t209 + t213 * t212) / 0.2e1 + t194 * (t213 * t209 - t219 * t212) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + ((t155 * t211 + t158 * t208) * t195 + (t154 * t211 + t157 * t208) * t194 + (t180 * t211 + t183 * t208 + Icges(2,3)) * t201) * t201 / 0.2e1 + ((-t209 * t181 + t184 * t212 + Icges(1,4)) * V_base(5) + (-t209 * t182 + t185 * t212 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t181 * t212 + t209 * t184 + Icges(1,2)) * V_base(5) + (t182 * t212 + t209 * t185 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t145 * t256 + t146 * t255 + t153 * t166 + t156 * t167 + t241 * t252) * t186 + (t124 * t166 + t126 * t167 + t145 * t259 + t146 * t257 + t241 * t253) * t165 + (t166 * t123 + t167 * t125 + t260 * t145 + t258 * t146 + t254 * t241) * t164) * t164 / 0.2e1 + ((t147 * t256 + t148 * t255 + t168 * t153 + t169 * t156 + t240 * t252) * t186 + (t168 * t124 + t169 * t126 + t259 * t147 + t257 * t148 + t253 * t240) * t165 + (t168 * t123 + t169 * t125 + t147 * t260 + t258 * t148 + t254 * t240) * t164) * t165 / 0.2e1 + ((-t164 * t254 - t165 * t253 - t186 * t252) * t211 + ((-t153 * t207 + t156 * t210 + t199 * t256 + t200 * t255) * t186 + (-t124 * t207 + t126 * t210 + t199 * t259 + t200 * t257) * t165 + (-t123 * t207 + t125 * t210 + t199 * t260 + t258 * t200) * t164) * t208) * t186 / 0.2e1 + V_base(4) * t201 * (Icges(2,5) * t212 - Icges(2,6) * t209) + V_base(5) * t201 * (Icges(2,5) * t209 + Icges(2,6) * t212) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
