% Calculate kinetic energy for
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:37
% EndTime: 2019-12-31 20:05:40
% DurationCPUTime: 2.74s
% Computational Cost: add. (1273->280), mult. (1871->403), div. (0->0), fcn. (1831->8), ass. (0->131)
t271 = Icges(5,1) + Icges(6,1);
t270 = -Icges(5,4) + Icges(6,5);
t269 = Icges(6,4) + Icges(5,5);
t268 = Icges(5,2) + Icges(6,3);
t267 = -Icges(6,6) + Icges(5,6);
t266 = -Icges(5,3) - Icges(6,2);
t265 = rSges(6,1) + pkin(4);
t264 = rSges(6,3) + qJ(5);
t205 = pkin(8) + qJ(4);
t199 = sin(t205);
t200 = cos(t205);
t212 = cos(qJ(1));
t210 = sin(qJ(1));
t211 = cos(qJ(2));
t240 = t210 * t211;
t150 = t199 * t240 + t200 * t212;
t151 = -t199 * t212 + t200 * t240;
t209 = sin(qJ(2));
t243 = t209 * t210;
t263 = t268 * t150 + t270 * t151 - t267 * t243;
t239 = t211 * t212;
t152 = t199 * t239 - t210 * t200;
t153 = t210 * t199 + t200 * t239;
t242 = t209 * t212;
t262 = t268 * t152 + t270 * t153 - t267 * t242;
t261 = -t267 * t150 + t269 * t151 - t266 * t243;
t260 = -t267 * t152 + t269 * t153 - t266 * t242;
t259 = t270 * t150 + t271 * t151 + t269 * t243;
t258 = t270 * t152 + t271 * t153 + t269 * t242;
t257 = t267 * t211 + (t268 * t199 + t270 * t200) * t209;
t256 = t266 * t211 + (-t267 * t199 + t269 * t200) * t209;
t255 = -t269 * t211 + (t270 * t199 + t271 * t200) * t209;
t207 = cos(pkin(8));
t248 = pkin(3) * t207;
t247 = Icges(2,4) * t210;
t246 = Icges(3,4) * t209;
t245 = Icges(3,4) * t211;
t206 = sin(pkin(8));
t244 = t206 * t212;
t241 = t210 * t206;
t237 = rSges(6,2) * t243 + t264 * t150 + t265 * t151;
t236 = rSges(6,2) * t242 + t264 * t152 + t265 * t153;
t235 = -rSges(6,2) * t211 + (t264 * t199 + t265 * t200) * t209;
t226 = pkin(2) * t211 + qJ(3) * t209;
t170 = t226 * t210;
t191 = t210 * pkin(1) - pkin(6) * t212;
t234 = -t170 - t191;
t233 = qJD(3) * t209;
t232 = qJD(4) * t209;
t231 = V_base(5) * pkin(5) + V_base(1);
t195 = qJD(2) * t210 + V_base(4);
t201 = V_base(6) + qJD(1);
t187 = pkin(2) * t209 - qJ(3) * t211;
t194 = -qJD(2) * t212 + V_base(5);
t228 = t194 * t187 + t212 * t233 + t231;
t227 = rSges(3,1) * t211 - rSges(3,2) * t209;
t225 = Icges(3,1) * t211 - t246;
t224 = -Icges(3,2) * t209 + t245;
t223 = Icges(3,5) * t211 - Icges(3,6) * t209;
t192 = pkin(1) * t212 + t210 * pkin(6);
t222 = -V_base(4) * pkin(5) + t201 * t192 + V_base(2);
t221 = V_base(4) * t191 - t192 * V_base(5) + V_base(3);
t220 = (-Icges(3,3) * t212 + t210 * t223) * t194 + (Icges(3,3) * t210 + t212 * t223) * t195 + (Icges(3,5) * t209 + Icges(3,6) * t211) * t201;
t219 = pkin(7) * t209 + t211 * t248;
t171 = t226 * t212;
t218 = t201 * t171 + t210 * t233 + t222;
t217 = -qJD(3) * t211 + t195 * t170 + t221;
t129 = -pkin(3) * t244 + t210 * t219;
t135 = -pkin(7) * t211 + t209 * t248;
t216 = t194 * t135 + (-t129 + t234) * t201 + t228;
t130 = pkin(3) * t241 + t212 * t219;
t215 = t201 * t130 + (-t135 - t187) * t195 + t218;
t214 = t195 * t129 + (-t130 - t171) * t194 + t217;
t156 = -Icges(3,6) * t212 + t210 * t224;
t157 = Icges(3,6) * t210 + t212 * t224;
t158 = -Icges(3,5) * t212 + t210 * t225;
t159 = Icges(3,5) * t210 + t212 * t225;
t180 = Icges(3,2) * t211 + t246;
t183 = Icges(3,1) * t209 + t245;
t213 = (-t157 * t209 + t159 * t211) * t195 + (-t156 * t209 + t158 * t211) * t194 + (-t180 * t209 + t183 * t211) * t201;
t203 = Icges(2,4) * t212;
t190 = rSges(2,1) * t212 - t210 * rSges(2,2);
t189 = t210 * rSges(2,1) + rSges(2,2) * t212;
t188 = rSges(3,1) * t209 + rSges(3,2) * t211;
t186 = -qJD(4) * t211 + t201;
t185 = Icges(2,1) * t212 - t247;
t184 = Icges(2,1) * t210 + t203;
t182 = -Icges(2,2) * t210 + t203;
t181 = Icges(2,2) * t212 + t247;
t176 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t175 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t174 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t169 = t212 * t232 + t195;
t168 = t210 * t232 + t194;
t167 = t207 * t239 + t241;
t166 = -t206 * t239 + t210 * t207;
t165 = t207 * t240 - t244;
t164 = -t206 * t240 - t207 * t212;
t162 = t210 * rSges(3,3) + t212 * t227;
t161 = -rSges(3,3) * t212 + t210 * t227;
t149 = -rSges(4,3) * t211 + (rSges(4,1) * t207 - rSges(4,2) * t206) * t209;
t147 = -Icges(4,5) * t211 + (Icges(4,1) * t207 - Icges(4,4) * t206) * t209;
t146 = -Icges(4,6) * t211 + (Icges(4,4) * t207 - Icges(4,2) * t206) * t209;
t145 = -Icges(4,3) * t211 + (Icges(4,5) * t207 - Icges(4,6) * t206) * t209;
t143 = -rSges(5,3) * t211 + (rSges(5,1) * t200 - rSges(5,2) * t199) * t209;
t134 = V_base(5) * rSges(2,3) - t189 * t201 + t231;
t133 = t190 * t201 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t132 = t189 * V_base(4) - t190 * V_base(5) + V_base(3);
t128 = t167 * rSges(4,1) + t166 * rSges(4,2) + rSges(4,3) * t242;
t127 = rSges(4,1) * t165 + rSges(4,2) * t164 + rSges(4,3) * t243;
t126 = Icges(4,1) * t167 + Icges(4,4) * t166 + Icges(4,5) * t242;
t125 = Icges(4,1) * t165 + Icges(4,4) * t164 + Icges(4,5) * t243;
t124 = Icges(4,4) * t167 + Icges(4,2) * t166 + Icges(4,6) * t242;
t123 = Icges(4,4) * t165 + Icges(4,2) * t164 + Icges(4,6) * t243;
t122 = Icges(4,5) * t167 + Icges(4,6) * t166 + Icges(4,3) * t242;
t121 = Icges(4,5) * t165 + Icges(4,6) * t164 + Icges(4,3) * t243;
t117 = t153 * rSges(5,1) - t152 * rSges(5,2) + rSges(5,3) * t242;
t115 = rSges(5,1) * t151 - rSges(5,2) * t150 + rSges(5,3) * t243;
t100 = t188 * t194 + (-t161 - t191) * t201 + t231;
t99 = t162 * t201 - t188 * t195 + t222;
t98 = t161 * t195 - t162 * t194 + t221;
t97 = t149 * t194 + (-t127 + t234) * t201 + t228;
t96 = t128 * t201 + (-t149 - t187) * t195 + t218;
t95 = t127 * t195 + (-t128 - t171) * t194 + t217;
t94 = -t115 * t186 + t143 * t168 + t216;
t93 = t117 * t186 - t143 * t169 + t215;
t92 = t115 * t169 - t117 * t168 + t214;
t91 = qJD(5) * t152 + t168 * t235 - t186 * t237 + t216;
t90 = qJD(5) * t150 - t169 * t235 + t186 * t236 + t215;
t89 = qJD(5) * t199 * t209 - t168 * t236 + t169 * t237 + t214;
t1 = m(1) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(2) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + ((t150 * t257 + t151 * t255 + t243 * t256) * t186 + (t150 * t262 + t151 * t258 + t243 * t260) * t169 + (t263 * t150 + t259 * t151 + t261 * t243) * t168) * t168 / 0.2e1 + ((t152 * t257 + t153 * t255 + t242 * t256) * t186 + (t262 * t152 + t258 * t153 + t260 * t242) * t169 + (t152 * t263 + t259 * t153 + t261 * t242) * t168) * t169 / 0.2e1 + ((-t168 * t261 - t169 * t260 - t186 * t256) * t211 + ((t199 * t257 + t200 * t255) * t186 + (t199 * t262 + t200 * t258) * t169 + (t199 * t263 + t259 * t200) * t168) * t209) * t186 / 0.2e1 + (t210 * t213 - t212 * t220 + (t122 * t243 + t124 * t164 + t126 * t165) * t195 + (t121 * t243 + t123 * t164 + t125 * t165) * t194 + (t145 * t243 + t146 * t164 + t147 * t165) * t201) * t194 / 0.2e1 + (t210 * t220 + t212 * t213 + (t122 * t242 + t166 * t124 + t167 * t126) * t195 + (t121 * t242 + t166 * t123 + t167 * t125) * t194 + (t145 * t242 + t166 * t146 + t167 * t147) * t201) * t195 / 0.2e1 + ((-t210 * t181 + t184 * t212 + Icges(1,4)) * V_base(5) + (-t210 * t182 + t185 * t212 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t181 * t212 + t210 * t184 + Icges(1,2)) * V_base(5) + (t182 * t212 + t210 * t185 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t157 * t211 + t159 * t209) * t195 + (t156 * t211 + t158 * t209) * t194 + (-t121 * t194 - t122 * t195) * t211 + ((-t124 * t206 + t126 * t207) * t195 + (-t123 * t206 + t125 * t207) * t194) * t209 + (Icges(2,3) + (t180 - t145) * t211 + (-t146 * t206 + t147 * t207 + t183) * t209) * t201) * t201 / 0.2e1 + t201 * V_base(4) * (Icges(2,5) * t212 - Icges(2,6) * t210) + t201 * V_base(5) * (Icges(2,5) * t210 + Icges(2,6) * t212) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
