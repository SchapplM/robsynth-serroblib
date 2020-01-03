% Calculate kinetic energy for
% S5RRRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:03
% EndTime: 2019-12-31 21:56:06
% DurationCPUTime: 2.18s
% Computational Cost: add. (1329->254), mult. (1648->372), div. (0->0), fcn. (1552->8), ass. (0->132)
t268 = Icges(5,1) + Icges(6,1);
t267 = -Icges(5,4) + Icges(6,5);
t266 = Icges(6,4) + Icges(5,5);
t265 = Icges(5,2) + Icges(6,3);
t264 = -Icges(6,6) + Icges(5,6);
t263 = -Icges(5,3) - Icges(6,2);
t262 = rSges(6,1) + pkin(4);
t261 = rSges(6,3) + qJ(5);
t197 = qJ(2) + qJ(3);
t194 = cos(t197);
t201 = cos(qJ(4));
t203 = cos(qJ(1));
t234 = t201 * t203;
t198 = sin(qJ(4));
t200 = sin(qJ(1));
t237 = t198 * t200;
t154 = t194 * t237 + t234;
t235 = t200 * t201;
t236 = t198 * t203;
t155 = t194 * t235 - t236;
t193 = sin(t197);
t239 = t193 * t200;
t260 = -t264 * t154 + t266 * t155 - t263 * t239;
t259 = t265 * t154 + t267 * t155 - t264 * t239;
t156 = t194 * t236 - t235;
t157 = t194 * t234 + t237;
t238 = t193 * t203;
t258 = t265 * t156 + t267 * t157 - t264 * t238;
t257 = -t264 * t156 + t266 * t157 - t263 * t238;
t256 = t267 * t154 + t268 * t155 + t266 * t239;
t255 = t267 * t156 + t268 * t157 + t266 * t238;
t254 = t264 * t194 + (t265 * t198 + t267 * t201) * t193;
t253 = t263 * t194 + (-t264 * t198 + t266 * t201) * t193;
t252 = -t266 * t194 + (t267 * t198 + t268 * t201) * t193;
t199 = sin(qJ(2));
t247 = pkin(2) * t199;
t202 = cos(qJ(2));
t246 = pkin(2) * t202;
t244 = Icges(2,4) * t200;
t243 = Icges(3,4) * t199;
t242 = Icges(3,4) * t202;
t241 = Icges(4,4) * t193;
t240 = Icges(4,4) * t194;
t233 = rSges(6,2) * t239 + t261 * t154 + t155 * t262;
t232 = rSges(6,2) * t238 + t261 * t156 + t157 * t262;
t231 = -rSges(6,2) * t194 + (t261 * t198 + t201 * t262) * t193;
t130 = -pkin(7) * t203 + t200 * t246;
t184 = t200 * pkin(1) - t203 * pkin(6);
t230 = -t130 - t184;
t229 = qJD(4) * t193;
t228 = V_base(5) * pkin(5) + V_base(1);
t187 = qJD(2) * t200 + V_base(4);
t190 = V_base(6) + qJD(1);
t186 = -qJD(2) * t203 + V_base(5);
t225 = t186 * t247 + t228;
t165 = qJD(3) * t200 + t187;
t224 = pkin(3) * t194 + pkin(8) * t193;
t223 = rSges(3,1) * t202 - rSges(3,2) * t199;
t222 = rSges(4,1) * t194 - rSges(4,2) * t193;
t221 = Icges(3,1) * t202 - t243;
t220 = Icges(4,1) * t194 - t241;
t219 = -Icges(3,2) * t199 + t242;
t218 = -Icges(4,2) * t193 + t240;
t217 = Icges(3,5) * t202 - Icges(3,6) * t199;
t216 = Icges(4,5) * t194 - Icges(4,6) * t193;
t185 = t203 * pkin(1) + t200 * pkin(6);
t215 = -V_base(4) * pkin(5) + t190 * t185 + V_base(2);
t214 = V_base(4) * t184 - t185 * V_base(5) + V_base(3);
t164 = V_base(5) + (-qJD(2) - qJD(3)) * t203;
t213 = (-Icges(4,3) * t203 + t200 * t216) * t164 + (Icges(4,3) * t200 + t203 * t216) * t165 + (Icges(4,5) * t193 + Icges(4,6) * t194) * t190;
t212 = (-Icges(3,3) * t203 + t200 * t217) * t186 + (Icges(3,3) * t200 + t203 * t217) * t187 + (Icges(3,5) * t199 + Icges(3,6) * t202) * t190;
t152 = t224 * t200;
t163 = pkin(3) * t193 - pkin(8) * t194;
t211 = t164 * t163 + (-t152 + t230) * t190 + t225;
t131 = pkin(7) * t200 + t203 * t246;
t210 = t187 * t130 - t131 * t186 + t214;
t209 = t190 * t131 - t187 * t247 + t215;
t153 = t224 * t203;
t208 = t165 * t152 - t153 * t164 + t210;
t207 = t190 * t153 - t163 * t165 + t209;
t135 = -Icges(4,6) * t203 + t200 * t218;
t136 = Icges(4,6) * t200 + t203 * t218;
t137 = -Icges(4,5) * t203 + t200 * t220;
t138 = Icges(4,5) * t200 + t203 * t220;
t160 = Icges(4,2) * t194 + t241;
t161 = Icges(4,1) * t193 + t240;
t206 = (-t136 * t193 + t138 * t194) * t165 + (-t135 * t193 + t137 * t194) * t164 + (-t160 * t193 + t161 * t194) * t190;
t145 = -Icges(3,6) * t203 + t200 * t219;
t146 = Icges(3,6) * t200 + t203 * t219;
t147 = -Icges(3,5) * t203 + t200 * t221;
t148 = Icges(3,5) * t200 + t203 * t221;
t175 = Icges(3,2) * t202 + t243;
t178 = Icges(3,1) * t199 + t242;
t205 = (-t146 * t199 + t148 * t202) * t187 + (-t145 * t199 + t147 * t202) * t186 + (-t175 * t199 + t178 * t202) * t190;
t195 = Icges(2,4) * t203;
t183 = rSges(2,1) * t203 - rSges(2,2) * t200;
t182 = rSges(2,1) * t200 + rSges(2,2) * t203;
t181 = rSges(3,1) * t199 + rSges(3,2) * t202;
t180 = Icges(2,1) * t203 - t244;
t179 = Icges(2,1) * t200 + t195;
t177 = -Icges(2,2) * t200 + t195;
t176 = Icges(2,2) * t203 + t244;
t171 = -qJD(4) * t194 + t190;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t162 = rSges(4,1) * t193 + rSges(4,2) * t194;
t150 = rSges(3,3) * t200 + t203 * t223;
t149 = -rSges(3,3) * t203 + t200 * t223;
t142 = t203 * t229 + t165;
t141 = t200 * t229 + t164;
t140 = rSges(4,3) * t200 + t203 * t222;
t139 = -rSges(4,3) * t203 + t200 * t222;
t129 = -rSges(5,3) * t194 + (rSges(5,1) * t201 - rSges(5,2) * t198) * t193;
t121 = V_base(5) * rSges(2,3) - t182 * t190 + t228;
t120 = t183 * t190 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t118 = t182 * V_base(4) - t183 * V_base(5) + V_base(3);
t112 = rSges(5,1) * t157 - rSges(5,2) * t156 + rSges(5,3) * t238;
t110 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t239;
t96 = t181 * t186 + (-t149 - t184) * t190 + t228;
t95 = t150 * t190 - t181 * t187 + t215;
t94 = t149 * t187 - t150 * t186 + t214;
t93 = t162 * t164 + (-t139 + t230) * t190 + t225;
t92 = t140 * t190 - t162 * t165 + t209;
t91 = t139 * t165 - t140 * t164 + t210;
t90 = -t110 * t171 + t129 * t141 + t211;
t89 = t112 * t171 - t129 * t142 + t207;
t88 = t110 * t142 - t112 * t141 + t208;
t87 = qJD(5) * t156 + t141 * t231 - t171 * t233 + t211;
t86 = qJD(5) * t154 - t142 * t231 + t171 * t232 + t207;
t85 = qJD(5) * t193 * t198 - t141 * t232 + t142 * t233 + t208;
t1 = m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t118 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t187 * (t212 * t200 + t205 * t203) / 0.2e1 + t186 * (t205 * t200 - t212 * t203) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t165 * (t213 * t200 + t206 * t203) / 0.2e1 + t164 * (t206 * t200 - t213 * t203) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t154 * t254 + t155 * t252 + t239 * t253) * t171 + (t154 * t258 + t155 * t255 + t239 * t257) * t142 + (t259 * t154 + t256 * t155 + t260 * t239) * t141) * t141 / 0.2e1 + ((t156 * t254 + t157 * t252 + t238 * t253) * t171 + (t258 * t156 + t255 * t157 + t257 * t238) * t142 + (t259 * t156 + t256 * t157 + t238 * t260) * t141) * t142 / 0.2e1 + ((-t141 * t260 - t257 * t142 - t253 * t171) * t194 + ((t198 * t254 + t201 * t252) * t171 + (t198 * t258 + t201 * t255) * t142 + (t198 * t259 + t201 * t256) * t141) * t193) * t171 / 0.2e1 + ((-t176 * t200 + t179 * t203 + Icges(1,4)) * V_base(5) + (-t200 * t177 + t203 * t180 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t203 * t176 + t200 * t179 + Icges(1,2)) * V_base(5) + (t177 * t203 + t180 * t200 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t146 * t202 + t148 * t199) * t187 + (t145 * t202 + t147 * t199) * t186 + (t136 * t194 + t138 * t193) * t165 + (t135 * t194 + t137 * t193) * t164 + (t194 * t160 + t193 * t161 + t202 * t175 + t199 * t178 + Icges(2,3)) * t190) * t190 / 0.2e1 + t190 * V_base(4) * (Icges(2,5) * t203 - Icges(2,6) * t200) + V_base(5) * t190 * (Icges(2,5) * t200 + Icges(2,6) * t203) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
