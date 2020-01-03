% Calculate kinetic energy for
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:21
% EndTime: 2019-12-31 21:32:23
% DurationCPUTime: 2.76s
% Computational Cost: add. (1108->288), mult. (2272->422), div. (0->0), fcn. (2386->8), ass. (0->131)
t255 = Icges(4,1) + Icges(5,1);
t254 = -Icges(4,4) + Icges(5,5);
t253 = Icges(5,4) + Icges(4,5);
t252 = Icges(4,2) + Icges(5,3);
t251 = -Icges(5,6) + Icges(4,6);
t250 = -Icges(4,3) - Icges(5,2);
t203 = sin(qJ(3));
t207 = cos(qJ(3));
t209 = cos(qJ(1));
t205 = sin(qJ(1));
t208 = cos(qJ(2));
t231 = t205 * t208;
t162 = t203 * t231 + t207 * t209;
t163 = -t203 * t209 + t207 * t231;
t204 = sin(qJ(2));
t233 = t204 * t205;
t249 = t252 * t162 + t254 * t163 - t251 * t233;
t230 = t208 * t209;
t164 = t203 * t230 - t205 * t207;
t165 = t205 * t203 + t207 * t230;
t232 = t204 * t209;
t248 = t252 * t164 + t254 * t165 - t251 * t232;
t247 = -t251 * t162 + t253 * t163 - t250 * t233;
t246 = -t251 * t164 + t253 * t165 - t250 * t232;
t245 = t254 * t162 + t255 * t163 + t253 * t233;
t244 = t254 * t164 + t255 * t165 + t253 * t232;
t243 = t251 * t208 + (t252 * t203 + t254 * t207) * t204;
t242 = t250 * t208 + (-t251 * t203 + t253 * t207) * t204;
t241 = -t253 * t208 + (t254 * t203 + t255 * t207) * t204;
t236 = Icges(2,4) * t205;
t235 = Icges(3,4) * t204;
t234 = Icges(3,4) * t208;
t229 = qJD(3) * t204;
t228 = qJD(5) * t204;
t227 = V_base(5) * pkin(5) + V_base(1);
t193 = qJD(2) * t205 + V_base(4);
t198 = V_base(6) + qJD(1);
t161 = t209 * t229 + t193;
t224 = pkin(2) * t208 + pkin(7) * t204;
t192 = -qJD(2) * t209 + V_base(5);
t223 = rSges(3,1) * t208 - rSges(3,2) * t204;
t222 = Icges(3,1) * t208 - t235;
t221 = -Icges(3,2) * t204 + t234;
t220 = Icges(3,5) * t208 - Icges(3,6) * t204;
t160 = t205 * t229 + t192;
t191 = pkin(1) * t209 + t205 * pkin(6);
t219 = -V_base(4) * pkin(5) + t198 * t191 + V_base(2);
t190 = t205 * pkin(1) - pkin(6) * t209;
t218 = V_base(4) * t190 - t191 * V_base(5) + V_base(3);
t168 = t224 * t205;
t189 = pkin(2) * t204 - pkin(7) * t208;
t217 = t192 * t189 + (-t168 - t190) * t198 + t227;
t216 = (-Icges(3,3) * t209 + t205 * t220) * t192 + (Icges(3,3) * t205 + t209 * t220) * t193 + (Icges(3,5) * t204 + Icges(3,6) * t208) * t198;
t169 = t224 * t209;
t215 = t198 * t169 - t189 * t193 + t219;
t166 = (pkin(3) * t207 + qJ(4) * t203) * t204;
t214 = qJD(4) * t164 + t160 * t166 + t217;
t213 = t193 * t168 - t169 * t192 + t218;
t128 = pkin(3) * t165 + qJ(4) * t164;
t185 = -qJD(3) * t208 + t198;
t212 = qJD(4) * t162 + t185 * t128 + t215;
t127 = pkin(3) * t163 + qJ(4) * t162;
t211 = qJD(4) * t204 * t203 + t161 * t127 + t213;
t145 = -Icges(3,6) * t209 + t205 * t221;
t146 = Icges(3,6) * t205 + t209 * t221;
t149 = -Icges(3,5) * t209 + t205 * t222;
t150 = Icges(3,5) * t205 + t209 * t222;
t179 = Icges(3,2) * t208 + t235;
t182 = Icges(3,1) * t204 + t234;
t210 = (-t146 * t204 + t150 * t208) * t193 + (-t145 * t204 + t149 * t208) * t192 + (-t179 * t204 + t182 * t208) * t198;
t206 = cos(qJ(5));
t202 = sin(qJ(5));
t200 = Icges(2,4) * t209;
t188 = rSges(2,1) * t209 - t205 * rSges(2,2);
t187 = t205 * rSges(2,1) + rSges(2,2) * t209;
t186 = rSges(3,1) * t204 + rSges(3,2) * t208;
t184 = Icges(2,1) * t209 - t236;
t183 = Icges(2,1) * t205 + t200;
t181 = -Icges(2,2) * t205 + t200;
t180 = Icges(2,2) * t209 + t236;
t175 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t174 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t173 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = pkin(4) * t204 * t207 + pkin(8) * t208;
t170 = (-qJD(3) + qJD(5)) * t208 + t198;
t156 = (t202 * t203 + t206 * t207) * t204;
t155 = (-t202 * t207 + t203 * t206) * t204;
t154 = t205 * rSges(3,3) + t209 * t223;
t153 = -rSges(3,3) * t209 + t205 * t223;
t152 = -rSges(4,3) * t208 + (rSges(4,1) * t207 - rSges(4,2) * t203) * t204;
t151 = -rSges(5,2) * t208 + (rSges(5,1) * t207 + rSges(5,3) * t203) * t204;
t137 = -t209 * t228 + t161;
t136 = -t205 * t228 + t160;
t134 = t165 * pkin(4) - pkin(8) * t232;
t133 = pkin(4) * t163 - pkin(8) * t233;
t132 = V_base(5) * rSges(2,3) - t187 * t198 + t227;
t131 = t188 * t198 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t130 = t187 * V_base(4) - t188 * V_base(5) + V_base(3);
t126 = t164 * t202 + t165 * t206;
t125 = t164 * t206 - t165 * t202;
t124 = t162 * t202 + t163 * t206;
t123 = t162 * t206 - t163 * t202;
t122 = t165 * rSges(4,1) - t164 * rSges(4,2) + rSges(4,3) * t232;
t121 = t165 * rSges(5,1) + rSges(5,2) * t232 + t164 * rSges(5,3);
t120 = rSges(4,1) * t163 - rSges(4,2) * t162 + rSges(4,3) * t233;
t119 = rSges(5,1) * t163 + rSges(5,2) * t233 + rSges(5,3) * t162;
t106 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t208;
t104 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t208;
t103 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t208;
t102 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t208;
t100 = t186 * t192 + (-t153 - t190) * t198 + t227;
t99 = t154 * t198 - t186 * t193 + t219;
t98 = t153 * t193 - t154 * t192 + t218;
t97 = t126 * rSges(6,1) + t125 * rSges(6,2) - rSges(6,3) * t232;
t96 = rSges(6,1) * t124 + rSges(6,2) * t123 - rSges(6,3) * t233;
t95 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t232;
t94 = Icges(6,1) * t124 + Icges(6,4) * t123 - Icges(6,5) * t233;
t93 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t232;
t92 = Icges(6,4) * t124 + Icges(6,2) * t123 - Icges(6,6) * t233;
t91 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t232;
t90 = Icges(6,5) * t124 + Icges(6,6) * t123 - Icges(6,3) * t233;
t89 = -t120 * t185 + t152 * t160 + t217;
t88 = t122 * t185 - t152 * t161 + t215;
t87 = t120 * t161 - t122 * t160 + t213;
t86 = t151 * t160 + (-t119 - t127) * t185 + t214;
t85 = t121 * t185 + (-t151 - t166) * t161 + t212;
t84 = t119 * t161 + (-t121 - t128) * t160 + t211;
t83 = t106 * t136 + t160 * t172 - t170 * t96 + (-t127 - t133) * t185 + t214;
t82 = -t106 * t137 + t134 * t185 + t170 * t97 + (-t166 - t172) * t161 + t212;
t81 = t133 * t161 - t136 * t97 + t137 * t96 + (-t128 - t134) * t160 + t211;
t1 = m(1) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(2) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t193 * (t216 * t205 + t210 * t209) / 0.2e1 + t192 * (t210 * t205 - t216 * t209) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t137 * ((t125 * t93 + t126 * t95 - t91 * t232) * t137 + (t125 * t92 + t126 * t94 - t232 * t90) * t136 + (-t102 * t232 + t125 * t103 + t126 * t104) * t170) / 0.2e1 + t136 * ((t123 * t93 + t124 * t95 - t233 * t91) * t137 + (t123 * t92 + t124 * t94 - t90 * t233) * t136 + (-t102 * t233 + t103 * t123 + t104 * t124) * t170) / 0.2e1 + t170 * ((t155 * t93 + t156 * t95 + t208 * t91) * t137 + (t155 * t92 + t156 * t94 + t208 * t90) * t136 + (t102 * t208 + t103 * t155 + t104 * t156) * t170) / 0.2e1 + ((t162 * t243 + t163 * t241 + t233 * t242) * t185 + (t162 * t248 + t163 * t244 + t233 * t246) * t161 + (t249 * t162 + t245 * t163 + t247 * t233) * t160) * t160 / 0.2e1 + ((t164 * t243 + t165 * t241 + t232 * t242) * t185 + (t248 * t164 + t244 * t165 + t246 * t232) * t161 + (t164 * t249 + t245 * t165 + t247 * t232) * t160) * t161 / 0.2e1 + ((-t160 * t247 - t161 * t246 - t185 * t242) * t208 + ((t203 * t243 + t207 * t241) * t185 + (t203 * t248 + t207 * t244) * t161 + (t203 * t249 + t245 * t207) * t160) * t204) * t185 / 0.2e1 + ((t146 * t208 + t150 * t204) * t193 + (t145 * t208 + t149 * t204) * t192 + (t179 * t208 + t182 * t204 + Icges(2,3)) * t198) * t198 / 0.2e1 + ((-t205 * t180 + t183 * t209 + Icges(1,4)) * V_base(5) + (-t205 * t181 + t184 * t209 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t180 * t209 + t205 * t183 + Icges(1,2)) * V_base(5) + (t181 * t209 + t205 * t184 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t198 * (Icges(2,5) * t209 - Icges(2,6) * t205) + V_base(5) * t198 * (Icges(2,5) * t205 + Icges(2,6) * t209) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
