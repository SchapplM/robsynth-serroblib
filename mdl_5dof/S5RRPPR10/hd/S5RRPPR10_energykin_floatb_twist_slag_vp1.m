% Calculate kinetic energy for
% S5RRPPR10
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:58
% EndTime: 2019-12-31 19:43:00
% DurationCPUTime: 2.50s
% Computational Cost: add. (1045->290), mult. (2167->405), div. (0->0), fcn. (2281->8), ass. (0->132)
t257 = Icges(4,1) + Icges(5,1);
t256 = -Icges(4,4) + Icges(5,5);
t255 = Icges(5,4) + Icges(4,5);
t254 = Icges(4,2) + Icges(5,3);
t253 = -Icges(5,6) + Icges(4,6);
t252 = -Icges(4,3) - Icges(5,2);
t201 = sin(pkin(8));
t202 = cos(pkin(8));
t208 = cos(qJ(1));
t205 = sin(qJ(1));
t207 = cos(qJ(2));
t234 = t205 * t207;
t160 = t201 * t234 + t202 * t208;
t161 = -t201 * t208 + t202 * t234;
t204 = sin(qJ(2));
t236 = t204 * t205;
t251 = t254 * t160 + t256 * t161 - t253 * t236;
t233 = t207 * t208;
t162 = t201 * t233 - t205 * t202;
t163 = t205 * t201 + t202 * t233;
t235 = t204 * t208;
t250 = t254 * t162 + t256 * t163 - t253 * t235;
t249 = -t253 * t160 + t255 * t161 - t252 * t236;
t248 = -t253 * t162 + t255 * t163 - t252 * t235;
t247 = t256 * t160 + t257 * t161 + t255 * t236;
t246 = t256 * t162 + t257 * t163 + t255 * t235;
t245 = t253 * t207 + (t254 * t201 + t256 * t202) * t204;
t244 = t252 * t207 + (-t253 * t201 + t255 * t202) * t204;
t243 = -t255 * t207 + (t256 * t201 + t257 * t202) * t204;
t239 = Icges(2,4) * t205;
t238 = Icges(3,4) * t204;
t237 = Icges(3,4) * t207;
t130 = pkin(3) * t163 + qJ(4) * t162;
t221 = pkin(2) * t207 + qJ(3) * t204;
t168 = t221 * t208;
t232 = -t130 - t168;
t166 = (pkin(3) * t202 + qJ(4) * t201) * t204;
t185 = pkin(2) * t204 - qJ(3) * t207;
t231 = -t166 - t185;
t167 = t221 * t205;
t189 = t205 * pkin(1) - pkin(6) * t208;
t230 = -t167 - t189;
t229 = qJD(3) * t204;
t228 = qJD(5) * t204;
t227 = V_base(5) * pkin(5) + V_base(1);
t129 = pkin(3) * t161 + qJ(4) * t160;
t224 = -t129 + t230;
t192 = qJD(2) * t205 + V_base(4);
t197 = V_base(6) + qJD(1);
t191 = -qJD(2) * t208 + V_base(5);
t223 = t191 * t185 + t208 * t229 + t227;
t222 = rSges(3,1) * t207 - rSges(3,2) * t204;
t220 = Icges(3,1) * t207 - t238;
t219 = -Icges(3,2) * t204 + t237;
t218 = Icges(3,5) * t207 - Icges(3,6) * t204;
t190 = pkin(1) * t208 + t205 * pkin(6);
t217 = -V_base(4) * pkin(5) + t197 * t190 + V_base(2);
t216 = V_base(4) * t189 - t190 * V_base(5) + V_base(3);
t215 = qJD(4) * t162 + t191 * t166 + t223;
t214 = (-Icges(3,3) * t208 + t205 * t218) * t191 + (Icges(3,3) * t205 + t208 * t218) * t192 + (Icges(3,5) * t204 + Icges(3,6) * t207) * t197;
t213 = t197 * t168 + t205 * t229 + t217;
t212 = -qJD(3) * t207 + t192 * t167 + t216;
t211 = qJD(4) * t160 + t197 * t130 + t213;
t210 = qJD(4) * t204 * t201 + t192 * t129 + t212;
t149 = -Icges(3,6) * t208 + t205 * t219;
t150 = Icges(3,6) * t205 + t208 * t219;
t151 = -Icges(3,5) * t208 + t205 * t220;
t152 = Icges(3,5) * t205 + t208 * t220;
t178 = Icges(3,2) * t207 + t238;
t181 = Icges(3,1) * t204 + t237;
t209 = (-t150 * t204 + t152 * t207) * t192 + (-t149 * t204 + t151 * t207) * t191 + (-t178 * t204 + t181 * t207) * t197;
t206 = cos(qJ(5));
t203 = sin(qJ(5));
t199 = Icges(2,4) * t208;
t188 = rSges(2,1) * t208 - t205 * rSges(2,2);
t187 = t205 * rSges(2,1) + rSges(2,2) * t208;
t186 = rSges(3,1) * t204 + rSges(3,2) * t207;
t184 = qJD(5) * t207 + t197;
t183 = Icges(2,1) * t208 - t239;
t182 = Icges(2,1) * t205 + t199;
t180 = -Icges(2,2) * t205 + t199;
t179 = Icges(2,2) * t208 + t239;
t174 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t173 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t172 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t171 = pkin(4) * t202 * t204 + pkin(7) * t207;
t165 = -t208 * t228 + t192;
t164 = -t205 * t228 + t191;
t156 = (t201 * t203 + t202 * t206) * t204;
t155 = (t201 * t206 - t202 * t203) * t204;
t154 = t205 * rSges(3,3) + t208 * t222;
t153 = -rSges(3,3) * t208 + t205 * t222;
t146 = -rSges(4,3) * t207 + (rSges(4,1) * t202 - rSges(4,2) * t201) * t204;
t145 = -rSges(5,2) * t207 + (rSges(5,1) * t202 + rSges(5,3) * t201) * t204;
t135 = t163 * pkin(4) - pkin(7) * t235;
t134 = pkin(4) * t161 - pkin(7) * t236;
t133 = V_base(5) * rSges(2,3) - t187 * t197 + t227;
t132 = t188 * t197 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t131 = t187 * V_base(4) - t188 * V_base(5) + V_base(3);
t128 = t162 * t203 + t163 * t206;
t127 = t162 * t206 - t163 * t203;
t126 = t160 * t203 + t161 * t206;
t125 = t160 * t206 - t161 * t203;
t123 = t163 * rSges(4,1) - t162 * rSges(4,2) + rSges(4,3) * t235;
t122 = t163 * rSges(5,1) + rSges(5,2) * t235 + t162 * rSges(5,3);
t121 = rSges(4,1) * t161 - rSges(4,2) * t160 + rSges(4,3) * t236;
t120 = rSges(5,1) * t161 + rSges(5,2) * t236 + rSges(5,3) * t160;
t106 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t207;
t105 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t207;
t104 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t207;
t103 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t207;
t102 = t186 * t191 + (-t153 - t189) * t197 + t227;
t101 = t154 * t197 - t186 * t192 + t217;
t100 = t153 * t192 - t154 * t191 + t216;
t99 = t128 * rSges(6,1) + t127 * rSges(6,2) - rSges(6,3) * t235;
t98 = rSges(6,1) * t126 + rSges(6,2) * t125 - rSges(6,3) * t236;
t97 = Icges(6,1) * t128 + Icges(6,4) * t127 - Icges(6,5) * t235;
t96 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t236;
t95 = Icges(6,4) * t128 + Icges(6,2) * t127 - Icges(6,6) * t235;
t94 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t236;
t93 = Icges(6,5) * t128 + Icges(6,6) * t127 - Icges(6,3) * t235;
t92 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t236;
t91 = t146 * t191 + (-t121 + t230) * t197 + t223;
t90 = t123 * t197 + (-t146 - t185) * t192 + t213;
t89 = t121 * t192 + (-t123 - t168) * t191 + t212;
t88 = t145 * t191 + (-t120 + t224) * t197 + t215;
t87 = t122 * t197 + (-t145 + t231) * t192 + t211;
t86 = t120 * t192 + (-t122 + t232) * t191 + t210;
t85 = t106 * t164 + t171 * t191 - t184 * t98 + (-t134 + t224) * t197 + t215;
t84 = -t106 * t165 + t135 * t197 + t184 * t99 + (-t171 + t231) * t192 + t211;
t83 = t134 * t192 - t164 * t99 + t165 * t98 + (-t135 + t232) * t191 + t210;
t1 = m(1) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(2) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t165 * ((t127 * t95 + t128 * t97 - t93 * t235) * t165 + (t127 * t94 + t128 * t96 - t235 * t92) * t164 + (-t103 * t235 + t127 * t104 + t128 * t105) * t184) / 0.2e1 + t164 * ((t125 * t95 + t126 * t97 - t236 * t93) * t165 + (t125 * t94 + t126 * t96 - t92 * t236) * t164 + (-t103 * t236 + t104 * t125 + t105 * t126) * t184) / 0.2e1 + t184 * ((t155 * t95 + t156 * t97 + t207 * t93) * t165 + (t155 * t94 + t156 * t96 + t207 * t92) * t164 + (t103 * t207 + t104 * t155 + t105 * t156) * t184) / 0.2e1 + ((-t205 * t179 + t182 * t208 + Icges(1,4)) * V_base(5) + (-t205 * t180 + t183 * t208 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t179 * t208 + t205 * t182 + Icges(1,2)) * V_base(5) + (t180 * t208 + t205 * t183 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t209 * t205 - t214 * t208 + (t160 * t245 + t161 * t243 + t236 * t244) * t197 + (t160 * t250 + t161 * t246 + t236 * t248) * t192 + (t251 * t160 + t247 * t161 + t249 * t236) * t191) * t191 / 0.2e1 + (t214 * t205 + t209 * t208 + (t162 * t245 + t163 * t243 + t235 * t244) * t197 + (t250 * t162 + t246 * t163 + t248 * t235) * t192 + (t162 * t251 + t247 * t163 + t249 * t235) * t191) * t192 / 0.2e1 + (((t150 - t248) * t192 + (t149 - t249) * t191) * t207 + ((t201 * t250 + t202 * t246 + t152) * t192 + (t201 * t251 + t247 * t202 + t151) * t191) * t204 + (Icges(2,3) + (t178 - t244) * t207 + (t201 * t245 + t202 * t243 + t181) * t204) * t197) * t197 / 0.2e1 + t197 * V_base(4) * (Icges(2,5) * t208 - Icges(2,6) * t205) + V_base(5) * t197 * (Icges(2,5) * t205 + Icges(2,6) * t208) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
