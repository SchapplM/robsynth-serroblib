% Calculate kinetic energy for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:43
% EndTime: 2022-01-23 09:31:46
% DurationCPUTime: 3.04s
% Computational Cost: add. (1371->308), mult. (1951->415), div. (0->0), fcn. (1901->8), ass. (0->148)
t284 = Icges(5,1) + Icges(6,1);
t283 = Icges(5,4) + Icges(6,4);
t282 = -Icges(6,5) - Icges(5,5);
t281 = Icges(5,2) + Icges(6,2);
t280 = -Icges(6,6) - Icges(5,6);
t279 = -Icges(6,3) - Icges(5,3);
t216 = qJ(3) + qJ(4);
t207 = sin(t216);
t208 = cos(t216);
t222 = cos(qJ(1));
t218 = cos(pkin(8));
t220 = sin(qJ(1));
t257 = t218 * t220;
t161 = -t207 * t257 - t208 * t222;
t162 = -t207 * t222 + t208 * t257;
t217 = sin(pkin(8));
t259 = t217 * t220;
t278 = -t280 * t161 - t282 * t162 - t279 * t259;
t256 = t218 * t222;
t163 = -t207 * t256 + t208 * t220;
t164 = t207 * t220 + t208 * t256;
t258 = t217 * t222;
t277 = -t280 * t163 - t282 * t164 - t279 * t258;
t276 = t281 * t161 + t283 * t162 - t280 * t259;
t275 = t281 * t163 + t283 * t164 - t280 * t258;
t274 = t283 * t161 + t284 * t162 - t282 * t259;
t273 = t283 * t163 + t284 * t164 - t282 * t258;
t272 = t279 * t218 + (t280 * t207 - t282 * t208) * t217;
t271 = t280 * t218 + (-t281 * t207 + t283 * t208) * t217;
t270 = t282 * t218 + (-t283 * t207 + t284 * t208) * t217;
t223 = pkin(7) + pkin(6);
t219 = sin(qJ(3));
t265 = pkin(3) * t219;
t221 = cos(qJ(3));
t264 = pkin(3) * t221;
t263 = pkin(2) * t218 + pkin(1);
t262 = Icges(2,4) * t220;
t261 = Icges(3,4) * t217;
t260 = Icges(3,4) * t218;
t255 = t219 * t220;
t254 = t219 * t222;
t253 = t220 * t221;
t252 = t221 * t222;
t212 = t220 * pkin(1);
t160 = t217 * t223 + t218 * t264 + t263;
t239 = -pkin(2) - t264;
t178 = pkin(4) * t208 - t239;
t215 = -qJ(5) - t223;
t232 = t178 * t218 - t215 * t217 - t160;
t204 = qJ(2) + t265;
t247 = pkin(4) * t207 - t204 + t265;
t251 = t212 + (-qJ(2) - t247) * t222 + t232 * t220 + rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t259;
t209 = t220 * qJ(2);
t213 = t222 * pkin(1);
t200 = t213 + t209;
t250 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t258 + t220 * t247 + t222 * t232 + t200;
t249 = (t215 + t223 - rSges(6,3)) * t218 + (rSges(6,1) * t208 - rSges(6,2) * t207 + t178 + t239) * t217;
t185 = pkin(6) * t217 + t263;
t248 = t160 - t185;
t246 = qJD(3) * t217;
t245 = qJD(4) * t217;
t244 = qJD(5) * t217;
t198 = -qJ(2) * t222 + t212;
t243 = V_base(4) * t198 + V_base(3);
t242 = V_base(5) * pkin(5) + V_base(1);
t184 = t222 * t246 + V_base(4);
t183 = t220 * t246 + V_base(5);
t205 = V_base(6) + qJD(1);
t238 = qJD(2) * t220 + t242;
t237 = rSges(3,1) * t218 - rSges(3,2) * t217;
t236 = Icges(3,1) * t218 - t261;
t235 = -Icges(3,2) * t217 + t260;
t234 = Icges(3,5) * t218 - Icges(3,6) * t217;
t233 = -qJD(2) * t222 + t205 * t200 + V_base(2);
t165 = t185 * t220 - t212;
t191 = t217 * pkin(2) - t218 * pkin(6);
t231 = V_base(5) * t191 + (-t165 - t198) * t205 + t238;
t166 = t185 * t222 - t213;
t230 = V_base(4) * t165 + (-t166 - t200) * V_base(5) + t243;
t229 = (-Icges(3,3) * t222 + t220 * t234) * V_base(5) + (Icges(3,3) * t220 + t222 * t234) * V_base(4) + (Icges(3,5) * t217 + Icges(3,6) * t218) * t205;
t228 = t205 * t166 + (-pkin(5) - t191) * V_base(4) + t233;
t107 = (qJ(2) - t204) * t222 + t248 * t220;
t169 = t217 * t264 + (pkin(6) - t223) * t218;
t190 = -qJD(3) * t218 + t205;
t227 = -t107 * t190 + t183 * t169 + t231;
t108 = t204 * t220 + t222 * t248 - t209;
t226 = t184 * t107 - t108 * t183 + t230;
t225 = t190 * t108 - t169 * t184 + t228;
t153 = -Icges(3,6) * t222 + t220 * t235;
t154 = Icges(3,6) * t220 + t222 * t235;
t155 = -Icges(3,5) * t222 + t220 * t236;
t156 = Icges(3,5) * t220 + t222 * t236;
t187 = Icges(3,2) * t218 + t261;
t188 = Icges(3,1) * t217 + t260;
t224 = (-t154 * t217 + t156 * t218) * V_base(4) + (-t153 * t217 + t155 * t218) * V_base(5) + (-t187 * t217 + t188 * t218) * t205;
t211 = Icges(2,4) * t222;
t201 = rSges(2,1) * t222 - rSges(2,2) * t220;
t199 = rSges(2,1) * t220 + rSges(2,2) * t222;
t197 = Icges(2,1) * t222 - t262;
t196 = Icges(2,1) * t220 + t211;
t195 = -Icges(2,2) * t220 + t211;
t194 = Icges(2,2) * t222 + t262;
t193 = Icges(2,5) * t222 - Icges(2,6) * t220;
t192 = Icges(2,5) * t220 + Icges(2,6) * t222;
t189 = rSges(3,1) * t217 + rSges(3,2) * t218;
t181 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t180 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t179 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t174 = (-qJD(3) - qJD(4)) * t218 + t205;
t173 = t218 * t252 + t255;
t172 = -t218 * t254 + t253;
t171 = t218 * t253 - t254;
t170 = -t218 * t255 - t252;
t168 = t222 * t245 + t184;
t167 = t220 * t245 + t183;
t159 = rSges(3,3) * t220 + t222 * t237;
t158 = -rSges(3,3) * t222 + t220 * t237;
t157 = -rSges(4,3) * t218 + (rSges(4,1) * t221 - rSges(4,2) * t219) * t217;
t150 = -Icges(4,5) * t218 + (Icges(4,1) * t221 - Icges(4,4) * t219) * t217;
t149 = -Icges(4,6) * t218 + (Icges(4,4) * t221 - Icges(4,2) * t219) * t217;
t148 = -Icges(4,3) * t218 + (Icges(4,5) * t221 - Icges(4,6) * t219) * t217;
t146 = -rSges(5,3) * t218 + (rSges(5,1) * t208 - rSges(5,2) * t207) * t217;
t137 = V_base(5) * rSges(2,3) - t199 * t205 + t242;
t136 = t201 * t205 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t135 = t199 * V_base(4) - t201 * V_base(5) + V_base(3);
t132 = rSges(4,1) * t173 + rSges(4,2) * t172 + rSges(4,3) * t258;
t131 = rSges(4,1) * t171 + rSges(4,2) * t170 + rSges(4,3) * t259;
t130 = Icges(4,1) * t173 + Icges(4,4) * t172 + Icges(4,5) * t258;
t129 = Icges(4,1) * t171 + Icges(4,4) * t170 + Icges(4,5) * t259;
t128 = Icges(4,4) * t173 + Icges(4,2) * t172 + Icges(4,6) * t258;
t127 = Icges(4,4) * t171 + Icges(4,2) * t170 + Icges(4,6) * t259;
t126 = Icges(4,5) * t173 + Icges(4,6) * t172 + Icges(4,3) * t258;
t125 = Icges(4,5) * t171 + Icges(4,6) * t170 + Icges(4,3) * t259;
t124 = rSges(5,1) * t164 + rSges(5,2) * t163 + rSges(5,3) * t258;
t122 = rSges(5,1) * t162 + rSges(5,2) * t161 + rSges(5,3) * t259;
t104 = t189 * V_base(5) + (-t158 - t198) * t205 + t238;
t103 = t159 * t205 + (-pkin(5) - t189) * V_base(4) + t233;
t102 = t158 * V_base(4) + (-t159 - t200) * V_base(5) + t243;
t99 = -t131 * t190 + t157 * t183 + t231;
t98 = t132 * t190 - t157 * t184 + t228;
t97 = t131 * t184 - t132 * t183 + t230;
t96 = -t122 * t174 + t146 * t167 + t227;
t95 = t124 * t174 - t146 * t168 + t225;
t94 = t122 * t168 - t124 * t167 + t226;
t93 = t167 * t249 - t174 * t251 + t222 * t244 + t227;
t92 = -t168 * t249 + t174 * t250 + t220 * t244 + t225;
t91 = -qJD(5) * t218 - t167 * t250 + t168 * t251 + t226;
t1 = m(1) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(2) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(3) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t184 * ((t126 * t258 + t172 * t128 + t173 * t130) * t184 + (t125 * t258 + t127 * t172 + t129 * t173) * t183 + (t148 * t258 + t149 * t172 + t150 * t173) * t190) / 0.2e1 + t183 * ((t126 * t259 + t128 * t170 + t130 * t171) * t184 + (t125 * t259 + t170 * t127 + t171 * t129) * t183 + (t148 * t259 + t149 * t170 + t150 * t171) * t190) / 0.2e1 + t190 * ((-t125 * t183 - t126 * t184 - t148 * t190) * t218 + ((-t128 * t219 + t130 * t221) * t184 + (-t127 * t219 + t129 * t221) * t183 + (-t149 * t219 + t150 * t221) * t190) * t217) / 0.2e1 + m(5) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + ((t271 * t161 + t270 * t162 + t272 * t259) * t174 + (t275 * t161 + t273 * t162 + t277 * t259) * t168 + (t276 * t161 + t274 * t162 + t278 * t259) * t167) * t167 / 0.2e1 + ((t271 * t163 + t270 * t164 + t272 * t258) * t174 + (t275 * t163 + t273 * t164 + t277 * t258) * t168 + (t276 * t163 + t274 * t164 + t278 * t258) * t167) * t168 / 0.2e1 + ((-t278 * t167 - t277 * t168 - t272 * t174) * t218 + ((-t271 * t207 + t270 * t208) * t174 + (-t275 * t207 + t273 * t208) * t168 + (-t276 * t207 + t274 * t208) * t167) * t217) * t174 / 0.2e1 + ((t153 * t218 + t155 * t217 + t192) * V_base(5) + (t154 * t218 + t156 * t217 + t193) * V_base(4) + (t218 * t187 + t217 * t188 + Icges(2,3)) * t205) * t205 / 0.2e1 + (t193 * t205 + t229 * t220 + t224 * t222 + (-t194 * t220 + t196 * t222 + Icges(1,4)) * V_base(5) + (-t195 * t220 + t197 * t222 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t192 * t205 + t224 * t220 - t229 * t222 + (t194 * t222 + t196 * t220 + Icges(1,2)) * V_base(5) + (t195 * t222 + t197 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
