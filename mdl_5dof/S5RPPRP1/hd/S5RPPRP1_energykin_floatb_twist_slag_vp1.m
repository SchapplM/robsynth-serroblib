% Calculate kinetic energy for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:09
% EndTime: 2022-01-23 09:12:10
% DurationCPUTime: 1.93s
% Computational Cost: add. (1213->247), mult. (1371->329), div. (0->0), fcn. (1263->8), ass. (0->123)
t253 = Icges(5,1) + Icges(6,1);
t252 = Icges(5,4) + Icges(6,4);
t251 = -Icges(6,5) - Icges(5,5);
t250 = Icges(5,2) + Icges(6,2);
t249 = -Icges(6,6) - Icges(5,6);
t248 = -Icges(6,3) - Icges(5,3);
t183 = qJ(1) + pkin(7);
t177 = sin(t183);
t178 = cos(t183);
t189 = cos(qJ(4));
t185 = cos(pkin(8));
t187 = sin(qJ(4));
t219 = t185 * t187;
t134 = -t177 * t219 - t178 * t189;
t218 = t185 * t189;
t220 = t178 * t187;
t135 = t177 * t218 - t220;
t184 = sin(pkin(8));
t223 = t177 * t184;
t247 = -t249 * t134 - t251 * t135 - t248 * t223;
t136 = t177 * t189 - t178 * t219;
t222 = t177 * t187;
t137 = t178 * t218 + t222;
t221 = t178 * t184;
t246 = -t249 * t136 - t251 * t137 - t248 * t221;
t245 = t250 * t134 + t252 * t135 - t249 * t223;
t244 = t250 * t136 + t252 * t137 - t249 * t221;
t243 = t252 * t134 + t253 * t135 - t251 * t223;
t242 = t252 * t136 + t253 * t137 - t251 * t221;
t241 = t248 * t185 + (t249 * t187 - t251 * t189) * t184;
t240 = t249 * t185 + (-t250 * t187 + t252 * t189) * t184;
t239 = t251 * t185 + (-t252 * t187 + t253 * t189) * t184;
t188 = sin(qJ(1));
t190 = cos(qJ(1));
t238 = Icges(2,5) * t188 + Icges(3,5) * t177 + Icges(2,6) * t190 + Icges(3,6) * t178;
t237 = Icges(2,5) * t190 + Icges(3,5) * t178 - Icges(2,6) * t188 - Icges(3,6) * t177;
t232 = pkin(1) * t188;
t231 = pkin(1) * t190;
t230 = pkin(4) * t189;
t229 = -pkin(5) - qJ(2);
t227 = Icges(2,4) * t188;
t226 = Icges(3,4) * t177;
t225 = Icges(4,4) * t184;
t224 = Icges(4,4) * t185;
t196 = qJ(5) * t184 + t185 * t230;
t217 = rSges(6,1) * t135 + rSges(6,2) * t134 + rSges(6,3) * t223 - pkin(4) * t220 + t177 * t196;
t216 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t221 + pkin(4) * t222 + t178 * t196;
t215 = (-qJ(5) - rSges(6,3)) * t185 + (rSges(6,1) * t189 - rSges(6,2) * t187 + t230) * t184;
t214 = qJD(4) * t184;
t213 = qJD(5) * t184;
t179 = V_base(6) + qJD(1);
t212 = t179 * t231 + V_base(2);
t211 = V_base(5) * pkin(5) + V_base(1);
t148 = pkin(2) * t177 - qJ(3) * t178;
t208 = -t148 - t232;
t150 = pkin(2) * t178 + qJ(3) * t177;
t207 = -t150 - t231;
t206 = V_base(5) * qJ(2) + t211;
t205 = V_base(4) * t232 + qJD(2) + V_base(3);
t204 = qJD(3) * t177 + t206;
t203 = pkin(3) * t185 + pkin(6) * t184;
t202 = V_base(4) * t148 + t205;
t201 = rSges(4,1) * t185 - rSges(4,2) * t184;
t200 = Icges(4,1) * t185 - t225;
t199 = -Icges(4,2) * t184 + t224;
t198 = Icges(4,5) * t185 - Icges(4,6) * t184;
t197 = -qJD(3) * t178 + t179 * t150 + t212;
t195 = (-Icges(4,3) * t178 + t177 * t198) * V_base(5) + (Icges(4,3) * t177 + t178 * t198) * V_base(4) + (Icges(4,5) * t184 + Icges(4,6) * t185) * t179;
t138 = t203 * t177;
t163 = pkin(3) * t184 - pkin(6) * t185;
t194 = V_base(5) * t163 + (-t138 + t208) * t179 + t204;
t139 = t203 * t178;
t193 = V_base(4) * t138 + (-t139 + t207) * V_base(5) + t202;
t192 = t179 * t139 + (-t163 + t229) * V_base(4) + t197;
t115 = -Icges(4,6) * t178 + t177 * t199;
t116 = Icges(4,6) * t177 + t178 * t199;
t117 = -Icges(4,5) * t178 + t177 * t200;
t118 = Icges(4,5) * t177 + t178 * t200;
t159 = Icges(4,2) * t185 + t225;
t160 = Icges(4,1) * t184 + t224;
t191 = (-t116 * t184 + t118 * t185) * V_base(4) + (-t115 * t184 + t117 * t185) * V_base(5) + (-t159 * t184 + t160 * t185) * t179;
t181 = Icges(2,4) * t190;
t175 = Icges(3,4) * t178;
t171 = rSges(2,1) * t190 - t188 * rSges(2,2);
t170 = t188 * rSges(2,1) + rSges(2,2) * t190;
t169 = Icges(2,1) * t190 - t227;
t168 = Icges(2,1) * t188 + t181;
t167 = -Icges(2,2) * t188 + t181;
t166 = Icges(2,2) * t190 + t227;
t162 = -qJD(4) * t185 + t179;
t161 = rSges(4,1) * t184 + rSges(4,2) * t185;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t153 = t178 * t214 + V_base(4);
t152 = t177 * t214 + V_base(5);
t151 = rSges(3,1) * t178 - rSges(3,2) * t177;
t149 = rSges(3,1) * t177 + rSges(3,2) * t178;
t147 = Icges(3,1) * t178 - t226;
t146 = Icges(3,1) * t177 + t175;
t145 = -Icges(3,2) * t177 + t175;
t144 = Icges(3,2) * t178 + t226;
t133 = -rSges(5,3) * t185 + (rSges(5,1) * t189 - rSges(5,2) * t187) * t184;
t122 = V_base(5) * rSges(2,3) - t170 * t179 + t211;
t121 = t171 * t179 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t120 = rSges(4,3) * t177 + t178 * t201;
t119 = -rSges(4,3) * t178 + t177 * t201;
t112 = t170 * V_base(4) - t171 * V_base(5) + V_base(3);
t111 = V_base(5) * rSges(3,3) + (-t149 - t232) * t179 + t206;
t110 = t151 * t179 + (-rSges(3,3) + t229) * V_base(4) + t212;
t109 = V_base(4) * t149 + (-t151 - t231) * V_base(5) + t205;
t108 = rSges(5,1) * t137 + rSges(5,2) * t136 + rSges(5,3) * t221;
t106 = rSges(5,1) * t135 + rSges(5,2) * t134 + rSges(5,3) * t223;
t90 = t161 * V_base(5) + (-t119 + t208) * t179 + t204;
t89 = t120 * t179 + (-t161 + t229) * V_base(4) + t197;
t88 = V_base(4) * t119 + (-t120 + t207) * V_base(5) + t202;
t87 = -t106 * t162 + t133 * t152 + t194;
t86 = t108 * t162 - t133 * t153 + t192;
t85 = t153 * t106 - t152 * t108 + t193;
t84 = t152 * t215 - t162 * t217 + t178 * t213 + t194;
t83 = -t153 * t215 + t162 * t216 + t177 * t213 + t192;
t82 = -qJD(5) * t185 - t152 * t216 + t153 * t217 + t193;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t112 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(3) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + ((t134 * t240 + t135 * t239 + t223 * t241) * t162 + (t134 * t244 + t135 * t242 + t223 * t246) * t153 + (t245 * t134 + t243 * t135 + t247 * t223) * t152) * t152 / 0.2e1 + ((t136 * t240 + t137 * t239 + t221 * t241) * t162 + (t244 * t136 + t242 * t137 + t246 * t221) * t153 + (t245 * t136 + t243 * t137 + t221 * t247) * t152) * t153 / 0.2e1 + ((-t152 * t247 - t246 * t153 - t241 * t162) * t185 + ((-t187 * t240 + t189 * t239) * t162 + (-t187 * t244 + t189 * t242) * t153 + (-t187 * t245 + t189 * t243) * t152) * t184) * t162 / 0.2e1 + ((t115 * t185 + t117 * t184 + t238) * V_base(5) + (t116 * t185 + t118 * t184 + t237) * V_base(4) + (t185 * t159 + t184 * t160 + Icges(2,3) + Icges(3,3)) * t179) * t179 / 0.2e1 + (t195 * t177 + t191 * t178 + t237 * t179 + (-t144 * t177 + t146 * t178 - t188 * t166 + t168 * t190 + Icges(1,4)) * V_base(5) + (-t177 * t145 + t178 * t147 - t188 * t167 + t190 * t169 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t191 * t177 - t195 * t178 + t238 * t179 + (t178 * t144 + t177 * t146 + t190 * t166 + t188 * t168 + Icges(1,2)) * V_base(5) + (t145 * t178 + t147 * t177 + t167 * t190 + t188 * t169 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
