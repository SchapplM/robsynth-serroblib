% Calculate kinetic energy for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:06
% EndTime: 2019-12-31 18:42:08
% DurationCPUTime: 2.07s
% Computational Cost: add. (1297->243), mult. (1413->337), div. (0->0), fcn. (1305->8), ass. (0->121)
t253 = Icges(5,1) + Icges(6,1);
t252 = Icges(5,4) + Icges(6,4);
t251 = -Icges(6,5) - Icges(5,5);
t250 = Icges(5,2) + Icges(6,2);
t249 = -Icges(6,6) - Icges(5,6);
t248 = -Icges(6,3) - Icges(5,3);
t185 = qJ(1) + pkin(8);
t179 = sin(t185);
t180 = cos(t185);
t190 = cos(qJ(4));
t187 = sin(qJ(4));
t191 = cos(qJ(3));
t219 = t187 * t191;
t129 = -t179 * t219 - t180 * t190;
t218 = t190 * t191;
t221 = t180 * t187;
t130 = t179 * t218 - t221;
t188 = sin(qJ(3));
t222 = t179 * t188;
t247 = -t249 * t129 - t251 * t130 - t248 * t222;
t131 = t179 * t190 - t180 * t219;
t223 = t179 * t187;
t132 = t180 * t218 + t223;
t220 = t180 * t188;
t246 = -t249 * t131 - t251 * t132 - t248 * t220;
t245 = t250 * t129 + t252 * t130 - t249 * t222;
t244 = t250 * t131 + t252 * t132 - t249 * t220;
t243 = t252 * t129 + t253 * t130 - t251 * t222;
t242 = t252 * t131 + t253 * t132 - t251 * t220;
t241 = t248 * t191 + (t249 * t187 - t251 * t190) * t188;
t240 = t249 * t191 + (-t250 * t187 + t252 * t190) * t188;
t239 = t251 * t191 + (-t252 * t187 + t253 * t190) * t188;
t189 = sin(qJ(1));
t232 = pkin(1) * t189;
t192 = cos(qJ(1));
t231 = pkin(1) * t192;
t230 = pkin(4) * t190;
t229 = -pkin(5) - qJ(2);
t227 = Icges(2,4) * t189;
t226 = Icges(3,4) * t179;
t225 = Icges(4,4) * t188;
t224 = Icges(4,4) * t191;
t200 = qJ(5) * t188 + t191 * t230;
t217 = rSges(6,1) * t130 + rSges(6,2) * t129 + rSges(6,3) * t222 - pkin(4) * t221 + t179 * t200;
t216 = rSges(6,1) * t132 + rSges(6,2) * t131 + rSges(6,3) * t220 + pkin(4) * t223 + t180 * t200;
t215 = (-qJ(5) - rSges(6,3)) * t191 + (rSges(6,1) * t190 - rSges(6,2) * t187 + t230) * t188;
t214 = qJD(4) * t188;
t213 = qJD(5) * t188;
t181 = V_base(6) + qJD(1);
t212 = t181 * t231 + V_base(2);
t211 = V_base(5) * pkin(5) + V_base(1);
t159 = qJD(3) * t179 + V_base(4);
t153 = pkin(2) * t179 - pkin(6) * t180;
t208 = -t153 - t232;
t207 = V_base(5) * qJ(2) + t211;
t206 = V_base(4) * t232 + qJD(2) + V_base(3);
t205 = pkin(3) * t191 + pkin(7) * t188;
t158 = -qJD(3) * t180 + V_base(5);
t204 = rSges(4,1) * t191 - rSges(4,2) * t188;
t203 = Icges(4,1) * t191 - t225;
t202 = -Icges(4,2) * t188 + t224;
t201 = Icges(4,5) * t191 - Icges(4,6) * t188;
t199 = (-Icges(4,3) * t180 + t179 * t201) * t158 + (Icges(4,3) * t179 + t180 * t201) * t159 + (Icges(4,5) * t188 + Icges(4,6) * t191) * t181;
t154 = pkin(2) * t180 + pkin(6) * t179;
t198 = t181 * t154 + t229 * V_base(4) + t212;
t142 = t205 * t179;
t173 = pkin(3) * t188 - pkin(7) * t191;
t197 = t158 * t173 + (-t142 + t208) * t181 + t207;
t196 = V_base(4) * t153 + (-t154 - t231) * V_base(5) + t206;
t143 = t205 * t180;
t195 = t181 * t143 - t159 * t173 + t198;
t194 = t159 * t142 - t158 * t143 + t196;
t116 = -Icges(4,6) * t180 + t179 * t202;
t117 = Icges(4,6) * t179 + t180 * t202;
t118 = -Icges(4,5) * t180 + t179 * t203;
t119 = Icges(4,5) * t179 + t180 * t203;
t163 = Icges(4,2) * t191 + t225;
t166 = Icges(4,1) * t188 + t224;
t193 = (-t117 * t188 + t119 * t191) * t159 + (-t116 * t188 + t118 * t191) * t158 + (-t163 * t188 + t166 * t191) * t181;
t183 = Icges(2,4) * t192;
t177 = Icges(3,4) * t180;
t172 = rSges(2,1) * t192 - t189 * rSges(2,2);
t171 = t189 * rSges(2,1) + rSges(2,2) * t192;
t170 = rSges(4,1) * t188 + rSges(4,2) * t191;
t169 = -qJD(4) * t191 + t181;
t168 = Icges(2,1) * t192 - t227;
t167 = Icges(2,1) * t189 + t183;
t165 = -Icges(2,2) * t189 + t183;
t164 = Icges(2,2) * t192 + t227;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t152 = rSges(3,1) * t180 - rSges(3,2) * t179;
t151 = rSges(3,1) * t179 + rSges(3,2) * t180;
t150 = Icges(3,1) * t180 - t226;
t149 = Icges(3,1) * t179 + t177;
t148 = -Icges(3,2) * t179 + t177;
t147 = Icges(3,2) * t180 + t226;
t140 = -rSges(5,3) * t191 + (rSges(5,1) * t190 - rSges(5,2) * t187) * t188;
t128 = t180 * t214 + t159;
t127 = t179 * t214 + t158;
t123 = rSges(4,3) * t179 + t180 * t204;
t122 = -rSges(4,3) * t180 + t179 * t204;
t121 = V_base(5) * rSges(2,3) - t171 * t181 + t211;
t120 = t172 * t181 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t113 = t171 * V_base(4) - t172 * V_base(5) + V_base(3);
t111 = V_base(5) * rSges(3,3) + (-t151 - t232) * t181 + t207;
t110 = t152 * t181 + (-rSges(3,3) + t229) * V_base(4) + t212;
t109 = V_base(4) * t151 + (-t152 - t231) * V_base(5) + t206;
t108 = rSges(5,1) * t132 + rSges(5,2) * t131 + rSges(5,3) * t220;
t106 = rSges(5,1) * t130 + rSges(5,2) * t129 + rSges(5,3) * t222;
t90 = t158 * t170 + (-t122 + t208) * t181 + t207;
t89 = t123 * t181 - t159 * t170 + t198;
t88 = t159 * t122 - t158 * t123 + t196;
t87 = -t106 * t169 + t127 * t140 + t197;
t86 = t108 * t169 - t128 * t140 + t195;
t85 = t128 * t106 - t127 * t108 + t194;
t84 = t127 * t215 - t169 * t217 + t180 * t213 + t197;
t83 = -t128 * t215 + t169 * t216 + t179 * t213 + t195;
t82 = -qJD(5) * t191 - t127 * t216 + t128 * t217 + t194;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t113 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + t159 * (t199 * t179 + t193 * t180) / 0.2e1 + t158 * (t193 * t179 - t199 * t180) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + ((t240 * t129 + t239 * t130 + t241 * t222) * t169 + (t244 * t129 + t242 * t130 + t246 * t222) * t128 + (t245 * t129 + t243 * t130 + t247 * t222) * t127) * t127 / 0.2e1 + ((t240 * t131 + t239 * t132 + t241 * t220) * t169 + (t244 * t131 + t242 * t132 + t246 * t220) * t128 + (t245 * t131 + t243 * t132 + t247 * t220) * t127) * t128 / 0.2e1 + ((-t247 * t127 - t246 * t128 - t241 * t169) * t191 + ((-t240 * t187 + t239 * t190) * t169 + (-t244 * t187 + t242 * t190) * t128 + (-t245 * t187 + t243 * t190) * t127) * t188) * t169 / 0.2e1 + ((t117 * t191 + t119 * t188) * t159 + (t116 * t191 + t118 * t188) * t158 + (t191 * t163 + t188 * t166 + Icges(2,3) + Icges(3,3)) * t181) * t181 / 0.2e1 + ((-t147 * t179 + t149 * t180 - t189 * t164 + t167 * t192 + Icges(1,4)) * V_base(5) + (-t179 * t148 + t180 * t150 - t189 * t165 + t192 * t168 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t180 * t147 + t179 * t149 + t192 * t164 + t189 * t167 + Icges(1,2)) * V_base(5) + (t148 * t180 + t150 * t179 + t165 * t192 + t189 * t168 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t181 * (Icges(2,5) * t189 + Icges(3,5) * t179 + Icges(2,6) * t192 + Icges(3,6) * t180) + V_base(4) * t181 * (Icges(2,5) * t192 + Icges(3,5) * t180 - Icges(2,6) * t189 - Icges(3,6) * t179) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
