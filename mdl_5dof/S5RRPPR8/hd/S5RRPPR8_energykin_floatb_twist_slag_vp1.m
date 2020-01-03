% Calculate kinetic energy for
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:02
% EndTime: 2019-12-31 19:38:05
% DurationCPUTime: 2.82s
% Computational Cost: add. (980->263), mult. (1571->362), div. (0->0), fcn. (1509->8), ass. (0->132)
t261 = Icges(3,4) - Icges(4,5);
t260 = Icges(3,1) + Icges(4,1);
t259 = Icges(3,2) + Icges(4,3);
t195 = cos(qJ(2));
t258 = t261 * t195;
t193 = sin(qJ(2));
t257 = t261 * t193;
t256 = Icges(4,4) + Icges(3,5);
t255 = Icges(3,6) - Icges(4,6);
t254 = t259 * t193 - t258;
t253 = t260 * t195 - t257;
t194 = sin(qJ(1));
t196 = cos(qJ(1));
t252 = t254 * t194 + t255 * t196;
t251 = -t255 * t194 + t254 * t196;
t250 = t253 * t194 - t256 * t196;
t249 = t256 * t194 + t253 * t196;
t248 = -t259 * t195 - t257;
t247 = t260 * t193 + t258;
t246 = t255 * t193 - t256 * t195;
t245 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t177 = -qJD(2) * t196 + V_base(5);
t178 = qJD(2) * t194 + V_base(4);
t184 = V_base(6) + qJD(1);
t244 = (t248 * t193 + t247 * t195) * t184 + (t251 * t193 + t249 * t195) * t178 + (t252 * t193 + t250 * t195) * t177;
t191 = cos(pkin(8));
t190 = sin(pkin(8));
t230 = t190 * t195;
t148 = t191 * t193 - t230;
t136 = t148 * t194;
t231 = t190 * t193;
t209 = t191 * t195 + t231;
t137 = t209 * t194;
t138 = t148 * t196;
t139 = t209 * t196;
t243 = (Icges(5,5) * t148 - Icges(5,6) * t209 - t256 * t193 - t255 * t195) * t184 + (Icges(5,5) * t139 + Icges(5,6) * t138 - t245 * t194 + t246 * t196) * t178 + (Icges(5,5) * t137 + Icges(5,6) * t136 + t246 * t194 + t245 * t196) * t177;
t239 = pkin(3) * t193;
t238 = pkin(3) * t195;
t237 = pkin(4) * t191;
t236 = Icges(2,4) * t194;
t217 = pkin(2) * t195 + qJ(3) * t193;
t144 = t217 * t194;
t175 = t194 * pkin(1) - pkin(6) * t196;
t228 = -t144 - t175;
t145 = t217 * t196;
t154 = -t194 * qJ(4) + t196 * t238;
t227 = -t145 - t154;
t226 = qJD(3) * t193;
t225 = V_base(5) * pkin(5) + V_base(1);
t153 = qJ(4) * t196 + t194 * t238;
t222 = -t153 + t228;
t170 = pkin(2) * t193 - qJ(3) * t195;
t221 = -t170 - t239;
t220 = t177 * t170 + t196 * t226 + t225;
t219 = rSges(3,1) * t195 - rSges(3,2) * t193;
t218 = rSges(4,1) * t195 + rSges(4,3) * t193;
t189 = pkin(8) + qJ(5);
t182 = sin(t189);
t183 = cos(t189);
t143 = -t182 * t195 + t183 * t193;
t210 = t182 * t193 + t183 * t195;
t176 = pkin(1) * t196 + t194 * pkin(6);
t208 = -V_base(4) * pkin(5) + t184 * t176 + V_base(2);
t207 = V_base(4) * t175 - t176 * V_base(5) + V_base(3);
t206 = pkin(4) * t231 + t237 * t195;
t203 = t184 * t145 + t194 * t226 + t208;
t202 = -qJD(4) * t194 + t177 * t239 + t220;
t201 = -qJD(3) * t195 + t178 * t144 + t207;
t200 = qJD(4) * t196 + t184 * t154 + t203;
t199 = t178 * t153 + t201;
t187 = Icges(2,4) * t196;
t174 = rSges(2,1) * t196 - t194 * rSges(2,2);
t173 = t194 * rSges(2,1) + rSges(2,2) * t196;
t172 = rSges(3,1) * t193 + rSges(3,2) * t195;
t171 = rSges(4,1) * t193 - rSges(4,3) * t195;
t169 = Icges(2,1) * t196 - t236;
t168 = Icges(2,1) * t194 + t187;
t165 = -Icges(2,2) * t194 + t187;
t164 = Icges(2,2) * t196 + t236;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = -qJD(5) * t194 + t178;
t149 = V_base(5) + (-qJD(2) + qJD(5)) * t196;
t135 = t194 * rSges(3,3) + t219 * t196;
t134 = t194 * rSges(4,2) + t218 * t196;
t133 = -rSges(3,3) * t196 + t219 * t194;
t132 = -rSges(4,2) * t196 + t218 * t194;
t117 = -pkin(4) * t230 + t237 * t193;
t116 = t210 * t196;
t115 = t143 * t196;
t114 = t210 * t194;
t113 = t143 * t194;
t111 = V_base(5) * rSges(2,3) - t173 * t184 + t225;
t110 = t174 * t184 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t109 = t173 * V_base(4) - t174 * V_base(5) + V_base(3);
t108 = rSges(5,1) * t148 - rSges(5,2) * t209;
t107 = Icges(5,1) * t148 - Icges(5,4) * t209;
t106 = Icges(5,4) * t148 - Icges(5,2) * t209;
t104 = rSges(6,1) * t143 - rSges(6,2) * t210;
t103 = Icges(6,1) * t143 - Icges(6,4) * t210;
t102 = Icges(6,4) * t143 - Icges(6,2) * t210;
t101 = Icges(6,5) * t143 - Icges(6,6) * t210;
t100 = -pkin(7) * t194 + t206 * t196;
t99 = pkin(7) * t196 + t206 * t194;
t98 = rSges(5,1) * t139 + rSges(5,2) * t138 - rSges(5,3) * t194;
t97 = t137 * rSges(5,1) + t136 * rSges(5,2) + rSges(5,3) * t196;
t96 = Icges(5,1) * t139 + Icges(5,4) * t138 - Icges(5,5) * t194;
t95 = Icges(5,1) * t137 + Icges(5,4) * t136 + Icges(5,5) * t196;
t94 = Icges(5,4) * t139 + Icges(5,2) * t138 - Icges(5,6) * t194;
t93 = Icges(5,4) * t137 + Icges(5,2) * t136 + Icges(5,6) * t196;
t90 = rSges(6,1) * t116 + rSges(6,2) * t115 - rSges(6,3) * t194;
t89 = t114 * rSges(6,1) + t113 * rSges(6,2) + rSges(6,3) * t196;
t88 = Icges(6,1) * t116 + Icges(6,4) * t115 - Icges(6,5) * t194;
t87 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t196;
t86 = Icges(6,4) * t116 + Icges(6,2) * t115 - Icges(6,6) * t194;
t85 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t196;
t84 = Icges(6,5) * t116 + Icges(6,6) * t115 - Icges(6,3) * t194;
t83 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t196;
t82 = t172 * t177 + (-t133 - t175) * t184 + t225;
t81 = t135 * t184 - t172 * t178 + t208;
t80 = t133 * t178 - t135 * t177 + t207;
t79 = t171 * t177 + (-t132 + t228) * t184 + t220;
t78 = t134 * t184 + (-t170 - t171) * t178 + t203;
t77 = t132 * t178 + (-t134 - t145) * t177 + t201;
t76 = t108 * t177 + (-t97 + t222) * t184 + t202;
t75 = t184 * t98 + (-t108 + t221) * t178 + t200;
t74 = t178 * t97 + (-t98 + t227) * t177 + t199;
t73 = t104 * t149 + t117 * t177 + (-t89 - t99 + t222) * t184 + t202;
t72 = -t104 * t150 + (t100 + t90) * t184 + (-t117 + t221) * t178 + t200;
t71 = -t149 * t90 + t150 * t89 + t178 * t99 + (-t100 + t227) * t177 + t199;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(3) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + t150 * ((t115 * t86 + t116 * t88 - t194 * t84) * t150 + (t115 * t85 + t116 * t87 - t194 * t83) * t149 + (-t101 * t194 + t102 * t115 + t103 * t116) * t184) / 0.2e1 + t149 * ((t113 * t86 + t114 * t88 + t196 * t84) * t150 + (t113 * t85 + t114 * t87 + t196 * t83) * t149 + (t101 * t196 + t113 * t102 + t114 * t103) * t184) / 0.2e1 + ((-t194 * t164 + t168 * t196 + Icges(1,4)) * V_base(5) + (-t194 * t165 + t196 * t169 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t196 * t164 + t194 * t168 + Icges(1,2)) * V_base(5) + (t165 * t196 + t194 * t169 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t94 + t137 * t96) * t178 + (t136 * t93 + t137 * t95) * t177 + (t136 * t106 + t137 * t107) * t184 + t243 * t196 + t244 * t194) * t177 / 0.2e1 + ((t138 * t94 + t139 * t96) * t178 + (t138 * t93 + t139 * t95) * t177 + (t106 * t138 + t107 * t139) * t184 + t244 * t196 - t243 * t194) * t178 / 0.2e1 + ((t143 * t88 - t210 * t86) * t150 + (t143 * t87 - t210 * t85) * t149 + (t148 * t96 + t249 * t193 - t251 * t195 - t209 * t94) * t178 + (t148 * t95 + t250 * t193 - t252 * t195 - t209 * t93) * t177 + (-t210 * t102 + t143 * t103 - t209 * t106 + t148 * t107 + t247 * t193 - t248 * t195 + Icges(2,3)) * t184) * t184 / 0.2e1 + t184 * V_base(4) * (Icges(2,5) * t196 - Icges(2,6) * t194) + V_base(5) * t184 * (Icges(2,5) * t194 + Icges(2,6) * t196) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
