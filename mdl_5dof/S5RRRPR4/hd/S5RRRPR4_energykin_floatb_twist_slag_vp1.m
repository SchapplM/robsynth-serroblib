% Calculate kinetic energy for
% S5RRRPR4
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:10:56
% EndTime: 2019-12-31 21:10:58
% DurationCPUTime: 2.38s
% Computational Cost: add. (1102->235), mult. (1140->319), div. (0->0), fcn. (998->8), ass. (0->119)
t245 = Icges(4,4) - Icges(5,5);
t244 = Icges(4,1) + Icges(5,1);
t243 = Icges(4,2) + Icges(5,3);
t181 = sin(qJ(3));
t242 = t245 * t181;
t184 = cos(qJ(3));
t241 = t245 * t184;
t240 = Icges(5,4) + Icges(4,5);
t239 = Icges(4,6) - Icges(5,6);
t238 = t243 * t181 - t241;
t237 = t244 * t184 - t242;
t236 = Icges(5,2) + Icges(4,3);
t179 = qJ(1) + qJ(2);
t174 = sin(t179);
t175 = cos(t179);
t235 = t238 * t174 + t239 * t175;
t234 = -t239 * t174 + t238 * t175;
t233 = t237 * t174 - t240 * t175;
t232 = t240 * t174 + t237 * t175;
t231 = -t243 * t184 - t242;
t230 = t244 * t181 + t241;
t229 = -t239 * t181 + t240 * t184;
t147 = -qJD(3) * t175 + V_base(5);
t148 = qJD(3) * t174 + V_base(4);
t173 = V_base(6) + qJD(1);
t172 = qJD(2) + t173;
t228 = (t231 * t181 + t230 * t184) * t172 + (t234 * t181 + t232 * t184) * t148 + (t235 * t181 + t233 * t184) * t147;
t227 = (t240 * t181 + t239 * t184) * t172 + (t236 * t174 + t229 * t175) * t148 + (t229 * t174 - t236 * t175) * t147;
t226 = -pkin(5) - pkin(6);
t182 = sin(qJ(1));
t222 = pkin(1) * t182;
t185 = cos(qJ(1));
t221 = pkin(1) * t185;
t220 = pkin(4) * t181;
t219 = pkin(4) * t184;
t218 = Icges(2,4) * t182;
t217 = Icges(3,4) * t174;
t203 = pkin(3) * t184 + qJ(4) * t181;
t125 = t203 * t174;
t140 = pkin(2) * t174 - pkin(7) * t175;
t212 = -t125 - t140;
t211 = qJD(4) * t181;
t210 = t173 * t221 + V_base(2);
t209 = V_base(4) * t222 + V_base(3);
t208 = V_base(5) * pkin(5) + V_base(1);
t205 = rSges(4,1) * t184 - rSges(4,2) * t181;
t204 = rSges(5,1) * t184 + rSges(5,3) * t181;
t180 = sin(qJ(5));
t183 = cos(qJ(5));
t143 = -t180 * t184 + t181 * t183;
t196 = t180 * t181 + t183 * t184;
t195 = V_base(5) * pkin(6) - t173 * t222 + t208;
t141 = pkin(2) * t175 + pkin(7) * t174;
t192 = t172 * t141 + t226 * V_base(4) + t210;
t161 = pkin(3) * t181 - qJ(4) * t184;
t191 = t147 * t161 + t175 * t211 + t195;
t190 = V_base(4) * t140 + (-t141 - t221) * V_base(5) + t209;
t126 = t203 * t175;
t189 = t172 * t126 + t174 * t211 + t192;
t188 = -qJD(4) * t184 + t148 * t125 + t190;
t176 = Icges(2,4) * t185;
t171 = Icges(3,4) * t175;
t165 = rSges(2,1) * t185 - t182 * rSges(2,2);
t164 = t182 * rSges(2,1) + rSges(2,2) * t185;
t163 = rSges(4,1) * t181 + rSges(4,2) * t184;
t162 = rSges(5,1) * t181 - rSges(5,3) * t184;
t160 = Icges(2,1) * t185 - t218;
t159 = Icges(2,1) * t182 + t176;
t156 = -Icges(2,2) * t182 + t176;
t155 = Icges(2,2) * t185 + t218;
t146 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t145 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t144 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t139 = rSges(3,1) * t175 - rSges(3,2) * t174;
t138 = rSges(3,1) * t174 + rSges(3,2) * t175;
t137 = Icges(3,1) * t175 - t217;
t136 = Icges(3,1) * t174 + t171;
t135 = -Icges(3,2) * t174 + t171;
t134 = Icges(3,2) * t175 + t217;
t131 = -pkin(8) * t174 + t175 * t219;
t130 = pkin(8) * t175 + t174 * t219;
t128 = -qJD(5) * t174 + t148;
t127 = V_base(5) + (-qJD(3) + qJD(5)) * t175;
t122 = t196 * t175;
t121 = t143 * t175;
t120 = t196 * t174;
t119 = t143 * t174;
t118 = rSges(4,3) * t174 + t175 * t205;
t117 = rSges(5,2) * t174 + t175 * t204;
t116 = -rSges(4,3) * t175 + t174 * t205;
t115 = -rSges(5,2) * t175 + t174 * t204;
t101 = V_base(5) * rSges(2,3) - t164 * t173 + t208;
t100 = t165 * t173 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t99 = t164 * V_base(4) - t165 * V_base(5) + V_base(3);
t97 = rSges(6,1) * t143 - rSges(6,2) * t196;
t96 = Icges(6,1) * t143 - Icges(6,4) * t196;
t95 = Icges(6,4) * t143 - Icges(6,2) * t196;
t94 = Icges(6,5) * t143 - Icges(6,6) * t196;
t93 = V_base(5) * rSges(3,3) - t138 * t172 + t195;
t92 = t139 * t172 + (-rSges(3,3) + t226) * V_base(4) + t210;
t91 = V_base(4) * t138 + (-t139 - t221) * V_base(5) + t209;
t90 = rSges(6,1) * t122 + rSges(6,2) * t121 - rSges(6,3) * t174;
t89 = rSges(6,1) * t120 + rSges(6,2) * t119 + rSges(6,3) * t175;
t88 = Icges(6,1) * t122 + Icges(6,4) * t121 - Icges(6,5) * t174;
t87 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t175;
t86 = Icges(6,4) * t122 + Icges(6,2) * t121 - Icges(6,6) * t174;
t85 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t175;
t84 = Icges(6,5) * t122 + Icges(6,6) * t121 - Icges(6,3) * t174;
t83 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t175;
t82 = t147 * t163 + (-t116 - t140) * t172 + t195;
t81 = t118 * t172 - t148 * t163 + t192;
t80 = t148 * t116 - t147 * t118 + t190;
t79 = t147 * t162 + (-t115 + t212) * t172 + t191;
t78 = t117 * t172 + (-t161 - t162) * t148 + t189;
t77 = t148 * t115 + (-t117 - t126) * t147 + t188;
t76 = t147 * t220 + t127 * t97 + (-t130 - t89 + t212) * t172 + t191;
t75 = -t128 * t97 + (t131 + t90) * t172 + (-t161 - t220) * t148 + t189;
t74 = -t127 * t90 + t128 * t89 + t148 * t130 + (-t126 - t131) * t147 + t188;
t1 = m(1) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t128 * ((t121 * t86 + t122 * t88 - t174 * t84) * t128 + (t121 * t85 + t122 * t87 - t174 * t83) * t127 + (t121 * t95 + t122 * t96 - t174 * t94) * t172) / 0.2e1 + t127 * ((t119 * t86 + t120 * t88 + t175 * t84) * t128 + (t119 * t85 + t120 * t87 + t175 * t83) * t127 + (t119 * t95 + t120 * t96 + t175 * t94) * t172) / 0.2e1 + (t228 * t174 - t227 * t175) * t147 / 0.2e1 + (t227 * t174 + t228 * t175) * t148 / 0.2e1 + ((-t134 * t174 + t136 * t175 - t182 * t155 + t159 * t185 + Icges(1,4)) * V_base(5) + (-t135 * t174 + t137 * t175 - t182 * t156 + t160 * t185 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t134 * t175 + t136 * t174 + t155 * t185 + t182 * t159 + Icges(1,2)) * V_base(5) + (t135 * t175 + t137 * t174 + t156 * t185 + t182 * t160 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t143 * t88 - t196 * t86) * t128 + (t143 * t87 - t196 * t85) * t127 + (t232 * t181 - t234 * t184) * t148 + (t233 * t181 - t235 * t184) * t147 + (t143 * t96 + t230 * t181 - t231 * t184 - t196 * t95 + Icges(3,3)) * t172) * t172 / 0.2e1 + V_base(4) * t172 * (Icges(3,5) * t175 - Icges(3,6) * t174) + V_base(5) * t172 * (Icges(3,5) * t174 + Icges(3,6) * t175) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t182 + Icges(2,6) * t185) * V_base(5) + (Icges(2,5) * t185 - Icges(2,6) * t182) * V_base(4) + Icges(2,3) * t173 / 0.2e1) * t173;
T = t1;
