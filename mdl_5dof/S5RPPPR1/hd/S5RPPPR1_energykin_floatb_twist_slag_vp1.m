% Calculate kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:02
% EndTime: 2020-01-03 11:20:04
% DurationCPUTime: 2.06s
% Computational Cost: add. (1318->299), mult. (1336->398), div. (0->0), fcn. (1228->10), ass. (0->140)
t178 = qJ(1) + pkin(7);
t170 = sin(t178);
t172 = cos(t178);
t184 = sin(qJ(1));
t185 = cos(qJ(1));
t235 = Icges(2,5) * t184 + Icges(3,5) * t170 + Icges(2,6) * t185 + Icges(3,6) * t172;
t234 = -Icges(2,5) * t185 - Icges(3,5) * t172 + Icges(2,6) * t184 + Icges(3,6) * t170;
t180 = sin(pkin(8));
t182 = cos(pkin(8));
t222 = Icges(4,4) * t182;
t196 = -Icges(4,2) * t180 + t222;
t103 = -Icges(4,6) * t172 + t170 * t196;
t223 = Icges(4,4) * t180;
t197 = Icges(4,1) * t182 - t223;
t105 = -Icges(4,5) * t172 + t170 * t197;
t224 = Icges(3,4) * t172;
t233 = Icges(3,1) * t170 - t103 * t180 + t105 * t182 + t224;
t104 = -Icges(4,6) * t170 - t172 * t196;
t106 = -Icges(4,5) * t170 - t172 * t197;
t167 = Icges(3,4) * t170;
t232 = -Icges(3,1) * t172 - t104 * t180 + t106 * t182 + t167;
t181 = cos(pkin(9));
t227 = pkin(4) * t181;
t231 = pkin(6) * t180 + t182 * t227;
t229 = pkin(1) * t184;
t228 = pkin(1) * t185;
t226 = -pkin(5) - qJ(2);
t225 = Icges(2,4) * t185;
t179 = sin(pkin(9));
t221 = t170 * t179;
t220 = t170 * t180;
t219 = t170 * t182;
t218 = t172 * t179;
t217 = t172 * t180;
t216 = t172 * t182;
t215 = t179 * t182;
t214 = t181 * t182;
t198 = pkin(3) * t182 + qJ(4) * t180;
t131 = t198 * t172;
t142 = -pkin(2) * t172 - qJ(3) * t170;
t212 = t131 - t142;
t211 = qJD(4) * t180;
t210 = qJD(5) * t180;
t173 = V_base(4) + qJD(1);
t209 = t173 * t229 + V_base(3);
t208 = V_base(6) * pkin(5) + V_base(2);
t205 = qJD(2) + V_base(1);
t153 = pkin(3) * t180 - qJ(4) * t182;
t204 = -t153 + t226;
t140 = pkin(2) * t170 - qJ(3) * t172;
t203 = -t140 - t229;
t202 = V_base(5) * t142 + t205;
t130 = t198 * t170;
t201 = -t130 + t203;
t200 = V_base(6) * qJ(2) + t173 * t228 + t208;
t199 = rSges(4,1) * t182 - rSges(4,2) * t180;
t195 = Icges(4,5) * t182 - Icges(4,6) * t180;
t151 = Icges(4,2) * t182 + t223;
t152 = Icges(4,1) * t180 + t222;
t192 = t151 * t180 - t152 * t182;
t191 = -qJD(3) * t170 + t173 * t140 + t209;
t190 = -qJD(3) * t172 + t200;
t189 = -qJD(4) * t182 - V_base(5) * t131 + t202;
t188 = V_base(6) * t153 + t170 * t211 + t190;
t187 = -(-Icges(4,3) * t172 + t170 * t195) * V_base(5) - (-Icges(4,3) * t170 - t172 * t195) * V_base(6) - (Icges(4,5) * t180 + Icges(4,6) * t182) * t173;
t186 = t173 * t130 - t172 * t211 + t191;
t177 = pkin(9) + qJ(5);
t175 = Icges(2,4) * t184;
t171 = cos(t177);
t169 = sin(t177);
t164 = -rSges(2,1) * t185 + t184 * rSges(2,2);
t163 = t184 * rSges(2,1) + rSges(2,2) * t185;
t162 = -Icges(2,1) * t185 + t175;
t161 = Icges(2,1) * t184 + t225;
t160 = Icges(2,2) * t184 - t225;
t159 = Icges(2,2) * t185 + t175;
t156 = -qJD(5) * t182 + t173;
t154 = rSges(4,1) * t180 + rSges(4,2) * t182;
t149 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t148 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t147 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t146 = t170 * t210 + V_base(5);
t145 = -t172 * t210 + V_base(6);
t143 = -rSges(3,1) * t172 + rSges(3,2) * t170;
t141 = rSges(3,1) * t170 + rSges(3,2) * t172;
t137 = Icges(3,2) * t170 - t224;
t136 = Icges(3,2) * t172 + t167;
t129 = -t172 * t214 - t221;
t128 = -t170 * t181 + t172 * t215;
t127 = t170 * t214 - t218;
t126 = -t170 * t215 - t172 * t181;
t125 = -rSges(5,3) * t182 + (rSges(5,1) * t181 - rSges(5,2) * t179) * t180;
t124 = -Icges(5,5) * t182 + (Icges(5,1) * t181 - Icges(5,4) * t179) * t180;
t123 = -Icges(5,6) * t182 + (Icges(5,4) * t181 - Icges(5,2) * t179) * t180;
t122 = -Icges(5,3) * t182 + (Icges(5,5) * t181 - Icges(5,6) * t179) * t180;
t120 = -t169 * t170 - t171 * t216;
t119 = t169 * t216 - t170 * t171;
t118 = -t169 * t172 + t171 * t219;
t117 = -t169 * t219 - t171 * t172;
t116 = -rSges(6,3) * t182 + (rSges(6,1) * t171 - rSges(6,2) * t169) * t180;
t115 = -Icges(6,5) * t182 + (Icges(6,1) * t171 - Icges(6,4) * t169) * t180;
t114 = -Icges(6,6) * t182 + (Icges(6,4) * t171 - Icges(6,2) * t169) * t180;
t113 = -Icges(6,3) * t182 + (Icges(6,5) * t171 - Icges(6,6) * t169) * t180;
t111 = V_base(6) * rSges(2,3) - t164 * t173 + t208;
t110 = t163 * t173 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t109 = -pkin(6) * t182 + t180 * t227;
t108 = -rSges(4,3) * t170 - t172 * t199;
t107 = -rSges(4,3) * t172 + t170 * t199;
t100 = -t163 * V_base(6) + t164 * V_base(5) + V_base(1);
t99 = V_base(6) * rSges(3,3) - t143 * t173 + t200;
t98 = t141 * t173 + (-rSges(3,3) + t226) * V_base(5) + t209;
t97 = -V_base(6) * t141 + V_base(5) * t143 + (-t184 * V_base(6) - t185 * V_base(5)) * pkin(1) + t205;
t96 = rSges(5,1) * t129 + rSges(5,2) * t128 - rSges(5,3) * t217;
t95 = rSges(5,1) * t127 + rSges(5,2) * t126 + rSges(5,3) * t220;
t94 = -pkin(4) * t221 - t172 * t231;
t93 = -pkin(4) * t218 + t170 * t231;
t92 = Icges(5,1) * t129 + Icges(5,4) * t128 - Icges(5,5) * t217;
t91 = Icges(5,1) * t127 + Icges(5,4) * t126 + Icges(5,5) * t220;
t90 = Icges(5,4) * t129 + Icges(5,2) * t128 - Icges(5,6) * t217;
t89 = Icges(5,4) * t127 + Icges(5,2) * t126 + Icges(5,6) * t220;
t88 = Icges(5,5) * t129 + Icges(5,6) * t128 - Icges(5,3) * t217;
t87 = Icges(5,5) * t127 + Icges(5,6) * t126 + Icges(5,3) * t220;
t86 = rSges(6,1) * t120 + rSges(6,2) * t119 - rSges(6,3) * t217;
t85 = rSges(6,1) * t118 + rSges(6,2) * t117 + rSges(6,3) * t220;
t84 = Icges(6,1) * t120 + Icges(6,4) * t119 - Icges(6,5) * t217;
t83 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t220;
t82 = Icges(6,4) * t120 + Icges(6,2) * t119 - Icges(6,6) * t217;
t81 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t220;
t80 = Icges(6,5) * t120 + Icges(6,6) * t119 - Icges(6,3) * t217;
t79 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t220;
t78 = t154 * V_base(6) + (-t108 - t142) * t173 + t190;
t77 = t107 * t173 + (-t154 + t226) * V_base(5) + t191;
t76 = (t108 - t228) * V_base(5) + (-t107 + t203) * V_base(6) + t202;
t75 = t125 * V_base(6) + (-t96 + t212) * t173 + t188;
t74 = t173 * t95 + (-t125 + t204) * V_base(5) + t186;
t73 = (t96 - t228) * V_base(5) + (t201 - t95) * V_base(6) + t189;
t72 = t109 * V_base(6) + t116 * t145 - t156 * t86 + (-t94 + t212) * t173 + t188;
t71 = -t116 * t146 + t156 * t85 + t173 * t93 + (-t109 + t204) * V_base(5) + t186;
t70 = -t145 * t85 + t146 * t86 + (t94 - t228) * V_base(5) + (t201 - t93) * V_base(6) + t189;
t1 = m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(3) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t156 * ((-t113 * t156 - t80 * t145 - t79 * t146) * t182 + ((-t114 * t169 + t115 * t171) * t156 + (-t169 * t81 + t171 * t83) * t146 + (-t169 * t82 + t171 * t84) * t145) * t180) / 0.2e1 + t146 * ((t113 * t220 + t114 * t117 + t115 * t118) * t156 + (t117 * t81 + t118 * t83 + t79 * t220) * t146 + (t117 * t82 + t118 * t84 + t220 * t80) * t145) / 0.2e1 + t145 * ((-t113 * t217 + t114 * t119 + t115 * t120) * t156 + (t119 * t81 + t120 * t83 - t217 * t79) * t146 + (t119 * t82 + t120 * t84 - t80 * t217) * t145) / 0.2e1 + (((t104 - t88) * t182 + (-t179 * t90 + t181 * t92 + t106) * t180 + t234) * V_base(6) + ((t103 - t87) * t182 + (-t179 * t89 + t181 * t91 + t105) * t180 + t235) * V_base(5) + (Icges(2,3) + Icges(3,3) + (t151 - t122) * t182 + (-t123 * t179 + t124 * t181 + t152) * t180) * t173) * t173 / 0.2e1 + (t187 * t172 + (t122 * t220 + t123 * t126 + t124 * t127 - t192 * t170 + t235) * t173 + (t126 * t90 + t127 * t92 + t137 * t172 + t160 * t185 + t184 * t162 + t170 * t232 + t220 * t88 + Icges(1,6)) * V_base(6) + (t126 * t89 + t127 * t91 + t136 * t172 + t159 * t185 + t184 * t161 + t170 * t233 + t220 * t87 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t187 * t170 + (-t122 * t217 + t123 * t128 + t124 * t129 + t192 * t172 + t234) * t173 + (t128 * t90 + t129 * t92 + t137 * t170 + t184 * t160 - t162 * t185 - t172 * t232 - t217 * t88 + Icges(1,3)) * V_base(6) + (t128 * t89 + t129 * t91 + t136 * t170 + t184 * t159 - t161 * t185 - t172 * t233 - t217 * t87 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
