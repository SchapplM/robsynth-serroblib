% Calculate kinetic energy for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:34
% EndTime: 2019-12-05 17:41:36
% DurationCPUTime: 1.98s
% Computational Cost: add. (1226->250), mult. (928->334), div. (0->0), fcn. (708->10), ass. (0->130)
t169 = qJ(1) + pkin(8);
t159 = sin(t169);
t161 = cos(t169);
t173 = sin(qJ(1));
t174 = cos(qJ(1));
t237 = -Icges(2,5) * t173 - Icges(3,5) * t159 - Icges(2,6) * t174 - Icges(3,6) * t161;
t236 = Icges(2,5) * t174 + Icges(3,5) * t161 - Icges(2,6) * t173 - Icges(3,6) * t159;
t170 = sin(pkin(9));
t171 = cos(pkin(9));
t219 = Icges(4,4) * t171;
t192 = -Icges(4,2) * t170 + t219;
t100 = Icges(4,6) * t161 - t159 * t192;
t220 = Icges(4,4) * t170;
t195 = Icges(4,1) * t171 - t220;
t102 = Icges(4,5) * t161 - t159 * t195;
t221 = Icges(3,4) * t161;
t235 = -Icges(3,1) * t159 - t100 * t170 + t102 * t171 - t221;
t101 = Icges(4,6) * t159 + t161 * t192;
t103 = Icges(4,5) * t159 + t161 * t195;
t222 = Icges(3,4) * t159;
t234 = Icges(3,1) * t161 - t101 * t170 + t103 * t171 - t222;
t137 = qJD(4) * t159 + V_base(6);
t109 = qJD(5) * t159 + t137;
t138 = qJD(4) * t161 + V_base(5);
t110 = qJD(5) * t161 + t138;
t168 = pkin(9) + qJ(4);
t162 = qJ(5) + t168;
t156 = cos(t162);
t155 = sin(t162);
t216 = Icges(6,4) * t155;
t113 = Icges(6,2) * t156 + t216;
t215 = Icges(6,4) * t156;
t114 = Icges(6,1) * t155 + t215;
t163 = V_base(4) + qJD(1);
t190 = -Icges(6,2) * t155 + t215;
t83 = Icges(6,6) * t161 - t159 * t190;
t84 = Icges(6,6) * t159 + t161 * t190;
t193 = Icges(6,1) * t156 - t216;
t85 = Icges(6,5) * t161 - t159 * t193;
t86 = Icges(6,5) * t159 + t161 * t193;
t233 = (t113 * t155 - t114 * t156) * t163 + (t155 * t83 - t156 * t85) * t110 + (t155 * t84 - t156 * t86) * t109;
t160 = cos(t168);
t158 = sin(t168);
t218 = Icges(5,4) * t158;
t119 = Icges(5,2) * t160 + t218;
t217 = Icges(5,4) * t160;
t122 = Icges(5,1) * t158 + t217;
t191 = -Icges(5,2) * t158 + t217;
t91 = Icges(5,6) * t161 - t159 * t191;
t92 = Icges(5,6) * t159 + t161 * t191;
t194 = Icges(5,1) * t160 - t218;
t93 = Icges(5,5) * t161 - t159 * t194;
t94 = Icges(5,5) * t159 + t161 * t194;
t232 = (t119 * t158 - t122 * t160) * t163 + (t158 * t91 - t160 * t93) * t138 + (t158 * t92 - t160 * t94) * t137;
t230 = pkin(1) * t173;
t229 = pkin(1) * t174;
t228 = pkin(3) * t170;
t227 = pkin(4) * t158;
t226 = t171 * pkin(3);
t225 = -pkin(5) - qJ(2);
t224 = Icges(2,4) * t173;
t223 = Icges(2,4) * t174;
t213 = pkin(4) * t160;
t211 = V_base(6) * pkin(5) + V_base(2);
t128 = pkin(2) * t161 + qJ(3) * t159;
t208 = -t128 - t229;
t126 = -pkin(2) * t159 + qJ(3) * t161;
t207 = qJD(3) * t159 + t163 * t126 + V_base(3);
t206 = V_base(6) * qJ(2) + t211;
t80 = pkin(6) * t159 + t161 * t226;
t205 = t208 - t80;
t204 = qJD(3) * t161 + t206;
t203 = V_base(5) * t229 + V_base(6) * t230 + qJD(2) + V_base(1);
t202 = rSges(4,1) * t171 - rSges(4,2) * t170;
t201 = rSges(5,1) * t160 - rSges(5,2) * t158;
t200 = rSges(6,1) * t156 - rSges(6,2) * t155;
t189 = Icges(4,5) * t171 - Icges(4,6) * t170;
t188 = Icges(5,5) * t160 - Icges(5,6) * t158;
t187 = Icges(6,5) * t156 - Icges(6,6) * t155;
t135 = Icges(4,2) * t171 + t220;
t136 = Icges(4,1) * t170 + t219;
t182 = t135 * t170 - t136 * t171;
t181 = V_base(6) * t228 + t204;
t180 = V_base(5) * t128 + t203;
t179 = (Icges(6,3) * t159 + t161 * t187) * t109 + (Icges(6,3) * t161 - t159 * t187) * t110 + (Icges(6,5) * t155 + Icges(6,6) * t156) * t163;
t178 = (Icges(5,5) * t158 + Icges(5,6) * t160) * t163 + (Icges(5,3) * t159 + t161 * t188) * t137 + (Icges(5,3) * t161 - t159 * t188) * t138;
t177 = (Icges(4,5) * t170 + Icges(4,6) * t171) * t163 + (Icges(4,3) * t161 - t159 * t189) * V_base(5) + (Icges(4,3) * t159 + t161 * t189) * V_base(6);
t79 = pkin(6) * t161 - t159 * t226;
t176 = t163 * t79 + (t225 - t228) * V_base(5) + t207;
t175 = V_base(5) * t80 + (-t126 - t79) * V_base(6) + t180;
t147 = rSges(2,1) * t174 - t173 * rSges(2,2);
t146 = -t173 * rSges(2,1) - rSges(2,2) * t174;
t145 = Icges(2,1) * t174 - t224;
t144 = -Icges(2,1) * t173 - t223;
t143 = -Icges(2,2) * t173 + t223;
t142 = -Icges(2,2) * t174 - t224;
t139 = rSges(4,1) * t170 + rSges(4,2) * t171;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t161 - rSges(3,2) * t159;
t127 = -rSges(3,1) * t159 - rSges(3,2) * t161;
t125 = rSges(5,1) * t158 + rSges(5,2) * t160;
t121 = -Icges(3,2) * t159 + t221;
t120 = -Icges(3,2) * t161 - t222;
t115 = rSges(6,1) * t155 + rSges(6,2) * t156;
t107 = V_base(6) * rSges(2,3) - t147 * t163 + t211;
t106 = t146 * t163 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t105 = rSges(4,3) * t159 + t161 * t202;
t104 = rSges(4,3) * t161 - t159 * t202;
t97 = -t146 * V_base(6) + t147 * V_base(5) + V_base(1);
t96 = rSges(5,3) * t159 + t161 * t201;
t95 = rSges(5,3) * t161 - t159 * t201;
t88 = rSges(6,3) * t159 + t161 * t200;
t87 = rSges(6,3) * t161 - t159 * t200;
t76 = V_base(6) * rSges(3,3) + (-t129 - t229) * t163 + t206;
t75 = V_base(3) + (t127 - t230) * t163 + (-rSges(3,3) + t225) * V_base(5);
t74 = pkin(7) * t159 + t161 * t213;
t73 = pkin(7) * t161 - t159 * t213;
t72 = -t127 * V_base(6) + t129 * V_base(5) + t203;
t71 = V_base(6) * t139 + (-t105 + t208) * t163 + t204;
t70 = (t104 - t230) * t163 + (-t139 + t225) * V_base(5) + t207;
t69 = t105 * V_base(5) + (-t104 - t126) * V_base(6) + t180;
t68 = t137 * t125 + (t205 - t96) * t163 + t181;
t67 = -t125 * t138 + (t95 - t230) * t163 + t176;
t66 = -t137 * t95 + t138 * t96 + t175;
t65 = t137 * t227 + t109 * t115 + (t205 - t74 - t88) * t163 + t181;
t64 = -t138 * t227 - t110 * t115 + (t73 + t87 - t230) * t163 + t176;
t63 = -t109 * t87 + t110 * t88 - t137 * t73 + t138 * t74 + t175;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t138 * (t232 * t159 + t178 * t161) / 0.2e1 + t137 * (t178 * t159 - t232 * t161) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t110 * (t233 * t159 + t179 * t161) / 0.2e1 + t109 * (t179 * t159 - t233 * t161) / 0.2e1 + (t177 * t161 + (t182 * t159 + t237) * t163 + (-t121 * t161 - t143 * t174 - t173 * t145 - t234 * t159 + Icges(1,6)) * V_base(6) + (-t161 * t120 - t174 * t142 - t173 * t144 - t235 * t159 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t177 * t159 + (-t182 * t161 + t236) * t163 + (-t159 * t121 - t173 * t143 + t174 * t145 + t234 * t161 + Icges(1,3)) * V_base(6) + (-t120 * t159 - t173 * t142 + t144 * t174 + t235 * t161 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t158 * t93 + t160 * t91) * t138 + (t158 * t94 + t160 * t92) * t137 + (t155 * t85 + t156 * t83) * t110 + (t155 * t86 + t156 * t84) * t109 + (t101 * t171 + t103 * t170 + t236) * V_base(6) + (t100 * t171 + t102 * t170 + t237) * V_base(5) + (t156 * t113 + t155 * t114 + t160 * t119 + t158 * t122 + t171 * t135 + t170 * t136 + Icges(2,3) + Icges(3,3)) * t163) * t163 / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
