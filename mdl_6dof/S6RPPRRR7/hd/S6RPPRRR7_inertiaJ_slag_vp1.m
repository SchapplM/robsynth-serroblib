% Calculate joint inertia matrix for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:13
% EndTime: 2019-03-09 02:33:17
% DurationCPUTime: 1.58s
% Computational Cost: add. (4803->288), mult. (4481->425), div. (0->0), fcn. (4611->10), ass. (0->147)
t137 = pkin(10) + qJ(4);
t130 = qJ(5) + t137;
t125 = sin(t130);
t126 = cos(t130);
t143 = sin(qJ(6));
t145 = cos(qJ(6));
t68 = rSges(7,3) * t125 + (rSges(7,1) * t145 - rSges(7,2) * t143) * t126;
t207 = pkin(5) * t126 + pkin(9) * t125 + t68;
t144 = sin(qJ(1));
t146 = cos(qJ(1));
t206 = t144 * t146;
t205 = (rSges(6,1) * t125 + rSges(6,2) * t126) * t146;
t129 = cos(t137);
t128 = sin(t137);
t194 = rSges(5,1) * t128;
t204 = (rSges(5,2) * t129 + t194) * t146;
t138 = t144 ^ 2;
t139 = t146 ^ 2;
t203 = t144 / 0.2e1;
t202 = t146 / 0.2e1;
t122 = t138 + t139;
t111 = m(5) * t122;
t110 = m(6) * t122;
t140 = sin(pkin(10));
t201 = pkin(3) * t140;
t200 = pkin(4) * t129;
t199 = pkin(5) * t125;
t142 = -pkin(7) - qJ(3);
t178 = t144 * t145;
t180 = t143 * t146;
t92 = t125 * t180 + t178;
t177 = t145 * t146;
t179 = t144 * t143;
t93 = -t125 * t177 + t179;
t166 = -t93 * rSges(7,1) - t92 * rSges(7,2);
t182 = t126 * t146;
t55 = rSges(7,3) * t182 - t166;
t198 = t146 * t55 + t139 * (pkin(9) * t126 - t199);
t184 = t125 * t144;
t116 = pkin(5) * t184;
t183 = t126 * t144;
t90 = -t125 * t179 + t177;
t91 = t125 * t178 + t180;
t195 = rSges(7,1) * t91 + rSges(7,2) * t90;
t54 = -rSges(7,3) * t183 + t195;
t197 = pkin(9) * t183 - t116 - t54;
t65 = Icges(7,3) * t125 + (Icges(7,5) * t145 - Icges(7,6) * t143) * t126;
t67 = Icges(7,5) * t125 + (Icges(7,1) * t145 - Icges(7,4) * t143) * t126;
t196 = t126 * t145 * t67 + t125 * t65;
t46 = t207 * t144;
t66 = Icges(7,6) * t125 + (Icges(7,4) * t145 - Icges(7,2) * t143) * t126;
t193 = t143 * t66;
t48 = Icges(7,5) * t91 + Icges(7,6) * t90 - Icges(7,3) * t183;
t50 = Icges(7,4) * t91 + Icges(7,2) * t90 - Icges(7,6) * t183;
t52 = Icges(7,1) * t91 + Icges(7,4) * t90 - Icges(7,5) * t183;
t21 = t125 * t48 + (-t143 * t50 + t145 * t52) * t126;
t192 = t21 * t146;
t49 = Icges(7,5) * t93 + Icges(7,6) * t92 + Icges(7,3) * t182;
t51 = Icges(7,4) * t93 + Icges(7,2) * t92 + Icges(7,6) * t182;
t53 = Icges(7,1) * t93 + Icges(7,4) * t92 + Icges(7,5) * t182;
t22 = t125 * t49 + (-t143 * t51 + t145 * t53) * t126;
t191 = t22 * t144;
t190 = rSges(4,3) + qJ(3);
t188 = Icges(5,4) * t128;
t187 = Icges(5,4) * t129;
t186 = Icges(6,4) * t125;
t185 = Icges(6,4) * t126;
t181 = t129 * t144;
t113 = pkin(4) * t128 + t201;
t136 = -pkin(8) + t142;
t176 = -t113 * t146 - t136 * t144;
t175 = t142 * t144 + t146 * t201;
t174 = pkin(1) * t146 + qJ(2) * t144;
t173 = t144 * t201;
t132 = t146 * qJ(2);
t172 = t132 - t176;
t78 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t146;
t171 = rSges(5,2) * t181 + rSges(5,3) * t146 + t144 * t194;
t170 = (-rSges(7,3) - pkin(9)) * t126;
t17 = t182 * t48 + t50 * t92 + t52 * t93;
t18 = t182 * t49 + t51 * t92 + t53 * t93;
t10 = t144 * t18 + t146 * t17;
t151 = Icges(6,5) * t125 + Icges(6,6) * t126;
t72 = Icges(6,3) * t146 + t144 * t151;
t73 = Icges(6,3) * t144 - t146 * t151;
t15 = -t183 * t48 + t50 * t90 + t52 * t91;
t16 = -t183 * t49 + t51 * t90 + t53 * t91;
t9 = t144 * t16 + t146 * t15;
t169 = (t139 * t72 + t9) * t146 + (t138 * t73 + t10 + (t144 * t72 + t146 * t73) * t146) * t144;
t26 = -t183 * t65 + t66 * t90 + t67 * t91;
t3 = t26 * t125 + (-t144 * t15 + t146 * t16) * t126;
t27 = t182 * t65 + t66 * t92 + t67 * t93;
t4 = t27 * t125 + (-t144 * t17 + t146 * t18) * t126;
t168 = t3 * t202 + t4 * t203 + t125 * (t191 + t192) / 0.2e1 - t9 * t183 / 0.2e1 + t10 * t182 / 0.2e1;
t167 = t110 + t111 + (m(4) + m(7)) * t122;
t141 = cos(pkin(10));
t165 = rSges(4,1) * t140 + rSges(4,2) * t141;
t98 = -Icges(6,2) * t125 + t185;
t99 = Icges(6,1) * t126 - t186;
t160 = t125 * t99 + t126 * t98;
t100 = rSges(6,1) * t126 - rSges(6,2) * t125;
t121 = pkin(4) * t181;
t69 = t100 * t144 + t121;
t70 = (-t100 - t200) * t146;
t157 = t144 * t69 - t146 * t70;
t156 = Icges(5,1) * t128 + t187;
t155 = Icges(6,1) * t125 + t185;
t154 = Icges(5,2) * t129 + t188;
t153 = Icges(6,2) * t126 + t186;
t152 = Icges(5,5) * t128 + Icges(5,6) * t129;
t102 = t144 * t113;
t150 = -t136 * t146 + t102 + t174;
t56 = t132 + t204 + (-rSges(5,3) - pkin(1)) * t144 + t175;
t57 = -t142 * t146 + t171 + t173 + t174;
t149 = m(5) * (t144 * t56 - t146 * t57);
t42 = t205 + (-rSges(6,3) - pkin(1)) * t144 + t172;
t43 = t150 + t78;
t148 = m(6) * (t144 * t42 - t146 * t43);
t97 = Icges(6,5) * t126 - Icges(6,6) * t125;
t147 = t192 / 0.2e1 + t191 / 0.2e1 + (-t125 * (Icges(6,6) * t144 - t146 * t153) + t126 * (Icges(6,5) * t144 - t146 * t155) + t144 * t97 - t146 * t160 + t27) * t203 + (-t125 * (Icges(6,6) * t146 + t144 * t153) + t126 * (Icges(6,5) * t146 + t144 * t155) + t160 * t144 + t146 * t97 + t26) * t202;
t118 = rSges(2,1) * t146 - rSges(2,2) * t144;
t117 = -rSges(2,1) * t144 - rSges(2,2) * t146;
t108 = rSges(5,1) * t129 - rSges(5,2) * t128;
t95 = -rSges(3,2) * t146 + rSges(3,3) * t144 + t174;
t94 = rSges(3,3) * t146 + t132 + (rSges(3,2) - pkin(1)) * t144;
t81 = Icges(5,3) * t144 - t146 * t152;
t80 = Icges(5,3) * t146 + t144 * t152;
t71 = t146 * (t144 * rSges(6,3) - t205);
t63 = t144 * t165 + t146 * t190 + t174;
t62 = t132 + t165 * t146 + (-pkin(1) - t190) * t144;
t61 = -t173 + t102 + (-t136 + t142) * t146;
t59 = t146 * (t175 + t176);
t47 = t207 * t146;
t45 = -t144 * t171 + (t144 * rSges(5,3) - t204) * t146;
t41 = -t144 * t78 + t71;
t40 = (-t207 - t200) * t146;
t39 = t121 + t46;
t34 = -t125 * t55 + t182 * t68;
t33 = t125 * t54 + t183 * t68;
t32 = t144 * t170 + t116 + t150 + t195;
t31 = -t144 * pkin(1) + (t170 + t199) * t146 + t166 + t172;
t30 = (-t126 * t193 + t196) * t125;
t29 = (-t144 * t55 - t146 * t54) * t126;
t28 = t59 + t71 + (-t61 - t78) * t144;
t23 = t144 * t197 + t198;
t12 = t59 + (-t61 + t197) * t144 + t198;
t1 = [Icges(4,1) * t141 ^ 2 - t128 * (-Icges(5,2) * t128 + t187) + t129 * (Icges(5,1) * t129 - t188) - t125 * t98 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t141 + Icges(4,2) * t140) * t140 + (t99 - t193) * t126 + m(7) * (t31 ^ 2 + t32 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t56 ^ 2 + t57 ^ 2) + m(3) * (t94 ^ 2 + t95 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(2) * (t117 ^ 2 + t118 ^ 2) + t196; m(7) * (t144 * t31 - t146 * t32) + t148 + t149 + m(3) * (t144 * t94 - t146 * t95) + m(4) * (t144 * t62 - t146 * t63); m(3) * t122 + t167; m(7) * (t144 * t32 + t146 * t31) + m(6) * (t144 * t43 + t146 * t42) + m(5) * (t144 * t57 + t146 * t56) + m(4) * (t144 * t63 + t146 * t62); 0; t167; (-t128 * (Icges(5,6) * t144 - t146 * t154) + t129 * (Icges(5,5) * t144 - t146 * t156)) * t203 + (-t128 * (Icges(5,6) * t146 + t144 * t154) + t129 * (Icges(5,5) * t146 + t144 * t156)) * t202 + m(7) * (t31 * t39 + t32 * t40) + m(6) * (t42 * t69 + t43 * t70) + t108 * t149 + (t139 / 0.2e1 + t138 / 0.2e1) * (Icges(5,5) * t129 - Icges(5,6) * t128) + t147; m(6) * t157 + m(7) * (t144 * t39 - t146 * t40) + t108 * t111; m(6) * (t144 * t70 + t146 * t69) + m(7) * (t144 * t40 + t146 * t39); m(7) * (t12 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t28 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t108 ^ 2 * t122 + t45 ^ 2) + t146 * (t139 * t80 + t206 * t81) + t144 * (t138 * t81 + t206 * t80) + t169; m(7) * (t31 * t46 - t32 * t47) + t100 * t148 + t147; m(7) * (t144 * t46 + t146 * t47) + t100 * t110; m(7) * (-t144 * t47 + t146 * t46); m(7) * (t12 * t23 + t39 * t46 - t40 * t47) + m(6) * (t100 * t157 + t41 * t28) + t169; m(6) * (t100 ^ 2 * t122 + t41 ^ 2) + m(7) * (t23 ^ 2 + t46 ^ 2 + t47 ^ 2) + t169; t30 + m(7) * (t31 * t34 + t32 * t33) + ((t22 / 0.2e1 + t27 / 0.2e1) * t146 + (-t26 / 0.2e1 - t21 / 0.2e1) * t144) * t126; m(7) * (t144 * t34 - t146 * t33); m(7) * (t144 * t33 + t146 * t34); m(7) * (t12 * t29 + t33 * t40 + t34 * t39) + t168; m(7) * (t23 * t29 - t33 * t47 + t34 * t46) + t168; t125 * t30 + m(7) * (t29 ^ 2 + t33 ^ 2 + t34 ^ 2) + (-t144 * t3 + t146 * t4 + t125 * (-t144 * t21 + t146 * t22)) * t126;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
