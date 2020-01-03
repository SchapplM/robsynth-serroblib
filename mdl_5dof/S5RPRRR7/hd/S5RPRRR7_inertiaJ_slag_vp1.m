% Calculate joint inertia matrix for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:03:17
% DurationCPUTime: 1.71s
% Computational Cost: add. (6043->297), mult. (5847->430), div. (0->0), fcn. (6348->10), ass. (0->162)
t143 = sin(qJ(3));
t204 = Icges(4,5) * t143;
t203 = t204 / 0.2e1;
t145 = cos(qJ(4));
t134 = t145 * pkin(4) + pkin(3);
t148 = -pkin(8) - pkin(7);
t140 = qJ(1) + pkin(9);
t136 = cos(t140);
t146 = cos(qJ(3));
t178 = t136 * t146;
t179 = t136 * t143;
t135 = sin(t140);
t142 = sin(qJ(4));
t182 = t135 * t142;
t141 = qJ(4) + qJ(5);
t138 = cos(t141);
t137 = sin(t141);
t176 = t137 * t146;
t97 = t135 * t138 - t136 * t176;
t175 = t138 * t146;
t98 = t135 * t137 + t136 * t175;
t66 = t98 * rSges(6,1) + t97 * rSges(6,2) + rSges(6,3) * t179;
t202 = pkin(4) * t182 + t134 * t178 - t148 * t179 + t66;
t201 = t135 ^ 2;
t200 = t136 ^ 2;
t181 = t135 * t143;
t95 = -t135 * t176 - t136 * t138;
t96 = t135 * t175 - t136 * t137;
t59 = Icges(6,5) * t96 + Icges(6,6) * t95 + Icges(6,3) * t181;
t61 = Icges(6,4) * t96 + Icges(6,2) * t95 + Icges(6,6) * t181;
t63 = Icges(6,1) * t96 + Icges(6,4) * t95 + Icges(6,5) * t181;
t19 = t59 * t181 + t95 * t61 + t96 * t63;
t60 = Icges(6,5) * t98 + Icges(6,6) * t97 + Icges(6,3) * t179;
t62 = Icges(6,4) * t98 + Icges(6,2) * t97 + Icges(6,6) * t179;
t64 = Icges(6,1) * t98 + Icges(6,4) * t97 + Icges(6,5) * t179;
t20 = t60 * t181 + t95 * t62 + t96 * t64;
t100 = -Icges(6,6) * t146 + (Icges(6,4) * t138 - Icges(6,2) * t137) * t143;
t101 = -Icges(6,5) * t146 + (Icges(6,1) * t138 - Icges(6,4) * t137) * t143;
t99 = -Icges(6,3) * t146 + (Icges(6,5) * t138 - Icges(6,6) * t137) * t143;
t39 = t95 * t100 + t96 * t101 + t99 * t181;
t5 = -t39 * t146 + (t135 * t19 + t136 * t20) * t143;
t21 = t59 * t179 + t97 * t61 + t98 * t63;
t22 = t60 * t179 + t97 * t62 + t98 * t64;
t40 = t97 * t100 + t98 * t101 + t99 * t179;
t6 = -t40 * t146 + (t135 * t21 + t136 * t22) * t143;
t199 = t6 * t179 + t5 * t181;
t198 = t135 / 0.2e1;
t197 = -t136 / 0.2e1;
t196 = t136 / 0.2e1;
t195 = -t146 / 0.2e1;
t120 = t143 * rSges(4,1) + t146 * rSges(4,2);
t194 = m(4) * t120;
t193 = pkin(3) * t146;
t144 = sin(qJ(1));
t192 = t144 * pkin(1);
t191 = -pkin(3) + t134;
t190 = pkin(7) + t148;
t170 = pkin(3) * t178 + pkin(7) * t179;
t189 = -t170 + t202;
t102 = -t146 * rSges(6,3) + (rSges(6,1) * t138 - rSges(6,2) * t137) * t143;
t158 = -t96 * rSges(6,1) - t95 * rSges(6,2);
t65 = rSges(6,3) * t181 - t158;
t46 = t102 * t181 + t146 * t65;
t188 = t201 * (pkin(7) * t143 + t193) + t136 * t170;
t187 = t136 * rSges(4,3);
t177 = t137 * t100;
t81 = t143 * t138 * t101;
t48 = -t143 * t177 - t146 * t99 + t81;
t186 = t48 * t146;
t94 = t191 * t143 + t190 * t146;
t185 = -t102 - t94;
t183 = Icges(4,4) * t146;
t180 = t136 * t142;
t110 = -Icges(5,6) * t146 + (Icges(5,4) * t145 - Icges(5,2) * t142) * t143;
t174 = t142 * t110;
t173 = t142 * t146;
t172 = t145 * t146;
t112 = -t146 * rSges(5,3) + (rSges(5,1) * t145 - rSges(5,2) * t142) * t143;
t126 = t143 * pkin(3) - t146 * pkin(7);
t171 = -t112 - t126;
t169 = pkin(4) * t180;
t105 = -t135 * t173 - t136 * t145;
t106 = t135 * t172 - t180;
t67 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t181;
t69 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t181;
t71 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t181;
t33 = -t146 * t67 + (-t142 * t69 + t145 * t71) * t143;
t109 = -Icges(5,3) * t146 + (Icges(5,5) * t145 - Icges(5,6) * t142) * t143;
t111 = -Icges(5,5) * t146 + (Icges(5,1) * t145 - Icges(5,4) * t142) * t143;
t44 = t105 * t110 + t106 * t111 + t109 * t181;
t168 = t44 / 0.2e1 + t33 / 0.2e1;
t107 = t135 * t145 - t136 * t173;
t108 = t136 * t172 + t182;
t68 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t179;
t70 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t179;
t72 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t179;
t34 = -t146 * t68 + (-t142 * t70 + t145 * t72) * t143;
t45 = t107 * t110 + t108 * t111 + t109 * t179;
t167 = t45 / 0.2e1 + t34 / 0.2e1;
t166 = -t126 + t185;
t74 = t108 * rSges(5,1) + t107 * rSges(5,2) + rSges(5,3) * t179;
t147 = cos(qJ(1));
t139 = t147 * pkin(1);
t165 = t136 * pkin(2) + t135 * pkin(6) + t139;
t164 = t181 / 0.2e1;
t163 = t179 / 0.2e1;
t162 = t136 * pkin(6) - t192;
t29 = -t146 * t59 + (-t137 * t61 + t138 * t63) * t143;
t30 = -t146 * t60 + (-t137 * t62 + t138 * t64) * t143;
t161 = (t29 + t39) * t164 + (t30 + t40) * t163;
t9 = -t186 + (t135 * t29 + t136 * t30) * t143;
t160 = -t146 * t9 + t199;
t12 = t20 * t135 - t19 * t136;
t13 = t22 * t135 - t21 * t136;
t159 = t12 * t164 + t13 * t163 + t5 * t197 + t6 * t198 + (t30 * t135 - t29 * t136) * t195;
t157 = rSges(4,1) * t146 - rSges(4,2) * t143;
t156 = -t106 * rSges(5,1) - t105 * rSges(5,2);
t152 = -Icges(4,2) * t143 + t183;
t151 = Icges(4,5) * t146 - Icges(4,6) * t143;
t150 = rSges(4,1) * t178 - rSges(4,2) * t179 + t135 * rSges(4,3);
t122 = t147 * rSges(2,1) - t144 * rSges(2,2);
t121 = -t144 * rSges(2,1) - t147 * rSges(2,2);
t117 = Icges(4,6) * t146 + t204;
t114 = t136 * rSges(3,1) - t135 * rSges(3,2) + t139;
t113 = -t135 * rSges(3,1) - t136 * rSges(3,2) - t192;
t89 = t143 * t145 * t111;
t84 = Icges(4,3) * t135 + t151 * t136;
t83 = -Icges(4,3) * t136 + t151 * t135;
t80 = t171 * t136;
t79 = t171 * t135;
t78 = t150 + t165;
t77 = t187 + (-pkin(2) - t157) * t135 + t162;
t75 = -t169 + (-t190 * t143 + t191 * t146) * t135;
t73 = rSges(5,3) * t181 - t156;
t58 = t136 * t150 + (t157 * t135 - t187) * t135;
t56 = t65 * t179;
t55 = t166 * t136;
t54 = t166 * t135;
t53 = -t146 * t109 - t143 * t174 + t89;
t52 = -t112 * t179 - t146 * t74;
t51 = t112 * t181 + t146 * t73;
t50 = t165 + t74 + t170;
t49 = (-t193 - pkin(2) + (-rSges(5,3) - pkin(7)) * t143) * t135 + t156 + t162;
t47 = -t102 * t179 - t146 * t66;
t43 = t165 + t202;
t42 = t169 + (-t134 * t146 - pkin(2) + (-rSges(6,3) + t148) * t143) * t135 + t158 + t162;
t41 = (-t135 * t74 + t136 * t73) * t143;
t36 = -t66 * t181 + t56;
t35 = t135 * t73 + t136 * t74 + t188;
t32 = -t189 * t146 + t185 * t179;
t31 = t146 * t75 + t94 * t181 + t46;
t28 = t107 * t70 + t108 * t72 + t68 * t179;
t27 = t107 * t69 + t108 * t71 + t67 * t179;
t26 = t105 * t70 + t106 * t72 + t68 * t181;
t25 = t105 * t69 + t106 * t71 + t67 * t181;
t18 = t56 + (-t189 * t135 + t136 * t75) * t143;
t17 = t189 * t136 + (t65 + t75) * t135 + t188;
t15 = t28 * t135 - t27 * t136;
t14 = t26 * t135 - t25 * t136;
t8 = -t45 * t146 + (t135 * t27 + t136 * t28) * t143;
t7 = -t44 * t146 + (t135 * t25 + t136 * t26) * t143;
t1 = [Icges(2,3) + Icges(3,3) + t81 + t89 + (Icges(4,4) * t143 + Icges(4,2) * t146 - t109 - t99) * t146 + (Icges(4,1) * t143 - t174 - t177 + t183) * t143 + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t49 ^ 2 + t50 ^ 2) + m(4) * (t77 ^ 2 + t78 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2) + m(2) * (t121 ^ 2 + t122 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(6) * (t42 * t55 + t43 * t54) + m(5) * (t49 * t80 + t50 * t79) + (-t39 / 0.2e1 - t29 / 0.2e1 + t136 * t203 + (-Icges(4,6) * t136 + t152 * t135) * t195 - t77 * t194 + t117 * t196 - t168) * t136 + (t40 / 0.2e1 + t30 / 0.2e1 + t135 * t203 + t146 * (Icges(4,6) * t135 + t152 * t136) / 0.2e1 - t78 * t194 + t117 * t198 + t167) * t135; m(4) * t58 + m(5) * t35 + m(6) * t17; m(6) * (t17 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t35 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(4) * (t58 ^ 2 + (t200 + t201) * t120 ^ 2) + (t201 * t84 + t13 + t15) * t135 + (-t200 * t83 - t12 - t14 + (-t135 * t83 + t136 * t84) * t135) * t136; (-t48 - t53) * t146 + m(6) * (t31 * t42 + t32 * t43) + m(5) * (t49 * t51 + t50 * t52) + (t168 * t135 + t167 * t136) * t143 + t161; m(5) * t41 + m(6) * t18; t8 * t198 + (t34 * t135 - t33 * t136) * t195 + t7 * t197 + (t14 * t198 + t15 * t196) * t143 + m(6) * (t17 * t18 + t31 * t55 + t32 * t54) + m(5) * (t35 * t41 + t51 * t80 + t52 * t79) + t159; (t53 * t146 - t9) * t146 + m(6) * (t18 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t41 ^ 2 + t51 ^ 2 + t52 ^ 2) + (t135 * t7 + t136 * t8 - t146 * (t135 * t33 + t136 * t34)) * t143 + t199; -t186 + m(6) * (t42 * t46 + t43 * t47) + t161; m(6) * t36; m(6) * (t17 * t36 + t46 * t55 + t47 * t54) + t159; m(6) * (t18 * t36 + t31 * t46 + t32 * t47) + t160; m(6) * (t36 ^ 2 + t46 ^ 2 + t47 ^ 2) + t160;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
