% Calculate joint inertia matrix for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR13_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:21
% EndTime: 2024-09-27 17:30:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (9099->253), mult. (5527->359), div. (0->0), fcn. (5133->16), ass. (0->131)
t136 = sin(pkin(5));
t137 = cos(pkin(5));
t138 = sin(qJ(4));
t140 = cos(qJ(4));
t135 = qJ(1) + qJ(2);
t132 = qJ(3) + t135;
t124 = cos(t132);
t158 = t124 * t136;
t123 = sin(t132);
t155 = t137 * t140;
t86 = -t123 * t138 + t124 * t155;
t156 = t137 * t138;
t87 = t123 * t140 + t124 * t156;
t54 = Icges(5,5) * t87 + Icges(5,6) * t86 - Icges(5,3) * t158;
t56 = Icges(5,4) * t87 + Icges(5,2) * t86 - Icges(5,6) * t158;
t58 = Icges(5,1) * t87 + Icges(5,4) * t86 - Icges(5,5) * t158;
t90 = Icges(5,3) * t137 + (Icges(5,5) * t138 + Icges(5,6) * t140) * t136;
t91 = Icges(5,6) * t137 + (Icges(5,4) * t138 + Icges(5,2) * t140) * t136;
t92 = Icges(5,5) * t137 + (Icges(5,1) * t138 + Icges(5,4) * t140) * t136;
t166 = t137 * t54 + (t138 * t58 + t140 * t56) * t136 - t90 * t158 + t86 * t91 + t87 * t92;
t160 = t123 * t136;
t88 = -t123 * t155 - t124 * t138;
t89 = -t123 * t156 + t124 * t140;
t55 = Icges(5,5) * t89 + Icges(5,6) * t88 + Icges(5,3) * t160;
t57 = Icges(5,4) * t89 + Icges(5,2) * t88 + Icges(5,6) * t160;
t59 = Icges(5,1) * t89 + Icges(5,4) * t88 + Icges(5,5) * t160;
t165 = t137 * t55 + (t138 * t59 + t140 * t57) * t136 + t90 * t160 + t88 * t91 + t89 * t92;
t134 = qJ(4) + qJ(5);
t127 = pkin(5) - t134;
t116 = cos(t127) / 0.2e1;
t126 = pkin(5) + t134;
t121 = cos(t126);
t102 = t116 - t121 / 0.2e1;
t115 = sin(t126) / 0.2e1;
t120 = sin(t127);
t99 = t115 + t120 / 0.2e1;
t67 = Icges(6,5) * t102 + Icges(6,6) * t99 + Icges(6,3) * t137;
t68 = Icges(6,4) * t102 + Icges(6,2) * t99 + Icges(6,6) * t137;
t69 = Icges(6,1) * t102 + Icges(6,4) * t99 + Icges(6,5) * t137;
t101 = t116 + t121 / 0.2e1;
t128 = sin(t134);
t75 = -t123 * t101 - t124 * t128;
t100 = t115 - t120 / 0.2e1;
t130 = cos(t134);
t76 = -t123 * t100 + t124 * t130;
t18 = t67 * t160 + t75 * t68 + t76 * t69;
t153 = t102 * t69 + t137 * t67 + t99 * t68;
t20 = t153 * t137;
t73 = t124 * t101 - t123 * t128;
t74 = t124 * t100 + t123 * t130;
t44 = Icges(6,5) * t74 + Icges(6,6) * t73 - Icges(6,3) * t158;
t45 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t160;
t46 = Icges(6,4) * t74 + Icges(6,2) * t73 - Icges(6,6) * t158;
t47 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t160;
t48 = Icges(6,1) * t74 + Icges(6,4) * t73 - Icges(6,5) * t158;
t49 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t160;
t7 = t102 * t48 + t137 * t44 + t99 * t46;
t8 = t102 * t49 + t137 * t45 + t99 * t47;
t164 = ((t45 * t160 + t75 * t47 + t76 * t49) * t160 - (t44 * t160 + t75 * t46 + t76 * t48) * t158 + t18 * t137) * t160 + t137 * (t20 + (t123 * t8 - t124 * t7) * t136);
t129 = sin(t135);
t163 = pkin(2) * t129;
t139 = sin(qJ(1));
t162 = t139 * pkin(1);
t50 = t74 * rSges(6,1) + t73 * rSges(6,2) - rSges(6,3) * t158;
t51 = t76 * rSges(6,1) + t75 * rSges(6,2) + rSges(6,3) * t160;
t19 = t51 * t158 + t50 * t160;
t161 = t123 * t124;
t142 = pkin(9) + pkin(10);
t103 = pkin(4) * t156 - t136 * t142;
t159 = t124 * t103;
t157 = t136 * t138;
t154 = t124 * pkin(3) + pkin(9) * t160;
t152 = t136 * t140 * t91 + t137 * t90 + t92 * t157;
t61 = t89 * rSges(5,1) + t88 * rSges(5,2) + rSges(5,3) * t160;
t151 = t160 / 0.2e1;
t150 = -t158 / 0.2e1;
t70 = t102 * rSges(6,1) + t99 * rSges(6,2) + t137 * rSges(6,3);
t149 = (-t70 - pkin(4) * t157 - (-pkin(9) + t142) * t137) * t136;
t131 = cos(t135);
t105 = t131 * rSges(3,1) - t129 * rSges(3,2);
t98 = t124 * rSges(4,1) - t123 * rSges(4,2);
t125 = t140 * pkin(4) + pkin(3);
t148 = -t123 * t103 + t124 * t125;
t17 = -t67 * t158 + t73 * t68 + t74 * t69;
t147 = t20 + (t18 + t8) * t151 + (t17 + t7) * t150;
t122 = pkin(2) * t131;
t85 = t122 + t98;
t53 = t61 + t154;
t2 = (-t45 * t158 + t73 * t47 + t74 * t49) * t160 - (-t44 * t158 + t73 * t46 + t74 * t48) * t158 + t17 * t137;
t146 = -t2 * t158 + t164;
t104 = -t129 * rSges(3,1) - t131 * rSges(3,2);
t97 = -t123 * rSges(4,1) - t124 * rSges(4,2);
t60 = t87 * rSges(5,1) + t86 * rSges(5,2) - rSges(5,3) * t158;
t43 = t122 + t53;
t145 = Icges(4,3) + t152 + t153;
t33 = t148 + t51;
t144 = Icges(3,3) + t145;
t84 = t97 - t163;
t29 = t122 + t33;
t41 = t152 * t137;
t143 = t166 * t150 + t165 * t151 + t147 + t41;
t112 = pkin(9) * t158;
t52 = -t123 * pkin(3) + t112 - t60;
t32 = -t123 * t125 - t159 - t50;
t42 = t52 - t163;
t28 = t32 - t163;
t141 = cos(qJ(1));
t133 = t141 * pkin(1);
t114 = t141 * rSges(2,1) - t139 * rSges(2,2);
t113 = -t139 * rSges(2,1) - t141 * rSges(2,2);
t96 = t105 + t133;
t95 = t104 - t162;
t93 = t137 * rSges(5,3) + (rSges(5,1) * t138 + rSges(5,2) * t140) * t136;
t80 = t133 + t85;
t79 = t84 - t162;
t65 = t148 - t154;
t64 = t159 + t112 + (-pkin(3) + t125) * t123;
t40 = t137 * t51;
t39 = t133 + t43;
t38 = t42 - t162;
t35 = -t137 * t60 - t93 * t158;
t34 = t137 * t61 - t93 * t160;
t27 = t133 + t29;
t26 = t28 - t162;
t23 = -t137 * t50 - t70 * t158;
t22 = -t70 * t160 + t40;
t21 = (t123 * t60 + t124 * t61) * t136;
t10 = (-t50 - t64) * t137 + t124 * t149;
t9 = t123 * t149 + t137 * t65 + t40;
t4 = (t123 * t64 + t124 * t65) * t136 + t19;
t1 = [Icges(2,3) + m(6) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t38 ^ 2 + t39 ^ 2) + m(4) * (t79 ^ 2 + t80 ^ 2) + m(3) * (t95 ^ 2 + t96 ^ 2) + m(2) * (t113 ^ 2 + t114 ^ 2) + t144; m(6) * (t28 * t26 + t27 * t29) + m(5) * (t42 * t38 + t43 * t39) + m(4) * (t84 * t79 + t85 * t80) + m(3) * (t104 * t95 + t105 * t96) + t144; m(6) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t84 ^ 2 + t85 ^ 2) + m(3) * (t104 ^ 2 + t105 ^ 2) + t144; m(6) * (t32 * t26 + t33 * t27) + m(5) * (t52 * t38 + t53 * t39) + m(4) * (t97 * t79 + t98 * t80) + t145; m(6) * (t32 * t28 + t33 * t29) + m(5) * (t52 * t42 + t53 * t43) + m(4) * (t97 * t84 + t98 * t85) + t145; m(6) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t97 ^ 2 + t98 ^ 2) + t145; m(6) * (t10 * t26 + t27 * t9) + m(5) * (t34 * t39 + t35 * t38) + t143; m(6) * (t10 * t28 + t29 * t9) + m(5) * (t34 * t43 + t35 * t42) + t143; m(6) * (t10 * t32 + t9 * t33) + m(5) * (t34 * t53 + t35 * t52) + t143; t137 * t41 + (-t124 * t2 + (t123 * ((t88 * t57 + t89 * t59) * t123 - (t88 * t56 + t89 * t58) * t124) - t124 * ((t86 * t57 + t87 * t59) * t123 - (t86 * t56 + t87 * t58) * t124) + (t123 * (t123 ^ 2 * t55 - t54 * t161) - t124 * (t124 ^ 2 * t54 - t55 * t161)) * t136) * t136 + (t165 * t123 - t166 * t124) * t137) * t136 + m(5) * (t21 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t10 ^ 2 + t4 ^ 2 + t9 ^ 2) + t164; m(6) * (t22 * t27 + t23 * t26) + t147; m(6) * (t22 * t29 + t23 * t28) + t147; m(6) * (t22 * t33 + t23 * t32) + t147; m(6) * (t10 * t23 + t19 * t4 + t22 * t9) + t146; m(6) * (t19 ^ 2 + t22 ^ 2 + t23 ^ 2) + t146;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
