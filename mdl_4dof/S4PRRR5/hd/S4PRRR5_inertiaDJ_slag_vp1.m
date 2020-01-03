% Calculate time derivative of joint inertia matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:39
% DurationCPUTime: 2.30s
% Computational Cost: add. (7276->252), mult. (10085->417), div. (0->0), fcn. (9487->8), ass. (0->153)
t135 = qJD(2) + qJD(3);
t137 = sin(pkin(7));
t133 = t137 ^ 2;
t138 = cos(pkin(7));
t134 = t138 ^ 2;
t234 = t133 + t134;
t229 = t135 * t234;
t140 = sin(qJ(2));
t142 = cos(qJ(2));
t235 = t137 * qJD(2) * (-Icges(3,5) * t140 - Icges(3,6) * t142);
t136 = qJ(2) + qJ(3);
t131 = sin(t136);
t132 = cos(t136);
t233 = -Icges(4,5) * t131 - Icges(4,6) * t132;
t223 = t233 * t229;
t141 = cos(qJ(4));
t193 = t138 * t141;
t139 = sin(qJ(4));
t196 = t137 * t139;
t122 = -t132 * t196 - t193;
t194 = t138 * t139;
t195 = t137 * t141;
t123 = t132 * t195 - t194;
t199 = t132 * t135;
t198 = t135 * t137;
t187 = t132 * t198;
t202 = t131 * t135;
t189 = t139 * t202;
t78 = -qJD(4) * t123 + t137 * t189;
t188 = t141 * t202;
t79 = qJD(4) * t122 - t137 * t188;
t46 = Icges(5,5) * t79 + Icges(5,6) * t78 + Icges(5,3) * t187;
t201 = t131 * t137;
t67 = Icges(5,5) * t123 + Icges(5,6) * t122 + Icges(5,3) * t201;
t155 = t131 * t46 + t199 * t67;
t48 = Icges(5,4) * t79 + Icges(5,2) * t78 + Icges(5,6) * t187;
t50 = Icges(5,1) * t79 + Icges(5,4) * t78 + Icges(5,5) * t187;
t69 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t201;
t71 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t201;
t13 = t122 * t48 + t123 * t50 + t137 * t155 + t69 * t78 + t71 * t79;
t197 = t135 * t138;
t186 = t132 * t197;
t125 = t132 * t193 + t196;
t80 = -qJD(4) * t125 + t138 * t189;
t124 = -t132 * t194 + t195;
t81 = qJD(4) * t124 - t138 * t188;
t47 = Icges(5,5) * t81 + Icges(5,6) * t80 + Icges(5,3) * t186;
t200 = t131 * t138;
t68 = Icges(5,5) * t125 + Icges(5,6) * t124 + Icges(5,3) * t200;
t154 = t131 * t47 + t199 * t68;
t49 = Icges(5,4) * t81 + Icges(5,2) * t80 + Icges(5,6) * t186;
t51 = Icges(5,1) * t81 + Icges(5,4) * t80 + Icges(5,5) * t186;
t70 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t200;
t72 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t200;
t14 = t122 * t49 + t123 * t51 + t137 * t154 + t70 * t78 + t72 * t79;
t9 = -t13 * t138 + t137 * t14;
t148 = t135 * t233;
t92 = t137 * t148;
t93 = t138 * t148;
t232 = -t134 * t92 - (-t138 * t93 + t223) * t137 - t9;
t231 = t138 * t235;
t226 = qJD(2) * (rSges(3,1) * t140 + rSges(3,2) * t142);
t222 = 2 * m(4);
t221 = 2 * m(5);
t15 = t124 * t48 + t125 * t50 + t138 * t155 + t69 * t80 + t71 * t81;
t16 = t124 * t49 + t125 * t51 + t138 * t154 + t70 * t80 + t72 * t81;
t10 = t137 * t16 - t138 * t15;
t220 = (t133 * t93 + t10 + (-t137 * t92 + t223) * t138) * t137;
t219 = pkin(2) * t140;
t126 = rSges(4,1) * t131 + rSges(4,2) * t132;
t57 = t126 * t229;
t216 = t234 * pkin(2) * t142;
t177 = rSges(4,1) * t132 - rSges(4,2) * t131;
t62 = t234 * t177;
t215 = pkin(2) * qJD(2);
t162 = Icges(5,5) * t141 - Icges(5,6) * t139;
t214 = t132 * (t162 * t199 + (Icges(5,3) * t135 + (-Icges(5,5) * t139 - Icges(5,6) * t141) * qJD(4)) * t131);
t98 = -Icges(5,3) * t132 + t131 * t162;
t213 = t132 * t98;
t30 = t124 * t69 + t125 * t71 + t200 * t67;
t212 = t137 * t30;
t29 = t122 * t70 + t123 * t72 + t201 * t68;
t211 = t138 * t29;
t179 = pkin(3) * t132 + pkin(6) * t131;
t176 = rSges(5,1) * t141 - rSges(5,2) * t139;
t61 = t176 * t199 + (rSges(5,3) * t135 + (-rSges(5,1) * t139 - rSges(5,2) * t141) * qJD(4)) * t131;
t210 = -t135 * t179 - t61;
t205 = Icges(5,4) * t139;
t204 = Icges(5,4) * t141;
t101 = -t132 * rSges(5,3) + t131 * t176;
t127 = pkin(3) * t131 - pkin(6) * t132;
t192 = -t101 - t127;
t190 = t142 * t215;
t184 = -t126 - t219;
t52 = rSges(5,1) * t79 + rSges(5,2) * t78 + rSges(5,3) * t187;
t53 = rSges(5,1) * t81 + rSges(5,2) * t80 + rSges(5,3) * t186;
t26 = -t127 * t229 + t137 * t52 + t138 * t53;
t73 = rSges(5,1) * t123 + rSges(5,2) * t122 + rSges(5,3) * t201;
t74 = rSges(5,1) * t125 + rSges(5,2) * t124 + rSges(5,3) * t200;
t34 = t137 * t73 + t138 * t74 + t179 * t234;
t181 = t192 - t219;
t114 = t177 * t135;
t180 = -t114 - t190;
t173 = -t139 * t69 + t141 * t71;
t32 = t131 * t173 - t132 * t67;
t172 = -t139 * t70 + t141 * t72;
t33 = t131 * t172 - t132 * t68;
t175 = t32 * t137 + t33 * t138;
t174 = -t137 * t74 + t138 * t73;
t166 = Icges(5,1) * t141 - t205;
t100 = -Icges(5,5) * t132 + t131 * t166;
t163 = -Icges(5,2) * t139 + t204;
t99 = -Icges(5,6) * t132 + t131 * t163;
t169 = t100 * t141 - t139 * t99;
t158 = t138 * t232 + t220;
t157 = -t190 + t210;
t156 = t234 * t140 * t215;
t11 = (t135 * t173 - t46) * t132 + (t135 * t67 - t139 * t48 + t141 * t50 + (-t139 * t71 - t141 * t69) * qJD(4)) * t131;
t12 = (t135 * t172 - t47) * t132 + (t135 * t68 - t139 * t49 + t141 * t51 + (-t139 * t72 - t141 * t70) * qJD(4)) * t131;
t28 = t122 * t69 + t123 * t71 + t201 * t67;
t35 = t100 * t123 + t122 * t99 + t201 * t98;
t59 = t163 * t199 + (Icges(5,6) * t135 + (-Icges(5,2) * t141 - t205) * qJD(4)) * t131;
t60 = t166 * t199 + (Icges(5,5) * t135 + (-Icges(5,1) * t139 - t204) * qJD(4)) * t131;
t3 = (-t79 * t100 - t122 * t59 - t123 * t60 - t78 * t99 + (t211 + (t28 - t213) * t137) * t135) * t132 + (t35 * t135 + t14 * t138 + (t13 - t214) * t137) * t131;
t31 = t124 * t70 + t125 * t72 + t200 * t68;
t36 = t100 * t125 + t124 * t99 + t200 * t98;
t4 = (-t81 * t100 - t124 * t59 - t125 * t60 - t80 * t99 + (t212 + (t31 - t213) * t138) * t135) * t132 + (t36 * t135 + t15 * t137 + (t16 - t214) * t138) * t131;
t152 = t137 * t4 / 0.2e1 + (t137 * t33 - t138 * t32) * t202 / 0.2e1 - t138 * t3 / 0.2e1 - t132 * (-t11 * t138 + t12 * t137) / 0.2e1 + t9 * t201 / 0.2e1 + t10 * t200 / 0.2e1 + (t137 * (t137 * t29 - t138 * t28) + t138 * (t137 * t31 - t138 * t30)) * t199 / 0.2e1;
t103 = t184 * t138;
t102 = t184 * t137;
t83 = t180 * t138;
t82 = t180 * t137;
t77 = t234 * t226;
t76 = t192 * t138;
t75 = t192 * t137;
t64 = t181 * t138;
t63 = t181 * t137;
t56 = -t156 - t57;
t55 = t210 * t138;
t54 = t210 * t137;
t45 = -t101 * t200 - t132 * t74;
t44 = t101 * t201 + t132 * t73;
t41 = t157 * t138;
t40 = t157 * t137;
t39 = t131 * t169 - t213;
t38 = t62 + t216;
t37 = t174 * t131;
t27 = t34 + t216;
t25 = (-t101 * t197 - t53) * t132 + (t135 * t74 - t138 * t61) * t131;
t24 = (t101 * t198 + t52) * t132 + (-t135 * t73 + t137 * t61) * t131;
t23 = -t156 + t26;
t22 = t174 * t199 + (-t137 * t53 + t138 * t52) * t131;
t1 = [0; -m(3) * t77 + m(4) * t56 + m(5) * t23; t133 * t231 + (t102 * t82 + t103 * t83 + t38 * t56) * t222 + (t27 * t23 + t40 * t63 + t41 * t64) * t221 + t220 + 0.2e1 * m(3) * (t226 - t77) * t234 * (rSges(3,1) * t142 - rSges(3,2) * t140) + (t138 * t231 - t234 * t235 + t232) * t138; -m(4) * t57 + m(5) * t26; m(4) * (-t38 * t57 + t56 * t62 + (-t137 * t82 - t138 * t83) * t126 + (-t102 * t137 - t103 * t138) * t114) + m(5) * (t34 * t23 + t26 * t27 + t40 * t75 + t41 * t76 + t54 * t63 + t55 * t64) + t158; (t114 * t126 * t234 - t57 * t62) * t222 + (t34 * t26 + t54 * t75 + t55 * t76) * t221 + t158; m(5) * t22; m(5) * (t22 * t27 + t37 * t23 + t24 * t64 + t25 * t63 + t45 * t40 + t44 * t41) + t152; m(5) * (t22 * t34 + t24 * t76 + t25 * t75 + t37 * t26 + t44 * t55 + t45 * t54) + t152; (t22 * t37 + t24 * t44 + t25 * t45) * t221 + (-t36 * t132 + (t138 * t31 + t212) * t131) * t186 + t4 * t200 + (-t35 * t132 + (t137 * t28 + t211) * t131) * t187 + t3 * t201 + (t131 * t175 - t39 * t132) * t202 - t132 * ((t214 + (-t132 * t169 + t175) * t135) * t132 + (t12 * t138 + t11 * t137 - (t135 * t98 - t139 * t59 + t141 * t60 + (-t100 * t139 - t141 * t99) * qJD(4)) * t132 + t39 * t135) * t131);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
