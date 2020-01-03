% Calculate time derivative of joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:25
% DurationCPUTime: 1.70s
% Computational Cost: add. (5693->277), mult. (4836->401), div. (0->0), fcn. (4330->10), ass. (0->135)
t139 = qJ(1) + qJ(2);
t134 = sin(t139);
t135 = cos(t139);
t99 = t134 * rSges(3,1) + t135 * rSges(3,2);
t142 = sin(qJ(5));
t144 = cos(qJ(5));
t140 = sin(pkin(9));
t141 = cos(pkin(9));
t174 = Icges(6,4) * t144;
t83 = -Icges(6,6) * t141 + (-Icges(6,2) * t142 + t174) * t140;
t175 = Icges(6,4) * t142;
t84 = -Icges(6,5) * t141 + (Icges(6,1) * t144 - t175) * t140;
t186 = -t142 * t84 - t144 * t83;
t133 = pkin(8) + t139;
t128 = sin(t133);
t129 = cos(t133);
t130 = pkin(2) * t134;
t76 = t128 * rSges(4,1) + t129 * rSges(4,2) + t130;
t185 = 2 * m(3);
t184 = 2 * m(4);
t183 = 2 * m(5);
t182 = 2 * m(6);
t181 = rSges(3,2) * t134;
t180 = pkin(1) * qJD(1);
t161 = qJD(5) * t140;
t89 = (-Icges(6,2) * t144 - t175) * t161;
t178 = t142 * t89;
t88 = (-Icges(6,5) * t142 - Icges(6,6) * t144) * t161;
t90 = (-Icges(6,1) * t142 - t174) * t161;
t151 = t140 * t144 * t90 - t141 * t88;
t26 = (t186 * qJD(5) - t178) * t140 + t151;
t176 = t26 * t141;
t173 = t128 * t140;
t172 = t128 * t141;
t171 = t129 * t140;
t170 = t129 * t141;
t138 = qJD(1) + qJD(2);
t169 = t135 * t138;
t168 = t138 * t128;
t167 = t138 * t129;
t166 = t138 * t140;
t165 = t141 * t142;
t164 = t141 * t144;
t163 = qJ(4) * t167 + qJD(4) * t128;
t162 = t128 * pkin(3) + t130;
t158 = t129 * t166;
t79 = t128 * t164 - t129 * t142;
t80 = -t128 * t144 + t129 * t165;
t61 = -t79 * qJD(5) - t80 * t138;
t147 = t128 * t142 + t129 * t164;
t78 = -t128 * t165 - t129 * t144;
t62 = t78 * qJD(5) + t147 * t138;
t34 = t62 * rSges(6,1) + t61 * rSges(6,2) + rSges(6,3) * t158;
t53 = t79 * rSges(6,1) + t78 * rSges(6,2) + rSges(6,3) * t173;
t143 = sin(qJ(1));
t160 = t143 * t180;
t159 = t128 * t166;
t157 = t141 * t167;
t116 = pkin(2) * t169;
t156 = pkin(3) * t167 + qJ(4) * t168 + t116;
t131 = pkin(2) * t135;
t155 = t129 * pkin(3) + t128 * qJ(4) + t131;
t100 = t135 * rSges(3,1) - t181;
t86 = rSges(3,1) * t169 - t138 * t181;
t77 = t129 * rSges(4,1) - rSges(4,2) * t128 + t131;
t59 = t147 * qJD(5) + t78 * t138;
t60 = t80 * qJD(5) + t79 * t138;
t148 = -rSges(6,1) * t60 - t59 * rSges(6,2);
t68 = rSges(4,1) * t167 - rSges(4,2) * t168 + t116;
t54 = -rSges(6,1) * t147 + rSges(6,2) * t80 - rSges(6,3) * t171;
t85 = t99 * t138;
t64 = rSges(5,1) * t170 - rSges(5,2) * t171 + t128 * rSges(5,3) + t155;
t67 = t76 * t138;
t47 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t173;
t49 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t173;
t51 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t173;
t18 = -t141 * t47 + (-t142 * t49 + t144 * t51) * t140;
t48 = -Icges(6,5) * t147 + Icges(6,6) * t80 - Icges(6,3) * t171;
t50 = -Icges(6,4) * t147 + Icges(6,2) * t80 - Icges(6,6) * t171;
t52 = -Icges(6,1) * t147 + Icges(6,4) * t80 - Icges(6,5) * t171;
t19 = -t141 * t48 + (-t142 * t50 + t144 * t52) * t140;
t28 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t158;
t30 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t158;
t32 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t158;
t3 = -t141 * t28 + (-t142 * t30 + t144 * t32 + (-t142 * t51 - t144 * t49) * qJD(5)) * t140;
t82 = -Icges(6,3) * t141 + (Icges(6,5) * t144 - Icges(6,6) * t142) * t140;
t35 = t82 * t173 + t78 * t83 + t79 * t84;
t36 = -t147 * t84 - t82 * t171 + t80 * t83;
t27 = Icges(6,5) * t60 + Icges(6,6) * t59 + Icges(6,3) * t159;
t29 = Icges(6,4) * t60 + Icges(6,2) * t59 + Icges(6,6) * t159;
t31 = Icges(6,1) * t60 + Icges(6,4) * t59 + Icges(6,5) * t159;
t4 = -t141 * t27 + (-t142 * t29 + t144 * t31 + (-t142 * t52 - t144 * t50) * qJD(5)) * t140;
t8 = t59 * t83 + t60 * t84 + t80 * t89 - t147 * t90 + (-t129 * t88 + t82 * t168) * t140;
t9 = t61 * t83 + t62 * t84 + t78 * t89 + t79 * t90 + (t128 * t88 + t82 * t167) * t140;
t146 = -t176 + (t3 + t9) * t173 / 0.2e1 - (t4 + t8) * t171 / 0.2e1 + ((t19 + t36) * t128 + (t18 + t35) * t129) * t166 / 0.2e1;
t39 = pkin(4) * t172 + pkin(7) * t173 - qJ(4) * t129 + t162 + t53;
t21 = pkin(4) * t157 + pkin(7) * t158 - qJD(4) * t129 + t156 + t34;
t63 = -rSges(5,2) * t173 + rSges(5,1) * t172 + (-rSges(5,3) - qJ(4)) * t129 + t162;
t40 = pkin(4) * t170 + pkin(7) * t171 + t155 - t54;
t46 = rSges(5,3) * t168 + rSges(5,1) * t157 + (-rSges(5,2) * t166 - qJD(4)) * t129 + t156;
t45 = rSges(5,3) * t167 + rSges(5,2) * t159 + (-t130 + (-rSges(5,1) * t141 - pkin(3)) * t128) * t138 + t163;
t20 = (-t130 + (-pkin(4) * t141 - pkin(3) + (-rSges(6,3) - pkin(7)) * t140) * t128) * t138 + t148 + t163;
t145 = cos(qJ(1));
t137 = t145 * pkin(1);
t136 = t143 * pkin(1);
t132 = t145 * t180;
t93 = t100 + t137;
t92 = t136 + t99;
t91 = (-rSges(6,1) * t142 - rSges(6,2) * t144) * t161;
t87 = -rSges(6,3) * t141 + (rSges(6,1) * t144 - rSges(6,2) * t142) * t140;
t72 = t132 + t86;
t71 = -t85 - t160;
t70 = t137 + t77;
t69 = t136 + t76;
t66 = t132 + t68;
t65 = -t67 - t160;
t56 = t137 + t64;
t55 = t136 + t63;
t44 = t132 + t46;
t43 = t45 - t160;
t42 = t141 * t54 - t87 * t171;
t41 = -t141 * t53 - t87 * t173;
t38 = t137 + t40;
t37 = t136 + t39;
t33 = rSges(6,3) * t159 - t148;
t23 = -t141 * t34 + (-t128 * t91 - t87 * t167) * t140;
t22 = t141 * t33 + (-t129 * t91 + t87 * t168) * t140;
t17 = t132 + t21;
t16 = t20 - t160;
t13 = -t147 * t52 - t48 * t171 + t50 * t80;
t12 = -t147 * t51 - t47 * t171 + t49 * t80;
t11 = t48 * t173 + t50 * t78 + t52 * t79;
t10 = t47 * t173 + t49 * t78 + t51 * t79;
t5 = ((t138 * t54 + t34) * t129 + (-t138 * t53 + t33) * t128) * t140;
t1 = [(t71 * t93 + t72 * t92) * t185 + (t65 * t70 + t66 * t69) * t184 + (t43 * t56 + t44 * t55) * t183 - t140 * t178 + (t16 * t38 + t17 * t37) * t182 + t151 + t186 * t161; m(3) * (t100 * t71 + t72 * t99 - t85 * t93 + t86 * t92) + m(4) * (t65 * t77 + t66 * t76 - t67 * t70 + t68 * t69) + m(5) * (t43 * t64 + t44 * t63 + t45 * t56 + t46 * t55) + m(6) * (t16 * t40 + t17 * t39 + t20 * t38 + t21 * t37) + t26; (t20 * t40 + t21 * t39) * t182 + (t45 * t64 + t46 * t63) * t183 + (-t67 * t77 + t68 * t76) * t184 + (-t100 * t85 + t86 * t99) * t185 + t26; 0; 0; 0; m(5) * ((-t138 * t55 - t43) * t129 + (t138 * t56 - t44) * t128) + m(6) * ((-t138 * t37 - t16) * t129 + (t138 * t38 - t17) * t128); m(6) * ((-t138 * t39 - t20) * t129 + (t138 * t40 - t21) * t128) + m(5) * ((-t138 * t63 - t45) * t129 + (t138 * t64 - t46) * t128); 0; 0; m(6) * (t16 * t42 + t17 * t41 + t22 * t38 + t23 * t37) + t146; m(6) * (t20 * t42 + t21 * t41 + t22 * t40 + t23 * t39) + t146; m(6) * t5; m(6) * ((-t138 * t41 - t22) * t129 + (t138 * t42 - t23) * t128); (t42 * t22 + t41 * t23 + (t128 * t54 + t129 * t53) * t5 * t140) * t182 - t141 * (-t176 + ((t138 * t18 - t4) * t129 + (t138 * t19 + t3) * t128) * t140) + (-t141 * t35 + (t10 * t128 - t11 * t129) * t140) * t158 + (-t9 * t141 + ((t78 * t30 + t79 * t32 + t61 * t49 + t62 * t51) * t128 + t10 * t167 - (t78 * t29 + t79 * t31 + t61 * t50 + t62 * t52) * t129 + t11 * t168 + ((t128 * t28 + t47 * t167) * t128 - (t128 * t27 + t48 * t167) * t129) * t140) * t140) * t173 + (-t141 * t36 + (t12 * t128 - t129 * t13) * t140) * t159 - (-t8 * t141 + ((-t147 * t32 + t80 * t30 + t59 * t49 + t60 * t51) * t128 + t12 * t167 - (-t147 * t31 + t80 * t29 + t59 * t50 + t60 * t52) * t129 + t13 * t168 + ((-t129 * t28 + t47 * t168) * t128 - (-t129 * t27 + t48 * t168) * t129) * t140) * t140) * t171;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
