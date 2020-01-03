% Calculate time derivative of joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:25
% DurationCPUTime: 1.71s
% Computational Cost: add. (4960->208), mult. (2858->291), div. (0->0), fcn. (1986->10), ass. (0->127)
t124 = sin(qJ(5));
t155 = qJD(5) * t124;
t126 = cos(qJ(5));
t154 = qJD(5) * t126;
t122 = qJD(1) + qJD(2);
t116 = qJD(4) + t122;
t138 = Icges(6,5) * t126 - Icges(6,6) * t124;
t167 = Icges(6,4) * t126;
t139 = -Icges(6,2) * t124 + t167;
t168 = Icges(6,4) * t124;
t140 = Icges(6,1) * t126 - t168;
t88 = Icges(6,2) * t126 + t168;
t89 = Icges(6,1) * t124 + t167;
t141 = t124 * t88 - t126 * t89;
t190 = -t138 * qJD(5) + (t124 * t140 + t126 * t139 - t141) * t116;
t123 = qJ(1) + qJ(2);
t117 = pkin(9) + t123;
t114 = qJ(4) + t117;
t105 = sin(t114);
t106 = cos(t114);
t43 = -Icges(6,6) * t105 - t106 * t139;
t45 = -Icges(6,5) * t105 - t140 * t106;
t143 = t124 * t43 - t126 * t45;
t187 = t143 * t105;
t110 = sin(t117);
t118 = sin(t123);
t112 = pkin(2) * t118;
t157 = pkin(3) * t110 + t112;
t119 = cos(t123);
t77 = t118 * rSges(3,1) + t119 * rSges(3,2);
t111 = cos(t117);
t58 = t110 * rSges(4,1) + t111 * rSges(4,2) + t112;
t87 = Icges(6,5) * t124 + Icges(6,6) * t126;
t186 = -Icges(6,3) * t116 + qJD(5) * t87;
t182 = 2 * m(3);
t181 = 2 * m(4);
t180 = 2 * m(5);
t179 = 2 * m(6);
t170 = rSges(6,2) * t124;
t153 = t105 * t170;
t161 = t106 * t116;
t176 = rSges(6,3) * t161 + t116 * t153;
t163 = t105 * t116;
t173 = rSges(6,1) * t126;
t84 = t106 * t173;
t175 = rSges(6,3) * t163 + t116 * t84;
t159 = t111 * t122;
t158 = t119 * t122;
t91 = pkin(2) * t158;
t174 = pkin(3) * t159 + t91;
t62 = t105 * rSges(5,1) + t106 * rSges(5,2);
t172 = rSges(3,2) * t118;
t171 = rSges(4,2) * t110;
t169 = pkin(1) * qJD(1);
t162 = t105 * t126;
t160 = t106 * t124;
t113 = pkin(2) * t119;
t156 = pkin(3) * t111 + t113;
t125 = sin(qJ(1));
t152 = t125 * t169;
t63 = t106 * rSges(5,1) - rSges(5,2) * t105;
t78 = t119 * rSges(3,1) - t172;
t50 = t157 + t62;
t61 = rSges(3,1) * t158 - t122 * t172;
t53 = rSges(5,1) * t161 - rSges(5,2) * t163;
t149 = rSges(6,1) * t162 - t153;
t59 = t111 * rSges(4,1) + t113 - t171;
t92 = rSges(6,1) * t124 + rSges(6,2) * t126;
t42 = -Icges(6,6) * t106 + t105 * t139;
t44 = -Icges(6,5) * t106 + t140 * t105;
t146 = t124 * t44 + t126 * t42;
t145 = t124 * t42 - t126 * t44;
t144 = t124 * t45 + t126 * t43;
t49 = rSges(4,1) * t159 - t122 * t171 + t91;
t47 = rSges(6,2) * t160 - t105 * rSges(6,3) - t84;
t51 = t63 + t156;
t35 = t53 + t174;
t137 = t157 * t122;
t136 = (t140 - t88) * t155 + (t139 + t89) * t154;
t135 = -(-t143 * qJD(5) + t190 * t105) * t105 / 0.2e1 - (-t145 * qJD(5) + t190 * t106) * t106 / 0.2e1 + (-t141 * t105 - t106 * t87 + t146) * t163 / 0.2e1 - (-t105 * t87 + t141 * t106 + t144) * t161 / 0.2e1;
t60 = t77 * t122;
t52 = t62 * t116;
t134 = t145 * t106;
t133 = qJD(5) * t92;
t130 = t138 * t116;
t33 = t106 * pkin(4) + t105 * pkin(8) - t47;
t48 = t58 * t122;
t32 = t105 * pkin(4) + (-rSges(6,3) - pkin(8)) * t106 + t149;
t29 = t33 + t156;
t28 = t32 + t157;
t34 = -t137 - t52;
t15 = -t105 * rSges(6,1) * t155 + pkin(8) * t163 + pkin(4) * t161 + (-t105 * t154 - t116 * t160) * rSges(6,2) + t175;
t13 = t15 + t174;
t14 = -rSges(6,2) * t106 * t154 - pkin(4) * t163 + pkin(8) * t161 + (-t106 * t155 - t116 * t162) * rSges(6,1) + t176;
t12 = -t137 + t14;
t127 = cos(qJ(1));
t121 = t127 * pkin(1);
t120 = t125 * pkin(1);
t115 = t127 * t169;
t82 = (-t170 + t173) * qJD(5);
t65 = t121 + t78;
t64 = t120 + t77;
t57 = t115 + t61;
t56 = -t60 - t152;
t55 = t121 + t59;
t54 = t120 + t58;
t46 = -rSges(6,3) * t106 + t149;
t41 = -Icges(6,3) * t105 - t106 * t138;
t40 = -Icges(6,3) * t106 + t105 * t138;
t39 = t121 + t51;
t38 = t120 + t50;
t37 = t115 + t49;
t36 = -t48 - t152;
t31 = t115 + t35;
t30 = t34 - t152;
t27 = t121 + t29;
t26 = t120 + t28;
t19 = -t186 * t105 + t106 * t130;
t18 = t105 * t130 + t186 * t106;
t11 = t115 + t13;
t10 = t12 - t152;
t9 = -t105 * t41 + t143 * t106;
t8 = -t105 * t40 + t134;
t7 = -t106 * t41 - t187;
t6 = -t145 * t105 - t106 * t40;
t1 = (-t106 * t133 + t116 * t46 + t176) * t106 + (-t105 * t133 + (t47 + (-t170 - t173) * t106) * t116 + t175) * t105;
t2 = [(t56 * t65 + t57 * t64) * t182 + (t36 * t55 + t37 * t54) * t181 + (t30 * t39 + t31 * t38) * t180 + (t10 * t27 + t11 * t26) * t179 + t136; m(3) * (t56 * t78 + t57 * t77 - t60 * t65 + t61 * t64) + m(4) * (t36 * t59 + t37 * t58 - t48 * t55 + t49 * t54) + m(5) * (t30 * t51 + t31 * t50 + t34 * t39 + t35 * t38) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + t136; (t12 * t29 + t13 * t28) * t179 + (t34 * t51 + t35 * t50) * t180 + (-t48 * t59 + t49 * t58) * t181 + (-t60 * t78 + t61 * t77) * t182 + t136; 0; 0; 0; m(5) * (t30 * t63 + t31 * t62 + t38 * t53 - t39 * t52) + m(6) * (t10 * t33 + t11 * t32 + t14 * t27 + t15 * t26) + t136; m(6) * (t12 * t33 + t13 * t32 + t14 * t29 + t15 * t28) + m(5) * (t34 * t63 + t35 * t62 + t50 * t53 - t51 * t52) + t136; 0; (-t52 * t63 + t53 * t62) * t180 + (t14 * t33 + t15 * t32) * t179 + t136; m(6) * ((-t105 * t27 + t106 * t26) * t82 + ((-t116 * t27 + t11) * t106 + (-t116 * t26 - t10) * t105) * t92) + t135; m(6) * ((-t105 * t29 + t106 * t28) * t82 + ((-t116 * t29 + t13) * t106 + (-t116 * t28 - t12) * t105) * t92) + t135; m(6) * t1; m(6) * ((-t105 * t33 + t106 * t32) * t82 + ((-t116 * t33 + t15) * t106 + (-t116 * t32 - t14) * t105) * t92) + t135; ((t105 * t46 - t106 * t47) * t1 + (t105 ^ 2 + t106 ^ 2) * t92 * t82) * t179 + (-t105 * t7 - t106 * t6) * t163 - t106 * ((t106 * t19 + (-t7 + t134) * t116) * t106 + (t6 * t116 + (t43 * t154 + t45 * t155) * t105 + (t146 * qJD(5) + t143 * t116 + t18) * t106) * t105) - (-t105 * t9 - t106 * t8) * t161 - t105 * ((t105 * t18 + (t8 + t187) * t116) * t105 + (-t9 * t116 + (-t42 * t154 - t44 * t155) * t106 + (-t144 * qJD(5) + t145 * t116 + t19) * t105) * t106);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
