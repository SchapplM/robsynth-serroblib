% Calculate time derivative of joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:06
% EndTime: 2020-01-03 11:20:13
% DurationCPUTime: 1.68s
% Computational Cost: add. (3787->245), mult. (3828->388), div. (0->0), fcn. (3531->10), ass. (0->118)
t87 = pkin(9) + qJ(5);
t83 = cos(t87);
t129 = Icges(6,4) * t83;
t81 = sin(t87);
t90 = sin(pkin(8));
t92 = cos(pkin(8));
t57 = -Icges(6,6) * t92 + (-Icges(6,2) * t81 + t129) * t90;
t130 = Icges(6,4) * t81;
t58 = -Icges(6,5) * t92 + (Icges(6,1) * t83 - t130) * t90;
t146 = -t57 * t83 - t58 * t81;
t91 = cos(pkin(9));
t139 = (pkin(4) * t91 + pkin(3)) * t92;
t93 = -pkin(6) - qJ(4);
t145 = (rSges(6,3) - t93) * t90 + t139;
t89 = sin(pkin(9));
t97 = (rSges(5,3) + qJ(4)) * t90 + (rSges(5,1) * t91 - rSges(5,2) * t89 + pkin(3)) * t92;
t144 = 2 * m(6);
t143 = pkin(4) * t89;
t85 = sin(qJ(1)) * pkin(1);
t86 = cos(qJ(1)) * pkin(1);
t120 = qJD(5) * t90;
t64 = (-Icges(6,5) * t81 - Icges(6,6) * t83) * t120;
t66 = (-Icges(6,1) * t81 - t129) * t120;
t110 = t90 * t83 * t66 - t92 * t64;
t65 = (-Icges(6,2) * t83 - t130) * t120;
t138 = t81 * t65;
t142 = ((t146 * qJD(5) - t138) * t90 + t110) * t92;
t88 = qJ(1) + pkin(7);
t82 = sin(t88);
t137 = t82 * t90;
t136 = t82 * t92;
t84 = cos(t88);
t135 = t84 * t90;
t134 = t84 * t92;
t119 = qJ(3) * qJD(1);
t132 = qJD(3) * t82 + t84 * t119;
t131 = t82 * pkin(2) + t85;
t128 = Icges(6,5) * t90;
t127 = Icges(6,6) * t90;
t126 = Icges(6,3) * t90;
t124 = qJD(1) * t82;
t123 = qJD(1) * t84;
t122 = qJD(1) * t90;
t121 = qJD(4) * t90;
t118 = t82 * t143;
t114 = t84 * t122;
t61 = t83 * t136 - t81 * t84;
t62 = t81 * t134 - t82 * t83;
t49 = -t62 * qJD(1) - t61 * qJD(5);
t100 = t83 * t134 + t81 * t82;
t60 = -t81 * t136 - t83 * t84;
t50 = t100 * qJD(1) + t60 * qJD(5);
t26 = t50 * rSges(6,1) + t49 * rSges(6,2) + rSges(6,3) * t114;
t41 = t61 * rSges(6,1) + t60 * rSges(6,2) + rSges(6,3) * t137;
t117 = t84 * t121 + t132;
t116 = t84 * pkin(2) + t82 * qJ(3) + t86;
t115 = t82 * t122;
t113 = Icges(6,5) * t122;
t112 = Icges(6,6) * t122;
t111 = Icges(6,3) * t122;
t107 = rSges(4,1) * t92 - rSges(4,2) * t90;
t106 = t89 * rSges(5,1) + t91 * rSges(5,2);
t47 = t60 * qJD(1) + t100 * qJD(5);
t48 = t61 * qJD(1) + t62 * qJD(5);
t105 = -rSges(6,1) * t48 - rSges(6,2) * t47;
t104 = rSges(6,1) * t100 - rSges(6,2) * t62;
t103 = -t90 * t93 + t139;
t102 = pkin(2) * t123 + qJD(1) * t86 - qJD(3) * t84 + t82 * t119;
t99 = t82 * t121 + t102;
t98 = rSges(4,3) * t82 + t107 * t84;
t96 = t106 * t82 + t97 * t84;
t67 = (-rSges(6,1) * t81 - rSges(6,2) * t83) * t120;
t59 = -rSges(6,3) * t92 + (rSges(6,1) * t83 - rSges(6,2) * t81) * t90;
t56 = -Icges(6,3) * t92 + (Icges(6,5) * t83 - Icges(6,6) * t81) * t90;
t52 = t98 + t116;
t51 = (-rSges(4,3) - qJ(3)) * t84 + t107 * t82 + t131;
t44 = t98 * qJD(1) + t102;
t43 = (rSges(4,3) * t84 - t85 + (-pkin(2) - t107) * t82) * qJD(1) + t132;
t42 = -rSges(6,3) * t135 - t104;
t40 = -Icges(6,1) * t100 + Icges(6,4) * t62 - t84 * t128;
t39 = Icges(6,1) * t61 + Icges(6,4) * t60 + t82 * t128;
t38 = -Icges(6,4) * t100 + Icges(6,2) * t62 - t84 * t127;
t37 = Icges(6,4) * t61 + Icges(6,2) * t60 + t82 * t127;
t36 = -Icges(6,5) * t100 + Icges(6,6) * t62 - t84 * t126;
t35 = Icges(6,5) * t61 + Icges(6,6) * t60 + t82 * t126;
t34 = t96 + t116;
t33 = (-qJ(3) - t106) * t84 + t97 * t82 + t131;
t32 = -t59 * t135 + t42 * t92;
t31 = -t59 * t137 - t41 * t92;
t30 = t145 * t84 + t104 + t116 + t118;
t29 = (-qJ(3) - t143) * t84 + t103 * t82 + t41 + t131;
t28 = t96 * qJD(1) + t99;
t27 = (-t85 + t106 * t84 + (-pkin(2) - t97) * t82) * qJD(1) + t117;
t25 = rSges(6,3) * t115 - t105;
t24 = Icges(6,1) * t50 + Icges(6,4) * t49 + t84 * t113;
t23 = Icges(6,1) * t48 + Icges(6,4) * t47 + t82 * t113;
t22 = Icges(6,4) * t50 + Icges(6,2) * t49 + t84 * t112;
t21 = Icges(6,4) * t48 + Icges(6,2) * t47 + t82 * t112;
t20 = Icges(6,5) * t50 + Icges(6,6) * t49 + t84 * t111;
t19 = Icges(6,5) * t48 + Icges(6,6) * t47 + t82 * t111;
t18 = -t100 * t58 - t56 * t135 + t57 * t62;
t17 = t56 * t137 + t57 * t60 + t58 * t61;
t15 = -t26 * t92 + (-t59 * t123 - t67 * t82) * t90;
t14 = t25 * t92 + (t59 * t124 - t67 * t84) * t90;
t13 = (t103 * t84 + t118) * qJD(1) + t99 + t26;
t12 = (t84 * t143 - t85 + (-pkin(2) - t145) * t82) * qJD(1) + t105 + t117;
t11 = -t36 * t92 + (-t38 * t81 + t40 * t83) * t90;
t10 = -t35 * t92 + (-t37 * t81 + t39 * t83) * t90;
t9 = -t100 * t40 - t36 * t135 + t38 * t62;
t8 = -t100 * t39 - t35 * t135 + t37 * t62;
t7 = t36 * t137 + t38 * t60 + t40 * t61;
t6 = t35 * t137 + t37 * t60 + t39 * t61;
t5 = t49 * t57 + t50 * t58 + t60 * t65 + t61 * t66 + (t56 * t123 + t64 * t82) * t90;
t4 = t47 * t57 + t48 * t58 + t62 * t65 - t100 * t66 + (t56 * t124 - t64 * t84) * t90;
t3 = (t25 * t82 + t26 * t84 + (-t41 * t82 + t42 * t84) * qJD(1)) * t90;
t2 = -t19 * t92 + (-t21 * t81 + t23 * t83 + (-t38 * t83 - t40 * t81) * qJD(5)) * t90;
t1 = -t20 * t92 + (-t22 * t81 + t24 * t83 + (-t37 * t83 - t39 * t81) * qJD(5)) * t90;
t16 = [0.2e1 * m(4) * (t43 * t52 + t44 * t51) + 0.2e1 * m(5) * (t27 * t34 + t28 * t33) - t90 * t138 + (t12 * t30 + t13 * t29) * t144 + t110 + t146 * t120; 0; 0; m(4) * (-t43 * t84 - t44 * t82 + (-t51 * t84 + t52 * t82) * qJD(1)) + m(5) * (-t27 * t84 - t28 * t82 + (-t33 * t84 + t34 * t82) * qJD(1)) + m(6) * (-t12 * t84 - t13 * t82 + (-t29 * t84 + t30 * t82) * qJD(1)); 0; 0; 0.2e1 * (m(5) * (t34 * t123 + t33 * t124 + t27 * t82 - t28 * t84) / 0.2e1 + m(6) * (t12 * t82 + t30 * t123 + t29 * t124 - t13 * t84) / 0.2e1) * t90; 0; 0; 0; -t142 + m(6) * (t12 * t32 + t13 * t31 + t14 * t30 + t15 * t29) + ((-t2 / 0.2e1 - t4 / 0.2e1) * t84 + (t1 / 0.2e1 + t5 / 0.2e1) * t82 + ((t10 / 0.2e1 + t17 / 0.2e1) * t84 + (t11 / 0.2e1 + t18 / 0.2e1) * t82) * qJD(1)) * t90; m(6) * t3; m(6) * (-t14 * t84 - t15 * t82 + (-t31 * t84 + t32 * t82) * qJD(1)); m(6) * (-t3 * t92 + (t14 * t82 - t15 * t84 + (t31 * t82 + t32 * t84) * qJD(1)) * t90); (t32 * t14 + t31 * t15 + (t41 * t84 + t42 * t82) * t3 * t90) * t144 - t92 * (-t142 + (t1 * t82 - t2 * t84 + (t10 * t84 + t11 * t82) * qJD(1)) * t90) + (-t17 * t92 + (t6 * t82 - t7 * t84) * t90) * t114 + (-t5 * t92 + ((t60 * t22 + t61 * t24 + t49 * t37 + t50 * t39) * t82 + t6 * t123 - (t60 * t21 + t61 * t23 + t49 * t38 + t50 * t40) * t84 + t7 * t124 + ((t35 * t123 + t20 * t82) * t82 - (t36 * t123 + t19 * t82) * t84) * t90) * t90) * t137 + (-t18 * t92 + (t8 * t82 - t84 * t9) * t90) * t115 - (-t4 * t92 + ((-t100 * t24 + t62 * t22 + t47 * t37 + t48 * t39) * t82 + t8 * t123 - (-t100 * t23 + t62 * t21 + t47 * t38 + t48 * t40) * t84 + t9 * t124 + ((t35 * t124 - t20 * t84) * t82 - (t36 * t124 - t19 * t84) * t84) * t90) * t90) * t135;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
