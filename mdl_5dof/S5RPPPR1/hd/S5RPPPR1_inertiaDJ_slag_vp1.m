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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:50
% DurationCPUTime: 1.70s
% Computational Cost: add. (3787->236), mult. (3828->381), div. (0->0), fcn. (3531->10), ass. (0->119)
t77 = pkin(9) + qJ(5);
t75 = cos(t77);
t121 = Icges(6,4) * t75;
t73 = sin(t77);
t80 = sin(pkin(8));
t82 = cos(pkin(8));
t57 = -Icges(6,6) * t82 + (-Icges(6,2) * t73 + t121) * t80;
t122 = Icges(6,4) * t73;
t58 = -Icges(6,5) * t82 + (Icges(6,1) * t75 - t122) * t80;
t141 = -t57 * t75 - t58 * t73;
t81 = cos(pkin(9));
t140 = (rSges(6,3) + pkin(6) + qJ(4)) * t80 + (pkin(4) * t81 + pkin(3)) * t82;
t79 = sin(pkin(9));
t139 = (rSges(5,3) + qJ(4)) * t80 + (rSges(5,1) * t81 - rSges(5,2) * t79 + pkin(3)) * t82;
t138 = 2 * m(6);
t137 = pkin(4) * t79;
t136 = sin(qJ(1)) * pkin(1);
t135 = cos(qJ(1)) * pkin(1);
t111 = qJD(5) * t80;
t64 = (-Icges(6,5) * t73 - Icges(6,6) * t75) * t111;
t66 = (-Icges(6,1) * t73 - t121) * t111;
t103 = t66 * t75 * t80 - t82 * t64;
t65 = (-Icges(6,2) * t75 - t122) * t111;
t130 = t73 * t65;
t134 = ((qJD(5) * t141 - t130) * t80 + t103) * t82;
t78 = qJ(1) + pkin(7);
t74 = sin(t78);
t129 = t74 * t80;
t128 = t74 * t82;
t76 = cos(t78);
t127 = t76 * t80;
t126 = t76 * t82;
t92 = t126 * t73 - t74 * t75;
t93 = t128 * t75 - t73 * t76;
t49 = qJD(1) * t92 + qJD(5) * t93;
t60 = t128 * t73 + t75 * t76;
t63 = t126 * t75 + t73 * t74;
t50 = -qJD(1) * t63 + qJD(5) * t60;
t124 = rSges(6,1) * t50 + rSges(6,2) * t49;
t123 = -rSges(6,1) * t93 + rSges(6,2) * t60;
t120 = Icges(6,5) * t80;
t119 = Icges(6,6) * t80;
t118 = Icges(6,3) * t80;
t117 = -rSges(4,3) - qJ(3);
t115 = qJD(1) * t74;
t114 = qJD(1) * t76;
t113 = qJD(1) * t80;
t112 = qJD(4) * t80;
t110 = t74 * t113;
t109 = t76 * t113;
t108 = Icges(6,5) * t113;
t107 = Icges(6,6) * t113;
t106 = Icges(6,3) * t113;
t105 = t76 * qJ(3) - t136;
t104 = -qJ(3) - t137;
t69 = qJD(3) * t76;
t102 = -t112 * t74 + t69;
t101 = pkin(2) * t115 + qJD(1) * t136 - qJD(3) * t74;
t99 = rSges(4,1) * t82 - rSges(4,2) * t80;
t98 = rSges(5,1) * t79 + rSges(5,2) * t81;
t47 = qJD(1) * t60 - qJD(5) * t63;
t48 = -qJD(1) * t93 - qJD(5) * t92;
t97 = -rSges(6,1) * t48 - rSges(6,2) * t47;
t96 = -t63 * rSges(6,1) + rSges(6,2) * t92;
t95 = -pkin(2) - t99;
t91 = -qJ(3) - t98;
t90 = -pkin(2) - t140;
t88 = -t112 * t76 + t101;
t87 = -pkin(2) - t139;
t52 = t117 * t74 + t76 * t95 - t135;
t86 = t104 * t74 + t76 * t90 - t135;
t34 = t74 * t91 + t76 * t87 - t135;
t67 = (-rSges(6,1) * t73 - rSges(6,2) * t75) * t111;
t59 = -rSges(6,3) * t82 + (rSges(6,1) * t75 - rSges(6,2) * t73) * t80;
t56 = -Icges(6,3) * t82 + (Icges(6,5) * t75 - Icges(6,6) * t73) * t80;
t51 = rSges(4,3) * t76 + t74 * t95 + t105;
t44 = qJD(1) * t52 + t69;
t43 = (t117 * t76 + t74 * t99) * qJD(1) + t101;
t42 = rSges(6,3) * t127 - t96;
t41 = -rSges(6,3) * t129 + t123;
t40 = Icges(6,1) * t63 - Icges(6,4) * t92 + t120 * t76;
t39 = -Icges(6,1) * t93 + Icges(6,4) * t60 - t120 * t74;
t38 = Icges(6,4) * t63 - Icges(6,2) * t92 + t119 * t76;
t37 = -Icges(6,4) * t93 + Icges(6,2) * t60 - t119 * t74;
t36 = Icges(6,5) * t63 - Icges(6,6) * t92 + t118 * t76;
t35 = -Icges(6,5) * t93 + Icges(6,6) * t60 - t118 * t74;
t33 = t74 * t87 + t76 * t98 + t105;
t32 = t127 * t59 + t42 * t82;
t31 = t129 * t59 - t41 * t82;
t30 = t86 + t96;
t29 = t137 * t76 + t74 * t90 + t105 + t123;
t28 = qJD(1) * t34 + t102;
t27 = (t139 * t74 + t91 * t76) * qJD(1) + t88;
t26 = -rSges(6,3) * t109 + t124;
t25 = -rSges(6,3) * t110 - t97;
t24 = Icges(6,1) * t50 + Icges(6,4) * t49 - t108 * t76;
t23 = Icges(6,1) * t48 + Icges(6,4) * t47 - t108 * t74;
t22 = Icges(6,4) * t50 + Icges(6,2) * t49 - t107 * t76;
t21 = Icges(6,4) * t48 + Icges(6,2) * t47 - t107 * t74;
t20 = Icges(6,5) * t50 + Icges(6,6) * t49 - t106 * t76;
t19 = Icges(6,5) * t48 + Icges(6,6) * t47 - t106 * t74;
t18 = t127 * t56 - t57 * t92 + t58 * t63;
t17 = -t129 * t56 + t57 * t60 - t58 * t93;
t15 = -t26 * t82 + (t114 * t59 + t67 * t74) * t80;
t14 = t25 * t82 + (-t115 * t59 + t67 * t76) * t80;
t13 = qJD(1) * t86 + t102 + t124;
t12 = (t104 * t76 + t140 * t74) * qJD(1) + t88 + t97;
t11 = -t36 * t82 + (-t38 * t73 + t40 * t75) * t80;
t10 = -t35 * t82 + (-t37 * t73 + t39 * t75) * t80;
t9 = t127 * t36 - t38 * t92 + t40 * t63;
t8 = t127 * t35 - t37 * t92 + t39 * t63;
t7 = -t129 * t36 + t38 * t60 - t40 * t93;
t6 = -t129 * t35 + t37 * t60 - t39 * t93;
t5 = t49 * t57 + t50 * t58 + t60 * t65 - t93 * t66 + (-t114 * t56 - t64 * t74) * t80;
t4 = t47 * t57 + t48 * t58 - t92 * t65 + t63 * t66 + (-t115 * t56 + t64 * t76) * t80;
t3 = (-t25 * t74 - t26 * t76 + (t41 * t74 - t42 * t76) * qJD(1)) * t80;
t2 = -t19 * t82 + (-t21 * t73 + t23 * t75 + (-t38 * t75 - t40 * t73) * qJD(5)) * t80;
t1 = -t20 * t82 + (-t22 * t73 + t24 * t75 + (-t37 * t75 - t39 * t73) * qJD(5)) * t80;
t16 = [0.2e1 * m(4) * (t43 * t52 + t44 * t51) + 0.2e1 * m(5) * (t27 * t34 + t28 * t33) - t80 * t130 + (t12 * t30 + t13 * t29) * t138 + t103 + t141 * t111; 0; 0; m(4) * (t43 * t76 + t44 * t74 + (t51 * t76 - t52 * t74) * qJD(1)) + m(5) * (t27 * t76 + t28 * t74 + (t33 * t76 - t34 * t74) * qJD(1)) + m(6) * (t12 * t76 + t13 * t74 + (t29 * t76 - t30 * t74) * qJD(1)); 0; 0; 0.2e1 * (m(5) * (-t114 * t34 - t115 * t33 - t27 * t74 + t28 * t76) / 0.2e1 + m(6) * (-t114 * t30 - t115 * t29 - t12 * t74 + t13 * t76) / 0.2e1) * t80; 0; 0; 0; -t134 + m(6) * (t12 * t32 + t13 * t31 + t14 * t30 + t15 * t29) + ((t2 / 0.2e1 + t4 / 0.2e1) * t76 + (-t1 / 0.2e1 - t5 / 0.2e1) * t74 + ((-t10 / 0.2e1 - t17 / 0.2e1) * t76 + (-t11 / 0.2e1 - t18 / 0.2e1) * t74) * qJD(1)) * t80; m(6) * t3; m(6) * (t14 * t76 + t15 * t74 + (t31 * t76 - t32 * t74) * qJD(1)); m(6) * (-t3 * t82 + (-t14 * t74 + t15 * t76 + (-t31 * t74 - t32 * t76) * qJD(1)) * t80); (t32 * t14 + t31 * t15 + (-t41 * t76 - t42 * t74) * t3 * t80) * t138 - t82 * (-t134 + (-t1 * t74 + t2 * t76 + (-t10 * t76 - t11 * t74) * qJD(1)) * t80) - (-t17 * t82 + (-t6 * t74 + t7 * t76) * t80) * t109 - (-t5 * t82 + (-(t60 * t22 - t24 * t93 + t49 * t37 + t50 * t39) * t74 - t6 * t114 + (t60 * t21 - t23 * t93 + t49 * t38 + t50 * t40) * t76 - t7 * t115 + (-(-t114 * t35 - t20 * t74) * t74 + (-t114 * t36 - t19 * t74) * t76) * t80) * t80) * t129 - (-t18 * t82 + (-t74 * t8 + t76 * t9) * t80) * t110 + (-t4 * t82 + (-(-t22 * t92 + t63 * t24 + t47 * t37 + t48 * t39) * t74 - t8 * t114 + (-t21 * t92 + t63 * t23 + t47 * t38 + t48 * t40) * t76 - t9 * t115 + (-(-t115 * t35 + t20 * t76) * t74 + (-t115 * t36 + t19 * t76) * t76) * t80) * t80) * t127;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
