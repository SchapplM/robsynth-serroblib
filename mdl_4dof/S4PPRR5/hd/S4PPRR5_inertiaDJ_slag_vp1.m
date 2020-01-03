% Calculate time derivative of joint inertia matrix for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:45
% DurationCPUTime: 1.55s
% Computational Cost: add. (1892->208), mult. (5502->362), div. (0->0), fcn. (5259->6), ass. (0->121)
t82 = cos(pkin(6));
t80 = t82 ^ 2;
t81 = sin(pkin(6));
t147 = t81 ^ 2 + t80;
t145 = t147 * qJD(3);
t117 = t81 * qJD(3);
t84 = sin(qJ(3));
t113 = t84 * t117;
t86 = cos(qJ(3));
t108 = -rSges(4,1) * t84 - rSges(4,2) * t86;
t146 = t108 * t145;
t143 = 2 * m(5);
t142 = t81 / 0.2e1;
t141 = t82 / 0.2e1;
t136 = t81 * t86;
t122 = Icges(5,3) * t86;
t83 = sin(qJ(4));
t134 = t83 * t84;
t85 = cos(qJ(4));
t73 = t134 * t82 + t81 * t85;
t131 = t84 * t85;
t89 = t131 * t82 - t81 * t83;
t40 = -Icges(5,5) * t89 + Icges(5,6) * t73 + t122 * t82;
t123 = Icges(5,6) * t86;
t42 = -Icges(5,4) * t89 + Icges(5,2) * t73 + t123 * t82;
t124 = Icges(5,5) * t86;
t44 = -Icges(5,1) * t89 + Icges(5,4) * t73 + t124 * t82;
t71 = -t134 * t81 + t82 * t85;
t72 = t131 * t81 + t82 * t83;
t16 = -t136 * t40 + t71 * t42 + t72 * t44;
t140 = t16 * t82;
t135 = t82 * t86;
t39 = Icges(5,5) * t72 + Icges(5,6) * t71 - t122 * t81;
t41 = Icges(5,4) * t72 + Icges(5,2) * t71 - t123 * t81;
t43 = Icges(5,1) * t72 + Icges(5,4) * t71 - t124 * t81;
t17 = t135 * t39 + t73 * t41 - t43 * t89;
t139 = t17 * t81;
t107 = rSges(5,1) * t85 - rSges(5,2) * t83;
t64 = t84 * rSges(5,3) + t107 * t86;
t138 = t64 * t84;
t118 = qJD(4) * t86;
t90 = Icges(5,5) * t85 - Icges(5,6) * t83;
t133 = t84 * ((-Icges(5,5) * t83 - Icges(5,6) * t85) * t118 + (-t84 * t90 + t122) * qJD(3));
t61 = Icges(5,3) * t84 + t86 * t90;
t132 = t84 * t61;
t109 = -pkin(3) * t84 + pkin(5) * t86;
t50 = (-rSges(5,1) * t83 - rSges(5,2) * t85) * t118 + (rSges(5,3) * t86 - t107 * t84) * qJD(3);
t130 = t109 * qJD(3) + t50;
t78 = pkin(3) * t86 + t84 * pkin(5);
t129 = t64 + t78;
t126 = Icges(5,4) * t83;
t125 = Icges(5,4) * t85;
t120 = qJD(3) * t84;
t119 = qJD(3) * t86;
t116 = t82 * t120;
t115 = t83 * t119;
t114 = t85 * t119;
t112 = Icges(5,5) * t120;
t111 = Icges(5,6) * t120;
t110 = Icges(5,3) * t120;
t77 = rSges(4,1) * t86 - t84 * rSges(4,2);
t105 = t41 * t83 - t43 * t85;
t19 = -t105 * t86 + t84 * t39;
t104 = t42 * t83 - t44 * t85;
t20 = -t104 * t86 + t84 * t40;
t106 = -t19 * t81 + t20 * t82;
t45 = t72 * rSges(5,1) + t71 * rSges(5,2) - rSges(5,3) * t136;
t46 = -rSges(5,1) * t89 + t73 * rSges(5,2) + rSges(5,3) * t135;
t103 = t45 * t82 + t46 * t81;
t92 = -Icges(5,2) * t83 + t125;
t62 = Icges(5,6) * t84 + t86 * t92;
t95 = Icges(5,1) * t85 - t126;
t63 = Icges(5,5) * t84 + t86 * t95;
t100 = t62 * t83 - t63 * t85;
t91 = Icges(4,5) * t86 - Icges(4,6) * t84;
t53 = -qJD(4) * t72 - t115 * t81;
t54 = qJD(4) * t71 + t114 * t81;
t28 = Icges(5,5) * t54 + Icges(5,6) * t53 + t110 * t81;
t88 = t120 * t39 - t28 * t86;
t55 = qJD(4) * t89 + t115 * t82;
t56 = qJD(4) * t73 - t114 * t82;
t29 = Icges(5,5) * t56 + Icges(5,6) * t55 - t110 * t82;
t87 = t120 * t40 - t29 * t86;
t65 = t91 * t117;
t52 = t129 * t82;
t51 = t129 * t81;
t49 = (-Icges(5,1) * t83 - t125) * t118 + (-t84 * t95 + t124) * qJD(3);
t48 = (-Icges(5,2) * t85 - t126) * t118 + (-t84 * t92 + t123) * qJD(3);
t38 = t77 * t145;
t37 = t130 * t82;
t36 = t130 * t81;
t35 = rSges(5,1) * t56 + rSges(5,2) * t55 - rSges(5,3) * t116;
t34 = rSges(5,1) * t54 + rSges(5,2) * t53 + rSges(5,3) * t113;
t33 = Icges(5,1) * t56 + Icges(5,4) * t55 - t112 * t82;
t32 = Icges(5,1) * t54 + Icges(5,4) * t53 + t112 * t81;
t31 = Icges(5,4) * t56 + Icges(5,2) * t55 - t111 * t82;
t30 = Icges(5,4) * t54 + Icges(5,2) * t53 + t111 * t81;
t27 = t135 * t64 - t84 * t46;
t26 = t136 * t64 + t84 * t45;
t25 = -t100 * t86 + t132;
t24 = t135 * t61 + t73 * t62 - t63 * t89;
t23 = -t136 * t61 + t71 * t62 + t72 * t63;
t22 = t103 * t86;
t21 = (t109 * t82 + t46) * t82 + (t109 * t81 - t45) * t81;
t18 = t135 * t40 + t73 * t42 - t44 * t89;
t15 = -t136 * t39 + t71 * t41 + t72 * t43;
t14 = -t78 * t145 - t81 * t34 + t82 * t35;
t13 = t50 * t135 - t84 * t35 + (-t138 * t82 - t46 * t86) * qJD(3);
t12 = t50 * t136 + t84 * t34 + (-t138 * t81 + t45 * t86) * qJD(3);
t11 = (-t34 * t82 - t35 * t81) * t86 + t103 * t120;
t10 = t73 * t31 - t33 * t89 + t55 * t42 + t56 * t44 - t82 * t87;
t9 = t73 * t30 - t32 * t89 + t55 * t41 + t56 * t43 - t82 * t88;
t8 = t71 * t31 + t72 * t33 + t53 * t42 + t54 * t44 + t81 * t87;
t7 = t71 * t30 + t72 * t32 + t53 * t41 + t54 * t43 + t81 * t88;
t6 = (qJD(3) * t104 + t29) * t84 + (qJD(3) * t40 - t31 * t83 + t33 * t85 + (-t42 * t85 - t44 * t83) * qJD(4)) * t86;
t5 = (qJD(3) * t105 + t28) * t84 + (qJD(3) * t39 - t30 * t83 + t32 * t85 + (-t41 * t85 - t43 * t83) * qJD(4)) * t86;
t4 = t10 * t81 + t82 * t9;
t3 = t7 * t82 + t8 * t81;
t2 = (t73 * t48 - t49 * t89 + t55 * t62 + t56 * t63) * t84 + (-t9 * t81 + (t10 + t133) * t82) * t86 + (t24 * t86 + (t139 + (-t18 - t132) * t82) * t84) * qJD(3);
t1 = (t71 * t48 + t72 * t49 + t53 * t62 + t54 * t63) * t84 + (t8 * t82 + (-t7 - t133) * t81) * t86 + (t23 * t86 + (-t140 + (t15 + t132) * t81) * t84) * qJD(3);
t47 = [0; 0; 0; -m(4) * t38 + m(5) * t14; m(5) * (t36 * t81 + t37 * t82) + m(4) * t146; 0.2e1 * m(4) * (-t147 * t38 * t108 + t77 * t146) + (t14 * t21 + t36 * t51 + t37 * t52) * t143 + (t80 * t65 + t3) * t82 + (t4 + (-t91 * t145 + t81 * t65) * t82) * t81; m(5) * t11; m(5) * (-t12 * t82 + t13 * t81); m(5) * (t11 * t21 - t12 * t52 + t13 * t51 - t14 * t22 - t26 * t37 + t27 * t36) + t1 * t141 + t2 * t142 + t84 * (t5 * t82 + t6 * t81) / 0.2e1 + (-t81 * t3 / 0.2e1 + t4 * t141) * t86 + (t86 * (t19 * t82 + t20 * t81) / 0.2e1 + ((t15 * t82 + t16 * t81) * t142 - t82 * (t17 * t82 + t18 * t81) / 0.2e1) * t84) * qJD(3); (-t11 * t22 + t12 * t26 + t13 * t27) * t143 + (t23 * t84 + (-t81 * t15 + t140) * t86) * t113 - t1 * t136 - (t24 * t84 + (t82 * t18 - t139) * t86) * t116 + t2 * t135 + (t106 * t86 + t25 * t84) * t119 + t84 * ((t133 + (t100 * t84 - t106) * qJD(3)) * t84 + (-t5 * t81 + t6 * t82 + (qJD(3) * t61 - t48 * t83 + t49 * t85 + (-t62 * t85 - t63 * t83) * qJD(4)) * t84 + t25 * qJD(3)) * t86);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t47(1), t47(2), t47(4), t47(7); t47(2), t47(3), t47(5), t47(8); t47(4), t47(5), t47(6), t47(9); t47(7), t47(8), t47(9), t47(10);];
Mq = res;
