% Calculate time derivative of joint inertia matrix for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:30
% DurationCPUTime: 1.00s
% Computational Cost: add. (2370->151), mult. (1874->221), div. (0->0), fcn. (1344->8), ass. (0->100)
t89 = sin(qJ(4));
t122 = qJD(4) * t89;
t91 = cos(qJ(4));
t121 = qJD(4) * t91;
t88 = qJ(1) + qJ(2);
t84 = sin(t88);
t85 = cos(t88);
t56 = t85 * rSges(3,1) - rSges(3,2) * t84;
t83 = pkin(7) + t88;
t81 = cos(t83);
t82 = pkin(2) * t85;
t142 = -t81 * rSges(4,1) - t82;
t124 = Icges(5,4) * t91;
t105 = -Icges(5,2) * t89 + t124;
t101 = t105 * t81;
t80 = sin(t83);
t35 = Icges(5,6) * t80 + t101;
t125 = Icges(5,4) * t89;
t106 = Icges(5,1) * t91 - t125;
t102 = t106 * t81;
t37 = Icges(5,5) * t80 + t102;
t108 = t35 * t89 - t37 * t91;
t141 = t108 * t80;
t34 = -Icges(5,6) * t81 + t105 * t80;
t36 = -Icges(5,5) * t81 + t106 * t80;
t110 = t34 * t89 - t36 * t91;
t140 = t110 * t81;
t104 = Icges(5,5) * t91 - Icges(5,6) * t89;
t69 = Icges(5,2) * t91 + t125;
t70 = Icges(5,1) * t89 + t124;
t107 = t69 * t89 - t70 * t91;
t87 = qJD(1) + qJD(2);
t139 = t104 * qJD(4) + t107 * t87;
t138 = 2 * m(3);
t137 = 2 * m(4);
t136 = 2 * m(5);
t133 = pkin(2) * t84;
t90 = sin(qJ(1));
t132 = t90 * pkin(1);
t131 = rSges(5,1) * t91;
t129 = rSges(5,2) * t89;
t73 = t80 * rSges(5,3);
t128 = t80 * t87;
t127 = t81 * t87;
t66 = t80 * t129;
t126 = t81 * rSges(5,3) + t66;
t123 = pkin(1) * qJD(1);
t120 = t81 * t129;
t119 = -t87 * t120 + (-rSges(5,1) * t122 - rSges(5,2) * t121) * t80;
t118 = t90 * t123;
t92 = cos(qJ(1));
t117 = t92 * t123;
t114 = -pkin(3) - t131;
t47 = t56 * t87;
t45 = -rSges(4,2) * t80 - t142;
t55 = -rSges(3,1) * t84 - rSges(3,2) * t85;
t72 = rSges(5,1) * t89 + rSges(5,2) * t91;
t111 = t34 * t91 + t36 * t89;
t109 = t35 * t91 + t37 * t89;
t68 = Icges(5,5) * t89 + Icges(5,6) * t91;
t39 = t81 * t131 - t120 + t73;
t46 = t55 * t87;
t103 = (t106 - t69) * t122 + (t105 + t70) * t121;
t100 = t104 * t81;
t99 = (-t108 * qJD(4) + t139 * t80 + (-t91 * t105 - t89 * t106) * t128) * t80 / 0.2e1 - (-t110 * qJD(4) - t139 * t81 + (t91 * t101 + t89 * t102) * t87) * t81 / 0.2e1 + (-t107 * t80 - t81 * t68 + t111) * t128 / 0.2e1 + (-t107 * t81 + t80 * t68 + t109) * t127 / 0.2e1;
t44 = -rSges(4,1) * t80 - rSges(4,2) * t81 - t133;
t98 = t114 * t80 - t133;
t31 = rSges(4,2) * t128 + t142 * t87;
t27 = t81 * pkin(3) + t80 * pkin(6) + t39 + t82;
t30 = t44 * t87;
t97 = -t72 * t81 * qJD(4) + rSges(5,3) * t127 + t87 * t66;
t94 = Icges(5,3) * t87 - t68 * qJD(4);
t26 = t81 * pkin(6) + t126 + t98;
t13 = (-t82 + t114 * t81 + (-rSges(5,3) - pkin(6)) * t80) * t87 - t119;
t12 = pkin(6) * t127 + t98 * t87 + t97;
t86 = t92 * pkin(1);
t60 = (-t129 + t131) * qJD(4);
t49 = t56 + t86;
t48 = t55 - t132;
t43 = -t47 - t117;
t42 = t46 - t118;
t41 = t45 + t86;
t40 = t44 - t132;
t38 = t80 * t131 - t126;
t33 = Icges(5,3) * t80 + t100;
t32 = -Icges(5,3) * t81 + t104 * t80;
t29 = t31 - t117;
t28 = t30 - t118;
t25 = t27 + t86;
t24 = t26 - t132;
t19 = t87 * t100 + t94 * t80;
t18 = -t104 * t128 + t94 * t81;
t11 = t13 - t117;
t10 = t12 - t118;
t9 = -t108 * t81 + t80 * t33;
t8 = t80 * t32 - t140;
t7 = -t81 * t33 - t141;
t6 = -t110 * t80 - t81 * t32;
t1 = ((-t39 + t73) * t87 + t119) * t80 + (t87 * t38 + t97) * t81;
t2 = [(t42 * t49 + t43 * t48) * t138 + (t28 * t41 + t29 * t40) * t137 + (t10 * t25 + t11 * t24) * t136 + t103; m(3) * (t42 * t56 + t43 * t55 + t46 * t49 - t47 * t48) + m(4) * (t28 * t45 + t29 * t44 + t30 * t41 + t31 * t40) + m(5) * (t10 * t27 + t11 * t26 + t12 * t25 + t13 * t24) + t103; (t46 * t56 - t47 * t55) * t138 + (t30 * t45 + t31 * t44) * t137 + (t12 * t27 + t13 * t26) * t136 + t103; 0; 0; 0; m(5) * ((-t24 * t81 - t25 * t80) * t60 + ((-t25 * t87 - t11) * t81 + (t24 * t87 - t10) * t80) * t72) + t99; m(5) * ((-t26 * t81 - t27 * t80) * t60 + ((-t27 * t87 - t13) * t81 + (t26 * t87 - t12) * t80) * t72) + t99; m(5) * t1; ((t38 * t80 + t39 * t81) * t1 + (t80 ^ 2 + t81 ^ 2) * t72 * t60) * t136 + (-t8 * t81 + t9 * t80) * t127 + t80 * ((t80 * t18 + (t8 + t141) * t87) * t80 + (t9 * t87 + (t34 * t121 + t36 * t122) * t81 + (-t109 * qJD(4) - t110 * t87 - t19) * t80) * t81) + (-t6 * t81 + t7 * t80) * t128 - t81 * ((t81 * t19 + (t7 + t140) * t87) * t81 + (t6 * t87 + (-t35 * t121 - t37 * t122) * t80 + (t111 * qJD(4) - t108 * t87 - t18) * t81) * t80);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
