% Calculate time derivative of joint inertia matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:11
% DurationCPUTime: 1.30s
% Computational Cost: add. (3585->185), mult. (2614->265), div. (0->0), fcn. (1838->8), ass. (0->112)
t97 = sin(qJ(4));
t131 = qJD(4) * t97;
t99 = cos(qJ(4));
t130 = qJD(4) * t99;
t133 = Icges(5,4) * t99;
t111 = -Icges(5,2) * t97 + t133;
t96 = qJ(1) + qJ(2);
t93 = qJ(3) + t96;
t89 = cos(t93);
t107 = t111 * t89;
t88 = sin(t93);
t38 = Icges(5,6) * t88 + t107;
t134 = Icges(5,4) * t97;
t112 = Icges(5,1) * t99 - t134;
t108 = t112 * t89;
t40 = Icges(5,5) * t88 + t108;
t114 = t38 * t97 - t40 * t99;
t153 = t114 * t88;
t37 = -Icges(5,6) * t89 + t111 * t88;
t39 = -Icges(5,5) * t89 + t112 * t88;
t116 = t37 * t97 - t39 * t99;
t152 = t116 * t89;
t110 = Icges(5,5) * t99 - Icges(5,6) * t97;
t76 = Icges(5,2) * t99 + t134;
t77 = Icges(5,1) * t97 + t133;
t113 = t76 * t97 - t77 * t99;
t95 = qJD(1) + qJD(2);
t90 = qJD(3) + t95;
t151 = t110 * qJD(4) + t113 * t90;
t150 = 2 * m(3);
t149 = 2 * m(4);
t148 = 2 * m(5);
t98 = sin(qJ(1));
t145 = pkin(1) * t98;
t91 = sin(t96);
t144 = pkin(2) * t91;
t143 = rSges(5,1) * t99;
t142 = rSges(5,2) * t97;
t80 = t88 * rSges(5,3);
t141 = t88 * t90;
t140 = t88 * t99;
t139 = t89 * t90;
t138 = t91 * t95;
t92 = cos(t96);
t137 = t92 * t95;
t73 = t88 * t142;
t136 = rSges(5,3) * t139 + t90 * t73;
t135 = t89 * rSges(5,3) + t73;
t132 = pkin(1) * qJD(1);
t129 = pkin(2) * t138;
t128 = pkin(2) * t137;
t127 = t89 * t142;
t124 = rSges(5,2) * t130;
t126 = -t90 * t127 + (-t131 * rSges(5,1) - t124) * t88;
t125 = t98 * t132;
t100 = cos(qJ(1));
t123 = t100 * t132;
t120 = -pkin(3) - t143;
t63 = t92 * rSges(3,1) - rSges(3,2) * t91;
t58 = t89 * rSges(4,1) - rSges(4,2) * t88;
t52 = -rSges(3,1) * t137 + rSges(3,2) * t138;
t48 = -rSges(4,1) * t139 + rSges(4,2) * t141;
t87 = pkin(2) * t92;
t50 = t58 + t87;
t62 = -rSges(3,1) * t91 - rSges(3,2) * t92;
t57 = -rSges(4,1) * t88 - rSges(4,2) * t89;
t79 = rSges(5,1) * t97 + rSges(5,2) * t99;
t117 = t37 * t99 + t39 * t97;
t115 = t38 * t99 + t40 * t97;
t75 = Icges(5,5) * t97 + Icges(5,6) * t99;
t42 = t89 * t143 - t127 + t80;
t51 = t62 * t95;
t47 = t57 * t90;
t109 = (t112 - t76) * t131 + (t111 + t77) * t130;
t106 = t110 * t89;
t105 = (-t114 * qJD(4) + t151 * t88 + (-t111 * t99 - t112 * t97) * t141) * t88 / 0.2e1 - (-t116 * qJD(4) - t151 * t89 + (t107 * t99 + t108 * t97) * t90) * t89 / 0.2e1 + (-t113 * t88 - t75 * t89 + t117) * t141 / 0.2e1 + (-t113 * t89 + t75 * t88 + t115) * t139 / 0.2e1;
t30 = t89 * pkin(3) + t88 * pkin(7) + t42;
t49 = t57 - t144;
t34 = t48 - t128;
t28 = t30 + t87;
t29 = t89 * pkin(7) + t120 * t88 + t135;
t102 = Icges(5,3) * t90 - t75 * qJD(4);
t33 = t47 - t129;
t27 = t29 - t144;
t14 = (t120 * t89 + (-rSges(5,3) - pkin(7)) * t88) * t90 - t126;
t12 = t14 - t128;
t13 = -t89 * t124 - pkin(3) * t141 + pkin(7) * t139 + (-t89 * t131 - t90 * t140) * rSges(5,1) + t136;
t11 = t13 - t129;
t94 = t100 * pkin(1);
t69 = (-t142 + t143) * qJD(4);
t54 = t63 + t94;
t53 = t62 - t145;
t46 = t52 - t123;
t45 = t51 - t125;
t44 = t50 + t94;
t43 = t49 - t145;
t41 = rSges(5,1) * t140 - t135;
t36 = Icges(5,3) * t88 + t106;
t35 = -Icges(5,3) * t89 + t110 * t88;
t32 = t34 - t123;
t31 = t33 - t125;
t26 = t28 + t94;
t25 = t27 - t145;
t20 = t102 * t88 + t90 * t106;
t19 = t102 * t89 - t110 * t141;
t10 = t12 - t123;
t9 = t11 - t125;
t8 = -t114 * t89 + t36 * t88;
t7 = t35 * t88 - t152;
t6 = -t36 * t89 - t153;
t5 = -t116 * t88 - t35 * t89;
t1 = [(t45 * t54 + t46 * t53) * t150 + (t31 * t44 + t32 * t43) * t149 + (t10 * t25 + t26 * t9) * t148 + t109; m(3) * (t45 * t63 + t46 * t62 + t51 * t54 + t52 * t53) + m(4) * (t31 * t50 + t32 * t49 + t33 * t44 + t34 * t43) + m(5) * (t10 * t27 + t11 * t26 + t12 * t25 + t28 * t9) + t109; (t51 * t63 + t52 * t62) * t150 + (t33 * t50 + t34 * t49) * t149 + (t11 * t28 + t12 * t27) * t148 + t109; m(4) * (t31 * t58 + t32 * t57 + t43 * t48 + t44 * t47) + m(5) * (t10 * t29 + t13 * t26 + t14 * t25 + t30 * t9) + t109; m(4) * (t33 * t58 + t34 * t57 + t47 * t50 + t48 * t49) + m(5) * (t11 * t30 + t12 * t29 + t13 * t28 + t14 * t27) + t109; (t47 * t58 + t48 * t57) * t149 + (t13 * t30 + t14 * t29) * t148 + t109; m(5) * ((-t25 * t89 - t26 * t88) * t69 + ((-t26 * t90 - t10) * t89 + (t25 * t90 - t9) * t88) * t79) + t105; m(5) * ((-t27 * t89 - t28 * t88) * t69 + ((-t28 * t90 - t12) * t89 + (t27 * t90 - t11) * t88) * t79) + t105; m(5) * ((-t29 * t89 - t30 * t88) * t69 + ((-t30 * t90 - t14) * t89 + (t29 * t90 - t13) * t88) * t79) + t105; ((t41 * t88 + t42 * t89) * (((-t42 + t80) * t90 + t126) * t88 + (-t79 * t89 * qJD(4) + t90 * t41 + t136) * t89) + (t88 ^ 2 + t89 ^ 2) * t79 * t69) * t148 + (-t7 * t89 + t8 * t88) * t139 + t88 * ((t88 * t19 + (t7 + t153) * t90) * t88 + (t8 * t90 + (t37 * t130 + t39 * t131) * t89 + (-t115 * qJD(4) - t116 * t90 - t20) * t88) * t89) + (-t5 * t89 + t6 * t88) * t141 - t89 * ((t89 * t20 + (t6 + t152) * t90) * t89 + (t5 * t90 + (-t38 * t130 - t40 * t131) * t88 + (t117 * qJD(4) - t114 * t90 - t19) * t89) * t88);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
