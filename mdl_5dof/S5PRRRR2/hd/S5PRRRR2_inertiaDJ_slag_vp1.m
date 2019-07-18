% Calculate time derivative of joint inertia matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:11
% EndTime: 2019-07-18 13:30:14
% DurationCPUTime: 1.28s
% Computational Cost: add. (3446->180), mult. (2586->262), div. (0->0), fcn. (1816->8), ass. (0->112)
t97 = sin(qJ(5));
t130 = qJD(5) * t97;
t99 = cos(qJ(5));
t129 = qJD(5) * t99;
t132 = Icges(6,4) * t99;
t111 = -Icges(6,2) * t97 + t132;
t96 = qJ(2) + qJ(3);
t93 = qJ(4) + t96;
t89 = cos(t93);
t107 = t111 * t89;
t88 = sin(t93);
t39 = Icges(6,6) * t88 + t107;
t133 = Icges(6,4) * t97;
t112 = Icges(6,1) * t99 - t133;
t108 = t112 * t89;
t41 = Icges(6,5) * t88 + t108;
t114 = t39 * t97 - t41 * t99;
t151 = t114 * t88;
t38 = -Icges(6,6) * t89 + t111 * t88;
t40 = -Icges(6,5) * t89 + t112 * t88;
t116 = t38 * t97 - t40 * t99;
t150 = t116 * t89;
t110 = Icges(6,5) * t99 - Icges(6,6) * t97;
t77 = Icges(6,2) * t99 + t133;
t78 = Icges(6,1) * t97 + t132;
t113 = t77 * t97 - t78 * t99;
t95 = qJD(2) + qJD(3);
t90 = qJD(4) + t95;
t149 = t110 * qJD(5) + t113 * t90;
t148 = 2 * m(4);
t147 = 2 * m(5);
t146 = 2 * m(6);
t98 = sin(qJ(2));
t143 = pkin(2) * t98;
t91 = sin(t96);
t142 = pkin(3) * t91;
t141 = rSges(6,1) * t99;
t140 = rSges(6,2) * t97;
t81 = t88 * rSges(6,3);
t139 = t88 * t90;
t138 = t88 * t99;
t137 = t89 * t90;
t136 = t91 * t95;
t92 = cos(t96);
t135 = t92 * t95;
t74 = t88 * t140;
t134 = rSges(6,3) * t137 + t90 * t74;
t131 = pkin(2) * qJD(2);
t128 = pkin(3) * t136;
t127 = pkin(3) * t135;
t126 = t89 * t140;
t75 = t89 * t141;
t123 = rSges(6,2) * t129;
t125 = -t90 * t126 + (-t130 * rSges(6,1) - t123) * t88;
t124 = t98 * t131;
t100 = cos(qJ(2));
t122 = t100 * t131;
t64 = t92 * rSges(4,1) - rSges(4,2) * t91;
t59 = t89 * rSges(5,1) - rSges(5,2) * t88;
t53 = -rSges(4,1) * t135 + rSges(4,2) * t136;
t49 = -rSges(5,1) * t137 + rSges(5,2) * t139;
t87 = pkin(3) * t92;
t51 = t59 + t87;
t63 = -rSges(4,1) * t91 - rSges(4,2) * t92;
t58 = -rSges(5,1) * t88 - rSges(5,2) * t89;
t80 = rSges(6,1) * t97 + rSges(6,2) * t99;
t117 = t38 * t99 + t40 * t97;
t115 = t39 * t99 + t41 * t97;
t76 = Icges(6,5) * t97 + Icges(6,6) * t99;
t42 = rSges(6,1) * t138 - t89 * rSges(6,3) - t74;
t43 = t75 + t81 - t126;
t33 = t89 * pkin(6) - t42;
t32 = t88 * pkin(6) + t43;
t52 = t63 * t95;
t48 = t58 * t90;
t109 = (t112 - t77) * t130 + (t111 + t78) * t129;
t106 = t110 * t89;
t105 = (-t114 * qJD(5) + t149 * t88 + (-t111 * t99 - t112 * t97) * t139) * t88 / 0.2e1 - (-t116 * qJD(5) - t149 * t89 + (t107 * t99 + t108 * t97) * t90) * t89 / 0.2e1 + (-t113 * t88 - t76 * t89 + t117) * t139 / 0.2e1 + (-t113 * t89 + t76 * t88 + t115) * t137 / 0.2e1;
t29 = t32 + t87;
t50 = t58 - t142;
t35 = t49 - t127;
t28 = t33 - t142;
t102 = Icges(6,3) * t90 - t76 * qJD(5);
t34 = t48 - t128;
t14 = (-t75 + (-rSges(6,3) - pkin(6)) * t88) * t90 - t125;
t13 = t14 - t127;
t15 = -t89 * t123 + pkin(6) * t137 + (-t89 * t130 - t90 * t138) * rSges(6,1) + t134;
t12 = t15 - t128;
t94 = t100 * pkin(2);
t70 = (-t140 + t141) * qJD(5);
t55 = t64 + t94;
t54 = t63 - t143;
t47 = t53 - t122;
t46 = t52 - t124;
t45 = t51 + t94;
t44 = t50 - t143;
t37 = Icges(6,3) * t88 + t106;
t36 = -Icges(6,3) * t89 + t110 * t88;
t31 = t35 - t122;
t30 = t34 - t124;
t27 = t29 + t94;
t26 = t28 - t143;
t21 = t102 * t88 + t90 * t106;
t20 = t102 * t89 - t110 * t139;
t11 = t13 - t122;
t10 = t12 - t124;
t9 = -t114 * t89 + t37 * t88;
t8 = t36 * t88 - t150;
t7 = -t37 * t89 - t151;
t6 = -t116 * t88 - t36 * t89;
t1 = ((-t43 + t81) * t90 + t125) * t88 + (-t80 * t89 * qJD(5) + t90 * t42 + t134) * t89;
t2 = [0; 0; (t10 * t27 + t11 * t26) * t146 + (t30 * t45 + t31 * t44) * t147 + (t46 * t55 + t47 * t54) * t148 + t109; 0; m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + m(5) * (t30 * t51 + t31 * t50 + t34 * t45 + t35 * t44) + m(4) * (t46 * t64 + t47 * t63 + t52 * t55 + t53 * t54) + t109; (t52 * t64 + t53 * t63) * t148 + (t34 * t51 + t35 * t50) * t147 + (t12 * t29 + t13 * t28) * t146 + t109; 0; m(6) * (t10 * t32 + t11 * t33 + t14 * t26 + t15 * t27) + m(5) * (t30 * t59 + t31 * t58 + t44 * t49 + t45 * t48) + t109; m(5) * (t34 * t59 + t35 * t58 + t48 * t51 + t49 * t50) + m(6) * (t12 * t32 + t13 * t33 + t14 * t28 + t15 * t29) + t109; (t48 * t59 + t49 * t58) * t147 + (t14 * t33 + t15 * t32) * t146 + t109; m(6) * t1; m(6) * ((-t26 * t89 - t27 * t88) * t70 + ((-t27 * t90 - t11) * t89 + (t26 * t90 - t10) * t88) * t80) + t105; m(6) * ((-t28 * t89 - t29 * t88) * t70 + ((-t29 * t90 - t13) * t89 + (t28 * t90 - t12) * t88) * t80) + t105; m(6) * ((-t32 * t88 - t33 * t89) * t70 + ((-t32 * t90 - t14) * t89 + (t33 * t90 - t15) * t88) * t80) + t105; ((t42 * t88 + t43 * t89) * t1 + (t88 ^ 2 + t89 ^ 2) * t80 * t70) * t146 + (-t8 * t89 + t9 * t88) * t137 + t88 * ((t88 * t20 + (t8 + t151) * t90) * t88 + (t9 * t90 + (t38 * t129 + t40 * t130) * t89 + (-t115 * qJD(5) - t116 * t90 - t21) * t88) * t89) + (-t6 * t89 + t7 * t88) * t139 - t89 * ((t89 * t21 + (t7 + t150) * t90) * t89 + (t6 * t90 + (-t39 * t129 - t41 * t130) * t88 + (t117 * qJD(5) - t114 * t90 - t20) * t89) * t88);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq  = res;
