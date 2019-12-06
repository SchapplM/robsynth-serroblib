% Calculate time derivative of joint inertia matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:09
% DurationCPUTime: 1.43s
% Computational Cost: add. (4738->190), mult. (2664->274), div. (0->0), fcn. (1870->8), ass. (0->112)
t101 = sin(qJ(5));
t129 = qJD(5) * t101;
t102 = cos(qJ(5));
t128 = qJD(5) * t102;
t130 = Icges(6,4) * t102;
t113 = -Icges(6,2) * t101 + t130;
t99 = pkin(9) + qJ(2);
t98 = qJ(3) + t99;
t94 = qJ(4) + t98;
t90 = cos(t94);
t108 = t113 * t90;
t89 = sin(t94);
t39 = Icges(6,6) * t89 + t108;
t131 = Icges(6,4) * t101;
t114 = Icges(6,1) * t102 - t131;
t109 = t114 * t90;
t41 = Icges(6,5) * t89 + t109;
t116 = t101 * t39 - t102 * t41;
t153 = t116 * t89;
t38 = -Icges(6,6) * t90 + t113 * t89;
t40 = -Icges(6,5) * t90 + t114 * t89;
t117 = t101 * t38 - t102 * t40;
t152 = t117 * t90;
t112 = Icges(6,5) * t102 - Icges(6,6) * t101;
t78 = Icges(6,2) * t102 + t131;
t79 = Icges(6,1) * t101 + t130;
t115 = t101 * t78 - t102 * t79;
t100 = qJD(2) + qJD(3);
t97 = qJD(4) + t100;
t151 = t112 * qJD(5) + t115 * t97;
t150 = 2 * m(4);
t149 = 2 * m(5);
t148 = 2 * m(6);
t95 = sin(t99);
t145 = pkin(2) * t95;
t92 = sin(t98);
t144 = pkin(3) * t92;
t81 = t89 * rSges(6,3);
t141 = t89 * t97;
t140 = t90 * t97;
t136 = rSges(6,2) * t101;
t74 = t89 * t136;
t139 = rSges(6,3) * t140 + t97 * t74;
t138 = t90 * rSges(6,3) + t74;
t137 = rSges(6,1) * t102;
t135 = pkin(2) * qJD(2);
t134 = t100 * t92;
t93 = cos(t98);
t133 = t100 * t93;
t132 = t102 * t89;
t127 = pkin(3) * t134;
t126 = pkin(3) * t133;
t119 = rSges(6,2) * t128;
t124 = t90 * t136;
t125 = -t97 * t124 + (-t129 * rSges(6,1) - t119) * t89;
t123 = t95 * t135;
t96 = cos(t99);
t122 = t96 * t135;
t118 = -pkin(4) - t137;
t61 = t93 * rSges(4,1) - rSges(4,2) * t92;
t57 = t90 * rSges(5,1) - rSges(5,2) * t89;
t49 = -rSges(5,1) * t140 + rSges(5,2) * t141;
t88 = pkin(3) * t93;
t51 = t57 + t88;
t53 = -rSges(4,1) * t133 + rSges(4,2) * t134;
t60 = -rSges(4,1) * t92 - rSges(4,2) * t93;
t56 = -rSges(5,1) * t89 - rSges(5,2) * t90;
t80 = t101 * rSges(6,1) + rSges(6,2) * t102;
t43 = t90 * t137 - t124 + t81;
t77 = Icges(6,5) * t101 + Icges(6,6) * t102;
t48 = t56 * t97;
t52 = t60 * t100;
t111 = (t114 - t78) * t129 + (t113 + t79) * t128;
t110 = (-t116 * qJD(5) + t151 * t89 + (-t101 * t114 - t102 * t113) * t141) * t89 / 0.2e1 - (-t117 * qJD(5) - t151 * t90 + (t101 * t109 + t102 * t108) * t97) * t90 / 0.2e1 + (t101 * t40 + t102 * t38 - t115 * t89 - t90 * t77) * t141 / 0.2e1 + (t101 * t41 + t102 * t39 - t115 * t90 + t89 * t77) * t140 / 0.2e1;
t107 = t112 * t90;
t50 = t56 - t144;
t31 = t90 * pkin(4) + t89 * pkin(8) + t43;
t35 = t49 - t126;
t29 = t31 + t88;
t30 = t90 * pkin(8) + t118 * t89 + t138;
t34 = t48 - t127;
t104 = Icges(6,3) * t97 - t77 * qJD(5);
t28 = t30 - t144;
t15 = (t118 * t90 + (-rSges(6,3) - pkin(8)) * t89) * t97 - t125;
t13 = t15 - t126;
t14 = -t90 * t119 - pkin(4) * t141 + pkin(8) * t140 + (-t90 * t129 - t97 * t132) * rSges(6,1) + t139;
t12 = t14 - t127;
t91 = pkin(2) * t96;
t73 = (-t136 + t137) * qJD(5);
t55 = t61 + t91;
t54 = t60 - t145;
t47 = t53 - t122;
t46 = t52 - t123;
t45 = t51 + t91;
t44 = t50 - t145;
t42 = rSges(6,1) * t132 - t138;
t37 = Icges(6,3) * t89 + t107;
t36 = -Icges(6,3) * t90 + t112 * t89;
t33 = t35 - t122;
t32 = t34 - t123;
t27 = t29 + t91;
t26 = t28 - t145;
t19 = t104 * t89 + t97 * t107;
t18 = t104 * t90 - t112 * t141;
t11 = t13 - t122;
t10 = t12 - t123;
t9 = -t116 * t90 + t89 * t37;
t8 = t89 * t36 - t152;
t7 = -t90 * t37 - t153;
t6 = -t117 * t89 - t90 * t36;
t1 = ((-t43 + t81) * t97 + t125) * t89 + (-t80 * t90 * qJD(5) + t97 * t42 + t139) * t90;
t2 = [0; 0; (t10 * t27 + t11 * t26) * t148 + (t32 * t45 + t33 * t44) * t149 + (t46 * t55 + t47 * t54) * t150 + t111; 0; m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + m(5) * (t32 * t51 + t33 * t50 + t34 * t45 + t35 * t44) + m(4) * (t46 * t61 + t47 * t60 + t52 * t55 + t53 * t54) + t111; (t52 * t61 + t53 * t60) * t150 + (t34 * t51 + t35 * t50) * t149 + (t12 * t29 + t13 * t28) * t148 + t111; 0; m(6) * (t10 * t31 + t11 * t30 + t14 * t27 + t15 * t26) + m(5) * (t32 * t57 + t33 * t56 + t44 * t49 + t45 * t48) + t111; m(5) * (t34 * t57 + t35 * t56 + t48 * t51 + t49 * t50) + m(6) * (t12 * t31 + t13 * t30 + t14 * t29 + t15 * t28) + t111; (t48 * t57 + t49 * t56) * t149 + (t14 * t31 + t15 * t30) * t148 + t111; m(6) * t1; m(6) * ((-t26 * t90 - t27 * t89) * t73 + ((-t27 * t97 - t11) * t90 + (t26 * t97 - t10) * t89) * t80) + t110; m(6) * ((-t28 * t90 - t29 * t89) * t73 + ((-t29 * t97 - t13) * t90 + (t28 * t97 - t12) * t89) * t80) + t110; m(6) * ((-t30 * t90 - t31 * t89) * t73 + ((-t31 * t97 - t15) * t90 + (t30 * t97 - t14) * t89) * t80) + t110; ((t42 * t89 + t43 * t90) * t1 + (t89 ^ 2 + t90 ^ 2) * t80 * t73) * t148 + (-t8 * t90 + t89 * t9) * t140 + t89 * ((t89 * t18 + (t8 + t153) * t97) * t89 + (t9 * t97 + (t38 * t128 + t40 * t129) * t90 + (-t19 + (-qJD(5) * t39 + t40 * t97) * t102 + (-qJD(5) * t41 - t38 * t97) * t101) * t89) * t90) + (-t6 * t90 + t7 * t89) * t141 - t90 * ((t90 * t19 + (t7 + t152) * t97) * t90 + (t6 * t97 + (-t39 * t128 - t41 * t129) * t89 + (-t18 + (qJD(5) * t38 + t41 * t97) * t102 + (qJD(5) * t40 - t39 * t97) * t101) * t90) * t89);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
