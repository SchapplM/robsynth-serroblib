% Calculate time derivative of joint inertia matrix for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:45
% DurationCPUTime: 1.14s
% Computational Cost: add. (3037->183), mult. (2242->266), div. (0->0), fcn. (1580->8), ass. (0->106)
t96 = sin(qJ(5));
t98 = cos(qJ(5));
t124 = rSges(6,1) * t96 + rSges(6,2) * t98;
t95 = qJ(1) + pkin(8);
t92 = qJ(3) + t95;
t89 = cos(t92);
t143 = t89 * t124;
t88 = sin(t92);
t162 = t88 / 0.2e1;
t161 = rSges(5,2) - pkin(3);
t141 = pkin(2) * cos(t95) + cos(qJ(1)) * pkin(1);
t140 = Icges(6,4) * t96;
t114 = Icges(6,2) * t98 + t140;
t40 = Icges(6,6) * t89 + t114 * t88;
t139 = Icges(6,4) * t98;
t115 = Icges(6,1) * t96 + t139;
t42 = Icges(6,5) * t89 + t115 * t88;
t121 = t40 * t98 + t42 * t96;
t160 = t121 * t89;
t70 = Icges(6,5) * t98 - Icges(6,6) * t96;
t94 = qJD(1) + qJD(3);
t159 = -Icges(6,3) * t94 + t70 * qJD(5);
t71 = -Icges(6,2) * t96 + t139;
t158 = -Icges(6,6) * t94 + qJD(5) * t71;
t72 = Icges(6,1) * t98 - t140;
t157 = -Icges(6,5) * t94 + qJD(5) * t72;
t56 = t114 * qJD(5);
t57 = t115 * qJD(5);
t156 = (t71 * t96 - t72 * t98) * qJD(5) + t56 * t98 + t57 * t96 + t70 * t94;
t155 = 2 * m(4);
t154 = 2 * m(5);
t153 = 2 * m(6);
t84 = t89 * pkin(3);
t147 = t88 * rSges(6,3);
t146 = t88 * t94;
t145 = t89 * t94;
t77 = t89 * qJ(4);
t144 = qJD(4) * t88 + t94 * t77;
t142 = t88 * qJ(4) + t84;
t135 = qJD(5) * t96;
t134 = qJD(5) * t98;
t133 = -rSges(6,3) - pkin(3) - pkin(7);
t131 = rSges(6,1) * t134;
t132 = -t88 * t131 - t143 * t94;
t44 = t89 * rSges(6,3) + t124 * t88;
t130 = rSges(6,2) * t135;
t58 = t124 * qJD(5);
t128 = (t88 ^ 2 + t89 ^ 2) * t58;
t51 = t89 * rSges(4,1) - rSges(4,2) * t88;
t127 = t96 * t56 - t98 * t57;
t49 = -rSges(4,1) * t145 + rSges(4,2) * t146;
t126 = -pkin(2) * sin(t95) - sin(qJ(1)) * pkin(1);
t50 = -rSges(4,1) * t88 - rSges(4,2) * t89;
t120 = -t40 * t96 + t42 * t98;
t41 = Icges(6,6) * t88 - t114 * t89;
t43 = Icges(6,5) * t88 - t115 * t89;
t119 = t41 * t98 + t43 * t96;
t118 = t41 * t96 - t43 * t98;
t117 = t71 * t98 + t72 * t96;
t36 = t89 * rSges(5,3) + t161 * t88 + t77;
t113 = Icges(6,5) * t96 + Icges(6,6) * t98;
t37 = -rSges(5,2) * t89 + t88 * rSges(5,3) + t142;
t31 = t89 * pkin(7) + t142 + t44;
t48 = t50 * t94;
t112 = t126 * qJD(1);
t111 = t141 * qJD(1);
t110 = t119 * t88;
t109 = t115 * t94;
t108 = t114 * t94;
t107 = t113 * t94;
t103 = -t113 * qJD(5) + t117 * t94;
t106 = (-t118 / 0.2e1 - t117 * t89 / 0.2e1 + t70 * t162) * t145 + (-t119 * qJD(5) + t103 * t88 + t156 * t89 - t96 * (t88 * t108 - t158 * t89) + t98 * (t88 * t109 - t157 * t89)) * t162 + (-t121 * qJD(5) + t103 * t89 - t156 * t88 - t96 * (t89 * t108 + t158 * t88) + t98 * (t89 * t109 + t157 * t88)) * t89 / 0.2e1 - (t117 * t88 + t89 * t70 + t120) * t146 / 0.2e1;
t28 = rSges(5,3) * t145 + t161 * t146 + t144;
t30 = t133 * t88 + t143 + t77;
t75 = qJD(4) * t89;
t29 = rSges(5,2) * t145 + t75 + (-t84 + (-rSges(5,3) - qJ(4)) * t88) * t94;
t26 = t126 + t30;
t27 = t31 + t141;
t12 = (t133 * t94 - t130) * t88 - t132 + t144;
t6 = t112 + t12;
t60 = t89 * t131;
t13 = -t89 * t130 + t60 + t75 + (t133 * t89 + (-qJ(4) - t124) * t88) * t94;
t7 = -t111 + t13;
t102 = (t26 * t94 - t6) * t89 + (t27 * t94 + t7) * t88;
t101 = -t117 * qJD(5) + t127;
t100 = (t30 * t94 - t12) * t89 + (t31 * t94 + t13) * t88;
t73 = rSges(6,1) * t98 - rSges(6,2) * t96;
t47 = t51 + t141;
t46 = t126 + t50;
t45 = -t143 + t147;
t39 = Icges(6,3) * t88 - t113 * t89;
t38 = Icges(6,3) * t89 + t113 * t88;
t35 = -t111 + t49;
t34 = t48 + t112;
t33 = t37 + t141;
t32 = t126 + t36;
t21 = t89 * t107 + t159 * t88;
t20 = t88 * t107 - t159 * t89;
t17 = -t111 + t29;
t16 = t112 + t28;
t11 = -t119 * t89 + t88 * t39;
t10 = t88 * t38 - t160;
t9 = t89 * t39 + t110;
t8 = t121 * t88 + t89 * t38;
t1 = (-t94 * t44 - t60 + (rSges(6,3) * t94 + t130) * t89) * t89 + (t88 * t130 + (t147 - t45 + t143) * t94 + t132) * t88;
t2 = [-t72 * t135 - t71 * t134 + (t16 * t33 + t17 * t32) * t154 + (t26 * t7 + t27 * t6) * t153 + (t34 * t47 + t35 * t46) * t155 + t127; 0; 0; m(5) * (t16 * t37 + t17 * t36 + t28 * t33 + t29 * t32) + m(6) * (t12 * t27 + t13 * t26 + t30 * t7 + t31 * t6) + m(4) * (t34 * t51 + t35 * t50 + t46 * t49 + t47 * t48) + t101; 0; (t48 * t51 + t49 * t50) * t155 + (t28 * t37 + t29 * t36) * t154 + (t12 * t31 + t13 * t30) * t153 + t101; m(5) * ((t32 * t94 - t16) * t89 + (t33 * t94 + t17) * t88) + m(6) * t102; 0; m(5) * ((t36 * t94 - t28) * t89 + (t37 * t94 + t29) * t88) + m(6) * t100; 0; m(6) * (-(t26 * t88 - t27 * t89) * t58 + t102 * t73) + t106; m(6) * t1; m(6) * (-(t30 * t88 - t31 * t89) * t58 + t100 * t73) + t106; -m(6) * t128; ((-t44 * t88 + t45 * t89) * t1 - t73 * t128) * t153 - (t8 * t89 + t88 * t9) * t146 + t89 * ((t89 * t21 + (t9 + t160) * t94) * t89 + (-t8 * t94 + (t43 * t134 - t41 * t135) * t88 + (t120 * qJD(5) + t119 * t94 + t20) * t89) * t88) + (t10 * t89 + t11 * t88) * t145 + t88 * ((t88 * t20 + (-t10 + t110) * t94) * t88 + (t11 * t94 + (-t42 * t134 + t40 * t135) * t89 + (t118 * qJD(5) + t121 * t94 + t21) * t88) * t89);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
