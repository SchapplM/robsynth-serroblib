% Calculate time derivative of joint inertia matrix for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:29
% DurationCPUTime: 1.51s
% Computational Cost: add. (2615->171), mult. (2164->247), div. (0->0), fcn. (1564->8), ass. (0->106)
t98 = cos(pkin(7));
t166 = rSges(4,1) * t98 + pkin(2) - rSges(4,2) * sin(pkin(7));
t96 = qJ(1) + qJ(2);
t91 = sin(t96);
t92 = cos(t96);
t34 = (rSges(4,3) + qJ(3)) * t91 + t166 * t92;
t95 = qJD(1) + qJD(2);
t165 = t34 * t95;
t94 = pkin(7) + qJ(4);
t89 = sin(t94);
t135 = qJD(4) * t89;
t90 = cos(t94);
t134 = qJD(4) * t90;
t61 = rSges(5,1) * t89 + rSges(5,2) * t90;
t110 = t61 * qJD(4);
t161 = t92 * qJ(3) - t166 * t91;
t137 = Icges(5,4) * t90;
t113 = -Icges(5,2) * t89 + t137;
t108 = t113 * t92;
t38 = Icges(5,6) * t91 + t108;
t138 = Icges(5,4) * t89;
t114 = Icges(5,1) * t90 - t138;
t109 = t114 * t92;
t40 = Icges(5,5) * t91 + t109;
t116 = t38 * t89 - t40 * t90;
t160 = t116 * t91;
t37 = -Icges(5,6) * t92 + t113 * t91;
t39 = -Icges(5,5) * t92 + t114 * t91;
t117 = t37 * t89 - t39 * t90;
t159 = t117 * t92;
t112 = Icges(5,5) * t90 - Icges(5,6) * t89;
t59 = Icges(5,2) * t90 + t138;
t60 = Icges(5,1) * t89 + t137;
t115 = t59 * t89 - t60 * t90;
t158 = t112 * qJD(4) + t115 * t95;
t157 = 2 * m(3);
t156 = 2 * m(4);
t155 = 2 * m(5);
t150 = rSges(5,1) * t90;
t148 = rSges(5,2) * t89;
t100 = sin(qJ(1));
t146 = pkin(1) * t100;
t82 = t91 * rSges(5,3);
t143 = t91 * t95;
t142 = t92 * t95;
t99 = -pkin(6) - qJ(3);
t141 = t95 * t99;
t70 = t91 * t148;
t140 = rSges(5,3) * t142 + t95 * t70;
t139 = t92 * rSges(5,3) + t70;
t136 = pkin(1) * qJD(1);
t131 = t92 * t148;
t130 = -t91 * t110 - t95 * t131;
t129 = t100 * t136;
t101 = cos(qJ(1));
t128 = t101 * t136;
t87 = pkin(3) * t98 + pkin(2);
t123 = -t87 - t150;
t118 = t123 * t91;
t29 = -t92 * t99 + t118 + t139;
t27 = t29 - t146;
t78 = qJD(3) * t91;
t11 = t78 + t95 * t118 + (-t110 - t141) * t92 + t140;
t9 = t11 - t129;
t124 = t27 * t95 - t9;
t63 = t92 * rSges(3,1) - rSges(3,2) * t91;
t79 = qJD(3) * t92;
t12 = t91 * t141 + t79 + (t123 * t92 - t82) * t95 - t130;
t10 = t12 - t128;
t42 = t92 * t150 - t131 + t82;
t30 = t92 * t87 - t91 * t99 + t42;
t93 = t101 * pkin(1);
t28 = t30 + t93;
t122 = t28 * t95 + t10;
t121 = t29 * t95 - t11;
t120 = t30 * t95 + t12;
t49 = -rSges(3,1) * t142 + rSges(3,2) * t143;
t62 = -rSges(3,1) * t91 - rSges(3,2) * t92;
t58 = Icges(5,5) * t89 + Icges(5,6) * t90;
t48 = t62 * t95;
t111 = (t114 - t59) * t135 + (t113 + t60) * t134;
t107 = t112 * t92;
t106 = (-qJD(4) * t116 + t158 * t91 + (-t113 * t90 - t114 * t89) * t143) * t91 / 0.2e1 - (-qJD(4) * t117 - t158 * t92 + (t108 * t90 + t109 * t89) * t95) * t92 / 0.2e1 + (-t115 * t91 + t37 * t90 + t39 * t89 - t58 * t92) * t143 / 0.2e1 + (-t115 * t92 + t38 * t90 + t40 * t89 + t58 * t91) * t142 / 0.2e1;
t33 = t92 * rSges(4,3) + t161;
t25 = rSges(4,3) * t142 + t161 * t95 + t78;
t103 = Icges(5,3) * t95 - t58 * qJD(4);
t26 = t79 - t165;
t55 = (-t148 + t150) * qJD(4);
t51 = t63 + t93;
t50 = t62 - t146;
t44 = t49 - t128;
t43 = t48 - t129;
t41 = t91 * t150 - t139;
t36 = Icges(5,3) * t91 + t107;
t35 = -Icges(5,3) * t92 + t112 * t91;
t32 = t34 + t93;
t31 = t33 - t146;
t20 = t103 * t91 + t107 * t95;
t19 = t103 * t92 - t112 * t143;
t18 = t26 - t128;
t17 = t25 - t129;
t8 = -t116 * t92 + t91 * t36;
t7 = t91 * t35 - t159;
t6 = -t36 * t92 - t160;
t5 = -t117 * t91 - t35 * t92;
t1 = [(t43 * t51 + t44 * t50) * t157 + (t17 * t32 + t18 * t31) * t156 + (t10 * t27 + t28 * t9) * t155 + t111; m(3) * (t43 * t63 + t44 * t62 + t48 * t51 + t49 * t50) + m(4) * (t17 * t34 + t18 * t33 + t25 * t32 + t26 * t31) + m(5) * (t10 * t29 + t11 * t28 + t12 * t27 + t30 * t9) + t111; (t48 * t63 + t49 * t62) * t157 + (t25 * t34 + t26 * t33) * t156 + (t11 * t30 + t12 * t29) * t155 + t111; m(4) * ((t31 * t95 - t17) * t92 + (t32 * t95 + t18) * t91) + m(5) * (t122 * t91 + t124 * t92); m(4) * ((t33 * t95 - t25) * t92 + (t26 + t165) * t91) + m(5) * (t120 * t91 + t121 * t92); 0; m(5) * ((-t27 * t92 - t28 * t91) * t55 + (-t122 * t92 + t124 * t91) * t61) + t106; m(5) * ((-t29 * t92 - t30 * t91) * t55 + (-t120 * t92 + t121 * t91) * t61) + t106; 0; ((t41 * t91 + t42 * t92) * (((-t42 + t82) * t95 + t130) * t91 + (-t110 * t92 + t95 * t41 + t140) * t92) + (t91 ^ 2 + t92 ^ 2) * t61 * t55) * t155 + (-t7 * t92 + t8 * t91) * t142 + t91 * ((t91 * t19 + (t7 + t160) * t95) * t91 + (t8 * t95 + (t134 * t37 + t135 * t39) * t92 + (-t20 + (-qJD(4) * t38 + t39 * t95) * t90 + (-qJD(4) * t40 - t37 * t95) * t89) * t91) * t92) + (-t5 * t92 + t6 * t91) * t143 - t92 * ((t92 * t20 + (t6 + t159) * t95) * t92 + (t5 * t95 + (-t134 * t38 - t135 * t40) * t91 + (-t19 + (qJD(4) * t37 + t40 * t95) * t90 + (qJD(4) * t39 - t38 * t95) * t89) * t92) * t91);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
