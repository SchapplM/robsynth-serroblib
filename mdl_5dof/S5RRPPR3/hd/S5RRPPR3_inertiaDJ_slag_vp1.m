% Calculate time derivative of joint inertia matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:31
% DurationCPUTime: 1.26s
% Computational Cost: add. (3211->197), mult. (2356->280), div. (0->0), fcn. (1646->8), ass. (0->121)
t106 = sin(qJ(5));
t108 = cos(qJ(5));
t131 = rSges(6,1) * t106 + rSges(6,2) * t108;
t105 = qJ(1) + qJ(2);
t100 = pkin(8) + t105;
t98 = cos(t100);
t158 = t98 * t131;
t97 = sin(t100);
t177 = t97 / 0.2e1;
t101 = sin(t105);
t102 = cos(t105);
t62 = t102 * rSges(3,1) - rSges(3,2) * t101;
t142 = qJD(5) * t108;
t143 = qJD(5) * t106;
t176 = -rSges(6,1) * t142 + rSges(6,2) * t143;
t99 = pkin(2) * t102;
t175 = -t98 * rSges(4,1) - t99;
t174 = -t98 * pkin(3) - t99;
t147 = Icges(6,4) * t106;
t121 = Icges(6,2) * t108 + t147;
t42 = Icges(6,6) * t98 + t121 * t97;
t146 = Icges(6,4) * t108;
t122 = Icges(6,1) * t106 + t146;
t44 = Icges(6,5) * t98 + t122 * t97;
t128 = t106 * t44 + t108 * t42;
t173 = t128 * t98;
t120 = Icges(6,5) * t106 + Icges(6,6) * t108;
t172 = -Icges(6,3) * t97 + t120 * t98;
t171 = -Icges(6,6) * t97 + t121 * t98;
t170 = -Icges(6,5) * t97 + t122 * t98;
t104 = qJD(1) + qJD(2);
t64 = t121 * qJD(5);
t65 = t122 * qJD(5);
t78 = Icges(6,5) * t108 - Icges(6,6) * t106;
t79 = -Icges(6,2) * t106 + t146;
t80 = Icges(6,1) * t108 - t147;
t169 = (t106 * t79 - t108 * t80) * qJD(5) + t104 * t78 + t106 * t65 + t108 * t64;
t168 = 2 * m(3);
t167 = 2 * m(4);
t166 = 2 * m(5);
t165 = 2 * m(6);
t162 = pkin(2) * t101;
t107 = sin(qJ(1));
t161 = t107 * pkin(1);
t160 = t97 * rSges(6,3);
t89 = t98 * rSges(6,3);
t86 = t98 * qJ(4);
t159 = qJD(4) * t97 + t104 * t86;
t151 = pkin(1) * qJD(1);
t40 = Icges(6,3) * t98 + t120 * t97;
t150 = t104 * t40;
t149 = t104 * t97;
t148 = t104 * t98;
t145 = qJD(5) * t97;
t144 = qJD(5) * t98;
t141 = -rSges(6,3) - pkin(3) - pkin(7);
t46 = t131 * t97 + t89;
t140 = t97 * qJ(4) - t174;
t139 = t107 * t151;
t109 = cos(qJ(1));
t138 = t109 * t151;
t66 = t131 * qJD(5);
t134 = (t97 ^ 2 + t98 ^ 2) * t66;
t133 = t106 * t64 - t108 * t65;
t53 = -rSges(4,2) * t97 - t175;
t55 = t62 * t104;
t132 = -pkin(3) * t97 - t162;
t61 = -rSges(3,1) * t101 - rSges(3,2) * t102;
t127 = -t106 * t42 + t108 * t44;
t126 = -t106 * t170 - t108 * t171;
t125 = -t106 * t171 + t108 * t170;
t124 = t106 * t80 + t108 * t79;
t119 = t176 * t98;
t37 = -rSges(5,2) * t98 + t97 * rSges(5,3) + t140;
t31 = t98 * pkin(7) + t140 + t46;
t118 = t126 * t97;
t113 = -t120 * qJD(5) + t124 * t104;
t117 = (-t125 / 0.2e1 - t124 * t98 / 0.2e1 + t78 * t177) * t148 + (-t126 * qJD(5) - t106 * (t42 * t104 - t79 * t144) + t108 * (t44 * t104 - t80 * t144) + t113 * t97 + t169 * t98) * t177 + (-t128 * qJD(5) - t106 * (t171 * t104 + t79 * t145) + t108 * (t170 * t104 + t80 * t145) + t113 * t98 - t169 * t97) * t98 / 0.2e1 - (t124 * t97 + t98 * t78 + t127) * t149 / 0.2e1;
t54 = t61 * t104;
t116 = t141 * t97 - t162;
t52 = -rSges(4,1) * t97 - rSges(4,2) * t98 - t162;
t115 = t158 * t104 - t176 * t97;
t36 = t97 * rSges(5,2) + t98 * rSges(5,3) + t132 + t86;
t39 = rSges(4,2) * t149 + t175 * t104;
t38 = t52 * t104;
t30 = t86 + t116 + t158;
t26 = rSges(5,2) * t149 + rSges(5,3) * t148 + t132 * t104 + t159;
t112 = -t124 * qJD(5) + t133;
t28 = t30 - t161;
t103 = t109 * pkin(1);
t29 = t103 + t31;
t12 = t116 * t104 + t115 + t159;
t6 = t12 - t139;
t84 = qJD(4) * t98;
t13 = t84 + (-t99 + t141 * t98 + (-qJ(4) - t131) * t97) * t104 - t119;
t7 = t13 - t138;
t111 = -t98 * t6 + t97 * t7 + (t28 * t98 + t29 * t97) * t104;
t110 = -t98 * t12 + t97 * t13 + (t30 * t98 + t31 * t97) * t104;
t27 = rSges(5,2) * t148 + t84 + ((-rSges(5,3) - qJ(4)) * t97 + t174) * t104;
t82 = rSges(6,1) * t108 - rSges(6,2) * t106;
t57 = t103 + t62;
t56 = t61 - t161;
t51 = -t55 - t138;
t50 = t54 - t139;
t49 = t103 + t53;
t48 = t52 - t161;
t47 = -t158 + t160;
t35 = t39 - t138;
t34 = t38 - t139;
t33 = t103 + t37;
t32 = t36 - t161;
t21 = t172 * t104 + t78 * t145;
t20 = -t78 * t144 + t150;
t17 = t27 - t138;
t16 = t26 - t139;
t11 = -t126 * t98 - t172 * t97;
t10 = t97 * t40 - t173;
t9 = -t172 * t98 + t118;
t8 = t128 * t97 + t98 * t40;
t1 = t98 * t119 - t97 * t115 + ((-t46 + t89) * t98 + (t160 - t47 + t158) * t97) * t104;
t2 = [(t28 * t7 + t29 * t6) * t165 - t80 * t143 - t79 * t142 + (t16 * t33 + t17 * t32) * t166 + (t50 * t57 + t51 * t56) * t168 + (t34 * t49 + t35 * t48) * t167 + t133; m(6) * (t12 * t29 + t13 * t28 + t30 * t7 + t31 * t6) + m(5) * (t16 * t37 + t17 * t36 + t26 * t33 + t27 * t32) + m(3) * (t50 * t62 + t51 * t61 + t54 * t57 - t55 * t56) + m(4) * (t34 * t53 + t35 * t52 + t38 * t49 + t39 * t48) + t112; (t26 * t37 + t27 * t36) * t166 + (t12 * t31 + t13 * t30) * t165 + (t38 * t53 + t39 * t52) * t167 + (t54 * t62 - t55 * t61) * t168 + t112; 0; 0; 0; m(6) * t111 + m(5) * (-t98 * t16 + t97 * t17 + (t32 * t98 + t33 * t97) * t104); m(5) * (-t98 * t26 + t97 * t27 + (t36 * t98 + t37 * t97) * t104) + m(6) * t110; 0; 0; m(6) * (-(t28 * t97 - t29 * t98) * t66 + t111 * t82) + t117; m(6) * (-(t30 * t97 - t31 * t98) * t66 + t110 * t82) + t117; m(6) * t1; -m(6) * t134; ((-t46 * t97 + t47 * t98) * t1 - t82 * t134) * t165 - (t8 * t98 + t9 * t97) * t149 + t98 * ((t98 * t21 + (t9 + t173) * t104) * t98 + (-t8 * t104 + (-t142 * t170 + t143 * t171) * t97 + (t20 + t127 * qJD(5) + (t126 - t40) * t104) * t98) * t97) + (t10 * t98 + t11 * t97) * t148 + t97 * ((t97 * t20 + (-t10 + t118) * t104) * t97 + (t11 * t104 + (-t44 * t142 + t42 * t143 + t150) * t98 + (t125 * qJD(5) + t128 * t104 + t21) * t97) * t98);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
