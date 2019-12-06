% Calculate time derivative of joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:23
% DurationCPUTime: 1.29s
% Computational Cost: add. (4786->194), mult. (2744->272), div. (0->0), fcn. (1920->10), ass. (0->113)
t104 = sin(qJ(5));
t140 = qJD(5) * t104;
t106 = cos(qJ(5));
t139 = qJD(5) * t106;
t103 = qJ(1) + pkin(9);
t172 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t103);
t171 = rSges(6,1) * t140 + rSges(6,2) * t139;
t143 = Icges(6,4) * t106;
t117 = -Icges(6,2) * t104 + t143;
t101 = qJ(3) + t103;
t96 = qJ(4) + t101;
t92 = sin(t96);
t93 = cos(t96);
t42 = Icges(6,6) * t93 - t117 * t92;
t144 = Icges(6,4) * t104;
t118 = Icges(6,1) * t106 - t144;
t44 = Icges(6,5) * t93 - t118 * t92;
t123 = t104 * t42 - t106 * t44;
t170 = t123 * t93;
t116 = Icges(6,5) * t106 - Icges(6,6) * t104;
t41 = Icges(6,3) * t92 + t116 * t93;
t43 = Icges(6,6) * t92 + t117 * t93;
t45 = Icges(6,5) * t92 + t118 * t93;
t168 = 2 * m(4);
t167 = 2 * m(5);
t166 = 2 * m(6);
t94 = sin(t101);
t162 = pkin(3) * t94;
t95 = cos(t101);
t161 = pkin(3) * t95;
t160 = -rSges(6,3) - pkin(8);
t157 = t92 * rSges(6,3);
t88 = t93 * rSges(6,3);
t102 = qJD(1) + qJD(3);
t100 = qJD(4) + t102;
t147 = t100 * t93;
t148 = t100 * t92;
t48 = rSges(5,1) * t148 + rSges(5,2) * t147;
t153 = rSges(6,2) * t104;
t77 = t92 * t153;
t156 = t77 + t88;
t145 = t102 * t95;
t146 = t102 * t94;
t54 = rSges(4,1) * t146 + rSges(4,2) * t145;
t155 = t172 * qJD(1);
t154 = rSges(6,1) * t106;
t40 = Icges(6,3) * t93 - t116 * t92;
t149 = t100 * t40;
t142 = qJD(5) * t92;
t141 = qJD(5) * t93;
t138 = pkin(3) * t145;
t135 = t92 * t154;
t137 = -t100 * t135 - t171 * t93;
t78 = t93 * t153;
t136 = t100 * t78 + t171 * t92;
t82 = pkin(3) * t146;
t34 = t82 + t48;
t61 = -rSges(4,1) * t95 + rSges(4,2) * t94;
t57 = -rSges(5,1) * t93 + rSges(5,2) * t92;
t130 = -pkin(4) - t154;
t55 = -rSges(4,1) * t145 + rSges(4,2) * t146;
t49 = -rSges(5,1) * t147 + rSges(5,2) * t148;
t128 = -pkin(2) * cos(t103) - cos(qJ(1)) * pkin(1);
t60 = -rSges(4,1) * t94 - rSges(4,2) * t95;
t56 = -rSges(5,1) * t92 - rSges(5,2) * t93;
t124 = t104 * t44 + t106 * t42;
t122 = t104 * t45 + t106 * t43;
t121 = t104 * t43 - t106 * t45;
t84 = Icges(6,2) * t106 + t144;
t85 = Icges(6,1) * t104 + t143;
t119 = t104 * t84 - t106 * t85;
t83 = Icges(6,5) * t104 + Icges(6,6) * t106;
t53 = t57 - t161;
t115 = t128 * qJD(1);
t114 = (t118 - t84) * t140 + (t117 + t85) * t139;
t113 = t121 * t92;
t109 = qJD(5) * t116 + t100 * t119;
t112 = (-qJD(5) * t121 + t104 * (t100 * t44 - t141 * t85) + t106 * (t100 * t42 - t141 * t84) + t109 * t92) * t92 / 0.2e1 + (-qJD(5) * t123 + t104 * (-t100 * t45 + t142 * t85) + t106 * (-t100 * t43 + t142 * t84) + t109 * t93) * t93 / 0.2e1 - (t119 * t92 + t83 * t93 + t124) * t148 / 0.2e1 + (-t119 * t93 + t83 * t92 + t122) * t147 / 0.2e1;
t52 = t56 - t162;
t35 = t49 - t138;
t32 = pkin(8) * t93 + t130 * t92 + t156;
t108 = t130 * t93 + t160 * t92;
t33 = t78 + t108;
t28 = t32 - t162;
t29 = t33 - t161;
t14 = pkin(4) * t148 + (t160 * t93 - t77) * t100 - t137;
t12 = t82 + t14;
t15 = t100 * t108 + t136;
t13 = t15 - t138;
t86 = rSges(6,1) * t104 + rSges(6,2) * t106;
t76 = (-t153 + t154) * qJD(5);
t51 = t128 + t61;
t50 = -t172 + t60;
t47 = t154 * t93 + t157 - t78;
t46 = -t135 + t156;
t39 = t115 + t55;
t38 = t155 + t54;
t37 = t128 + t53;
t36 = -t172 + t52;
t31 = t115 + t35;
t30 = t34 + t155;
t27 = t128 + t29;
t26 = -t172 + t28;
t19 = -t100 * t41 + t142 * t83;
t18 = -t141 * t83 + t149;
t11 = t115 + t13;
t10 = t12 + t155;
t9 = -t121 * t93 + t41 * t92;
t8 = t40 * t92 - t170;
t7 = t41 * t93 + t113;
t6 = t123 * t92 + t40 * t93;
t1 = t93 * t137 - t92 * t136 + ((-t46 + t88) * t93 + (t157 - t47 + (t153 + t154) * t93) * t92) * t100;
t2 = [(t38 * t51 + t39 * t50) * t168 + (t30 * t37 + t31 * t36) * t167 + (t10 * t27 + t11 * t26) * t166 + t114; 0; 0; m(4) * (t38 * t61 + t39 * t60 + t50 * t55 + t51 * t54) + m(5) * (t30 * t53 + t31 * t52 + t34 * t37 + t35 * t36) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + t114; 0; (t54 * t61 + t55 * t60) * t168 + (t34 * t53 + t35 * t52) * t167 + (t12 * t29 + t13 * t28) * t166 + t114; m(5) * (t30 * t57 + t31 * t56 + t36 * t49 + t37 * t48) + m(6) * (t10 * t33 + t11 * t32 + t14 * t27 + t15 * t26) + t114; 0; m(5) * (t34 * t57 + t35 * t56 + t48 * t53 + t49 * t52) + m(6) * (t12 * t33 + t13 * t32 + t14 * t29 + t15 * t28) + t114; (t48 * t57 + t49 * t56) * t167 + (t14 * t33 + t15 * t32) * t166 + t114; m(6) * ((-t26 * t93 + t27 * t92) * t76 + (t10 * t92 - t11 * t93 + (t26 * t92 + t27 * t93) * t100) * t86) + t112; m(6) * t1; m(6) * ((-t28 * t93 + t29 * t92) * t76 + (t12 * t92 - t13 * t93 + (t28 * t92 + t29 * t93) * t100) * t86) + t112; m(6) * ((-t32 * t93 + t33 * t92) * t76 + (t14 * t92 - t15 * t93 + (t32 * t92 + t33 * t93) * t100) * t86) + t112; ((-t46 * t92 + t47 * t93) * t1 + (t92 ^ 2 + t93 ^ 2) * t86 * t76) * t166 - (t6 * t93 + t7 * t92) * t148 + t93 * ((t93 * t19 + (t7 + t170) * t100) * t93 + (-t6 * t100 + (t139 * t43 + t140 * t45) * t92 + (t18 + t124 * qJD(5) + (t121 - t40) * t100) * t93) * t92) + (t8 * t93 + t9 * t92) * t147 + t92 * ((t92 * t18 + (-t8 + t113) * t100) * t92 + (t9 * t100 + (-t139 * t42 - t140 * t44 + t149) * t93 + (-t122 * qJD(5) + t100 * t123 + t19) * t92) * t93);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
