% Calculate time derivative of joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:18
% EndTime: 2019-12-05 18:30:21
% DurationCPUTime: 1.54s
% Computational Cost: add. (4960->207), mult. (2858->283), div. (0->0), fcn. (1986->10), ass. (0->132)
t116 = sin(qJ(5));
t154 = qJD(5) * t116;
t118 = cos(qJ(5));
t153 = qJD(5) * t118;
t114 = qJD(1) + qJD(2);
t110 = qJD(4) + t114;
t130 = Icges(6,5) * t118 - Icges(6,6) * t116;
t163 = Icges(6,4) * t118;
t131 = -Icges(6,2) * t116 + t163;
t164 = Icges(6,4) * t116;
t132 = Icges(6,1) * t118 - t164;
t92 = Icges(6,2) * t118 + t164;
t93 = Icges(6,1) * t116 + t163;
t133 = t116 * t92 - t118 * t93;
t191 = t130 * qJD(5) + (-t116 * t132 - t118 * t131 + t133) * t110;
t188 = rSges(6,1) * t154 + rSges(6,2) * t153;
t115 = qJ(1) + qJ(2);
t111 = pkin(9) + t115;
t108 = qJ(4) + t111;
t104 = cos(t108);
t103 = sin(t108);
t42 = Icges(6,6) * t104 - t131 * t103;
t44 = Icges(6,5) * t104 - t132 * t103;
t137 = t116 * t42 - t118 * t44;
t187 = t137 * t104;
t91 = Icges(6,5) * t116 + Icges(6,6) * t118;
t186 = -Icges(6,3) * t110 + t91 * qJD(5);
t182 = 2 * m(3);
t181 = 2 * m(4);
t180 = 2 * m(5);
t179 = 2 * m(6);
t176 = -rSges(6,3) - pkin(8);
t117 = sin(qJ(1));
t175 = pkin(1) * t117;
t119 = cos(qJ(1));
t174 = pkin(1) * t119;
t112 = sin(t115);
t173 = pkin(2) * t112;
t113 = cos(t115);
t172 = pkin(2) * t113;
t158 = t104 * t110;
t159 = t103 * t110;
t52 = rSges(5,1) * t159 + rSges(5,2) * t158;
t106 = sin(t111);
t157 = t106 * t114;
t156 = t112 * t114;
t97 = pkin(2) * t156;
t171 = pkin(3) * t157 + t97;
t155 = t113 * t114;
t60 = rSges(3,1) * t156 + rSges(3,2) * t155;
t170 = rSges(6,1) * t118;
t107 = cos(t111);
t169 = rSges(4,2) * t107;
t168 = rSges(6,2) * t116;
t167 = pkin(1) * qJD(1);
t166 = t103 * rSges(6,3);
t100 = t104 * rSges(6,3);
t85 = t103 * t168;
t165 = t100 + t85;
t150 = t103 * t170;
t152 = -t188 * t104 - t110 * t150;
t86 = t104 * t168;
t151 = t188 * t103 + t110 * t86;
t48 = rSges(4,1) * t157 + t114 * t169 + t97;
t149 = t119 * t167;
t144 = -pkin(4) - t170;
t63 = -t104 * rSges(5,1) + t103 * rSges(5,2);
t34 = t171 + t52;
t76 = -t113 * rSges(3,1) + t112 * rSges(3,2);
t61 = -rSges(3,1) * t155 + rSges(3,2) * t156;
t53 = -rSges(5,1) * t158 + rSges(5,2) * t159;
t143 = -pkin(3) * t106 - t173;
t142 = -pkin(3) * t107 - t172;
t141 = -t107 * rSges(4,1) - t172;
t75 = -rSges(3,1) * t112 - rSges(3,2) * t113;
t62 = -rSges(5,1) * t103 - rSges(5,2) * t104;
t138 = t116 * t44 + t118 * t42;
t43 = Icges(6,6) * t103 + t131 * t104;
t45 = Icges(6,5) * t103 + t132 * t104;
t136 = t116 * t45 + t118 * t43;
t135 = t116 * t43 - t118 * t45;
t59 = t106 * rSges(4,2) + t141;
t129 = t142 * t114;
t128 = (t132 - t92) * t154 + (t131 + t93) * t153;
t127 = (-t135 * qJD(5) + t191 * t103) * t103 / 0.2e1 + (-t137 * qJD(5) + t191 * t104) * t104 / 0.2e1 - (t133 * t103 + t104 * t91 + t138) * t159 / 0.2e1 + (t103 * t91 - t133 * t104 + t136) * t158 / 0.2e1;
t126 = t135 * t103;
t123 = t130 * t110;
t58 = -rSges(4,1) * t106 - t169 - t173;
t49 = rSges(4,2) * t157 + t141 * t114;
t51 = t142 + t63;
t32 = t104 * pkin(8) + t144 * t103 + t165;
t120 = t176 * t103 + t144 * t104;
t50 = t143 + t62;
t33 = t86 + t120;
t35 = t129 + t53;
t14 = pkin(4) * t159 + (t176 * t104 - t85) * t110 - t152;
t28 = t143 + t32;
t15 = t120 * t110 + t151;
t29 = t142 + t33;
t12 = t14 + t171;
t13 = t129 + t15;
t109 = t117 * t167;
t98 = rSges(6,1) * t116 + rSges(6,2) * t118;
t84 = (-t168 + t170) * qJD(5);
t65 = t76 - t174;
t64 = t75 - t175;
t57 = t61 - t149;
t56 = t109 + t60;
t55 = t59 - t174;
t54 = t58 - t175;
t47 = t104 * t170 + t166 - t86;
t46 = -t150 + t165;
t41 = Icges(6,3) * t103 + t130 * t104;
t40 = Icges(6,3) * t104 - t130 * t103;
t39 = t51 - t174;
t38 = t50 - t175;
t37 = t49 - t149;
t36 = t109 + t48;
t31 = t35 - t149;
t30 = t109 + t34;
t27 = t29 - t174;
t26 = t28 - t175;
t19 = t186 * t103 - t104 * t123;
t18 = -t103 * t123 - t186 * t104;
t11 = t13 - t149;
t10 = t109 + t12;
t9 = t103 * t41 - t135 * t104;
t8 = t103 * t40 - t187;
t7 = t104 * t41 + t126;
t6 = t137 * t103 + t104 * t40;
t1 = t104 * t152 - t103 * t151 + ((-t46 + t100) * t104 + (t166 - t47 + (t168 + t170) * t104) * t103) * t110;
t2 = [(t56 * t65 + t57 * t64) * t182 + (t36 * t55 + t37 * t54) * t181 + (t30 * t39 + t31 * t38) * t180 + (t10 * t27 + t11 * t26) * t179 + t128; m(3) * (t56 * t76 + t57 * t75 + t60 * t65 + t61 * t64) + m(4) * (t36 * t59 + t37 * t58 + t48 * t55 + t49 * t54) + m(5) * (t30 * t51 + t31 * t50 + t34 * t39 + t35 * t38) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + t128; (t12 * t29 + t13 * t28) * t179 + (t34 * t51 + t35 * t50) * t180 + (t48 * t59 + t49 * t58) * t181 + (t60 * t76 + t61 * t75) * t182 + t128; 0; 0; 0; m(5) * (t30 * t63 + t31 * t62 + t38 * t53 + t39 * t52) + m(6) * (t10 * t33 + t11 * t32 + t14 * t27 + t15 * t26) + t128; m(6) * (t12 * t33 + t13 * t32 + t14 * t29 + t15 * t28) + m(5) * (t34 * t63 + t35 * t62 + t50 * t53 + t51 * t52) + t128; 0; (t52 * t63 + t53 * t62) * t180 + (t14 * t33 + t15 * t32) * t179 + t128; m(6) * ((t103 * t27 - t104 * t26) * t84 + ((t110 * t27 - t11) * t104 + (t110 * t26 + t10) * t103) * t98) + t127; m(6) * ((t103 * t29 - t104 * t28) * t84 + ((t110 * t29 - t13) * t104 + (t110 * t28 + t12) * t103) * t98) + t127; m(6) * t1; m(6) * ((t103 * t33 - t104 * t32) * t84 + ((t110 * t33 - t15) * t104 + (t110 * t32 + t14) * t103) * t98) + t127; ((-t103 * t46 + t104 * t47) * t1 + (t103 ^ 2 + t104 ^ 2) * t98 * t84) * t179 - (t103 * t7 + t104 * t6) * t159 + t104 * ((t104 * t19 + (t7 + t187) * t110) * t104 + (-t6 * t110 + (t43 * t153 + t45 * t154) * t103 + (t138 * qJD(5) + t135 * t110 + t18) * t104) * t103) + (t103 * t9 + t104 * t8) * t158 + t103 * ((t103 * t18 + (-t8 + t126) * t110) * t103 + (t9 * t110 + (-t42 * t153 - t44 * t154) * t104 + (-t136 * qJD(5) + t137 * t110 + t19) * t103) * t104);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
