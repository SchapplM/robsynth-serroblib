% Calculate time derivative of joint inertia matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:24
% DurationCPUTime: 2.33s
% Computational Cost: add. (4087->270), mult. (6757->455), div. (0->0), fcn. (6351->8), ass. (0->147)
t121 = sin(qJ(2));
t122 = cos(qJ(2));
t117 = sin(pkin(6));
t113 = t117 ^ 2;
t119 = cos(pkin(6));
t114 = t119 ^ 2;
t196 = t113 + t114;
t194 = qJD(2) * t196;
t200 = m(3) * (t121 * rSges(3,1) + rSges(3,2) * t122) * t194;
t199 = 0.1e1 - t196;
t116 = sin(pkin(7));
t118 = cos(pkin(7));
t155 = rSges(4,1) * t118 - rSges(4,2) * t116;
t195 = rSges(4,3) * t122 - t121 * t155;
t120 = -pkin(5) - qJ(3);
t166 = qJ(3) + t120;
t110 = pkin(3) * t118 + pkin(2);
t185 = pkin(2) - t110;
t193 = t121 * t185 - t122 * t166;
t190 = 2 * m(5);
t189 = m(4) / 0.2e1;
t188 = m(5) / 0.2e1;
t187 = t117 / 0.2e1;
t186 = t119 / 0.2e1;
t108 = t121 * pkin(2) - qJ(3) * t122;
t184 = t196 * (-qJD(2) * t108 + qJD(3) * t121);
t153 = pkin(2) * t122 + qJ(3) * t121;
t183 = t196 * t153;
t168 = t119 * t121;
t170 = t117 * t121;
t115 = pkin(7) + qJ(4);
t111 = sin(t115);
t112 = cos(t115);
t169 = t117 * t122;
t89 = -t111 * t169 - t119 * t112;
t90 = -t119 * t111 + t112 * t169;
t42 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t170;
t44 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t170;
t46 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t170;
t167 = t119 * t122;
t91 = -t111 * t167 + t117 * t112;
t92 = t117 * t111 + t112 * t167;
t18 = t168 * t42 + t44 * t91 + t46 * t92;
t181 = t117 * t18;
t43 = Icges(5,5) * t92 + Icges(5,6) * t91 + Icges(5,3) * t168;
t45 = Icges(5,4) * t92 + Icges(5,2) * t91 + Icges(5,6) * t168;
t47 = Icges(5,1) * t92 + Icges(5,4) * t91 + Icges(5,5) * t168;
t17 = t170 * t43 + t45 * t89 + t47 * t90;
t180 = t119 * t17;
t139 = Icges(5,5) * t112 - Icges(5,6) * t111;
t163 = qJD(4) * t121;
t179 = t122 * ((-Icges(5,5) * t111 - Icges(5,6) * t112) * t163 + (Icges(5,3) * t121 + t122 * t139) * qJD(2));
t80 = -Icges(5,3) * t122 + t121 * t139;
t178 = t122 * t80;
t102 = qJD(2) * t153 - qJD(3) * t122;
t177 = -t102 - (rSges(4,3) * t121 + t122 * t155) * qJD(2);
t176 = -t108 + t195;
t173 = Icges(5,4) * t111;
t172 = Icges(5,4) * t112;
t165 = qJD(2) * t121;
t164 = qJD(2) * t122;
t154 = rSges(5,1) * t112 - rSges(5,2) * t111;
t53 = (-rSges(5,1) * t111 - rSges(5,2) * t112) * t163 + (rSges(5,3) * t121 + t122 * t154) * qJD(2);
t162 = -t102 - t53 - (-t121 * t166 - t122 * t185) * qJD(2);
t84 = -rSges(5,3) * t122 + t121 * t154;
t161 = -t108 + t193 - t84;
t160 = t117 * t165;
t159 = t117 * t164;
t158 = t119 * t165;
t157 = t119 * t164;
t152 = -t111 * t44 + t112 * t46;
t151 = -t111 * t45 + t112 * t47;
t140 = -Icges(5,2) * t111 + t172;
t81 = -Icges(5,6) * t122 + t121 * t140;
t142 = Icges(5,1) * t112 - t173;
t82 = -Icges(5,5) * t122 + t121 * t142;
t150 = t111 * t81 - t112 * t82;
t20 = t121 * t152 - t122 * t42;
t21 = t121 * t151 - t122 * t43;
t148 = t20 * t117 + t21 * t119;
t48 = rSges(5,1) * t90 + rSges(5,2) * t89 + rSges(5,3) * t170;
t49 = rSges(5,1) * t92 + rSges(5,2) * t91 + rSges(5,3) * t168;
t147 = -t117 * t49 + t119 * t48;
t65 = -qJD(4) * t90 + t111 * t160;
t66 = qJD(4) * t89 - t112 * t160;
t32 = Icges(5,5) * t66 + Icges(5,6) * t65 + Icges(5,3) * t159;
t138 = t121 * t32 + t164 * t42;
t67 = -qJD(4) * t92 + t111 * t158;
t68 = qJD(4) * t91 - t112 * t158;
t33 = Icges(5,5) * t68 + Icges(5,6) * t67 + Icges(5,3) * t157;
t137 = t121 * t33 + t164 * t43;
t130 = qJD(2) * (-Icges(3,5) * t121 - Icges(3,6) * t122);
t126 = t110 * t122 - t120 * t121 - t153;
t125 = qJD(2) * (Icges(4,5) * t122 + (-Icges(4,1) * t118 + Icges(4,4) * t116) * t121);
t124 = qJD(2) * (Icges(4,6) * t122 + (-Icges(4,4) * t118 + Icges(4,2) * t116) * t121);
t106 = t117 * t116 + t118 * t167;
t105 = -t116 * t167 + t117 * t118;
t104 = -t119 * t116 + t118 * t169;
t103 = -t116 * t169 - t119 * t118;
t97 = t119 * t130;
t96 = t117 * t130;
t74 = t119 * t125;
t73 = t117 * t125;
t72 = t119 * t124;
t71 = t117 * t124;
t64 = t176 * t119;
t63 = t176 * t117;
t55 = t177 * t119;
t54 = t177 * t117;
t52 = (-Icges(5,1) * t111 - t172) * t163 + (Icges(5,5) * t121 + t122 * t142) * qJD(2);
t51 = (-Icges(5,2) * t112 - t173) * t163 + (Icges(5,6) * t121 + t122 * t140) * qJD(2);
t41 = t161 * t119;
t40 = t161 * t117;
t39 = t68 * rSges(5,1) + t67 * rSges(5,2) + rSges(5,3) * t157;
t38 = t66 * rSges(5,1) + t65 * rSges(5,2) + rSges(5,3) * t159;
t37 = Icges(5,1) * t68 + Icges(5,4) * t67 + Icges(5,5) * t157;
t36 = Icges(5,1) * t66 + Icges(5,4) * t65 + Icges(5,5) * t159;
t35 = Icges(5,4) * t68 + Icges(5,2) * t67 + Icges(5,6) * t157;
t34 = Icges(5,4) * t66 + Icges(5,2) * t65 + Icges(5,6) * t159;
t31 = -t122 * t49 - t168 * t84;
t30 = t122 * t48 + t170 * t84;
t29 = -t121 * t150 - t178;
t28 = t162 * t119;
t27 = t162 * t117;
t26 = t147 * t121;
t25 = t168 * t80 + t81 * t91 + t82 * t92;
t24 = t170 * t80 + t81 * t89 + t82 * t90;
t23 = t194 * t195 + t184;
t22 = t117 * (rSges(4,1) * t104 + rSges(4,2) * t103 + rSges(4,3) * t170) + t119 * (rSges(4,1) * t106 + rSges(4,2) * t105 + rSges(4,3) * t168) + t183;
t19 = t168 * t43 + t45 * t91 + t47 * t92;
t16 = t170 * t42 + t44 * t89 + t46 * t90;
t15 = -t53 * t168 - t122 * t39 + (t121 * t49 - t167 * t84) * qJD(2);
t14 = t53 * t170 + t122 * t38 + (-t121 * t48 + t169 * t84) * qJD(2);
t13 = (t119 * t126 + t49) * t119 + (t117 * t126 + t48) * t117 + t183;
t12 = (-t117 * t39 + t119 * t38) * t121 + t147 * t164;
t11 = t117 * t38 + t119 * t39 + t193 * t194 + t184;
t10 = t119 * t137 + t91 * t35 + t92 * t37 + t67 * t45 + t68 * t47;
t9 = t119 * t138 + t91 * t34 + t92 * t36 + t67 * t44 + t68 * t46;
t8 = t117 * t137 + t89 * t35 + t90 * t37 + t65 * t45 + t66 * t47;
t7 = t117 * t138 + t89 * t34 + t90 * t36 + t65 * t44 + t66 * t46;
t6 = (qJD(2) * t151 - t33) * t122 + (qJD(2) * t43 - t111 * t35 + t112 * t37 + (-t111 * t47 - t112 * t45) * qJD(4)) * t121;
t5 = (qJD(2) * t152 - t32) * t122 + (qJD(2) * t42 - t111 * t34 + t112 * t36 + (-t111 * t46 - t112 * t44) * qJD(4)) * t121;
t4 = t10 * t117 - t119 * t9;
t3 = t117 * t8 - t119 * t7;
t2 = -(t91 * t51 + t92 * t52 + t67 * t81 + t68 * t82) * t122 + (t9 * t117 + (t10 - t179) * t119) * t121 + (t25 * t121 + (t181 + (t19 - t178) * t119) * t122) * qJD(2);
t1 = -(t89 * t51 + t90 * t52 + t65 * t81 + t66 * t82) * t122 + (t8 * t119 + (t7 - t179) * t117) * t121 + (t24 * t121 + (t180 + (t16 - t178) * t117) * t122) * qJD(2);
t50 = [0; m(4) * t23 + m(5) * t11 - t200; 0.2e1 * m(4) * (t22 * t23 + t54 * t63 + t55 * t64) + (t11 * t13 + t27 * t40 + t28 * t41) * t190 + (-t114 * t96 - t3 + (t103 * t71 + t104 * t73) * t119) * t119 + (t113 * t97 + t4 + (t105 * t72 + t106 * t74) * t117 + (-t103 * t72 - t104 * t74 - t105 * t71 - t106 * t73 - t117 * t96 + t119 * t97) * t119) * t117 + 0.2e1 * t199 * (rSges(3,1) * t122 - rSges(3,2) * t121) * t200; (m(4) + m(5)) * t165; m(4) * (-t122 * t23 + t168 * t55 + t170 * t54) + m(5) * (-t11 * t122 + t168 * t28 + t170 * t27) + 0.2e1 * ((t121 * t22 + t167 * t64 + t169 * t63) * t189 + (t121 * t13 + t167 * t41 + t169 * t40) * t188) * qJD(2); -0.4e1 * (t189 + t188) * t199 * t121 * t164; m(5) * t12; -t122 * (t117 * t6 - t119 * t5) / 0.2e1 + t2 * t187 - t119 * t1 / 0.2e1 + m(5) * (t11 * t26 + t12 * t13 + t14 * t41 + t15 * t40 + t27 * t31 + t28 * t30) + (t4 * t186 + t3 * t187) * t121 + (t121 * (t117 * t21 - t119 * t20) / 0.2e1 + ((t117 * t19 - t119 * t18) * t186 + (t117 * t17 - t119 * t16) * t187) * t122) * qJD(2); m(5) * (-t12 * t122 + (t117 * t15 + t119 * t14) * t121 + (t121 * t26 + (t117 * t31 + t119 * t30) * t122) * qJD(2)); (t12 * t26 + t14 * t30 + t15 * t31) * t190 + (-t122 * t25 + (t119 * t19 + t181) * t121) * t157 + t2 * t168 + (-t122 * t24 + (t117 * t16 + t180) * t121) * t159 + t1 * t170 + (t121 * t148 - t122 * t29) * t165 - t122 * ((t179 + (t122 * t150 + t148) * qJD(2)) * t122 + (t6 * t119 + t5 * t117 - (qJD(2) * t80 - t111 * t51 + t112 * t52 + (-t111 * t82 - t112 * t81) * qJD(4)) * t122 + t29 * qJD(2)) * t121);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t50(1), t50(2), t50(4), t50(7); t50(2), t50(3), t50(5), t50(8); t50(4), t50(5), t50(6), t50(9); t50(7), t50(8), t50(9), t50(10);];
Mq = res;
