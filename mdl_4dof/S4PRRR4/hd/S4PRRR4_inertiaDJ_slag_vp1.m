% Calculate time derivative of joint inertia matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:34
% DurationCPUTime: 2.19s
% Computational Cost: add. (4501->262), mult. (4080->391), div. (0->0), fcn. (3158->6), ass. (0->156)
t106 = pkin(7) + qJ(2);
t102 = sin(t106);
t103 = cos(t106);
t108 = qJ(3) + qJ(4);
t104 = sin(t108);
t183 = rSges(5,2) * t104;
t105 = cos(t108);
t185 = rSges(5,1) * t105;
t140 = -t183 + t185;
t62 = t102 * rSges(5,3) + t140 * t103;
t109 = sin(qJ(3));
t110 = cos(qJ(3));
t92 = rSges(4,1) * t109 + rSges(4,2) * t110;
t120 = qJD(3) * t92;
t203 = t102 * t120;
t107 = qJD(3) + qJD(4);
t80 = rSges(5,1) * t104 + rSges(5,2) * t105;
t121 = t80 * t107;
t156 = qJD(3) * t109;
t152 = pkin(3) * t156;
t201 = t152 + t121;
t125 = Icges(5,5) * t105 - Icges(5,6) * t104;
t55 = -Icges(5,3) * t103 + t102 * t125;
t200 = qJD(2) * t55;
t173 = Icges(4,4) * t110;
t129 = -Icges(4,2) * t109 + t173;
t68 = Icges(4,6) * t102 + t103 * t129;
t174 = Icges(4,4) * t109;
t132 = Icges(4,1) * t110 - t174;
t70 = Icges(4,5) * t102 + t103 * t132;
t133 = t109 * t68 - t110 * t70;
t199 = t102 * t133;
t67 = -Icges(4,6) * t103 + t102 * t129;
t69 = -Icges(4,5) * t103 + t102 * t132;
t134 = t109 * t67 - t110 * t69;
t198 = t103 * t134;
t171 = Icges(5,4) * t105;
t127 = -Icges(5,2) * t104 + t171;
t58 = Icges(5,6) * t102 + t103 * t127;
t172 = Icges(5,4) * t104;
t130 = Icges(5,1) * t105 - t172;
t60 = Icges(5,5) * t102 + t103 * t130;
t138 = t104 * t58 - t105 * t60;
t197 = t138 * t102;
t57 = -Icges(5,6) * t103 + t102 * t127;
t59 = -Icges(5,5) * t103 + t102 * t130;
t139 = t104 * t57 - t105 * t59;
t196 = t139 * t103;
t100 = t102 ^ 2;
t101 = t103 ^ 2;
t160 = t100 + t101;
t78 = Icges(5,2) * t105 + t172;
t79 = Icges(5,1) * t104 + t171;
t137 = t104 * t78 - t105 * t79;
t195 = qJD(2) * t137 + t125 * t107;
t126 = Icges(4,5) * t110 - Icges(4,6) * t109;
t65 = -Icges(4,3) * t103 + t102 * t126;
t194 = 2 * m(4);
t193 = 2 * m(5);
t192 = m(4) * t92;
t191 = t102 / 0.2e1;
t190 = -t103 / 0.2e1;
t98 = t103 * pkin(5);
t111 = -pkin(6) - pkin(5);
t189 = -pkin(5) - t111;
t61 = -t103 * rSges(5,3) + t140 * t102;
t27 = t102 * t61 + t103 * t62;
t158 = qJD(2) * t103;
t159 = qJD(2) * t102;
t188 = rSges(5,3) * t158 + t159 * t183;
t186 = rSges(4,1) * t110;
t97 = t102 * rSges(4,3);
t187 = t103 * t186 + t97;
t184 = rSges(4,2) * t109;
t182 = rSges(4,3) * t103;
t181 = t107 * t78;
t180 = t107 * t79;
t179 = t109 * t69;
t178 = t109 * t70;
t177 = t110 * t67;
t176 = t110 * t68;
t175 = rSges(5,3) - t111;
t56 = Icges(5,3) * t102 + t103 * t125;
t164 = qJD(2) * t56;
t66 = Icges(4,3) * t102 + t103 * t126;
t163 = qJD(2) * t66;
t162 = t104 * t107;
t161 = t105 * t107;
t157 = qJD(2) * t109;
t155 = qJD(3) * t110;
t154 = t102 * (t62 * qJD(2) - t102 * t121) + t103 * (-t103 * rSges(5,2) * t161 + (-t103 * t162 - t105 * t159) * rSges(5,1) + t188) + t61 * t158;
t153 = t103 * t184;
t151 = t102 * t157;
t150 = -pkin(3) * t109 - t80;
t149 = t189 * t102;
t36 = -t57 * qJD(2) - t103 * t181;
t146 = t107 * t60 + t36;
t37 = qJD(2) * t58 - t102 * t181;
t145 = t107 * t59 + t37;
t38 = -t59 * qJD(2) - t103 * t180;
t144 = -t107 * t58 + t38;
t39 = qJD(2) * t60 - t102 * t180;
t143 = t107 * t57 - t39;
t13 = -t139 * t102 - t103 * t55;
t14 = -t103 * t56 - t197;
t15 = t102 * t55 - t196;
t16 = t102 * t56 - t138 * t103;
t77 = Icges(5,5) * t104 + Icges(5,6) * t105;
t117 = t77 * t107;
t34 = -t103 * t117 - t200;
t35 = -t102 * t117 + t164;
t142 = -t103 * ((t103 * t35 + (t14 + t196) * qJD(2)) * t103 + (t13 * qJD(2) + (-t104 * t36 + t105 * t38 - t58 * t161 - t60 * t162 + t164) * t102 + (-t34 + t143 * t105 + t145 * t104 + (-t138 - t55) * qJD(2)) * t103) * t102) + t102 * ((t102 * t34 + (t15 + t197) * qJD(2)) * t102 + (t16 * qJD(2) + (t104 * t37 - t105 * t39 + t57 * t161 + t59 * t162 - t200) * t103 + (-t35 + t144 * t105 - t146 * t104 + (-t139 + t56) * qJD(2)) * t102) * t103) + (t102 * t14 - t103 * t13) * t159 + (t102 * t16 - t103 * t15) * t158;
t141 = -t184 + t186;
t131 = Icges(4,1) * t109 + t173;
t128 = Icges(4,2) * t110 + t174;
t124 = -pkin(2) - t141;
t99 = pkin(3) * t110 + pkin(2);
t123 = -t140 - t99;
t74 = t127 * t107;
t75 = t130 * t107;
t112 = qJD(2) * t77 + (t75 - t181) * t105 + (-t74 - t180) * t104;
t122 = (t195 * t102 + t112 * t103 + t144 * t104 + t146 * t105) * t191 + (t112 * t102 - t195 * t103 - t143 * t104 + t145 * t105) * t190 + (-t102 * t137 - t103 * t77 + t104 * t59 + t105 * t57) * t159 / 0.2e1 + (t102 * t77 - t103 * t137 + t104 * t60 + t105 * t58) * t158 / 0.2e1;
t116 = qJD(3) * t131;
t115 = qJD(3) * t128;
t114 = qJD(3) * (-Icges(4,5) * t109 - Icges(4,6) * t110);
t113 = rSges(4,2) * t151 + rSges(4,3) * t158 - t103 * t120;
t86 = t103 * t99;
t85 = t141 * qJD(3);
t76 = t140 * t107;
t72 = -t153 + t187;
t71 = t141 * t102 - t182;
t64 = t150 * t103;
t63 = t150 * t102;
t54 = -t103 * pkin(2) + t149 + t86;
t53 = t103 * t111 + t98 + t102 * (-pkin(2) + t99);
t49 = t102 * pkin(5) + (pkin(2) - t184) * t103 + t187;
t48 = t102 * t124 + t182 + t98;
t47 = -t102 * t111 + t62 + t86;
t46 = t123 * t102 + t175 * t103;
t41 = t102 * t114 + t163;
t40 = -qJD(2) * t65 + t103 * t114;
t33 = t203 + ((-rSges(4,3) - pkin(5)) * t102 + t124 * t103) * qJD(2);
t32 = (t98 + (-pkin(2) - t186) * t102) * qJD(2) + t113;
t31 = -t80 * t158 - t102 * t76 + (-t102 * t155 - t103 * t157) * pkin(3);
t30 = t80 * t159 - t103 * t76 + (-t103 * t155 + t151) * pkin(3);
t22 = t201 * t102 + (-t175 * t102 + t123 * t103) * qJD(2);
t21 = (-t99 - t185) * t159 + (-qJD(2) * t111 - t201) * t103 + t188;
t20 = t66 * t102 - t103 * t133;
t19 = t102 * t65 - t198;
t18 = -t103 * t66 - t199;
t17 = -t102 * t134 - t103 * t65;
t12 = t102 * t53 + t103 * t54 + t27;
t11 = (qJD(2) * t71 + t113) * t103 + (-t203 + (-t153 - t72 + t97) * qJD(2)) * t102;
t10 = -t62 * t159 + t154;
t3 = -t160 * t152 + ((t189 * t103 + t53) * t103 + (-t54 - t62 + t149) * t102) * qJD(2) + t154;
t1 = [0; 0; t79 * t161 + t104 * t75 - t78 * t162 + t105 * t74 + (t32 * t49 + t33 * t48) * t194 + (t21 * t47 + t22 * t46) * t193 + (t132 - t128) * t156 + (t131 + t129) * t155; m(4) * t11 + m(5) * t3; (-qJD(3) * t133 + t109 * (-t69 * qJD(2) - t103 * t116) + t110 * (-t67 * qJD(2) - t103 * t115)) * t191 + (-qJD(3) * t134 + t109 * (qJD(2) * t70 - t102 * t116) + t110 * (qJD(2) * t68 - t102 * t115)) * t190 + m(4) * ((-t33 * t92 - t48 * t85) * t103 + (-t32 * t92 - t49 * t85) * t102) + m(5) * (t21 * t63 + t22 * t64 + t30 * t46 + t31 * t47) + (t100 / 0.2e1 + t101 / 0.2e1) * t126 * qJD(3) + ((t178 / 0.2e1 + t176 / 0.2e1 - t49 * t192) * t103 + (t179 / 0.2e1 + t177 / 0.2e1 + t48 * t192) * t102) * qJD(2) + t122; ((t102 * t71 + t103 * t72) * t11 + t160 * t92 * t85) * t194 + (t102 * t20 - t103 * t19) * t158 + t102 * ((t102 * t40 + (t19 + t199) * qJD(2)) * t102 + (t20 * qJD(2) + (t67 * t155 + t69 * t156) * t103 + (-t41 + (-t176 - t178) * qJD(3) + (-t134 + t66) * qJD(2)) * t102) * t103) + (t102 * t18 - t103 * t17) * t159 - t103 * ((t103 * t41 + (t18 + t198) * qJD(2)) * t103 + (t17 * qJD(2) + (-t68 * t155 - t70 * t156 + t163) * t102 + (-t40 + (t177 + t179) * qJD(3) - t133 * qJD(2)) * t103) * t102) + (t12 * t3 + t30 * t64 + t31 * t63) * t193 + t142; m(5) * t10; m(5) * ((-t102 * t47 - t103 * t46) * t76 + (-t102 * t21 - t103 * t22 + (t102 * t46 - t103 * t47) * qJD(2)) * t80) + t122; m(5) * (t10 * t12 + t27 * t3 + (-t102 * t63 - t103 * t64) * t76 + (-t102 * t31 - t103 * t30 + (t102 * t64 - t103 * t63) * qJD(2)) * t80) + t142; (t160 * t80 * t76 + t10 * t27) * t193 + t142;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
