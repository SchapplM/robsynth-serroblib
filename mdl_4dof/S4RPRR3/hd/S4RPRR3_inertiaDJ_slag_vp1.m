% Calculate time derivative of joint inertia matrix for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:09
% DurationCPUTime: 2.29s
% Computational Cost: add. (4535->269), mult. (4138->393), div. (0->0), fcn. (3194->8), ass. (0->160)
t108 = qJ(1) + pkin(7);
t102 = sin(t108);
t103 = cos(t108);
t109 = qJ(3) + qJ(4);
t104 = sin(t109);
t189 = rSges(5,2) * t104;
t105 = cos(t109);
t191 = rSges(5,1) * t105;
t144 = -t189 + t191;
t62 = t102 * rSges(5,3) + t144 * t103;
t110 = sin(qJ(3));
t112 = cos(qJ(3));
t92 = rSges(4,1) * t110 + rSges(4,2) * t112;
t124 = qJD(3) * t92;
t211 = t102 * t124;
t107 = qJD(3) + qJD(4);
t80 = rSges(5,1) * t104 + rSges(5,2) * t105;
t125 = t80 * t107;
t161 = qJD(3) * t110;
t157 = pkin(3) * t161;
t209 = t125 + t157;
t196 = sin(qJ(1)) * pkin(1);
t98 = t103 * pkin(5);
t208 = t98 - t196;
t129 = Icges(5,5) * t105 - Icges(5,6) * t104;
t55 = -Icges(5,3) * t103 + t129 * t102;
t207 = qJD(1) * t55;
t179 = Icges(4,4) * t112;
t133 = -Icges(4,2) * t110 + t179;
t68 = Icges(4,6) * t102 + t133 * t103;
t180 = Icges(4,4) * t110;
t136 = Icges(4,1) * t112 - t180;
t70 = Icges(4,5) * t102 + t136 * t103;
t137 = t110 * t68 - t112 * t70;
t206 = t137 * t102;
t67 = -Icges(4,6) * t103 + t133 * t102;
t69 = -Icges(4,5) * t103 + t136 * t102;
t138 = t110 * t67 - t112 * t69;
t205 = t138 * t103;
t177 = Icges(5,4) * t105;
t131 = -Icges(5,2) * t104 + t177;
t58 = Icges(5,6) * t102 + t131 * t103;
t178 = Icges(5,4) * t104;
t134 = Icges(5,1) * t105 - t178;
t60 = Icges(5,5) * t102 + t134 * t103;
t142 = t104 * t58 - t105 * t60;
t204 = t142 * t102;
t57 = -Icges(5,6) * t103 + t131 * t102;
t59 = -Icges(5,5) * t103 + t134 * t102;
t143 = t104 * t57 - t105 * t59;
t203 = t143 * t103;
t100 = t102 ^ 2;
t101 = t103 ^ 2;
t165 = t100 + t101;
t78 = Icges(5,2) * t105 + t178;
t79 = Icges(5,1) * t104 + t177;
t141 = t104 * t78 - t105 * t79;
t202 = t141 * qJD(1) + t129 * t107;
t130 = Icges(4,5) * t112 - Icges(4,6) * t110;
t65 = -Icges(4,3) * t103 + t130 * t102;
t201 = 2 * m(4);
t200 = 2 * m(5);
t199 = m(4) * t92;
t198 = t102 / 0.2e1;
t197 = -t103 / 0.2e1;
t106 = cos(qJ(1)) * pkin(1);
t114 = -pkin(6) - pkin(5);
t195 = -pkin(5) - t114;
t61 = -t103 * rSges(5,3) + t144 * t102;
t29 = t102 * t61 + t103 * t62;
t163 = qJD(1) * t103;
t164 = qJD(1) * t102;
t194 = rSges(5,3) * t163 + t164 * t189;
t192 = rSges(4,1) * t112;
t97 = t102 * rSges(4,3);
t193 = t103 * t192 + t97;
t190 = rSges(4,2) * t110;
t188 = rSges(4,3) * t103;
t187 = t107 * t78;
t186 = t107 * t79;
t185 = t110 * t69;
t184 = t110 * t70;
t183 = t112 * t67;
t182 = t112 * t68;
t181 = rSges(5,3) - t114;
t56 = Icges(5,3) * t102 + t129 * t103;
t170 = qJD(1) * t56;
t66 = Icges(4,3) * t102 + t130 * t103;
t169 = qJD(1) * t66;
t168 = t103 * t114;
t167 = t104 * t107;
t166 = t105 * t107;
t162 = qJD(1) * t110;
t160 = qJD(3) * t112;
t159 = t102 * (t62 * qJD(1) - t102 * t125) + t103 * (-t103 * rSges(5,2) * t166 + (-t103 * t167 - t105 * t164) * rSges(5,1) + t194) + t61 * t163;
t158 = t103 * t190;
t156 = t102 * t162;
t155 = -pkin(3) * t110 - t80;
t154 = t195 * t102;
t36 = -t57 * qJD(1) - t103 * t187;
t151 = t107 * t60 + t36;
t37 = t58 * qJD(1) - t102 * t187;
t150 = t107 * t59 + t37;
t38 = -t59 * qJD(1) - t103 * t186;
t149 = -t107 * t58 + t38;
t39 = t60 * qJD(1) - t102 * t186;
t148 = t107 * t57 - t39;
t13 = -t143 * t102 - t103 * t55;
t14 = -t103 * t56 - t204;
t15 = t102 * t55 - t203;
t16 = t102 * t56 - t142 * t103;
t77 = Icges(5,5) * t104 + Icges(5,6) * t105;
t121 = t77 * t107;
t34 = -t103 * t121 - t207;
t35 = -t102 * t121 + t170;
t147 = -t103 * ((t103 * t35 + (t14 + t203) * qJD(1)) * t103 + (t13 * qJD(1) + (-t104 * t36 + t105 * t38 - t58 * t166 - t60 * t167 + t170) * t102 + (-t34 + t148 * t105 + t150 * t104 + (-t142 - t55) * qJD(1)) * t103) * t102) + t102 * ((t102 * t34 + (t15 + t204) * qJD(1)) * t102 + (t16 * qJD(1) + (t104 * t37 - t105 * t39 + t57 * t166 + t59 * t167 - t207) * t103 + (-t35 + t149 * t105 - t151 * t104 + (-t143 + t56) * qJD(1)) * t102) * t103) + (t102 * t14 - t103 * t13) * t164 + (t102 * t16 - t103 * t15) * t163;
t145 = -t190 + t192;
t135 = Icges(4,1) * t110 + t179;
t132 = Icges(4,2) * t112 + t180;
t128 = -pkin(2) - t145;
t99 = pkin(3) * t112 + pkin(2);
t127 = -t144 - t99;
t74 = t131 * t107;
t75 = t134 * t107;
t115 = qJD(1) * t77 + (t75 - t187) * t105 + (-t74 - t186) * t104;
t126 = (t202 * t102 + t115 * t103 + t149 * t104 + t151 * t105) * t198 + (t115 * t102 - t202 * t103 - t148 * t104 + t150 * t105) * t197 + (-t141 * t102 - t103 * t77 + t104 * t59 + t105 * t57) * t164 / 0.2e1 + (t102 * t77 - t141 * t103 + t104 * t60 + t105 * t58) * t163 / 0.2e1;
t120 = qJD(3) * t135;
t119 = qJD(3) * t132;
t118 = qJD(3) * (-Icges(4,5) * t110 - Icges(4,6) * t112);
t116 = rSges(4,2) * t156 + rSges(4,3) * t163 - t103 * t124;
t86 = t103 * t99;
t85 = t145 * qJD(3);
t76 = t144 * t107;
t72 = -t158 + t193;
t71 = t145 * t102 - t188;
t64 = t155 * t103;
t63 = t155 * t102;
t54 = -t103 * pkin(2) + t154 + t86;
t53 = t168 + t98 + t102 * (-pkin(2) + t99);
t49 = t102 * pkin(5) + t106 + (pkin(2) - t190) * t103 + t193;
t48 = t128 * t102 + t188 + t208;
t43 = t102 * t118 + t169;
t42 = -qJD(1) * t65 + t103 * t118;
t41 = -t102 * t114 + t106 + t62 + t86;
t40 = t127 * t102 + t181 * t103 - t196;
t33 = -t80 * t163 - t102 * t76 + (-t102 * t160 - t103 * t162) * pkin(3);
t32 = t80 * t164 - t103 * t76 + (-t103 * t160 + t156) * pkin(3);
t28 = t211 + (-t106 + (-rSges(4,3) - pkin(5)) * t102 + t128 * t103) * qJD(1);
t27 = ((-pkin(2) - t192) * t102 + t208) * qJD(1) + t116;
t22 = t209 * t102 + (-t181 * t102 + t127 * t103 - t106) * qJD(1);
t21 = -t209 * t103 + (-t196 - t168 + (-t99 - t191) * t102) * qJD(1) + t194;
t20 = t102 * t66 - t137 * t103;
t19 = t102 * t65 - t205;
t18 = -t103 * t66 - t206;
t17 = -t138 * t102 - t103 * t65;
t12 = t102 * t53 + t103 * t54 + t29;
t11 = (qJD(1) * t71 + t116) * t103 + (-t211 + (-t158 - t72 + t97) * qJD(1)) * t102;
t10 = -t62 * t164 + t159;
t3 = -t165 * t157 + ((t195 * t103 + t53) * t103 + (-t54 - t62 + t154) * t102) * qJD(1) + t159;
t1 = [(t27 * t49 + t28 * t48) * t201 + t79 * t166 + t104 * t75 - t78 * t167 + t105 * t74 + (t21 * t41 + t22 * t40) * t200 + (t136 - t132) * t161 + (t135 + t133) * t160; 0; 0; (-t137 * qJD(3) + t110 * (-t69 * qJD(1) - t103 * t120) + t112 * (-t67 * qJD(1) - t103 * t119)) * t198 + (-t138 * qJD(3) + t110 * (t70 * qJD(1) - t102 * t120) + t112 * (t68 * qJD(1) - t102 * t119)) * t197 + m(4) * ((-t28 * t92 - t48 * t85) * t103 + (-t27 * t92 - t49 * t85) * t102) + m(5) * (t21 * t63 + t22 * t64 + t32 * t40 + t33 * t41) + (t100 / 0.2e1 + t101 / 0.2e1) * t130 * qJD(3) + ((t184 / 0.2e1 + t182 / 0.2e1 - t49 * t199) * t103 + (t185 / 0.2e1 + t183 / 0.2e1 + t48 * t199) * t102) * qJD(1) + t126; m(4) * t11 + m(5) * t3; ((t102 * t71 + t103 * t72) * t11 + t165 * t92 * t85) * t201 + (t102 * t20 - t103 * t19) * t163 + t102 * ((t102 * t42 + (t19 + t206) * qJD(1)) * t102 + (t20 * qJD(1) + (t67 * t160 + t69 * t161) * t103 + (-t43 + (-t182 - t184) * qJD(3) + (-t138 + t66) * qJD(1)) * t102) * t103) + (t102 * t18 - t103 * t17) * t164 - t103 * ((t103 * t43 + (t18 + t205) * qJD(1)) * t103 + (t17 * qJD(1) + (-t68 * t160 - t70 * t161 + t169) * t102 + (-t42 + (t183 + t185) * qJD(3) - t137 * qJD(1)) * t103) * t102) + (t12 * t3 + t32 * t64 + t33 * t63) * t200 + t147; m(5) * ((-t102 * t41 - t103 * t40) * t76 + (-t102 * t21 - t103 * t22 + (t102 * t40 - t103 * t41) * qJD(1)) * t80) + t126; m(5) * t10; m(5) * (t10 * t12 + t29 * t3 + (-t102 * t63 - t103 * t64) * t76 + (-t102 * t33 - t103 * t32 + (t102 * t64 - t103 * t63) * qJD(1)) * t80) + t147; (t165 * t80 * t76 + t10 * t29) * t200 + t147;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
