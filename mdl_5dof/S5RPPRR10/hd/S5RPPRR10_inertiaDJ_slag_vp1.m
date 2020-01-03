% Calculate time derivative of joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:06
% DurationCPUTime: 3.17s
% Computational Cost: add. (5754->405), mult. (9372->569), div. (0->0), fcn. (9077->8), ass. (0->187)
t151 = sin(pkin(8));
t155 = cos(qJ(4));
t152 = cos(pkin(8));
t153 = sin(qJ(4));
t203 = t152 * t153;
t128 = t151 * t155 - t203;
t156 = cos(qJ(1));
t120 = t128 * t156;
t137 = pkin(4) * t155 + pkin(3);
t224 = -pkin(3) + t137;
t150 = qJ(4) + qJ(5);
t142 = sin(t150);
t143 = cos(t150);
t172 = t142 * t151 + t143 * t152;
t173 = t142 * t152 - t143 * t151;
t96 = -rSges(6,1) * t173 - rSges(6,2) * t172;
t208 = -pkin(4) * t203 + t151 * t224 + t96;
t51 = t208 * t156;
t238 = qJD(1) * t51;
t154 = sin(qJ(1));
t237 = t172 * t154;
t198 = t154 ^ 2 + t156 ^ 2;
t145 = t156 * qJ(2);
t206 = qJ(3) * t151;
t166 = -t206 - pkin(1) + (-pkin(2) - pkin(3)) * t152;
t226 = -rSges(5,3) - pkin(6);
t159 = t154 * t166 + t156 * t226;
t118 = t128 * t154;
t205 = t151 * t153;
t171 = t152 * t155 + t205;
t119 = t171 * t154;
t179 = -rSges(5,1) * t119 - rSges(5,2) * t118;
t48 = t145 + t159 + t179;
t199 = t156 * pkin(1) + t154 * qJ(2);
t202 = t152 * t156;
t185 = pkin(2) * t202 + t156 * t206 + t199;
t193 = pkin(3) * t202;
t121 = t171 * t156;
t201 = rSges(5,1) * t121 + rSges(5,2) * t120;
t49 = t154 * t226 + t185 + t193 + t201;
t236 = t154 * t49 + t156 * t48;
t195 = qJD(3) * t151;
t196 = qJD(1) * t156;
t200 = qJ(2) * t196 + qJD(2) * t154;
t190 = t156 * t195 + t200;
t125 = t171 * qJD(4);
t197 = qJD(1) * t154;
t85 = -t125 * t156 - t128 * t197;
t234 = qJD(4) * t128;
t86 = t156 * t234 - t171 * t197;
t221 = rSges(5,1) * t86 + rSges(5,2) * t85;
t27 = qJD(1) * t159 + t190 + t221;
t139 = pkin(6) * t197;
t141 = qJD(2) * t156;
t181 = -t154 * t195 + t141;
t87 = qJD(1) * t120 - t125 * t154;
t88 = qJD(1) * t121 + t154 * t234;
t183 = rSges(5,1) * t88 + rSges(5,2) * t87;
t28 = t139 + ((rSges(5,3) - qJ(2)) * t154 + t166 * t156) * qJD(1) + t181 - t183;
t235 = t154 * t27 + t156 * t28;
t233 = 2 * m(5);
t232 = 2 * m(6);
t231 = m(5) / 0.2e1;
t230 = m(6) / 0.2e1;
t229 = -t154 / 0.2e1;
t227 = t156 / 0.2e1;
t113 = t173 * t156;
t114 = t172 * t156;
t111 = t173 * t154;
t61 = Icges(6,5) * t237 - Icges(6,6) * t111 + Icges(6,3) * t156;
t63 = Icges(6,4) * t237 - Icges(6,2) * t111 + Icges(6,6) * t156;
t65 = Icges(6,1) * t237 - Icges(6,4) * t111 + Icges(6,5) * t156;
t16 = -t113 * t63 + t114 * t65 - t154 * t61;
t62 = Icges(6,5) * t114 - Icges(6,6) * t113 - Icges(6,3) * t154;
t218 = t154 * t62;
t64 = Icges(6,4) * t114 - Icges(6,2) * t113 - Icges(6,6) * t154;
t66 = Icges(6,1) * t114 - Icges(6,4) * t113 - Icges(6,5) * t154;
t17 = -t113 * t64 + t114 * t66 - t218;
t211 = t156 * t61;
t175 = -t211 + t218;
t147 = qJD(4) + qJD(5);
t107 = t172 * t147;
t57 = -t107 * t156 + t173 * t197;
t58 = -qJD(1) * t237 - t113 * t147;
t31 = Icges(6,5) * t58 + Icges(6,6) * t57 - Icges(6,3) * t196;
t59 = -qJD(1) * t113 - t147 * t237;
t60 = qJD(1) * t114 - t111 * t147;
t32 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t197;
t33 = Icges(6,4) * t58 + Icges(6,2) * t57 - Icges(6,6) * t196;
t34 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t197;
t35 = Icges(6,1) * t58 + Icges(6,4) * t57 - Icges(6,5) * t196;
t36 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t197;
t2 = -(-t113 * t33 + t114 * t35 - t154 * t31 + t57 * t64 + t58 * t66) * t154 + (-t113 * t34 + t114 * t36 - t154 * t32 + t57 * t63 + t58 * t65) * t156 + (-t16 * t154 + (-t17 + t175) * t156) * qJD(1);
t225 = t154 * t2;
t223 = rSges(6,1) * t58 + rSges(6,2) * t57;
t157 = -pkin(7) - pkin(6);
t194 = pkin(4) * t205;
t178 = -rSges(6,1) * t237 + t111 * rSges(6,2);
t216 = t156 * rSges(6,3);
t67 = -t178 + t216;
t222 = -t67 - (-pkin(6) - t157) * t156 - t154 * (t152 * t224 + t194);
t108 = t173 * t147;
t72 = -rSges(6,1) * t107 + rSges(6,2) * t108;
t207 = -pkin(4) * t125 + t72;
t29 = t156 * t207 - t197 * t208;
t214 = t156 * t29;
t212 = t156 * t51;
t209 = -rSges(6,3) + t157;
t138 = t154 * t157;
t191 = t137 * t202 + t156 * t194 + t138;
t187 = t198 * t72;
t100 = -rSges(5,1) * t125 - rSges(5,2) * t234;
t186 = t100 * t198;
t14 = -t111 * t63 + t237 * t65 + t211;
t15 = -t111 * t64 + t156 * t62 + t237 * t66;
t1 = t156 * (-(-t111 * t33 + t156 * t31 + t237 * t35 + t59 * t64 + t60 * t66) * t154 + (-t111 * t34 + t156 * t32 + t237 * t36 + t59 * t63 + t60 * t65) * t156 + (-t15 * t156 + (-t14 + t175) * t154) * qJD(1));
t9 = -t154 * t17 + t156 * t16;
t184 = -t196 * t9 + t1;
t68 = rSges(6,1) * t114 - rSges(6,2) * t113 - rSges(6,3) * t154;
t182 = -t60 * rSges(6,1) - t59 * rSges(6,2);
t180 = rSges(3,1) * t152 - rSges(3,2) * t151;
t162 = -pkin(1) + (-pkin(2) - t137) * t152 + (-pkin(4) * t153 - qJ(3)) * t151;
t160 = t162 * t154;
t165 = pkin(4) * qJD(4) * t120 + t157 * t196;
t12 = (t160 - t216) * qJD(1) + t165 + t190 + t223;
t164 = pkin(4) * t234;
t13 = t141 + (-t164 - t195) * t154 + ((-qJ(2) - t209) * t154 + t162 * t156) * qJD(1) + t182;
t177 = t12 * t154 + t13 * t156;
t40 = t156 * t209 + t145 + t160 + t178;
t41 = t68 + t185 + t191;
t176 = t154 * t41 + t156 * t40;
t170 = -pkin(1) - t180;
t69 = -Icges(6,5) * t107 + Icges(6,6) * t108;
t70 = -Icges(6,4) * t107 + Icges(6,2) * t108;
t71 = -Icges(6,1) * t107 + Icges(6,4) * t108;
t93 = -Icges(6,5) * t173 - Icges(6,6) * t172;
t94 = -Icges(6,4) * t173 - Icges(6,2) * t172;
t95 = -Icges(6,1) * t173 - Icges(6,4) * t172;
t169 = (-t107 * t66 + t108 * t64 - t113 * t70 + t114 * t71 - t154 * t69 - t172 * t33 - t173 * t35 - t93 * t196 + t57 * t94 + t58 * t95) * t229 + (-t107 * t65 + t108 * t63 - t111 * t70 + t156 * t69 - t172 * t34 - t173 * t36 - t93 * t197 + t237 * t71 + t59 * t94 + t60 * t95) * t227 - (-t111 * t94 + t156 * t93 - t172 * t63 - t173 * t65 + t237 * t95) * t197 / 0.2e1 - (-t113 * t94 + t114 * t95 - t154 * t93 - t172 * t64 - t173 * t66) * t196 / 0.2e1;
t163 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t152 + (-rSges(4,3) - qJ(3)) * t151;
t161 = rSges(3,3) * t156 + t154 * t170;
t158 = rSges(4,2) * t156 + t154 * t163;
t106 = rSges(3,3) * t154 + t156 * t180 + t199;
t105 = t145 + t161;
t104 = rSges(5,1) * t128 - rSges(5,2) * t171;
t103 = Icges(5,1) * t128 - Icges(5,4) * t171;
t102 = Icges(5,4) * t128 - Icges(5,2) * t171;
t99 = -Icges(5,1) * t125 - Icges(5,4) * t234;
t98 = -Icges(5,4) * t125 - Icges(5,2) * t234;
t97 = -Icges(5,5) * t125 - Icges(5,6) * t234;
t92 = pkin(6) * t154 + t191 - t193;
t90 = t141 + ((-rSges(3,3) - qJ(2)) * t154 + t170 * t156) * qJD(1);
t89 = qJD(1) * t161 + t200;
t84 = rSges(4,2) * t154 + (rSges(4,1) * t152 + rSges(4,3) * t151) * t156 + t185;
t83 = t145 + t158;
t80 = -rSges(5,3) * t154 + t201;
t79 = rSges(5,3) * t156 - t179;
t78 = Icges(5,1) * t121 + Icges(5,4) * t120 - Icges(5,5) * t154;
t77 = Icges(5,1) * t119 + Icges(5,4) * t118 + Icges(5,5) * t156;
t76 = Icges(5,4) * t121 + Icges(5,2) * t120 - Icges(5,6) * t154;
t75 = Icges(5,4) * t119 + Icges(5,2) * t118 + Icges(5,6) * t156;
t74 = Icges(5,5) * t121 + Icges(5,6) * t120 - Icges(5,3) * t154;
t73 = Icges(5,5) * t119 + Icges(5,6) * t118 + Icges(5,3) * t156;
t54 = t68 * t197;
t53 = ((-rSges(4,2) - qJ(2)) * t154 + t163 * t156) * qJD(1) + t181;
t52 = qJD(1) * t158 + t190;
t50 = t208 * t154;
t47 = Icges(5,1) * t88 + Icges(5,4) * t87 - Icges(5,5) * t197;
t46 = Icges(5,1) * t86 + Icges(5,4) * t85 - Icges(5,5) * t196;
t45 = Icges(5,4) * t88 + Icges(5,2) * t87 - Icges(5,6) * t197;
t44 = Icges(5,4) * t86 + Icges(5,2) * t85 - Icges(5,6) * t196;
t43 = Icges(5,5) * t88 + Icges(5,6) * t87 - Icges(5,3) * t197;
t42 = Icges(5,5) * t86 + Icges(5,6) * t85 - Icges(5,3) * t196;
t39 = -t154 * t67 - t156 * t68;
t38 = -rSges(6,3) * t197 - t182;
t37 = -rSges(6,3) * t196 + t223;
t30 = t154 * t207 + t238;
t22 = t120 * t76 + t121 * t78 - t154 * t74;
t21 = t120 * t75 + t121 * t77 - t154 * t73;
t20 = t118 * t76 + t119 * t78 + t156 * t74;
t19 = t118 * t75 + t119 * t77 + t156 * t73;
t18 = (-t68 - t92) * t156 + t222 * t154;
t11 = -t154 * t183 - t156 * t221 + (rSges(5,3) * t198 + t154 * t80 - t156 * t79) * qJD(1);
t10 = -t154 * t38 + t54 + (-qJD(1) * t67 - t37) * t156;
t8 = t14 * t156 - t15 * t154;
t3 = t54 + (-t165 - t37) * t156 + (-t154 * t164 - t139 - t38) * t154 + ((t92 - t138) * t154 + (-t156 * pkin(6) + t222) * t156) * qJD(1);
t4 = [-t107 * t95 - t173 * t71 + t108 * t94 - t172 * t70 + (t12 * t41 + t13 * t40) * t232 + (t27 * t49 + t28 * t48) * t233 - t125 * t103 + t128 * t99 - t234 * t102 - t171 * t98 + 0.2e1 * m(4) * (t52 * t84 + t53 * t83) + 0.2e1 * m(3) * (t105 * t90 + t106 * t89); m(6) * (qJD(1) * t176 - t12 * t156 + t13 * t154) + m(5) * (qJD(1) * t236 + t154 * t28 - t156 * t27) + m(4) * (t154 * t53 - t156 * t52 + (t154 * t84 + t156 * t83) * qJD(1)) + m(3) * (t154 * t90 - t156 * t89 + (t105 * t156 + t106 * t154) * qJD(1)); 0; 0.2e1 * ((t196 * t41 - t197 * t40 + t177) * t230 + (t196 * t49 - t197 * t48 + t235) * t231 + m(4) * (t154 * t52 + t156 * t53 + t196 * t84 - t197 * t83) / 0.2e1) * t151; 0; 0; m(6) * (t12 * t50 + t13 * t51 + t29 * t40 + t30 * t41) + m(5) * (t236 * t100 + t235 * t104) + (m(5) * (-t154 * t48 + t156 * t49) * t104 - (t102 * t120 + t103 * t121 + t128 * t78 - t171 * t76) * t156 / 0.2e1) * qJD(1) + t169 + (t85 * t102 + t86 * t103 + t120 * t98 + t121 * t99 - t154 * t97 - t125 * t78 - t234 * t76 - t171 * t44 + t128 * t46 + (t102 * t118 + t103 * t119 + t128 * t77 - t171 * t75) * qJD(1)) * t229 + (t102 * t87 + t103 * t88 + t118 * t98 + t119 * t99 - t125 * t77 + t128 * t47 + t156 * t97 - t171 * t45 - t234 * t75) * t227; m(6) * (t154 * t29 - t156 * t30 + (t154 * t50 + t212) * qJD(1)); 0.2e1 * (-m(5) * t11 / 0.2e1 - m(6) * t3 / 0.2e1) * t152 + 0.2e1 * ((t154 * t30 + t214 + (-t154 * t51 + t156 * t50) * qJD(1)) * t230 + t186 * t231) * t151; ((-t154 * t79 - t156 * t80) * t11 + t104 * t186) * t233 - (-t154 * t22 + t156 * t21) * t196 - t154 * (-(t120 * t44 + t121 * t46 - t154 * t42 + t76 * t85 + t78 * t86) * t154 + (t120 * t45 + t121 * t47 - t154 * t43 + t85 * t75 + t86 * t77) * t156 + (-t21 * t154 - t22 * t156) * qJD(1)) + t156 * (-(t118 * t44 + t119 * t46 + t156 * t42 + t87 * t76 + t88 * t78) * t154 + (t118 * t45 + t119 * t47 + t156 * t43 + t75 * t87 + t77 * t88) * t156 + (-t154 * t19 - t20 * t156) * qJD(1)) + (t18 * t3 + t29 * t51 + t30 * t50) * t232 - t225 + t184 + (t154 * t20 - t19 * t156 - t8) * t197; m(6) * (t176 * t72 + ((-t154 * t40 + t156 * t41) * qJD(1) + t177) * t96) + t169; 0; m(6) * (-t10 * t152 + t151 * t187); m(6) * (t10 * t18 + t72 * t212 + t3 * t39 + (t196 * t50 + t214) * t96) + (m(6) * (t50 * t72 + (t30 - t238) * t96) - t2 - qJD(1) * t8) * t154 + t184; (t10 * t39 + t187 * t96) * t232 - t225 + t1 + (-t154 * t8 - t156 * t9) * qJD(1);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
