% Calculate time derivative of joint inertia matrix for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:43
% DurationCPUTime: 3.15s
% Computational Cost: add. (4580->328), mult. (4750->469), div. (0->0), fcn. (3626->8), ass. (0->189)
t134 = sin(qJ(1));
t135 = cos(qJ(1));
t127 = pkin(8) + qJ(4);
t117 = qJ(5) + t127;
t111 = sin(t117);
t112 = cos(t117);
t167 = rSges(6,1) * t111 + rSges(6,2) * t112;
t66 = t135 * rSges(6,3) + t167 * t134;
t116 = cos(t127);
t115 = sin(t127);
t237 = rSges(5,1) * t115;
t168 = rSges(5,2) * t116 + t237;
t133 = -pkin(6) - qJ(3);
t126 = -pkin(7) + t133;
t131 = sin(pkin(8));
t239 = pkin(3) * t131;
t99 = pkin(4) * t115 + t239;
t91 = t134 * t99;
t259 = -t126 * t135 + t91;
t144 = t135 * t168;
t221 = Icges(5,4) * t115;
t153 = Icges(5,2) * t116 + t221;
t70 = Icges(5,6) * t135 + t153 * t134;
t220 = Icges(5,4) * t116;
t156 = Icges(5,1) * t115 + t220;
t72 = Icges(5,5) * t135 + t156 * t134;
t161 = t115 * t72 + t116 * t70;
t258 = t161 * t135;
t219 = Icges(6,4) * t111;
t152 = Icges(6,2) * t112 + t219;
t62 = Icges(6,6) * t135 + t152 * t134;
t218 = Icges(6,4) * t112;
t155 = Icges(6,1) * t111 + t218;
t64 = Icges(6,5) * t135 + t155 * t134;
t166 = t111 * t64 + t112 * t62;
t257 = t166 * t135;
t122 = t135 * qJ(2);
t110 = t135 * t239;
t113 = t134 * t133;
t203 = t110 + t113;
t241 = -rSges(5,3) - pkin(1);
t49 = t241 * t134 + t122 + t144 + t203;
t190 = t134 * t239;
t201 = t135 * pkin(1) + t134 * qJ(2);
t124 = t135 * rSges(5,3);
t77 = t168 * t134 + t124;
t50 = -t133 * t135 + t190 + t201 + t77;
t256 = -t134 * t49 + t135 * t50;
t195 = qJD(4) * t134;
t183 = t115 * t195;
t198 = qJD(1) * t135;
t202 = qJ(2) * t198 + qJD(2) * t134;
t185 = qJD(3) * t135 + t202;
t182 = t116 * t195;
t184 = t116 * t198;
t186 = rSges(5,1) * t182 + rSges(5,2) * t184 + t198 * t237;
t22 = -rSges(5,2) * t183 + (t110 + (t133 + t241) * t134) * qJD(1) + t185 + t186;
t109 = t133 * t198;
t120 = qJD(2) * t135;
t170 = -qJD(3) * t134 + t120;
t194 = qJD(4) * t135;
t234 = rSges(5,2) * t115;
t95 = rSges(5,1) * t116 - t234;
t23 = t109 + t95 * t194 + (t241 * t135 + (-qJ(2) - t168 - t239) * t134) * qJD(1) + t170;
t255 = t134 * t23 - t135 * t22;
t149 = Icges(6,5) * t111 + Icges(6,6) * t112;
t254 = -Icges(6,3) * t134 + t149 * t135;
t150 = Icges(5,5) * t115 + Icges(5,6) * t116;
t253 = -Icges(5,3) * t134 + t150 * t135;
t252 = -Icges(6,6) * t134 + t152 * t135;
t251 = -Icges(5,6) * t134 + t153 * t135;
t250 = -Icges(6,5) * t134 + t155 * t135;
t249 = -Icges(5,5) * t134 + t156 * t135;
t128 = qJD(4) + qJD(5);
t75 = t152 * t128;
t76 = t155 * t128;
t86 = Icges(6,5) * t112 - Icges(6,6) * t111;
t87 = -Icges(6,2) * t111 + t218;
t88 = Icges(6,1) * t112 - t219;
t248 = (t128 * t87 + t76) * t111 - (t128 * t88 - t75) * t112 + qJD(1) * t86;
t247 = 2 * m(5);
t246 = 2 * m(6);
t129 = t134 ^ 2;
t130 = t135 ^ 2;
t245 = m(5) * t95;
t244 = t134 / 0.2e1;
t243 = t135 / 0.2e1;
t242 = rSges(3,2) - pkin(1);
t240 = -rSges(6,3) - pkin(1);
t238 = t190 - t91 - (-t126 + t133) * t135 - t66;
t235 = rSges(6,1) * t112;
t231 = t115 * t70;
t230 = t115 * t251;
t229 = t116 * t72;
t228 = t116 * t249;
t227 = t134 * rSges(5,3);
t222 = rSges(4,3) + qJ(3);
t60 = Icges(6,3) * t135 + t149 * t134;
t211 = qJD(1) * t60;
t68 = Icges(5,3) * t135 + t150 * t134;
t210 = qJD(1) * t68;
t208 = t111 * t128;
t207 = t112 * t128;
t205 = t128 * t134;
t204 = t128 * t135;
t200 = t129 + t130;
t199 = qJD(1) * t134;
t197 = qJD(4) * t115;
t196 = qJD(4) * t116;
t193 = -pkin(1) - t222;
t192 = t126 + t240;
t191 = t167 * t198 + t205 * t235;
t189 = rSges(6,2) * t208;
t188 = pkin(4) * t196;
t187 = -pkin(4) * t182 - t126 * t199 - t99 * t198;
t89 = -rSges(6,2) * t111 + t235;
t181 = pkin(4) * t116 + t89;
t34 = t62 * qJD(1) - t87 * t204;
t178 = t128 * t250 - t34;
t35 = t252 * qJD(1) + t87 * t205;
t177 = t128 * t64 + t35;
t36 = t64 * qJD(1) - t88 * t204;
t176 = -t128 * t252 - t36;
t37 = t250 * qJD(1) + t88 * t205;
t175 = -t128 * t62 + t37;
t79 = t167 * t128;
t172 = t200 * t79;
t85 = t168 * qJD(4);
t171 = t200 * t85;
t169 = rSges(4,1) * t131 + rSges(4,2) * cos(pkin(8));
t165 = -t111 * t250 - t112 * t252;
t164 = t111 * t88 + t112 * t87;
t160 = -t115 * t249 - t116 * t251;
t147 = t167 + t99;
t39 = t192 * t134 + t147 * t135 + t122;
t40 = t201 + t66 + t259;
t159 = t134 * t39 - t135 * t40;
t55 = t181 * t134;
t56 = t181 * t135;
t158 = t134 * t55 + t135 * t56;
t157 = Icges(5,1) * t116 - t221;
t154 = -Icges(5,2) * t115 + t220;
t151 = Icges(5,5) * t116 - Icges(5,6) * t115;
t14 = t166 * t134 + t135 * t60;
t143 = t165 * t134;
t15 = -t135 * t254 + t143;
t16 = t134 * t60 - t257;
t17 = -t134 * t254 - t165 * t135;
t32 = -t86 * t204 + t211;
t33 = t254 * qJD(1) + t86 * t205;
t148 = -(t15 * t134 + t135 * t14) * t199 + t134 * ((t134 * t32 + (-t16 + t143) * qJD(1)) * t134 + (t17 * qJD(1) + (-t111 * t37 - t112 * t35 - t64 * t207 + t62 * t208 + t211) * t135 + (t166 * qJD(1) + t176 * t111 + t178 * t112 + t33) * t134) * t135) + t135 * ((t135 * t33 + (t15 + t257) * qJD(1)) * t135 + (-t14 * qJD(1) + (t111 * t36 + t112 * t34 - t207 * t250 + t208 * t252) * t134 + (t32 + t177 * t112 + t175 * t111 + (t165 - t60) * qJD(1)) * t135) * t134) + (t17 * t134 + t135 * t16) * t198;
t146 = rSges(3,3) * t135 + t242 * t134;
t139 = t164 * qJD(1) - t149 * t128;
t145 = (t178 * t111 - t176 * t112 + t139 * t134 + t248 * t135) * t244 + (-t177 * t111 + t175 * t112 - t248 * t134 + t139 * t135) * t243 - (-t111 * t62 + t112 * t64 + t164 * t134 + t135 * t86) * t199 / 0.2e1 + (t111 * t252 - t112 * t250 + t134 * t86 - t164 * t135) * t198 / 0.2e1;
t142 = t160 * t134;
t141 = qJD(4) * t157;
t140 = qJD(4) * t154;
t138 = t193 * t134 + t169 * t135;
t12 = (t240 * qJD(1) - t189) * t134 + t185 - t187 + t191;
t13 = (t89 * t128 + t188) * t135 + (t192 * t135 + (-qJ(2) - t147) * t134) * qJD(1) + t170;
t137 = -t12 * t135 + t134 * t13 + (t134 * t40 + t135 * t39) * qJD(1);
t28 = t89 * t199 + t135 * t79 + (t115 * t194 + t116 * t199) * pkin(4);
t29 = t89 * t198 - t134 * t79 + (-t183 + t184) * pkin(4);
t136 = t29 * t134 - t28 * t135 + (-t134 * t56 + t135 * t55) * qJD(1);
t81 = -rSges(3,2) * t135 + t134 * rSges(3,3) + t201;
t80 = t122 + t146;
t78 = t227 - t144;
t67 = t134 * rSges(6,3) - t167 * t135;
t59 = t120 + (t242 * t135 + (-rSges(3,3) - qJ(2)) * t134) * qJD(1);
t58 = t146 * qJD(1) + t202;
t57 = t135 * t67;
t54 = t169 * t134 + t222 * t135 + t201;
t53 = t122 + t138;
t51 = -t134 * t126 - t135 * t99 + t203;
t44 = t253 * qJD(1) + t151 * t195;
t43 = -t151 * t194 + t210;
t42 = (t193 * t135 + (-qJ(2) - t169) * t134) * qJD(1) + t170;
t41 = t138 * qJD(1) + t185;
t38 = (-rSges(6,3) * qJD(1) - t189) * t134 + t191;
t31 = -t134 * t66 + t57;
t30 = t135 * (qJD(1) * t66 - t89 * t204);
t21 = -t134 * t253 - t160 * t135;
t20 = t134 * t68 - t258;
t19 = -t135 * t253 + t142;
t18 = t161 * t134 + t135 * t68;
t11 = t238 * t134 + t135 * t51 + t57;
t10 = -t134 * t38 + t30 + (-t134 * t67 - t135 * t66) * qJD(1);
t3 = t30 + (-t38 + (-t51 - t67 + t113) * qJD(1) + t187) * t134 + (-t135 * t188 + t109 + (t238 + t259) * qJD(1)) * t135;
t1 = [(t12 * t40 + t13 * t39) * t246 + 0.2e1 * m(3) * (t58 * t81 + t59 * t80) + 0.2e1 * m(4) * (t41 * t54 + t42 * t53) - t115 * t141 - t156 * t196 - t116 * t140 + t153 * t197 + (t22 * t50 + t23 * t49) * t247 - t88 * t208 - t112 * t76 - t87 * t207 + t111 * t75; m(6) * t137 + m(3) * (t134 * t59 - t135 * t58 + (t134 * t81 + t135 * t80) * qJD(1)) + m(4) * (t134 * t42 - t135 * t41 + (t134 * t54 + t135 * t53) * qJD(1)) + m(5) * ((t134 * t50 + t135 * t49) * qJD(1) + t255); 0; m(6) * (-t159 * qJD(1) + t134 * t12 + t13 * t135) + m(4) * (t134 * t41 + t135 * t42 + (-t134 * t53 + t135 * t54) * qJD(1)) + m(5) * (t256 * qJD(1) + t134 * t22 + t135 * t23); 0; 0; (-t161 * qJD(4) - t115 * (t251 * qJD(1) + t134 * t140) + t116 * (t249 * qJD(1) + t134 * t141)) * t243 + (-t160 * qJD(4) - t115 * (t70 * qJD(1) - t154 * t194) + t116 * (t72 * qJD(1) - t157 * t194)) * t244 + m(6) * (-t12 * t56 + t13 * t55 + t28 * t40 + t29 * t39) + m(5) * (t255 * t95 + t256 * t85) - (t130 / 0.2e1 + t129 / 0.2e1) * t150 * qJD(4) + ((t231 / 0.2e1 - t229 / 0.2e1 + t50 * t245) * t134 + (t230 / 0.2e1 - t228 / 0.2e1 + t49 * t245) * t135) * qJD(1) + t145; -m(5) * t171 + m(6) * t136; m(6) * (-t158 * qJD(1) + t28 * t134 + t135 * t29); ((-t134 * t77 + t135 * t78) * (-t134 * t186 + (t129 * t234 - t95 * t130) * qJD(4) + ((-t77 + t124) * t135 + (t144 - t78 + t227) * t134) * qJD(1)) - t95 * t171) * t247 - (t19 * t134 + t135 * t18) * t199 + t135 * ((t135 * t44 + (t19 + t258) * qJD(1)) * t135 + (-t18 * qJD(1) + (-t196 * t249 + t197 * t251) * t134 + (t43 + (t229 - t231) * qJD(4) + (t160 - t68) * qJD(1)) * t135) * t134) + (t21 * t134 + t135 * t20) * t198 + t134 * ((t134 * t43 + (-t20 + t142) * qJD(1)) * t134 + (t21 * qJD(1) + (-t72 * t196 + t70 * t197 + t210) * t135 + (t44 + (t228 - t230) * qJD(4) + t161 * qJD(1)) * t134) * t135) + (t11 * t3 - t28 * t56 + t29 * t55) * t246 + t148; m(6) * (t137 * t89 - t159 * t79) + t145; -m(6) * t172; 0; m(6) * (t10 * t11 + t136 * t89 - t158 * t79 + t31 * t3) + t148; (t10 * t31 - t89 * t172) * t246 + t148;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
