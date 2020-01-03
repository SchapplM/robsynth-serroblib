% Calculate time derivative of joint inertia matrix for
% S4RPRR7
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:43
% EndTime: 2019-12-31 16:53:48
% DurationCPUTime: 3.67s
% Computational Cost: add. (8274->407), mult. (12273->623), div. (0->0), fcn. (11743->8), ass. (0->213)
t142 = sin(qJ(1));
t144 = cos(qJ(1));
t137 = pkin(7) + qJ(3);
t133 = cos(t137);
t132 = sin(t137);
t221 = Icges(4,4) * t132;
t165 = Icges(4,1) * t133 - t221;
t100 = Icges(4,5) * t142 + t165 * t144;
t220 = Icges(4,4) * t133;
t162 = -Icges(4,2) * t132 + t220;
t98 = Icges(4,6) * t142 + t162 * t144;
t166 = t100 * t133 - t132 * t98;
t256 = t166 * t142;
t97 = -Icges(4,6) * t144 + t162 * t142;
t99 = -Icges(4,5) * t144 + t165 * t142;
t174 = t132 * t97 - t133 * t99;
t255 = t174 * t144;
t136 = t142 * rSges(4,3);
t211 = t132 * t144;
t254 = -rSges(4,2) * t211 + t136;
t140 = -pkin(5) - qJ(2);
t139 = cos(pkin(7));
t130 = pkin(2) * t139 + pkin(1);
t236 = rSges(4,1) * t133;
t179 = -rSges(4,2) * t132 + t236;
t154 = -t130 - t179;
t85 = (rSges(4,3) - t140) * t144 + t154 * t142;
t210 = t133 * t144;
t102 = rSges(4,1) * t210 + t254;
t184 = t144 * t130 - t142 * t140;
t86 = t102 + t184;
t253 = t142 * t86 + t144 * t85;
t141 = sin(qJ(4));
t143 = cos(qJ(4));
t218 = Icges(5,4) * t143;
t160 = -Icges(5,2) * t141 + t218;
t92 = -Icges(5,6) * t133 + t160 * t132;
t219 = Icges(5,4) * t141;
t163 = Icges(5,1) * t143 - t219;
t93 = -Icges(5,5) * t133 + t163 * t132;
t252 = -t141 * t93 - t143 * t92;
t159 = Icges(4,5) * t133 - Icges(4,6) * t132;
t95 = -Icges(4,3) * t144 + t159 * t142;
t182 = qJD(1) * t133 - qJD(4);
t197 = qJD(3) * t144;
t188 = t132 * t197;
t251 = t182 * t142 + t188;
t198 = qJD(3) * t143;
t201 = qJD(3) * t132;
t231 = t141 * t92;
t158 = Icges(5,5) * t143 - Icges(5,6) * t141;
t196 = qJD(4) * t132;
t59 = (-Icges(5,5) * t141 - Icges(5,6) * t143) * t196 + (Icges(5,3) * t132 + t158 * t133) * qJD(3);
t61 = (-Icges(5,1) * t141 - t218) * t196 + (Icges(5,5) * t132 + t163 * t133) * qJD(3);
t91 = -Icges(5,3) * t133 + t158 * t132;
t250 = t132 * t143 * t61 + t91 * t201 + (-qJD(3) * t231 + t93 * t198 - t59) * t133;
t155 = rSges(3,1) * t139 - rSges(3,2) * sin(pkin(7)) + pkin(1);
t226 = rSges(3,3) + qJ(2);
t90 = t226 * t142 + t155 * t144;
t249 = 2 * m(4);
t248 = 2 * m(5);
t247 = t142 ^ 2;
t246 = t144 ^ 2;
t245 = -t133 / 0.2e1;
t244 = t142 / 0.2e1;
t242 = t144 / 0.2e1;
t241 = -rSges(5,3) - pkin(6);
t119 = rSges(4,1) * t132 + rSges(4,2) * t133;
t240 = m(4) * t119;
t239 = pkin(3) * t132;
t238 = pkin(3) * t133;
t237 = qJD(1) / 0.2e1;
t235 = t132 * t99;
t234 = t133 * t97;
t233 = t133 * t98;
t60 = (-Icges(5,2) * t143 - t219) * t196 + (Icges(5,6) * t132 + t160 * t133) * qJD(3);
t232 = t141 * t60;
t180 = pkin(6) * t132 + t238;
t205 = t143 * t144;
t207 = t142 * t141;
t107 = -t133 * t207 - t205;
t206 = t142 * t143;
t208 = t141 * t144;
t108 = t133 * t206 - t208;
t178 = -t108 * rSges(5,1) - t107 * rSges(5,2);
t212 = t132 * t142;
t79 = rSges(5,3) * t212 - t178;
t225 = t180 * t142 + t79;
t109 = -t133 * t208 + t206;
t110 = t133 * t205 + t207;
t80 = t110 * rSges(5,1) + t109 * rSges(5,2) + rSges(5,3) * t211;
t224 = pkin(3) * t210 + pkin(6) * t211 + t80;
t177 = rSges(5,1) * t143 - rSges(5,2) * t141;
t66 = (-rSges(5,1) * t141 - rSges(5,2) * t143) * t196 + (rSges(5,3) * t132 + t177 * t133) * qJD(3);
t223 = -t180 * qJD(3) - t66;
t120 = -pkin(6) * t133 + t239;
t94 = -rSges(5,3) * t133 + t177 * t132;
t222 = t120 + t94;
t96 = Icges(4,3) * t142 + t159 * t144;
t214 = qJD(1) * t96;
t213 = t100 * t132;
t209 = t140 * t144;
t135 = qJD(2) * t144;
t203 = qJD(1) * t142;
t204 = t140 * t203 + t135;
t202 = qJD(1) * t144;
t200 = qJD(3) * t133;
t199 = qJD(3) * t142;
t186 = t133 * t197;
t183 = -qJD(4) * t133 + qJD(1);
t157 = t183 * t144;
t62 = t251 * t141 + t143 * t157;
t63 = t141 * t157 - t251 * t143;
t194 = t63 * rSges(5,1) + t62 * rSges(5,2) + rSges(5,3) * t186;
t75 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t212;
t77 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t212;
t173 = -t141 * t75 + t143 * t77;
t73 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t212;
t28 = t173 * t132 - t133 * t73;
t39 = t107 * t92 + t108 * t93 + t91 * t212;
t192 = t28 / 0.2e1 + t39 / 0.2e1;
t76 = Icges(5,4) * t110 + Icges(5,2) * t109 + Icges(5,6) * t211;
t78 = Icges(5,1) * t110 + Icges(5,4) * t109 + Icges(5,5) * t211;
t172 = -t141 * t76 + t143 * t78;
t74 = Icges(5,5) * t110 + Icges(5,6) * t109 + Icges(5,3) * t211;
t29 = t172 * t132 - t133 * t74;
t40 = t109 * t92 + t110 * t93 + t91 * t211;
t191 = t29 / 0.2e1 + t40 / 0.2e1;
t189 = t132 * t203;
t187 = t133 * t199;
t185 = t200 / 0.2e1;
t82 = t222 * t144;
t64 = t183 * t206 + (t132 * t199 - t182 * t144) * t141;
t65 = t182 * t205 + (-t132 * t198 + t183 * t141) * t142;
t181 = t65 * rSges(5,1) + t64 * rSges(5,2);
t24 = t107 * t75 + t108 * t77 + t73 * t212;
t25 = t107 * t76 + t108 * t78 + t74 * t212;
t17 = t25 * t142 - t144 * t24;
t171 = t142 * t24 + t144 * t25;
t26 = t109 * t75 + t110 * t77 + t73 * t211;
t27 = t109 * t76 + t110 * t78 + t74 * t211;
t18 = t27 * t142 - t144 * t26;
t170 = t142 * t26 + t144 * t27;
t169 = t29 * t142 - t144 * t28;
t168 = t142 * t28 + t144 * t29;
t167 = -t142 * t80 + t144 * t79;
t164 = Icges(4,1) * t132 + t220;
t161 = Icges(4,2) * t133 + t221;
t156 = -pkin(3) * t188 + pkin(6) * t186;
t153 = qJD(3) * t119;
t152 = qJD(3) * t164;
t151 = qJD(3) * t161;
t150 = qJD(3) * (-Icges(4,5) * t132 - Icges(4,6) * t133);
t149 = t241 * t132 - t130 - t238;
t148 = t132 * t202 + t187;
t147 = t186 - t189;
t146 = rSges(4,2) * t189 + rSges(4,3) * t202 - t144 * t153;
t145 = t149 * t142 - t209;
t89 = -t155 * t142 + t226 * t144;
t134 = qJD(2) * t142;
t114 = t179 * qJD(3);
t101 = -rSges(4,3) * t144 + t142 * t179;
t84 = -t90 * qJD(1) + t135;
t83 = t89 * qJD(1) + t134;
t81 = t222 * t142;
t68 = t142 * t150 + t214;
t67 = -qJD(1) * t95 + t144 * t150;
t55 = t119 * t199 + (t154 * t144 - t136) * qJD(1) + t204;
t54 = t134 + (-t209 + (-t130 - t236) * t142) * qJD(1) + t146;
t53 = t184 + t224;
t52 = t145 + t178;
t51 = -t133 * t80 - t94 * t211;
t50 = t133 * t79 + t94 * t212;
t49 = t142 * t96 + t166 * t144;
t48 = t142 * t95 - t255;
t47 = -t144 * t96 + t256;
t46 = -t174 * t142 - t144 * t95;
t45 = -t133 * t91 + (t143 * t93 - t231) * t132;
t44 = t45 * t201;
t43 = t167 * t132;
t42 = -qJD(1) * t82 + t223 * t142;
t41 = t223 * t144 + t222 * t203;
t38 = t148 * rSges(5,3) + t181;
t37 = -rSges(5,3) * t189 + t194;
t36 = Icges(5,1) * t65 + Icges(5,4) * t64 + t148 * Icges(5,5);
t35 = Icges(5,1) * t63 + Icges(5,4) * t62 + t147 * Icges(5,5);
t34 = Icges(5,4) * t65 + Icges(5,2) * t64 + t148 * Icges(5,6);
t33 = Icges(5,4) * t63 + Icges(5,2) * t62 + t147 * Icges(5,6);
t32 = Icges(5,5) * t65 + Icges(5,6) * t64 + t148 * Icges(5,3);
t31 = Icges(5,5) * t63 + Icges(5,6) * t62 + t147 * Icges(5,3);
t30 = t225 * t142 + t224 * t144;
t23 = (t241 * t133 + t239) * t199 + t149 * t202 - t181 + t204;
t22 = t145 * qJD(1) + t134 + t156 + t194;
t21 = (t94 * t199 + t38) * t133 + (-qJD(3) * t79 + t142 * t66 + t94 * t202) * t132;
t20 = (-t94 * t197 - t37) * t133 + (qJD(3) * t80 - t144 * t66 + t94 * t203) * t132;
t19 = (t252 * qJD(4) - t232) * t132 + t250;
t16 = t91 * t187 + t107 * t60 + t108 * t61 + t64 * t92 + t65 * t93 + (t142 * t59 + t91 * t202) * t132;
t15 = t91 * t186 + t109 * t60 + t110 * t61 + t62 * t92 + t63 * t93 + (t144 * t59 - t91 * t203) * t132;
t14 = t167 * t200 + (-t142 * t37 + t144 * t38 + (-t142 * t79 - t144 * t80) * qJD(1)) * t132;
t13 = (t156 + t37) * t144 + (-t120 * t199 + t38) * t142 + (-t224 * t142 + t225 * t144) * qJD(1);
t12 = t170 * t132 - t40 * t133;
t11 = t171 * t132 - t39 * t133;
t10 = (t172 * qJD(3) - t31) * t133 + (qJD(3) * t74 - t141 * t33 + t143 * t35 + (-t141 * t78 - t143 * t76) * qJD(4)) * t132;
t9 = (t173 * qJD(3) - t32) * t133 + (qJD(3) * t73 - t141 * t34 + t143 * t36 + (-t141 * t77 - t143 * t75) * qJD(4)) * t132;
t8 = t74 * t187 + t107 * t33 + t108 * t35 + t64 * t76 + t65 * t78 + (t142 * t31 + t74 * t202) * t132;
t7 = t73 * t187 + t107 * t34 + t108 * t36 + t64 * t75 + t65 * t77 + (t142 * t32 + t73 * t202) * t132;
t6 = t74 * t186 + t109 * t33 + t110 * t35 + t62 * t76 + t63 * t78 + (t144 * t31 - t74 * t203) * t132;
t5 = t73 * t186 + t109 * t34 + t110 * t36 + t62 * t75 + t63 * t77 + (t144 * t32 - t73 * t203) * t132;
t4 = t171 * qJD(1) + t8 * t142 - t144 * t7;
t3 = t170 * qJD(1) + t6 * t142 - t144 * t5;
t2 = (t171 * qJD(3) - t16) * t133 + (-t17 * qJD(1) + qJD(3) * t39 + t142 * t7 + t144 * t8) * t132;
t1 = (t170 * qJD(3) - t15) * t133 + (-t18 * qJD(1) + qJD(3) * t40 + t142 * t5 + t144 * t6) * t132;
t56 = [0.2e1 * m(3) * (t83 * t90 + t84 * t89) + (t54 * t86 + t55 * t85) * t249 - t132 * t232 + (t22 * t53 + t23 * t52) * t248 + (t165 - t161) * t201 + (t164 + t162) * t200 + t252 * t196 + t250; m(3) * (t142 * t84 - t144 * t83 + (t142 * t90 + t144 * t89) * qJD(1)) + m(4) * (t253 * qJD(1) + t142 * t55 - t144 * t54) + m(5) * (t142 * t23 - t144 * t22 + (t142 * t53 + t144 * t52) * qJD(1)); 0; m(4) * ((-t142 * t54 - t144 * t55) * t119 - t253 * t114) + m(5) * (-t22 * t81 - t23 * t82 + t41 * t52 + t42 * t53) + (t247 / 0.2e1 + t246 / 0.2e1) * t159 * qJD(3) + ((t213 / 0.2e1 + t233 / 0.2e1 - t86 * t240 + t191) * t144 + (t235 / 0.2e1 + t234 / 0.2e1 + t85 * t240 + t192) * t142) * qJD(1) + (t166 * qJD(3) + t132 * (-t99 * qJD(1) - t144 * t152) + t133 * (-t97 * qJD(1) - t144 * t151) + t10 + t15) * t244 - (-t174 * qJD(3) + t132 * (t100 * qJD(1) - t142 * t152) + t133 * (t98 * qJD(1) - t142 * t151) + t16 + t9) * t144 / 0.2e1; m(5) * (t41 * t142 - t144 * t42 + (-t142 * t81 - t144 * t82) * qJD(1)); ((t142 * t101 + t102 * t144) * ((qJD(1) * t101 + t146) * t144 + (-t142 * t153 + (-t102 + t254) * qJD(1)) * t142) + (t246 + t247) * t119 * t114) * t249 + t142 * ((t142 * t67 + (t48 - t256) * qJD(1)) * t142 + (t49 * qJD(1) + (t97 * t200 + t99 * t201) * t144 + (-t68 + (-t213 - t233) * qJD(3) + (-t174 + t96) * qJD(1)) * t142) * t144) - t144 * ((t144 * t68 + (t47 + t255) * qJD(1)) * t144 + (t46 * qJD(1) + (-t100 * t201 - t98 * t200 + t214) * t142 + (-t67 + (t234 + t235) * qJD(3) + t166 * qJD(1)) * t144) * t142) + (t30 * t13 - t41 * t82 - t42 * t81) * t248 + t142 * t3 - t144 * t4 + (t47 * t142 - t144 * t46 + t17) * t203 + (t49 * t142 - t144 * t48 + t18) * t202; t44 + m(5) * (t20 * t53 + t21 * t52 + t22 * t51 + t23 * t50) + (-t19 + (t192 * t142 + t191 * t144) * qJD(3)) * t133 + ((t10 / 0.2e1 + t15 / 0.2e1) * t144 + (t9 / 0.2e1 + t16 / 0.2e1) * t142 + (-t191 * t142 + t192 * t144) * qJD(1)) * t132; m(5) * (t21 * t142 - t144 * t20 + (t142 * t51 + t144 * t50) * qJD(1)); m(5) * (t43 * t13 + t14 * t30 - t20 * t81 - t21 * t82 + t50 * t41 + t51 * t42) + (t12 * t237 + t18 * t185 - t2 / 0.2e1 + (qJD(1) * t29 - t9) * t245) * t144 + (t1 / 0.2e1 + t11 * t237 + t17 * t185 + (qJD(1) * t28 + t10) * t245) * t142 + (t3 * t242 + t4 * t244 + qJD(3) * t169 / 0.2e1 + (-t142 * t18 / 0.2e1 + t17 * t242) * qJD(1)) * t132; (t14 * t43 + t20 * t51 + t21 * t50) * t248 + (t19 * t133 - t44 + (t142 * t11 + t144 * t12 - t168 * t133) * qJD(3)) * t133 + (t144 * t1 + t142 * t2 - t133 * (t10 * t144 + t9 * t142) + (t168 * t132 - t45 * t133) * qJD(3) + (t144 * t11 - t142 * t12 + t133 * t169) * qJD(1)) * t132;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t56(1), t56(2), t56(4), t56(7); t56(2), t56(3), t56(5), t56(8); t56(4), t56(5), t56(6), t56(9); t56(7), t56(8), t56(9), t56(10);];
Mq = res;
