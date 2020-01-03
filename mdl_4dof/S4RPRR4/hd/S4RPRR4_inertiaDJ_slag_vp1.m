% Calculate time derivative of joint inertia matrix for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:15
% EndTime: 2019-12-31 16:50:20
% DurationCPUTime: 2.77s
% Computational Cost: add. (8764->379), mult. (12007->578), div. (0->0), fcn. (11529->8), ass. (0->206)
t134 = qJ(1) + pkin(7);
t131 = sin(t134);
t136 = sin(qJ(3));
t139 = cos(qJ(3));
t132 = cos(t134);
t219 = Icges(4,4) * t139;
t162 = -Icges(4,2) * t136 + t219;
t90 = Icges(4,6) * t131 + t162 * t132;
t220 = Icges(4,4) * t136;
t165 = Icges(4,1) * t139 - t220;
t92 = Icges(4,5) * t131 + t165 * t132;
t166 = t136 * t90 - t139 * t92;
t247 = t166 * t131;
t89 = -Icges(4,6) * t132 + t162 * t131;
t91 = -Icges(4,5) * t132 + t165 * t131;
t167 = t136 * t89 - t139 * t91;
t246 = t167 * t132;
t207 = t132 * t136;
t245 = -rSges(4,2) * t207 + t131 * rSges(4,3);
t133 = cos(qJ(1)) * pkin(1);
t244 = t131 * pkin(5) + t133;
t159 = Icges(4,5) * t139 - Icges(4,6) * t136;
t87 = -Icges(4,3) * t132 + t159 * t131;
t243 = 2 * m(4);
t242 = 2 * m(5);
t241 = t131 ^ 2;
t240 = t132 ^ 2;
t239 = t131 / 0.2e1;
t237 = t132 / 0.2e1;
t236 = -t139 / 0.2e1;
t235 = -rSges(5,3) - pkin(6);
t119 = rSges(4,1) * t136 + rSges(4,2) * t139;
t234 = m(4) * t119;
t233 = sin(qJ(1)) * pkin(1);
t232 = pkin(3) * t136;
t231 = pkin(3) * t139;
t230 = qJD(1) / 0.2e1;
t229 = rSges(4,1) * t139;
t228 = rSges(4,3) * t132;
t227 = t136 * t91;
t226 = t136 * t92;
t225 = t139 * t89;
t224 = t139 * t90;
t180 = pkin(6) * t136 + t231;
t135 = sin(qJ(4));
t204 = t135 * t139;
t138 = cos(qJ(4));
t206 = t132 * t138;
t98 = -t131 * t204 - t206;
t203 = t138 * t139;
t208 = t132 * t135;
t99 = t131 * t203 - t208;
t181 = -rSges(5,1) * t99 - rSges(5,2) * t98;
t210 = t131 * t136;
t69 = rSges(5,3) * t210 - t181;
t223 = t131 * t180 + t69;
t205 = t132 * t139;
t209 = t131 * t138;
t100 = -t132 * t204 + t209;
t211 = t131 * t135;
t101 = t132 * t203 + t211;
t70 = t101 * rSges(5,1) + t100 * rSges(5,2) + rSges(5,3) * t207;
t222 = pkin(3) * t205 + pkin(6) * t207 + t70;
t177 = rSges(5,1) * t138 - rSges(5,2) * t135;
t195 = qJD(4) * t136;
t83 = (-rSges(5,1) * t135 - rSges(5,2) * t138) * t195 + (rSges(5,3) * t136 + t177 * t139) * qJD(3);
t221 = -t180 * qJD(3) - t83;
t218 = Icges(5,4) * t135;
t217 = Icges(5,4) * t138;
t88 = Icges(4,3) * t131 + t159 * t132;
t213 = qJD(1) * t88;
t163 = Icges(5,1) * t138 - t218;
t104 = -Icges(5,5) * t139 + t163 * t136;
t212 = t104 * t138;
t105 = -rSges(5,3) * t139 + t177 * t136;
t122 = -pkin(6) * t139 + t232;
t202 = t105 + t122;
t201 = qJD(1) * t131;
t200 = qJD(1) * t132;
t199 = qJD(1) * t136;
t198 = qJD(3) * t131;
t197 = qJD(3) * t136;
t196 = qJD(3) * t139;
t186 = t132 * t196;
t184 = -qJD(4) * t139 + qJD(1);
t145 = t135 * t197 + t184 * t138;
t183 = qJD(1) * t139 - qJD(4);
t59 = t145 * t132 + t183 * t211;
t144 = t184 * t135 - t138 * t197;
t60 = t144 * t132 - t183 * t209;
t194 = t60 * rSges(5,1) + t59 * rSges(5,2) + rSges(5,3) * t186;
t65 = Icges(5,4) * t99 + Icges(5,2) * t98 + Icges(5,6) * t210;
t67 = Icges(5,1) * t99 + Icges(5,4) * t98 + Icges(5,5) * t210;
t171 = -t135 * t65 + t138 * t67;
t63 = Icges(5,5) * t99 + Icges(5,6) * t98 + Icges(5,3) * t210;
t29 = t171 * t136 - t139 * t63;
t158 = Icges(5,5) * t138 - Icges(5,6) * t135;
t102 = -Icges(5,3) * t139 + t158 * t136;
t160 = -Icges(5,2) * t135 + t217;
t103 = -Icges(5,6) * t139 + t160 * t136;
t41 = t102 * t210 + t103 * t98 + t104 * t99;
t192 = t29 / 0.2e1 + t41 / 0.2e1;
t66 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t207;
t68 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t207;
t170 = -t135 * t66 + t138 * t68;
t64 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t207;
t30 = t170 * t136 - t139 * t64;
t42 = t100 * t103 + t101 * t104 + t102 * t207;
t191 = t30 / 0.2e1 + t42 / 0.2e1;
t190 = t132 * pkin(2) + t244;
t189 = t131 * t199;
t188 = t132 * t199;
t187 = t103 * t196;
t185 = t196 / 0.2e1;
t85 = t202 * t132;
t61 = t145 * t131 - t183 * t208;
t62 = t144 * t131 + t183 * t206;
t182 = t62 * rSges(5,1) + t61 * rSges(5,2);
t178 = -rSges(4,2) * t136 + t229;
t22 = t63 * t210 + t65 * t98 + t67 * t99;
t23 = t64 * t210 + t66 * t98 + t68 * t99;
t15 = t131 * t23 - t132 * t22;
t176 = t131 * t22 + t132 * t23;
t24 = t100 * t65 + t101 * t67 + t63 * t207;
t25 = t100 * t66 + t101 * t68 + t64 * t207;
t16 = t131 * t25 - t132 * t24;
t175 = t131 * t24 + t132 * t25;
t174 = t131 * t30 - t132 * t29;
t173 = t131 * t29 + t132 * t30;
t172 = -t131 * t70 + t132 * t69;
t164 = Icges(4,1) * t136 + t219;
t161 = Icges(4,2) * t139 + t220;
t157 = -t132 * pkin(3) * t197 + pkin(6) * t186;
t94 = rSges(4,1) * t205 + t245;
t156 = -pkin(2) - t178;
t146 = t186 - t189;
t32 = Icges(5,5) * t60 + Icges(5,6) * t59 + t146 * Icges(5,3);
t155 = t136 * t32 + t64 * t196;
t147 = t131 * t196 + t188;
t33 = Icges(5,5) * t62 + Icges(5,6) * t61 + t147 * Icges(5,3);
t154 = t136 * t33 + t63 * t196;
t153 = qJD(3) * t119;
t152 = t235 * t136 - pkin(2) - t231;
t151 = qJD(3) * t164;
t150 = qJD(3) * t161;
t149 = qJD(3) * (-Icges(4,5) * t136 - Icges(4,6) * t139);
t80 = (-Icges(5,5) * t135 - Icges(5,6) * t138) * t195 + (Icges(5,3) * t136 + t158 * t139) * qJD(3);
t82 = (-Icges(5,1) * t135 - t217) * t195 + (Icges(5,5) * t136 + t163 * t139) * qJD(3);
t143 = t102 * t197 - t139 * t80 + t196 * t212 + (-t103 * t195 + t136 * t82) * t138;
t142 = rSges(4,2) * t189 + rSges(4,3) * t200 - t132 * t153;
t141 = t152 * t131 - t233;
t129 = t132 * pkin(5);
t126 = pkin(5) * t200;
t111 = t178 * qJD(3);
t93 = t178 * t131 - t228;
t84 = t202 * t131;
t81 = (-Icges(5,2) * t138 - t218) * t195 + (Icges(5,6) * t136 + t160 * t139) * qJD(3);
t79 = t94 + t190;
t78 = t156 * t131 + t129 + t228 - t233;
t72 = t131 * t149 + t213;
t71 = -qJD(1) * t87 + t132 * t149;
t56 = t119 * t198 + (-t133 + (-rSges(4,3) - pkin(5)) * t131 + t156 * t132) * qJD(1);
t55 = t126 + (-t233 + (-pkin(2) - t229) * t131) * qJD(1) + t142;
t54 = -t102 * t139 + (-t103 * t135 + t212) * t136;
t53 = t54 * t197;
t52 = -t105 * t207 - t139 * t70;
t51 = t105 * t210 + t139 * t69;
t50 = t190 + t222;
t49 = t129 + t141 + t181;
t48 = t131 * t88 - t166 * t132;
t47 = t131 * t87 - t246;
t46 = -t132 * t88 - t247;
t45 = -t167 * t131 - t132 * t87;
t44 = -qJD(1) * t85 + t221 * t131;
t43 = t221 * t132 + t202 * t201;
t40 = t172 * t136;
t39 = t147 * rSges(5,3) + t182;
t38 = -rSges(5,3) * t189 + t194;
t37 = Icges(5,1) * t62 + Icges(5,4) * t61 + t147 * Icges(5,5);
t36 = Icges(5,1) * t60 + Icges(5,4) * t59 + t146 * Icges(5,5);
t35 = Icges(5,4) * t62 + Icges(5,2) * t61 + t147 * Icges(5,6);
t34 = Icges(5,4) * t60 + Icges(5,2) * t59 + t146 * Icges(5,6);
t31 = t223 * t131 + t222 * t132;
t28 = (qJD(1) * t93 + t142) * t132 + (-t131 * t153 + (-t94 + t245) * qJD(1)) * t131;
t27 = (t235 * t139 + t232) * t198 + (t152 * t132 - t244) * qJD(1) - t182;
t26 = t141 * qJD(1) + t126 + t157 + t194;
t21 = (-t187 + (-qJD(4) * t104 - t81) * t136) * t135 + t143;
t20 = (t105 * t198 + t39) * t139 + (-qJD(3) * t69 + t105 * t200 + t131 * t83) * t136;
t19 = (-qJD(3) * t105 * t132 - t38) * t139 + (qJD(3) * t70 + t105 * t201 - t132 * t83) * t136;
t18 = t147 * t102 + t103 * t61 + t104 * t62 + t80 * t210 + t81 * t98 + t82 * t99;
t17 = t100 * t81 + t101 * t82 + t146 * t102 + t103 * t59 + t104 * t60 + t80 * t207;
t14 = t172 * t196 + (-t131 * t38 + t132 * t39 + (-t131 * t69 - t132 * t70) * qJD(1)) * t136;
t13 = t175 * t136 - t139 * t42;
t12 = t176 * t136 - t139 * t41;
t11 = (t157 + t38) * t132 + (-t122 * t198 + t39) * t131 + (-t222 * t131 + t223 * t132) * qJD(1);
t10 = (t170 * qJD(3) - t32) * t139 + (qJD(3) * t64 - t135 * t34 + t138 * t36 + (-t135 * t68 - t138 * t66) * qJD(4)) * t136;
t9 = (t171 * qJD(3) - t33) * t139 + (qJD(3) * t63 - t135 * t35 + t138 * t37 + (-t135 * t67 - t138 * t65) * qJD(4)) * t136;
t8 = t155 * t131 + t64 * t188 + t34 * t98 + t36 * t99 + t61 * t66 + t62 * t68;
t7 = t154 * t131 + t63 * t188 + t35 * t98 + t37 * t99 + t61 * t65 + t62 * t67;
t6 = t100 * t34 + t101 * t36 + t155 * t132 - t64 * t189 + t59 * t66 + t60 * t68;
t5 = t100 * t35 + t101 * t37 + t154 * t132 - t63 * t189 + t59 * t65 + t60 * t67;
t4 = t176 * qJD(1) + t131 * t8 - t132 * t7;
t3 = t175 * qJD(1) + t131 * t6 - t132 * t5;
t2 = (t176 * qJD(3) - t18) * t139 + (-t15 * qJD(1) + qJD(3) * t41 + t131 * t7 + t132 * t8) * t136;
t1 = (t175 * qJD(3) - t17) * t139 + (-t16 * qJD(1) + qJD(3) * t42 + t131 * t5 + t132 * t6) * t136;
t57 = [(t55 * t79 + t56 * t78) * t243 + (t26 * t50 + t27 * t49) * t242 + t143 + (t165 - t161) * t197 + (t164 + t162) * t196 + (-t104 * t195 - t136 * t81 - t187) * t135; 0; 0; m(4) * ((-t131 * t55 - t132 * t56) * t119 + (-t131 * t79 - t132 * t78) * t111) + m(5) * (-t26 * t84 - t27 * t85 + t43 * t49 + t44 * t50) + (t241 / 0.2e1 + t240 / 0.2e1) * t159 * qJD(3) + ((t226 / 0.2e1 + t224 / 0.2e1 - t79 * t234 + t191) * t132 + (t227 / 0.2e1 + t225 / 0.2e1 + t78 * t234 + t192) * t131) * qJD(1) + (-t166 * qJD(3) + t136 * (-t91 * qJD(1) - t132 * t151) + t139 * (-t89 * qJD(1) - t132 * t150) + t10 + t17) * t239 - (-t167 * qJD(3) + t136 * (t92 * qJD(1) - t131 * t151) + t139 * (t90 * qJD(1) - t131 * t150) + t18 + t9) * t132 / 0.2e1; m(4) * t28 + m(5) * t11; ((t131 * t93 + t132 * t94) * t28 + (t240 + t241) * t119 * t111) * t243 + t131 * ((t131 * t71 + (t47 + t247) * qJD(1)) * t131 + (t48 * qJD(1) + (t89 * t196 + t91 * t197) * t132 + (-t72 + (-t224 - t226) * qJD(3) + (-t167 + t88) * qJD(1)) * t131) * t132) - t132 * ((t132 * t72 + (t46 + t246) * qJD(1)) * t132 + (t45 * qJD(1) + (-t90 * t196 - t92 * t197 + t213) * t131 + (-t71 + (t225 + t227) * qJD(3) - t166 * qJD(1)) * t132) * t131) + (t11 * t31 - t43 * t85 - t44 * t84) * t242 + t131 * t3 - t132 * t4 + (t131 * t46 - t132 * t45 + t15) * t201 + (t131 * t48 - t132 * t47 + t16) * t200; t53 + m(5) * (t19 * t50 + t20 * t49 + t26 * t52 + t27 * t51) + (-t21 + (t192 * t131 + t191 * t132) * qJD(3)) * t139 + ((t10 / 0.2e1 + t17 / 0.2e1) * t132 + (t9 / 0.2e1 + t18 / 0.2e1) * t131 + (-t191 * t131 + t192 * t132) * qJD(1)) * t136; m(5) * t14; m(5) * (t11 * t40 + t14 * t31 - t19 * t84 - t20 * t85 + t43 * t51 + t44 * t52) + (t13 * t230 + t16 * t185 - t2 / 0.2e1 + (qJD(1) * t30 - t9) * t236) * t132 + (t1 / 0.2e1 + t12 * t230 + t15 * t185 + (qJD(1) * t29 + t10) * t236) * t131 + (t3 * t237 + t4 * t239 + qJD(3) * t174 / 0.2e1 + (-t131 * t16 / 0.2e1 + t15 * t237) * qJD(1)) * t136; (t14 * t40 + t19 * t52 + t20 * t51) * t242 + (t21 * t139 - t53 + (t131 * t12 + t132 * t13 - t173 * t139) * qJD(3)) * t139 + (t132 * t1 + t131 * t2 - t139 * (t10 * t132 + t9 * t131) + (t173 * t136 - t54 * t139) * qJD(3) + (t132 * t12 - t131 * t13 + t139 * t174) * qJD(1)) * t136;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t57(1), t57(2), t57(4), t57(7); t57(2), t57(3), t57(5), t57(8); t57(4), t57(5), t57(6), t57(9); t57(7), t57(8), t57(9), t57(10);];
Mq = res;
