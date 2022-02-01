% Calculate time derivative of joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:20
% DurationCPUTime: 4.18s
% Computational Cost: add. (5212->383), mult. (7260->542), div. (0->0), fcn. (6751->8), ass. (0->181)
t139 = qJ(1) + pkin(7);
t136 = sin(t139);
t141 = cos(pkin(8));
t145 = cos(qJ(4));
t187 = t141 * t145;
t137 = cos(t139);
t143 = sin(qJ(4));
t190 = t137 * t143;
t109 = t136 * t187 - t190;
t188 = t141 * t143;
t155 = t136 * t188 + t137 * t145;
t140 = sin(pkin(8));
t193 = t136 * t140;
t65 = Icges(6,5) * t109 - Icges(6,6) * t155 + Icges(6,3) * t193;
t67 = Icges(5,5) * t109 - Icges(5,6) * t155 + Icges(5,3) * t193;
t226 = t65 + t67;
t174 = t137 * t188;
t110 = t136 * t145 - t174;
t192 = t136 * t143;
t111 = t137 * t187 + t192;
t191 = t137 * t140;
t66 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t191;
t68 = Icges(5,5) * t111 + Icges(5,6) * t110 + Icges(5,3) * t191;
t225 = t66 + t68;
t70 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t191;
t72 = Icges(5,4) * t111 + Icges(5,2) * t110 + Icges(5,6) * t191;
t224 = t70 + t72;
t69 = Icges(6,4) * t109 - Icges(6,2) * t155 + Icges(6,6) * t193;
t71 = Icges(5,4) * t109 - Icges(5,2) * t155 + Icges(5,6) * t193;
t223 = t71 + t69;
t74 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t191;
t76 = Icges(5,1) * t111 + Icges(5,4) * t110 + Icges(5,5) * t191;
t222 = t74 + t76;
t73 = Icges(6,1) * t109 - Icges(6,4) * t155 + Icges(6,5) * t193;
t75 = Icges(5,1) * t109 - Icges(5,4) * t155 + Icges(5,5) * t193;
t221 = t75 + t73;
t195 = Icges(6,4) * t145;
t102 = -Icges(6,6) * t141 + (-Icges(6,2) * t143 + t195) * t140;
t197 = Icges(5,4) * t145;
t103 = -Icges(5,6) * t141 + (-Icges(5,2) * t143 + t197) * t140;
t241 = -t102 - t103;
t196 = Icges(6,4) * t143;
t104 = -Icges(6,5) * t141 + (Icges(6,1) * t145 - t196) * t140;
t198 = Icges(5,4) * t143;
t105 = -Icges(5,5) * t141 + (Icges(5,1) * t145 - t198) * t140;
t240 = -t105 - t104;
t178 = qJD(4) * t140;
t114 = (-Icges(6,2) * t145 - t196) * t178;
t115 = (-Icges(5,2) * t145 - t198) * t178;
t239 = (-t114 - t115) * t143;
t112 = (-Icges(6,5) * t143 - Icges(6,6) * t145) * t178;
t113 = (-Icges(5,5) * t143 - Icges(5,6) * t145) * t178;
t116 = (-Icges(6,1) * t143 - t195) * t178;
t117 = (-Icges(5,1) * t143 - t197) * t178;
t237 = (t116 + t117) * t140 * t145 + (-t112 - t113) * t141;
t238 = (t237 + ((t240 * t143 + t241 * t145) * qJD(4) + t239) * t140) * t141;
t236 = t221 * t109 - t223 * t155 + t226 * t193;
t235 = -t222 * t109 + t224 * t155 - t225 * t193;
t234 = -t224 * t110 - t222 * t111 - t225 * t191;
t233 = t223 * t110 + t221 * t111 + t226 * t191;
t179 = qJD(1) * t140;
t171 = t137 * t179;
t89 = qJD(1) * t110 - qJD(4) * t109;
t150 = t155 * qJD(4);
t90 = qJD(1) * t111 - t150;
t42 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t171;
t44 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t171;
t231 = t44 + t42;
t172 = t136 * t179;
t87 = qJD(1) * t155 - qJD(4) * t111;
t88 = -qJD(1) * t109 + qJD(4) * t110;
t45 = Icges(6,4) * t88 + Icges(6,2) * t87 - Icges(6,6) * t172;
t47 = Icges(5,4) * t88 + Icges(5,2) * t87 - Icges(5,6) * t172;
t230 = t47 + t45;
t46 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t171;
t48 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t171;
t229 = t48 + t46;
t49 = Icges(6,1) * t88 + Icges(6,4) * t87 - Icges(6,5) * t172;
t51 = Icges(5,1) * t88 + Icges(5,4) * t87 - Icges(5,5) * t172;
t228 = t51 + t49;
t50 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t171;
t52 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t171;
t227 = t52 + t50;
t128 = pkin(4) * t192;
t135 = pkin(4) * t145 + pkin(3);
t194 = t135 * t141;
t220 = -t111 * rSges(6,1) - t110 * rSges(6,2) - rSges(6,3) * t191 - t137 * t194 - t128;
t208 = pkin(3) - t135;
t219 = t141 * t208;
t138 = cos(qJ(1)) * pkin(1);
t218 = t136 * qJ(3) + t138;
t142 = -qJ(5) - pkin(6);
t175 = pkin(4) * t190;
t176 = qJD(5) * t140;
t177 = qJD(4) * t145;
t200 = pkin(4) * qJD(4);
t217 = -t136 * pkin(4) * t177 - t88 * rSges(6,1) - t87 * rSges(6,2) - qJD(1) * t175 - t137 * t176 - t142 * t172 + t174 * t200;
t216 = rSges(6,1) * t109 - rSges(6,2) * t155 - t175;
t215 = -t90 * rSges(6,1) - t89 * rSges(6,2) - t136 * t176;
t41 = Icges(6,5) * t88 + Icges(6,6) * t87 - Icges(6,3) * t172;
t43 = Icges(5,5) * t88 + Icges(5,6) * t87 - Icges(5,3) * t172;
t214 = (t43 + t41) * t137;
t213 = 2 * m(5);
t212 = 2 * m(6);
t211 = sin(qJ(1)) * pkin(1);
t210 = pkin(3) * t141;
t209 = pkin(6) * t140;
t207 = pkin(6) + t142;
t181 = qJD(1) * t136;
t206 = rSges(6,3) * t172 - (t209 + t219) * t181 + t217;
t167 = t207 * t140;
t151 = -t167 - t219;
t205 = rSges(6,3) * t171 - pkin(4) * t150 + (t137 * t151 + t128) * qJD(1) - t215;
t204 = rSges(6,3) * t193 + t136 * t151 + t216;
t203 = -(-t167 - t210) * t137 + t220;
t201 = t88 * rSges(5,1) + t87 * rSges(5,2);
t199 = (t207 - rSges(6,3)) * t141 + (rSges(6,1) * t145 - rSges(6,2) * t143 - t208) * t140;
t170 = t143 * t178;
t184 = (-rSges(6,1) * t143 - rSges(6,2) * t145) * t178 - pkin(4) * t170 - qJD(5) * t141;
t180 = qJD(1) * t137;
t182 = qJ(3) * t180 + qJD(3) * t136;
t82 = t111 * rSges(5,1) + t110 * rSges(5,2) + rSges(5,3) * t191;
t173 = t137 * pkin(2) + t218;
t166 = -pkin(2) - t194;
t163 = qJD(1) * t199;
t162 = -t90 * rSges(5,1) - t89 * rSges(5,2);
t159 = rSges(4,1) * t141 - rSges(4,2) * t140;
t158 = -rSges(5,1) * t109 + rSges(5,2) * t155;
t156 = -pkin(2) - t159;
t154 = -t210 - pkin(2) + (-rSges(5,3) - pkin(6)) * t140;
t152 = (-rSges(6,3) + t142) * t140 + t166;
t148 = t136 * t154 - t211;
t147 = rSges(4,3) * t137 + t136 * t156 - t211;
t133 = t137 * qJ(3);
t131 = qJD(3) * t137;
t119 = (-rSges(5,1) * t143 - rSges(5,2) * t145) * t178;
t107 = -rSges(5,3) * t141 + (rSges(5,1) * t145 - rSges(5,2) * t143) * t140;
t101 = -Icges(5,3) * t141 + (Icges(5,5) * t145 - Icges(5,6) * t143) * t140;
t100 = -Icges(6,3) * t141 + (Icges(6,5) * t145 - Icges(6,6) * t143) * t140;
t92 = rSges(4,3) * t136 + t137 * t159 + t173;
t91 = t133 + t147;
t80 = rSges(5,3) * t193 - t158;
t64 = t131 + (-t138 + (-rSges(4,3) - qJ(3)) * t136 + t156 * t137) * qJD(1);
t63 = qJD(1) * t147 + t182;
t62 = -t107 * t191 - t141 * t82;
t61 = t107 * t193 + t141 * t80;
t60 = (t209 + t210) * t137 + t173 + t82;
t59 = t133 + t148 + t158;
t56 = rSges(5,3) * t171 - t162;
t54 = -rSges(5,3) * t172 + t201;
t40 = -t142 * t191 + t173 - t220;
t39 = t136 * t152 + t133 - t211 - t216;
t38 = t101 * t191 + t103 * t110 + t105 * t111;
t37 = t100 * t191 + t102 * t110 + t104 * t111;
t36 = t101 * t193 - t103 * t155 + t105 * t109;
t35 = t100 * t193 - t102 * t155 + t104 * t109;
t32 = t131 + (t137 * t154 - t218) * qJD(1) + t162;
t31 = qJD(1) * t148 + t182 + t201;
t30 = t141 * t56 + (t107 * t180 + t119 * t136) * t140;
t29 = -t141 * t54 + (t107 * t181 - t119 * t137) * t140;
t28 = t141 * t203 - t191 * t199;
t27 = t141 * t204 + t193 * t199;
t26 = -t141 * t68 + (-t143 * t72 + t145 * t76) * t140;
t25 = -t141 * t67 + (-t143 * t71 + t145 * t75) * t140;
t24 = -t141 * t66 + (-t143 * t70 + t145 * t74) * t140;
t23 = -t141 * t65 + (-t143 * t69 + t145 * t73) * t140;
t22 = t131 + t155 * t200 + (-t138 + (-pkin(4) * t143 - qJ(3)) * t136 + t152 * t137) * qJD(1) + t215;
t21 = (-t211 + (-rSges(6,3) * t140 + t166) * t136) * qJD(1) + t182 - t217;
t12 = t103 * t89 + t105 * t90 - t155 * t115 + t109 * t117 + (t101 * t180 + t113 * t136) * t140;
t11 = t102 * t89 + t104 * t90 - t155 * t114 + t109 * t116 + (t100 * t180 + t112 * t136) * t140;
t10 = t103 * t87 + t105 * t88 + t110 * t115 + t111 * t117 + (-t101 * t181 + t113 * t137) * t140;
t9 = t102 * t87 + t104 * t88 + t110 * t114 + t111 * t116 + (-t100 * t181 + t112 * t137) * t140;
t8 = t205 * t141 + (t136 * t184 + t137 * t163) * t140;
t7 = t206 * t141 + (t136 * t163 - t137 * t184) * t140;
t6 = (-t136 * t54 + t137 * t56 + (-t136 * t80 - t137 * t82) * qJD(1)) * t140;
t5 = -t141 * t43 + (-t143 * t47 + t145 * t51 + (-t143 * t76 - t145 * t72) * qJD(4)) * t140;
t4 = -t141 * t44 + (-t143 * t48 + t145 * t52 + (-t143 * t75 - t145 * t71) * qJD(4)) * t140;
t3 = -t141 * t41 + (-t143 * t45 + t145 * t49 + (-t143 * t74 - t145 * t70) * qJD(4)) * t140;
t2 = -t141 * t42 + (-t143 * t46 + t145 * t50 + (-t143 * t73 - t145 * t69) * qJD(4)) * t140;
t1 = (t205 * t137 + t206 * t136 + (-t136 * t204 + t137 * t203) * qJD(1)) * t140;
t13 = [(t21 * t40 + t22 * t39) * t212 + (t31 * t60 + t32 * t59) * t213 + 0.2e1 * m(4) * (t63 * t92 + t64 * t91) + t240 * t170 + t237 + (t241 * t177 + t239) * t140; 0; 0; m(6) * (t136 * t22 - t137 * t21 + (t136 * t40 + t137 * t39) * qJD(1)) + m(5) * (t136 * t32 - t137 * t31 + (t136 * t60 + t137 * t59) * qJD(1)) + m(4) * (t136 * t64 - t137 * t63 + (t136 * t92 + t137 * t91) * qJD(1)); 0; 0; -t238 + m(6) * (t21 * t28 + t22 * t27 + t39 * t8 + t40 * t7) + m(5) * (t29 * t60 + t30 * t59 + t31 * t62 + t32 * t61) + ((t10 / 0.2e1 + t9 / 0.2e1 + t5 / 0.2e1 + t3 / 0.2e1) * t137 + (t4 / 0.2e1 + t2 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t136 + ((t25 / 0.2e1 + t23 / 0.2e1 + t36 / 0.2e1 + t35 / 0.2e1) * t137 + (-t38 / 0.2e1 - t37 / 0.2e1 - t26 / 0.2e1 - t24 / 0.2e1) * t136) * qJD(1)) * t140; m(5) * t6 + m(6) * t1; m(5) * (t136 * t30 - t137 * t29 + (t136 * t62 + t137 * t61) * qJD(1)) + m(6) * (t136 * t8 - t137 * t7 + (t136 * t28 + t137 * t27) * qJD(1)); (t62 * t29 + t61 * t30) * t213 + (t27 * t8 + t28 * t7) * t212 + ((-t136 * t82 + t137 * t80) * t6 * t213 + (t136 * t203 + t204 * t137) * t1 * t212 + (-t233 * t136 + t234 * t137) * t172 + (t236 * t136 - t235 * t137) * t171 + (t235 * t181 + t236 * t180 + (t225 * t140 * t180 + t228 * t109 - t155 * t230 + t222 * t90 + t224 * t89) * t137 + (t221 * t90 + t223 * t89 + (t231 * t136 + t226 * t180 + t214) * t140 + t227 * t109 - t229 * t155) * t136) * t193 + (t234 * t181 + t233 * t180 + (t222 * t88 + t224 * t87 + (-t225 * t181 + t214) * t140 + t228 * t111 + t230 * t110) * t137 + ((t231 * t137 - t226 * t181) * t140 + t221 * t88 + t223 * t87 + t227 * t111 + t229 * t110) * t136) * t191) * t140 + (t238 + (-t12 - t11) * t193 + (-t10 - t9) * t191 + (t38 + t37) * t172 + (-t36 - t35) * t171 + ((-t3 - t5) * t137 + (-t2 - t4) * t136 + ((-t23 - t25) * t137 + (t24 + t26) * t136) * qJD(1)) * t140) * t141; m(6) * (t136 * t21 + t137 * t22 + (-t136 * t39 + t137 * t40) * qJD(1)) * t140; 0; 0; m(6) * (-t1 * t141 + (t136 * t7 + t137 * t8 + (-t136 * t27 + t137 * t28) * qJD(1)) * t140); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
