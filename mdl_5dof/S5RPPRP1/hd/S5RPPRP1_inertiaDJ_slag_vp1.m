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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:35:37
% EndTime: 2019-12-05 17:35:48
% DurationCPUTime: 4.91s
% Computational Cost: add. (5212->384), mult. (7260->545), div. (0->0), fcn. (6751->8), ass. (0->184)
t135 = qJ(1) + pkin(7);
t134 = cos(t135);
t141 = cos(qJ(4));
t133 = sin(t135);
t137 = cos(pkin(8));
t139 = sin(qJ(4));
t187 = t137 * t139;
t172 = t133 * t187;
t108 = t134 * t141 + t172;
t186 = t137 * t141;
t190 = t134 * t139;
t149 = t133 * t186 - t190;
t136 = sin(pkin(8));
t193 = t133 * t136;
t65 = -Icges(6,5) * t149 + Icges(6,6) * t108 - Icges(6,3) * t193;
t67 = -Icges(5,5) * t149 + Icges(5,6) * t108 - Icges(5,3) * t193;
t229 = t65 + t67;
t192 = t133 * t139;
t111 = t134 * t186 + t192;
t150 = -t133 * t141 + t134 * t187;
t191 = t134 * t136;
t66 = Icges(6,5) * t111 - Icges(6,6) * t150 + Icges(6,3) * t191;
t68 = Icges(5,5) * t111 - Icges(5,6) * t150 + Icges(5,3) * t191;
t228 = t66 + t68;
t69 = -Icges(6,4) * t149 + Icges(6,2) * t108 - Icges(6,6) * t193;
t71 = -Icges(5,4) * t149 + Icges(5,2) * t108 - Icges(5,6) * t193;
t227 = t71 + t69;
t70 = Icges(6,4) * t111 - Icges(6,2) * t150 + Icges(6,6) * t191;
t72 = Icges(5,4) * t111 - Icges(5,2) * t150 + Icges(5,6) * t191;
t226 = t72 + t70;
t74 = Icges(6,1) * t111 - Icges(6,4) * t150 + Icges(6,5) * t191;
t76 = Icges(5,1) * t111 - Icges(5,4) * t150 + Icges(5,5) * t191;
t225 = t74 + t76;
t73 = -Icges(6,1) * t149 + Icges(6,4) * t108 - Icges(6,5) * t193;
t75 = -Icges(5,1) * t149 + Icges(5,4) * t108 - Icges(5,5) * t193;
t224 = t75 + t73;
t195 = Icges(6,4) * t141;
t102 = -Icges(6,6) * t137 + (-Icges(6,2) * t139 + t195) * t136;
t197 = Icges(5,4) * t141;
t103 = -Icges(5,6) * t137 + (-Icges(5,2) * t139 + t197) * t136;
t244 = -t102 - t103;
t196 = Icges(6,4) * t139;
t104 = -Icges(6,5) * t137 + (Icges(6,1) * t141 - t196) * t136;
t198 = Icges(5,4) * t139;
t105 = -Icges(5,5) * t137 + (Icges(5,1) * t141 - t198) * t136;
t243 = -t105 - t104;
t177 = qJD(4) * t136;
t114 = (-Icges(6,2) * t141 - t196) * t177;
t115 = (-Icges(5,2) * t141 - t198) * t177;
t242 = (-t114 - t115) * t139;
t112 = (-Icges(6,5) * t139 - Icges(6,6) * t141) * t177;
t113 = (-Icges(5,5) * t139 - Icges(5,6) * t141) * t177;
t116 = (-Icges(6,1) * t139 - t195) * t177;
t117 = (-Icges(5,1) * t139 - t197) * t177;
t240 = (t116 + t117) * t136 * t141 + (-t113 - t112) * t137;
t241 = (t240 + ((t243 * t139 + t244 * t141) * qJD(4) + t242) * t136) * t137;
t239 = t227 * t108 - t224 * t149 - t229 * t193;
t238 = t226 * t108 - t225 * t149 - t228 * t193;
t237 = -t225 * t111 + t226 * t150 - t228 * t191;
t236 = -t224 * t111 + t227 * t150 - t229 * t191;
t178 = qJD(1) * t136;
t170 = t134 * t178;
t89 = qJD(1) * t150 + qJD(4) * t149;
t90 = -qJD(1) * t111 + qJD(4) * t108;
t42 = Icges(6,5) * t90 + Icges(6,6) * t89 - Icges(6,3) * t170;
t44 = Icges(5,5) * t90 + Icges(5,6) * t89 - Icges(5,3) * t170;
t234 = -t44 - t42;
t171 = t133 * t178;
t87 = qJD(1) * t108 - qJD(4) * t111;
t144 = t150 * qJD(4);
t88 = -qJD(1) * t149 - t144;
t45 = Icges(6,4) * t88 + Icges(6,2) * t87 - Icges(6,6) * t171;
t47 = Icges(5,4) * t88 + Icges(5,2) * t87 - Icges(5,6) * t171;
t233 = t45 + t47;
t46 = Icges(6,4) * t90 + Icges(6,2) * t89 - Icges(6,6) * t170;
t48 = Icges(5,4) * t90 + Icges(5,2) * t89 - Icges(5,6) * t170;
t232 = t48 + t46;
t49 = Icges(6,1) * t88 + Icges(6,4) * t87 - Icges(6,5) * t171;
t51 = Icges(5,1) * t88 + Icges(5,4) * t87 - Icges(5,5) * t171;
t231 = t49 + t51;
t50 = Icges(6,1) * t90 + Icges(6,4) * t89 - Icges(6,5) * t170;
t52 = Icges(5,1) * t90 + Icges(5,4) * t89 - Icges(5,5) * t170;
t230 = t52 + t50;
t138 = -qJ(5) - pkin(6);
t223 = (rSges(6,3) - t138) * t136;
t131 = pkin(4) * t141 + pkin(3);
t212 = pkin(3) - t131;
t222 = t137 * t212;
t175 = qJD(5) * t136;
t176 = qJD(4) * t141;
t202 = pkin(4) * qJD(4);
t221 = t134 * pkin(4) * t176 + t90 * rSges(6,1) + t89 * rSges(6,2) - t133 * t175 + t138 * t170 + t172 * t202;
t220 = -rSges(6,1) * t88 - rSges(6,2) * t87 - t134 * t175;
t189 = t136 * t138;
t219 = -rSges(6,1) * t149 + t108 * rSges(6,2) + pkin(4) * t190 + t133 * t189;
t41 = Icges(6,5) * t88 + Icges(6,6) * t87 - Icges(6,3) * t171;
t43 = Icges(5,5) * t88 + Icges(5,6) * t87 - Icges(5,3) * t171;
t218 = (t43 + t41) * t134;
t217 = 2 * m(5);
t216 = 2 * m(6);
t215 = sin(qJ(1)) * pkin(1);
t214 = cos(qJ(1)) * pkin(1);
t213 = pkin(3) * t137;
t211 = pkin(6) + t138;
t179 = qJD(1) * t134;
t180 = qJD(1) * t133;
t181 = pkin(6) * t171 + t180 * t213;
t194 = t131 * t137;
t210 = -rSges(6,3) * t171 + (t189 - t194) * t180 + (t139 * t179 - t144) * pkin(4) + t181 - t220;
t151 = pkin(6) * t136 + t222;
t174 = pkin(4) * t192;
t209 = rSges(6,3) * t170 - (t134 * t151 - t174) * qJD(1) - t221;
t208 = -rSges(6,3) * t193 + t133 * t151 + t219;
t154 = -t111 * rSges(6,1) + rSges(6,2) * t150;
t207 = t174 + (-t136 * t211 - t222) * t134 + rSges(6,3) * t191 - t154;
t205 = t90 * rSges(5,1) + t89 * rSges(5,2);
t203 = -rSges(5,1) * t149 + t108 * rSges(5,2);
t201 = -rSges(4,3) - qJ(3);
t199 = (t211 - rSges(6,3)) * t137 + (rSges(6,1) * t141 - rSges(6,2) * t139 - t212) * t136;
t130 = t134 * qJ(3);
t169 = t139 * t177;
t183 = (-rSges(6,1) * t139 - rSges(6,2) * t141) * t177 - pkin(4) * t169 - qJD(5) * t137;
t173 = rSges(5,3) * t193;
t166 = t130 - t215;
t165 = -pkin(2) - t194;
t164 = -pkin(4) * t139 - qJ(3);
t163 = t199 * t136;
t160 = -rSges(5,1) * t88 - rSges(5,2) * t87;
t158 = pkin(2) * t180 + qJD(1) * t215 - qJD(3) * t133;
t156 = rSges(4,1) * t137 - rSges(4,2) * t136;
t155 = -t111 * rSges(5,1) + rSges(5,2) * t150;
t153 = -pkin(2) - t156;
t152 = -rSges(6,3) * t136 + t165;
t148 = -t213 - pkin(2) + (-rSges(5,3) - pkin(6)) * t136;
t145 = t133 * t164 - t214;
t143 = -t133 * qJ(3) + t134 * t148 - t214;
t92 = t133 * t201 + t134 * t153 - t214;
t129 = qJD(3) * t134;
t119 = (-rSges(5,1) * t139 - rSges(5,2) * t141) * t177;
t107 = -rSges(5,3) * t137 + (rSges(5,1) * t141 - rSges(5,2) * t139) * t136;
t101 = -Icges(5,3) * t137 + (Icges(5,5) * t141 - Icges(5,6) * t139) * t136;
t100 = -Icges(6,3) * t137 + (Icges(6,5) * t141 - Icges(6,6) * t139) * t136;
t91 = rSges(4,3) * t134 + t133 * t153 + t166;
t82 = rSges(5,3) * t191 - t155;
t80 = -t173 + t203;
t64 = qJD(1) * t92 + t129;
t63 = (t133 * t156 + t134 * t201) * qJD(1) + t158;
t62 = t107 * t191 + t137 * t82;
t61 = t107 * t193 - t137 * t80;
t60 = t143 + t155;
t59 = t133 * t148 + t166 + t203;
t56 = -rSges(5,3) * t170 + t205;
t54 = -rSges(5,3) * t171 - t160;
t40 = (t165 - t223) * t134 + t145 + t154;
t39 = t133 * t152 + t166 + t219;
t38 = t101 * t191 - t103 * t150 + t105 * t111;
t37 = t100 * t191 - t102 * t150 + t104 * t111;
t36 = -t101 * t193 + t103 * t108 - t105 * t149;
t35 = -t100 * t193 + t102 * t108 - t104 * t149;
t32 = qJD(1) * t143 + t129 + t205;
t31 = (t173 - t130) * qJD(1) + t158 + t160 + t181;
t30 = -t137 * t56 + (t107 * t179 + t119 * t133) * t136;
t29 = t137 * t54 + (-t107 * t180 + t119 * t134) * t136;
t28 = t134 * t163 + t137 * t207;
t27 = t133 * t163 - t137 * t208;
t26 = -t137 * t68 + (-t139 * t72 + t141 * t76) * t136;
t25 = -t137 * t67 + (-t139 * t71 + t141 * t75) * t136;
t24 = -t137 * t66 + (-t139 * t70 + t141 * t74) * t136;
t23 = -t137 * t65 + (-t139 * t69 + t141 * t73) * t136;
t22 = t129 + (t134 * t152 + t145) * qJD(1) + t221;
t21 = t150 * t202 + (t164 * t134 + (t194 + t223) * t133) * qJD(1) + t158 + t220;
t12 = t103 * t89 + t105 * t90 + t108 * t115 - t149 * t117 + (-t101 * t179 - t113 * t133) * t136;
t11 = t102 * t89 + t104 * t90 + t108 * t114 - t149 * t116 + (-t100 * t179 - t112 * t133) * t136;
t10 = t103 * t87 + t105 * t88 - t150 * t115 + t111 * t117 + (-t101 * t180 + t113 * t134) * t136;
t9 = t102 * t87 + t104 * t88 - t150 * t114 + t111 * t116 + (-t100 * t180 + t112 * t134) * t136;
t8 = t209 * t137 + (t133 * t183 + t179 * t199) * t136;
t7 = t210 * t137 + (t134 * t183 - t180 * t199) * t136;
t6 = (-t133 * t54 - t134 * t56 + (t133 * t80 - t134 * t82) * qJD(1)) * t136;
t5 = -t137 * t43 + (-t139 * t47 + t141 * t51 + (-t139 * t76 - t141 * t72) * qJD(4)) * t136;
t4 = -t137 * t44 + (-t139 * t48 + t141 * t52 + (-t139 * t75 - t141 * t71) * qJD(4)) * t136;
t3 = -t137 * t41 + (-t139 * t45 + t141 * t49 + (-t139 * t74 - t141 * t70) * qJD(4)) * t136;
t2 = -t137 * t42 + (-t139 * t46 + t141 * t50 + (-t139 * t73 - t141 * t69) * qJD(4)) * t136;
t1 = (t209 * t134 - t210 * t133 + (t133 * t208 - t134 * t207) * qJD(1)) * t136;
t13 = [0.2e1 * m(4) * (t63 * t92 + t64 * t91) + (t31 * t60 + t32 * t59) * t217 + (t21 * t40 + t22 * t39) * t216 + t243 * t169 + t240 + (t244 * t176 + t242) * t136; 0; 0; m(4) * (t133 * t64 + t134 * t63 + (-t133 * t92 + t134 * t91) * qJD(1)) + m(5) * (t133 * t32 + t134 * t31 + (-t133 * t60 + t134 * t59) * qJD(1)) + m(6) * (t133 * t22 + t134 * t21 + (-t133 * t40 + t134 * t39) * qJD(1)); 0; 0; -t241 + m(5) * (t29 * t60 + t30 * t59 + t31 * t62 + t32 * t61) + m(6) * (t21 * t28 + t22 * t27 + t39 * t8 + t40 * t7) + ((t5 / 0.2e1 + t3 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t134 + (-t4 / 0.2e1 - t2 / 0.2e1 - t12 / 0.2e1 - t11 / 0.2e1) * t133 + ((-t25 / 0.2e1 - t23 / 0.2e1 - t36 / 0.2e1 - t35 / 0.2e1) * t134 + (-t26 / 0.2e1 - t24 / 0.2e1 - t38 / 0.2e1 - t37 / 0.2e1) * t133) * qJD(1)) * t136; m(5) * t6 + m(6) * t1; m(5) * (t133 * t30 + t134 * t29 + (-t133 * t62 + t134 * t61) * qJD(1)) + m(6) * (t133 * t8 + t134 * t7 + (-t133 * t28 + t134 * t27) * qJD(1)); (t62 * t29 + t61 * t30) * t217 + (t27 * t8 + t28 * t7) * t216 + ((-t133 * t82 - t134 * t80) * t6 * t217 + (-t133 * t207 - t134 * t208) * t1 * t216 + (-t236 * t133 + t237 * t134) * t171 + (t239 * t133 - t238 * t134) * t170 + (t238 * t180 + t239 * t179 + (t228 * t136 * t179 - t233 * t108 + t149 * t231 - t225 * t90 - t226 * t89) * t134 + (t224 * t90 + t227 * t89 + (t234 * t133 - t229 * t179 + t218) * t136 - t230 * t149 + t232 * t108) * t133) * t193 + (t237 * t180 + t236 * t179 + (t225 * t88 + t226 * t87 + (-t228 * t180 + t218) * t136 + t231 * t111 - t233 * t150) * t134 + ((t234 * t134 + t229 * t180) * t136 - t224 * t88 - t227 * t87 - t230 * t111 + t232 * t150) * t133) * t191) * t136 + (t241 + (t12 + t11) * t193 + (-t10 - t9) * t191 + (t38 + t37) * t171 + (t36 + t35) * t170 + ((-t3 - t5) * t134 + (t2 + t4) * t133 + ((t23 + t25) * t134 + (t24 + t26) * t133) * qJD(1)) * t136) * t137; m(6) * (-t133 * t21 + t134 * t22 + (-t133 * t39 - t134 * t40) * qJD(1)) * t136; 0; 0; m(6) * (-t1 * t137 + (-t133 * t7 + t134 * t8 + (-t133 * t27 - t134 * t28) * qJD(1)) * t136); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
