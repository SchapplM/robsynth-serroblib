% Calculate time derivative of joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:43
% DurationCPUTime: 4.32s
% Computational Cost: add. (5166->376), mult. (7182->539), div. (0->0), fcn. (6703->6), ass. (0->179)
t138 = pkin(7) + qJ(2);
t136 = sin(t138);
t140 = cos(pkin(8));
t143 = cos(qJ(4));
t182 = t140 * t143;
t137 = cos(t138);
t142 = sin(qJ(4));
t185 = t137 * t142;
t109 = t136 * t182 - t185;
t183 = t140 * t142;
t151 = t136 * t183 + t137 * t143;
t139 = sin(pkin(8));
t188 = t136 * t139;
t63 = Icges(6,5) * t109 - Icges(6,6) * t151 + Icges(6,3) * t188;
t65 = Icges(5,5) * t109 - Icges(5,6) * t151 + Icges(5,3) * t188;
t219 = t63 + t65;
t168 = t137 * t183;
t110 = t136 * t143 - t168;
t187 = t136 * t142;
t111 = t137 * t182 + t187;
t186 = t137 * t139;
t64 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t186;
t66 = Icges(5,5) * t111 + Icges(5,6) * t110 + Icges(5,3) * t186;
t218 = t64 + t66;
t68 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t186;
t70 = Icges(5,4) * t111 + Icges(5,2) * t110 + Icges(5,6) * t186;
t217 = t68 + t70;
t67 = Icges(6,4) * t109 - Icges(6,2) * t151 + Icges(6,6) * t188;
t69 = Icges(5,4) * t109 - Icges(5,2) * t151 + Icges(5,6) * t188;
t216 = t69 + t67;
t72 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t186;
t74 = Icges(5,1) * t111 + Icges(5,4) * t110 + Icges(5,5) * t186;
t215 = t72 + t74;
t71 = Icges(6,1) * t109 - Icges(6,4) * t151 + Icges(6,5) * t188;
t73 = Icges(5,1) * t109 - Icges(5,4) * t151 + Icges(5,5) * t188;
t214 = t73 + t71;
t190 = Icges(6,4) * t143;
t102 = -Icges(6,6) * t140 + (-Icges(6,2) * t142 + t190) * t139;
t192 = Icges(5,4) * t143;
t103 = -Icges(5,6) * t140 + (-Icges(5,2) * t142 + t192) * t139;
t234 = -t102 - t103;
t191 = Icges(6,4) * t142;
t104 = -Icges(6,5) * t140 + (Icges(6,1) * t143 - t191) * t139;
t193 = Icges(5,4) * t142;
t105 = -Icges(5,5) * t140 + (Icges(5,1) * t143 - t193) * t139;
t233 = -t104 - t105;
t172 = qJD(4) * t139;
t114 = (-Icges(6,2) * t143 - t191) * t172;
t115 = (-Icges(5,2) * t143 - t193) * t172;
t232 = (-t114 - t115) * t142;
t112 = (-Icges(6,5) * t142 - Icges(6,6) * t143) * t172;
t113 = (-Icges(5,5) * t142 - Icges(5,6) * t143) * t172;
t116 = (-Icges(6,1) * t142 - t190) * t172;
t117 = (-Icges(5,1) * t142 - t192) * t172;
t230 = (t116 + t117) * t139 * t143 + (-t112 - t113) * t140;
t231 = (t230 + ((t233 * t142 + t234 * t143) * qJD(4) + t232) * t139) * t140;
t229 = t214 * t109 - t216 * t151 + t219 * t188;
t228 = -t215 * t109 + t217 * t151 - t218 * t188;
t227 = -t217 * t110 - t215 * t111 - t218 * t186;
t226 = t216 * t110 + t214 * t111 + t219 * t186;
t173 = qJD(2) * t139;
t166 = t137 * t173;
t89 = qJD(2) * t110 - qJD(4) * t109;
t147 = t151 * qJD(4);
t90 = qJD(2) * t111 - t147;
t40 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t166;
t42 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t166;
t224 = t42 + t40;
t167 = t136 * t173;
t87 = qJD(2) * t151 - qJD(4) * t111;
t88 = -qJD(2) * t109 + qJD(4) * t110;
t43 = Icges(6,4) * t88 + Icges(6,2) * t87 - Icges(6,6) * t167;
t45 = Icges(5,4) * t88 + Icges(5,2) * t87 - Icges(5,6) * t167;
t223 = t43 + t45;
t44 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t166;
t46 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t166;
t222 = t46 + t44;
t47 = Icges(6,1) * t88 + Icges(6,4) * t87 - Icges(6,5) * t167;
t49 = Icges(5,1) * t88 + Icges(5,4) * t87 - Icges(5,5) * t167;
t221 = t47 + t49;
t48 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t166;
t50 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t166;
t220 = t50 + t48;
t128 = pkin(4) * t187;
t135 = pkin(4) * t143 + pkin(3);
t189 = t135 * t140;
t213 = -t111 * rSges(6,1) - t110 * rSges(6,2) - rSges(6,3) * t186 - t137 * t189 - t128;
t203 = pkin(3) - t135;
t212 = t140 * t203;
t141 = -qJ(5) - pkin(6);
t169 = pkin(4) * t185;
t170 = qJD(5) * t139;
t171 = qJD(4) * t143;
t195 = pkin(4) * qJD(4);
t211 = -t136 * pkin(4) * t171 - t88 * rSges(6,1) - t87 * rSges(6,2) - qJD(2) * t169 - t137 * t170 - t141 * t167 + t168 * t195;
t210 = -rSges(6,1) * t109 + rSges(6,2) * t151 + t169;
t209 = -t90 * rSges(6,1) - t89 * rSges(6,2) - t136 * t170;
t39 = Icges(6,5) * t88 + Icges(6,6) * t87 - Icges(6,3) * t167;
t41 = Icges(5,5) * t88 + Icges(5,6) * t87 - Icges(5,3) * t167;
t208 = (t41 + t39) * t137;
t207 = 2 * m(5);
t206 = 2 * m(6);
t205 = pkin(3) * t140;
t204 = pkin(6) * t139;
t202 = pkin(6) + t141;
t175 = qJD(2) * t136;
t201 = rSges(6,3) * t167 - (t204 + t212) * t175 + t211;
t162 = t202 * t139;
t148 = -t162 - t212;
t200 = rSges(6,3) * t166 - pkin(4) * t147 + (t137 * t148 + t128) * qJD(2) - t209;
t199 = rSges(6,3) * t188 + t136 * t148 - t210;
t198 = -(-t162 - t205) * t137 + t213;
t196 = t88 * rSges(5,1) + t87 * rSges(5,2);
t194 = (t202 - rSges(6,3)) * t140 + (rSges(6,1) * t143 - rSges(6,2) * t142 - t203) * t139;
t132 = t136 * qJ(3);
t165 = t142 * t172;
t179 = (-rSges(6,1) * t142 - rSges(6,2) * t143) * t172 - pkin(4) * t165 - qJD(5) * t140;
t174 = qJD(2) * t137;
t177 = qJ(3) * t174 + qJD(3) * t136;
t176 = t137 * pkin(2) + t132;
t80 = t111 * rSges(5,1) + t110 * rSges(5,2) + rSges(5,3) * t186;
t161 = -pkin(2) - t189;
t158 = qJD(2) * t194;
t157 = -rSges(5,1) * t90 - rSges(5,2) * t89;
t155 = rSges(4,1) * t140 - rSges(4,2) * t139;
t154 = -rSges(5,1) * t109 + rSges(5,2) * t151;
t152 = -pkin(2) - t155;
t150 = -t205 - pkin(2) + (-rSges(5,3) - pkin(6)) * t139;
t149 = (-rSges(6,3) + t141) * t139 + t161;
t146 = t150 * t136;
t144 = rSges(4,3) * t137 + t136 * t152;
t133 = t137 * qJ(3);
t131 = qJD(3) * t137;
t119 = (-rSges(5,1) * t142 - rSges(5,2) * t143) * t172;
t107 = -t140 * rSges(5,3) + (rSges(5,1) * t143 - rSges(5,2) * t142) * t139;
t101 = -Icges(5,3) * t140 + (Icges(5,5) * t143 - Icges(5,6) * t142) * t139;
t100 = -Icges(6,3) * t140 + (Icges(6,5) * t143 - Icges(6,6) * t142) * t139;
t92 = rSges(4,3) * t136 + t137 * t155 + t176;
t91 = t133 + t144;
t82 = t131 + ((-rSges(4,3) - qJ(3)) * t136 + t152 * t137) * qJD(2);
t81 = qJD(2) * t144 + t177;
t78 = rSges(5,3) * t188 - t154;
t62 = (t204 + t205) * t137 + t80 + t176;
t61 = t133 + t146 + t154;
t60 = -t107 * t186 - t140 * t80;
t59 = t107 * t188 + t140 * t78;
t56 = -t141 * t186 + t176 - t213;
t55 = t136 * t149 + t133 + t210;
t54 = rSges(5,3) * t166 - t157;
t52 = -rSges(5,3) * t167 + t196;
t38 = t101 * t186 + t103 * t110 + t105 * t111;
t37 = t100 * t186 + t102 * t110 + t104 * t111;
t36 = t101 * t188 - t103 * t151 + t105 * t109;
t35 = t100 * t188 - t102 * t151 + t104 * t109;
t32 = t131 + (t137 * t150 - t132) * qJD(2) + t157;
t31 = qJD(2) * t146 + t177 + t196;
t30 = t140 * t54 + (t107 * t174 + t119 * t136) * t139;
t29 = -t140 * t52 + (t107 * t175 - t119 * t137) * t139;
t28 = t140 * t198 - t186 * t194;
t27 = t140 * t199 + t188 * t194;
t26 = -t140 * t66 + (-t142 * t70 + t143 * t74) * t139;
t25 = -t140 * t65 + (-t142 * t69 + t143 * t73) * t139;
t24 = -t140 * t64 + (-t142 * t68 + t143 * t72) * t139;
t23 = -t140 * t63 + (-t142 * t67 + t143 * t71) * t139;
t22 = t131 + t151 * t195 + ((-pkin(4) * t142 - qJ(3)) * t136 + t149 * t137) * qJD(2) + t209;
t21 = (-rSges(6,3) * t139 + t161) * t175 + t177 - t211;
t12 = t103 * t89 + t105 * t90 - t151 * t115 + t109 * t117 + (t101 * t174 + t113 * t136) * t139;
t11 = t102 * t89 + t104 * t90 - t151 * t114 + t109 * t116 + (t100 * t174 + t112 * t136) * t139;
t10 = t103 * t87 + t105 * t88 + t110 * t115 + t111 * t117 + (-t101 * t175 + t113 * t137) * t139;
t9 = t102 * t87 + t104 * t88 + t110 * t114 + t111 * t116 + (-t100 * t175 + t112 * t137) * t139;
t8 = t200 * t140 + (t136 * t179 + t137 * t158) * t139;
t7 = t201 * t140 + (t136 * t158 - t137 * t179) * t139;
t6 = (-t136 * t52 + t137 * t54 + (-t136 * t78 - t137 * t80) * qJD(2)) * t139;
t5 = -t140 * t41 + (-t142 * t45 + t143 * t49 + (-t142 * t74 - t143 * t70) * qJD(4)) * t139;
t4 = -t140 * t42 + (-t142 * t46 + t143 * t50 + (-t142 * t73 - t143 * t69) * qJD(4)) * t139;
t3 = -t140 * t39 + (-t142 * t43 + t143 * t47 + (-t142 * t72 - t143 * t68) * qJD(4)) * t139;
t2 = -t140 * t40 + (-t142 * t44 + t143 * t48 + (-t142 * t71 - t143 * t67) * qJD(4)) * t139;
t1 = (t200 * t137 + t201 * t136 + (-t136 * t199 + t137 * t198) * qJD(2)) * t139;
t13 = [0; 0; 0.2e1 * m(4) * (t81 * t92 + t82 * t91) + (t31 * t62 + t32 * t61) * t207 + (t21 * t56 + t22 * t55) * t206 + t233 * t165 + (t234 * t171 + t232) * t139 + t230; 0; m(4) * (t136 * t82 - t137 * t81 + (t136 * t92 + t137 * t91) * qJD(2)) + m(5) * (t136 * t32 - t137 * t31 + (t136 * t62 + t137 * t61) * qJD(2)) + m(6) * (t136 * t22 - t137 * t21 + (t136 * t56 + t137 * t55) * qJD(2)); 0; m(5) * t6 + m(6) * t1; -t231 + m(5) * (t29 * t62 + t30 * t61 + t31 * t60 + t32 * t59) + m(6) * (t21 * t28 + t22 * t27 + t55 * t8 + t56 * t7) + ((t5 / 0.2e1 + t3 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t137 + (t4 / 0.2e1 + t2 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t136 + ((t25 / 0.2e1 + t23 / 0.2e1 + t36 / 0.2e1 + t35 / 0.2e1) * t137 + (-t26 / 0.2e1 - t24 / 0.2e1 - t38 / 0.2e1 - t37 / 0.2e1) * t136) * qJD(2)) * t139; m(5) * (t136 * t30 - t137 * t29 + (t136 * t60 + t137 * t59) * qJD(2)) + m(6) * (t136 * t8 - t137 * t7 + (t136 * t28 + t137 * t27) * qJD(2)); (t60 * t29 + t59 * t30) * t207 + (t27 * t8 + t28 * t7) * t206 + ((-t136 * t80 + t137 * t78) * t6 * t207 + (t136 * t198 + t137 * t199) * t1 * t206 + (-t136 * t226 + t137 * t227) * t167 + (t136 * t229 - t228 * t137) * t166 + (t228 * t175 + t229 * t174 + (t139 * t174 * t218 + t109 * t221 - t151 * t223 + t215 * t90 + t217 * t89) * t137 + (t214 * t90 + t216 * t89 + (t136 * t224 + t174 * t219 + t208) * t139 + t220 * t109 - t222 * t151) * t136) * t188 + (t227 * t175 + t226 * t174 + (t215 * t88 + t217 * t87 + (-t175 * t218 + t208) * t139 + t221 * t111 + t223 * t110) * t137 + ((t137 * t224 - t175 * t219) * t139 + t214 * t88 + t216 * t87 + t220 * t111 + t222 * t110) * t136) * t186) * t139 + (t231 + (-t12 - t11) * t188 + (-t10 - t9) * t186 + (t38 + t37) * t167 + (-t36 - t35) * t166 + ((-t3 - t5) * t137 + (-t2 - t4) * t136 + ((-t23 - t25) * t137 + (t24 + t26) * t136) * qJD(2)) * t139) * t140; 0; m(6) * (t136 * t21 + t137 * t22 + (-t136 * t55 + t137 * t56) * qJD(2)) * t139; 0; m(6) * (-t1 * t140 + (t136 * t7 + t137 * t8 + (-t136 * t27 + t137 * t28) * qJD(2)) * t139); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
