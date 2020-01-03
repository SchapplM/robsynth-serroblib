% Calculate time derivative of joint inertia matrix for
% S5RPPRR7
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:38
% DurationCPUTime: 3.71s
% Computational Cost: add. (9287->419), mult. (12679->629), div. (0->0), fcn. (12043->8), ass. (0->220)
t146 = qJ(1) + pkin(8);
t143 = sin(t146);
t144 = cos(t146);
t151 = cos(qJ(4));
t215 = qJD(4) * t151;
t148 = sin(qJ(4));
t220 = qJD(1) * t148;
t279 = t143 * t215 + t144 * t220;
t260 = rSges(6,3) + pkin(7);
t206 = t260 * t151;
t257 = pkin(4) * t148;
t278 = -t206 + t257;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t239 = Icges(6,4) * t150;
t171 = -Icges(6,2) * t147 + t239;
t108 = Icges(6,6) * t148 + t151 * t171;
t240 = Icges(6,4) * t147;
t174 = Icges(6,1) * t150 - t240;
t109 = Icges(6,5) * t148 + t151 * t174;
t277 = -t108 * t147 + t109 * t150;
t276 = -t108 * t150 - t109 * t147;
t242 = Icges(5,4) * t148;
t172 = Icges(5,2) * t151 + t242;
t94 = Icges(5,6) * t144 + t143 * t172;
t241 = Icges(5,4) * t151;
t175 = Icges(5,1) * t148 + t241;
t96 = Icges(5,5) * t144 + t143 * t175;
t178 = t148 * t96 + t151 * t94;
t275 = t144 * t178;
t190 = rSges(5,1) * t148 + rSges(5,2) * t151;
t162 = t144 * t190;
t266 = -pkin(2) - pkin(6);
t213 = -rSges(5,3) + t266;
t258 = sin(qJ(1)) * pkin(1);
t160 = t143 * t213 - t258;
t216 = qJD(4) * t148;
t199 = t143 * t216;
t219 = qJD(1) * t151;
t203 = t144 * t219;
t209 = rSges(5,1) * t279 + rSges(5,2) * t203;
t221 = qJD(1) * t144;
t223 = qJ(3) * t221 + qJD(3) * t143;
t51 = -rSges(5,2) * t199 + qJD(1) * t160 + t209 + t223;
t254 = rSges(5,2) * t148;
t128 = rSges(5,1) * t151 - t254;
t135 = qJD(3) * t144;
t145 = cos(qJ(1)) * pkin(1);
t217 = qJD(4) * t144;
t52 = t135 + t128 * t217 + (-t145 + t213 * t144 + (-qJ(3) - t190) * t143) * qJD(1);
t274 = t143 * t52 - t144 * t51;
t169 = Icges(5,5) * t148 + Icges(5,6) * t151;
t273 = -Icges(5,3) * t143 + t144 * t169;
t272 = -Icges(5,6) * t143 + t144 * t172;
t271 = -Icges(5,5) * t143 + t144 * t175;
t196 = qJD(5) * t148 + qJD(1);
t270 = t147 * t196 - t150 * t215;
t269 = t147 * t215 + t150 * t196;
t268 = 2 * m(5);
t267 = 2 * m(6);
t141 = t143 ^ 2;
t142 = t144 ^ 2;
t265 = -t143 / 0.2e1;
t263 = t144 / 0.2e1;
t262 = t148 / 0.2e1;
t261 = rSges(4,2) - pkin(2);
t259 = m(5) * t128;
t256 = pkin(4) * t151;
t168 = Icges(6,5) * t150 - Icges(6,6) * t147;
t107 = Icges(6,3) * t148 + t151 * t168;
t214 = qJD(5) * t151;
t83 = (-Icges(6,5) * t147 - Icges(6,6) * t150) * t214 + (Icges(6,3) * t151 - t148 * t168) * qJD(4);
t85 = (-Icges(6,1) * t147 - t239) * t214 + (Icges(6,5) * t151 - t148 * t174) * qJD(4);
t154 = t151 * t150 * t85 + t107 * t215 + t148 * t83 - t216 * t277;
t84 = (-Icges(6,2) * t150 - t240) * t214 + (Icges(6,6) * t151 - t148 * t171) * qJD(4);
t250 = t147 * t84;
t56 = t107 * t148 + t151 * t277;
t255 = ((qJD(5) * t276 - t250) * t151 + t154) * t148 + t56 * t215;
t253 = rSges(5,3) * t143;
t138 = t144 * rSges(5,3);
t249 = t148 * t94;
t248 = t148 * t272;
t247 = t151 * t96;
t246 = t151 * t271;
t230 = t143 * t148;
t132 = pkin(4) * t230;
t229 = t143 * t151;
t227 = t147 * t148;
t103 = -t143 * t227 + t144 * t150;
t226 = t148 * t150;
t104 = t143 * t226 + t144 * t147;
t225 = t104 * rSges(6,1) + t103 * rSges(6,2);
t69 = -rSges(6,3) * t229 + t225;
t245 = pkin(7) * t229 - t132 - t69;
t192 = pkin(7) * t151 - t257;
t105 = t143 * t150 + t144 * t227;
t106 = t143 * t147 - t144 * t226;
t189 = -rSges(6,1) * t106 - rSges(6,2) * t105;
t228 = t144 * t151;
t70 = rSges(6,3) * t228 - t189;
t244 = t192 * t144 + t70;
t188 = rSges(6,1) * t150 - rSges(6,2) * t147;
t86 = (-rSges(6,1) * t147 - rSges(6,2) * t150) * t214 + (rSges(6,3) * t151 - t148 * t188) * qJD(4);
t243 = t192 * qJD(4) + t86;
t92 = Icges(5,3) * t144 + t143 * t169;
t235 = qJD(1) * t92;
t110 = rSges(6,3) * t148 + t151 * t188;
t131 = pkin(7) * t148 + t256;
t224 = t110 + t131;
t222 = qJD(1) * t143;
t218 = qJD(4) * t143;
t195 = qJD(5) + t220;
t167 = t195 * t147;
t61 = -t143 * t269 - t144 * t167;
t166 = t195 * t150;
t62 = -t143 * t270 + t144 * t166;
t212 = t62 * rSges(6,1) + t61 * rSges(6,2) + rSges(6,3) * t199;
t65 = Icges(6,4) * t104 + Icges(6,2) * t103 - Icges(6,6) * t229;
t67 = Icges(6,1) * t104 + Icges(6,4) * t103 - Icges(6,5) * t229;
t182 = t147 * t65 - t150 * t67;
t63 = Icges(6,5) * t104 + Icges(6,6) * t103 - Icges(6,3) * t229;
t29 = t148 * t63 - t151 * t182;
t41 = t103 * t108 + t104 * t109 - t107 * t229;
t211 = t29 / 0.2e1 + t41 / 0.2e1;
t66 = Icges(6,4) * t106 + Icges(6,2) * t105 + Icges(6,6) * t228;
t68 = Icges(6,1) * t106 + Icges(6,4) * t105 + Icges(6,5) * t228;
t181 = t147 * t66 - t150 * t68;
t64 = Icges(6,5) * t106 + Icges(6,6) * t105 + Icges(6,3) * t228;
t30 = t148 * t64 - t151 * t181;
t42 = t105 * t108 + t106 * t109 + t107 * t228;
t210 = -t42 / 0.2e1 - t30 / 0.2e1;
t208 = -pkin(4) * t279 - pkin(7) * t199;
t98 = rSges(5,1) * t230 + rSges(5,2) * t229 + t138;
t207 = t144 * pkin(2) + t143 * qJ(3) + t145;
t205 = t143 * t219;
t116 = t190 * qJD(4);
t198 = (t141 + t142) * t116;
t197 = qJD(1) * t224;
t194 = t144 * pkin(6) + t207;
t59 = -t143 * t167 + t144 * t269;
t60 = t143 * t166 + t144 * t270;
t193 = t60 * rSges(6,1) + t59 * rSges(6,2);
t24 = t103 * t65 + t104 * t67 - t229 * t63;
t25 = t103 * t66 + t104 * t68 - t229 * t64;
t15 = t143 * t25 + t144 * t24;
t187 = t143 * t24 - t144 * t25;
t26 = t105 * t65 + t106 * t67 + t228 * t63;
t27 = t105 * t66 + t106 * t68 + t228 * t64;
t16 = t143 * t27 + t144 * t26;
t186 = t143 * t26 - t144 * t27;
t185 = t143 * t30 + t144 * t29;
t184 = t143 * t29 - t144 * t30;
t183 = t143 * t70 + t144 * t69;
t177 = -t148 * t271 - t151 * t272;
t176 = Icges(5,1) * t151 - t242;
t173 = -Icges(5,2) * t148 + t241;
t170 = Icges(5,5) * t151 - Icges(5,6) * t148;
t165 = t143 * t266 - t258;
t155 = -t144 * t216 - t205;
t32 = Icges(6,5) * t60 + Icges(6,6) * t59 + Icges(6,3) * t155;
t164 = -t151 * t32 + t216 * t64;
t156 = t199 - t203;
t33 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t156;
t163 = -t151 * t33 + t216 * t63;
t161 = t177 * t143;
t159 = qJD(4) * t176;
t158 = qJD(4) * t173;
t153 = rSges(4,3) * t144 + t143 * t261 - t258;
t137 = t144 * qJ(3);
t99 = t253 - t162;
t90 = -rSges(4,2) * t144 + rSges(4,3) * t143 + t207;
t89 = t137 + t153;
t88 = t224 * t144;
t87 = t224 * t143;
t81 = t135 + (-t145 + t261 * t144 + (-rSges(4,3) - qJ(3)) * t143) * qJD(1);
t80 = qJD(1) * t153 + t223;
t78 = t194 + t98;
t77 = t137 + t162 + t160;
t72 = qJD(1) * t273 + t170 * t218;
t71 = -t170 * t217 + t235;
t54 = t110 * t228 - t148 * t70;
t53 = t110 * t229 + t148 * t69;
t50 = -t143 * t206 + t132 + t194 + t225;
t49 = t144 * t278 + t137 + t165 + t189;
t48 = -t143 * t273 - t144 * t177;
t47 = t143 * t92 - t275;
t46 = -t144 * t273 + t161;
t45 = t178 * t143 + t144 * t92;
t44 = t143 * t243 + t144 * t197;
t43 = t143 * t197 - t144 * t243;
t40 = t183 * t151;
t39 = -rSges(6,3) * t203 + t212;
t38 = rSges(6,3) * t155 + t193;
t37 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t156;
t36 = Icges(6,1) * t60 + Icges(6,4) * t59 + Icges(6,5) * t155;
t35 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t156;
t34 = Icges(6,4) * t60 + Icges(6,2) * t59 + Icges(6,6) * t155;
t31 = t143 * t245 + t144 * t244;
t28 = -t143 * t209 + (-t128 * t142 + t141 * t254) * qJD(4) + ((-t98 + t138) * t144 + (t162 - t99 + t253) * t143) * qJD(1);
t23 = t135 + (t148 * t260 + t256) * t217 + (-t145 + t266 * t144 + (-qJ(3) - t278) * t143) * qJD(1) - t193;
t22 = (-t144 * t206 + t165) * qJD(1) - t208 + t212 + t223;
t20 = (-t110 * t218 + t39) * t148 + (qJD(4) * t69 + t110 * t221 + t143 * t86) * t151;
t19 = (-t110 * t217 - t38) * t148 + (-qJD(4) * t70 - t110 * t222 + t144 * t86) * t151;
t18 = t103 * t84 + t104 * t85 + t107 * t156 + t108 * t61 + t109 * t62 - t229 * t83;
t17 = t105 * t84 + t106 * t85 + t107 * t155 + t108 * t59 + t109 * t60 + t228 * t83;
t14 = t183 * t216 + (-t143 * t38 - t144 * t39 + (t143 * t69 - t144 * t70) * qJD(1)) * t151;
t13 = t148 * t42 - t151 * t186;
t12 = t148 * t41 - t151 * t187;
t11 = (-qJD(1) * t244 + t208 - t39) * t143 + (t38 - t131 * t217 + (t132 + t245) * qJD(1)) * t144;
t10 = (qJD(4) * t181 + t32) * t148 + (qJD(4) * t64 - t147 * t34 + t150 * t36 + (-t147 * t68 - t150 * t66) * qJD(5)) * t151;
t9 = (qJD(4) * t182 + t33) * t148 + (qJD(4) * t63 - t147 * t35 + t150 * t37 + (-t147 * t67 - t150 * t65) * qJD(5)) * t151;
t8 = t103 * t34 + t104 * t36 + t143 * t164 - t203 * t64 + t61 * t66 + t62 * t68;
t7 = t103 * t35 + t104 * t37 + t143 * t163 - t203 * t63 + t61 * t65 + t62 * t67;
t6 = t105 * t34 + t106 * t36 - t144 * t164 - t205 * t64 + t59 * t66 + t60 * t68;
t5 = t105 * t35 + t106 * t37 - t144 * t163 - t205 * t63 + t59 * t65 + t60 * t67;
t4 = -qJD(1) * t187 + t143 * t8 + t144 * t7;
t3 = -qJD(1) * t186 + t143 * t6 + t144 * t5;
t2 = (qJD(4) * t187 + t18) * t148 + (-qJD(1) * t15 + qJD(4) * t41 - t143 * t7 + t144 * t8) * t151;
t1 = (qJD(4) * t186 + t17) * t148 + (-qJD(1) * t16 + qJD(4) * t42 - t143 * t5 + t144 * t6) * t151;
t21 = [(t22 * t50 + t23 * t49) * t267 - t148 * t159 - t175 * t215 + t172 * t216 + (t51 * t78 + t52 * t77) * t268 + 0.2e1 * m(4) * (t80 * t90 + t81 * t89) + t154 + t276 * t214 + (-t250 - t158) * t151; 0; 0; m(6) * (t143 * t23 - t144 * t22 + (t143 * t50 + t144 * t49) * qJD(1)) + m(5) * ((t143 * t78 + t144 * t77) * qJD(1) + t274) + m(4) * (t143 * t81 - t144 * t80 + (t143 * t90 + t144 * t89) * qJD(1)); 0; 0; m(6) * (-t22 * t88 + t23 * t87 + t43 * t50 + t44 * t49) + m(5) * (t274 * t128 - (t143 * t77 - t144 * t78) * t116) - (t141 / 0.2e1 + t142 / 0.2e1) * t169 * qJD(4) + ((t249 / 0.2e1 - t247 / 0.2e1 + t78 * t259 - t211) * t143 + (t248 / 0.2e1 - t246 / 0.2e1 + t77 * t259 - t210) * t144) * qJD(1) + (-qJD(4) * t177 - t148 * (qJD(1) * t94 - t173 * t217) + t151 * (qJD(1) * t96 - t176 * t217) + t10 + t17) * t143 / 0.2e1 + (-qJD(4) * t178 - t148 * (qJD(1) * t272 + t143 * t158) + t151 * (qJD(1) * t271 + t143 * t159) + t18 + t9) * t263; m(5) * t28 + m(6) * t11; m(6) * (t44 * t143 - t43 * t144 + (-t143 * t88 + t144 * t87) * qJD(1)) - m(5) * t198; ((-t143 * t98 + t144 * t99) * t28 - t128 * t198) * t268 + t144 * ((t144 * t72 + (t46 + t275) * qJD(1)) * t144 + (-t45 * qJD(1) + (-t215 * t271 + t216 * t272) * t143 + (t71 + (t247 - t249) * qJD(4) + (t177 - t92) * qJD(1)) * t144) * t143) + t143 * ((t143 * t71 + (-t47 + t161) * qJD(1)) * t143 + (t48 * qJD(1) + (-t215 * t96 + t216 * t94 + t235) * t144 + (t72 + (t246 - t248) * qJD(4) + t178 * qJD(1)) * t143) * t144) + (t31 * t11 - t43 * t88 + t44 * t87) * t267 + t144 * t4 + t143 * t3 + (-t143 * t46 - t144 * t45 - t15) * t222 + (t143 * t48 + t144 * t47 + t16) * t221; m(6) * (t19 * t49 + t20 * t50 + t22 * t53 + t23 * t54) + (t143 * t211 + t144 * t210) * t216 + ((t17 / 0.2e1 + t10 / 0.2e1) * t144 + (-t9 / 0.2e1 - t18 / 0.2e1) * t143 + (t143 * t210 - t144 * t211) * qJD(1)) * t151 + t255; m(6) * t14; m(6) * (t143 * t19 - t144 * t20 + (t143 * t53 + t144 * t54) * qJD(1)); m(6) * (-t40 * t11 + t14 * t31 + t19 * t87 - t20 * t88 + t53 * t43 + t54 * t44) + (t2 / 0.2e1 + qJD(1) * t13 / 0.2e1 - t16 * t216 / 0.2e1 + (qJD(1) * t30 + t9) * t262) * t144 + (-qJD(1) * t12 / 0.2e1 + t15 * t216 / 0.2e1 + t1 / 0.2e1 + (-qJD(1) * t29 + t10) * t262) * t143 + (t4 * t265 + t3 * t263 + qJD(4) * t185 / 0.2e1 + (-t144 * t15 / 0.2e1 + t16 * t265) * qJD(1)) * t151; (-t14 * t40 + t19 * t54 + t20 * t53) * t267 + ((t143 * t12 - t144 * t13 + t148 * t184) * qJD(4) + t255) * t148 + (-t143 * t2 + t144 * t1 + t148 * (t10 * t144 - t143 * t9) + (t148 * t56 - t151 * t184) * qJD(4) + (-t144 * t12 - t143 * t13 - t148 * t185) * qJD(1)) * t151;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
