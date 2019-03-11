% Calculate time derivative of joint inertia matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:51
% EndTime: 2019-03-09 01:29:58
% DurationCPUTime: 4.11s
% Computational Cost: add. (9932->453), mult. (13435->675), div. (0->0), fcn. (12619->8), ass. (0->229)
t150 = qJ(1) + pkin(9);
t148 = cos(t150);
t152 = sin(qJ(5));
t220 = qJD(5) * t152;
t208 = t148 * t220;
t147 = sin(t150);
t155 = cos(qJ(5));
t222 = qJD(1) * t155;
t210 = t147 * t222;
t164 = t208 + t210;
t151 = sin(qJ(6));
t154 = cos(qJ(6));
t249 = Icges(7,4) * t154;
t177 = -Icges(7,2) * t151 + t249;
t112 = Icges(7,6) * t152 + t155 * t177;
t250 = Icges(7,4) * t151;
t179 = Icges(7,1) * t154 - t250;
t113 = Icges(7,5) * t152 + t155 * t179;
t286 = -t112 * t151 + t113 * t154;
t285 = -t112 * t154 - t113 * t151;
t277 = -t147 / 0.2e1;
t219 = qJD(5) * t155;
t207 = t148 * t219;
t204 = -qJD(6) * t152 - qJD(1);
t159 = -t151 * t219 + t154 * t204;
t203 = qJD(1) * t152 + qJD(6);
t240 = t147 * t151;
t63 = t148 * t159 + t203 * t240;
t158 = t151 * t204 + t154 * t219;
t239 = t147 * t154;
t64 = t148 * t158 - t203 * t239;
t38 = t64 * rSges(7,1) + t63 * rSges(7,2) + rSges(7,3) * t164;
t284 = pkin(5) * t207 + pkin(8) * t164 + t38;
t149 = cos(qJ(1)) * pkin(1);
t283 = t147 * qJ(3) + t149;
t145 = t147 ^ 2;
t146 = t148 ^ 2;
t282 = t145 + t146;
t143 = t148 * qJ(3);
t196 = rSges(6,1) * t152 + rSges(6,2) * t155;
t265 = -pkin(2) - qJ(4);
t170 = -t196 + t265;
t269 = sin(qJ(1)) * pkin(1);
t272 = -rSges(6,3) - pkin(7);
t157 = t147 * t170 + t148 * t272 - t269;
t61 = t143 + t157;
t211 = t148 * pkin(2) + t283;
t202 = t148 * qJ(4) + t211;
t234 = t148 * t155;
t236 = t148 * t152;
t227 = rSges(6,1) * t236 + rSges(6,2) * t234;
t62 = t147 * t272 + t202 + t227;
t281 = t147 * t62 + t148 * t61;
t123 = rSges(6,1) * t207;
t223 = qJD(1) * t148;
t226 = qJ(3) * t223 + qJD(3) * t147;
t212 = qJD(4) * t148 + t226;
t45 = -rSges(6,2) * t208 + qJD(1) * t157 + t123 + t212;
t263 = rSges(6,2) * t152;
t131 = rSges(6,1) * t155 - t263;
t140 = qJD(3) * t148;
t224 = qJD(1) * t147;
t225 = pkin(7) * t224 + t140;
t46 = (-t131 * qJD(5) - qJD(4)) * t147 + (-t149 + (rSges(6,3) - qJ(3)) * t147 + t170 * t148) * qJD(1) + t225;
t280 = t147 * t45 + t148 * t46;
t176 = Icges(6,5) * t152 + Icges(6,6) * t155;
t96 = Icges(6,3) * t148 + t147 * t176;
t252 = Icges(6,4) * t152;
t178 = Icges(6,2) * t155 + t252;
t98 = Icges(6,6) * t148 + t147 * t178;
t251 = Icges(6,4) * t155;
t180 = Icges(6,1) * t152 + t251;
t100 = Icges(6,5) * t148 + t147 * t180;
t279 = 2 * m(6);
t278 = 2 * m(7);
t276 = -t148 / 0.2e1;
t274 = t152 / 0.2e1;
t273 = rSges(4,2) - pkin(2);
t271 = rSges(7,3) + pkin(8);
t270 = m(6) * t131;
t268 = pkin(5) * t152;
t267 = pkin(5) * t155;
t266 = -qJD(1) / 0.2e1;
t175 = Icges(7,5) * t154 - Icges(7,6) * t151;
t111 = Icges(7,3) * t152 + t155 * t175;
t218 = qJD(6) * t155;
t85 = (-Icges(7,5) * t151 - Icges(7,6) * t154) * t218 + (Icges(7,3) * t155 - t152 * t175) * qJD(5);
t87 = (-Icges(7,1) * t151 - t249) * t218 + (Icges(7,5) * t155 - t152 * t179) * qJD(5);
t162 = t155 * t154 * t87 + t111 * t219 + t152 * t85 - t220 * t286;
t86 = (-Icges(7,2) * t154 - t250) * t218 + (Icges(7,6) * t155 - t152 * t177) * qJD(5);
t258 = t151 * t86;
t56 = t111 * t152 + t155 * t286;
t264 = ((qJD(6) * t285 - t258) * t155 + t162) * t152 + t56 * t219;
t257 = t152 * t98;
t99 = -Icges(6,6) * t147 + t148 * t178;
t256 = t152 * t99;
t198 = -pkin(8) * t155 + t268;
t233 = t151 * t152;
t235 = t148 * t154;
t107 = -t147 * t233 + t235;
t232 = t152 * t154;
t237 = t148 * t151;
t108 = t147 * t232 + t237;
t195 = -rSges(7,1) * t108 - rSges(7,2) * t107;
t238 = t147 * t155;
t73 = -rSges(7,3) * t238 - t195;
t255 = -t198 * t147 - t73;
t135 = pkin(5) * t236;
t216 = pkin(8) * t234;
t109 = -t148 * t233 - t239;
t110 = t148 * t232 - t240;
t229 = t110 * rSges(7,1) + t109 * rSges(7,2);
t74 = -rSges(7,3) * t234 + t229;
t254 = -t135 + t216 - t74;
t194 = rSges(7,1) * t154 - rSges(7,2) * t151;
t88 = (-rSges(7,1) * t151 - rSges(7,2) * t154) * t218 + (rSges(7,3) * t155 - t152 * t194) * qJD(5);
t253 = -t198 * qJD(5) + t88;
t97 = -Icges(6,3) * t147 + t148 * t176;
t245 = qJD(1) * t97;
t231 = t155 * t100;
t101 = -Icges(6,5) * t147 + t148 * t180;
t230 = t155 * t101;
t114 = rSges(7,3) * t152 + t155 * t194;
t134 = pkin(8) * t152 + t267;
t228 = t114 + t134;
t221 = qJD(5) * t147;
t217 = -rSges(5,3) + t265;
t69 = Icges(7,4) * t108 + Icges(7,2) * t107 - Icges(7,6) * t238;
t71 = Icges(7,1) * t108 + Icges(7,4) * t107 - Icges(7,5) * t238;
t186 = t151 * t69 - t154 * t71;
t67 = Icges(7,5) * t108 + Icges(7,6) * t107 - Icges(7,3) * t238;
t29 = t152 * t67 - t155 * t186;
t41 = t107 * t112 + t108 * t113 - t111 * t238;
t215 = t29 / 0.2e1 + t41 / 0.2e1;
t70 = Icges(7,4) * t110 + Icges(7,2) * t109 - Icges(7,6) * t234;
t72 = Icges(7,1) * t110 + Icges(7,4) * t109 - Icges(7,5) * t234;
t185 = t151 * t70 - t154 * t72;
t68 = Icges(7,5) * t110 + Icges(7,6) * t109 - Icges(7,3) * t234;
t30 = t152 * t68 - t155 * t185;
t42 = t109 * t112 + t110 * t113 - t111 * t234;
t214 = t30 / 0.2e1 + t42 / 0.2e1;
t209 = t148 * t222;
t206 = t220 / 0.2e1;
t120 = t196 * qJD(5);
t205 = t120 * t282;
t90 = t228 * t148;
t201 = t265 - t268;
t65 = t147 * t159 - t203 * t237;
t66 = t147 * t158 + t203 * t235;
t200 = t66 * rSges(7,1) + t65 * rSges(7,2);
t199 = -pkin(7) * t148 - t269;
t24 = t107 * t69 + t108 * t71 - t238 * t67;
t25 = t107 * t70 + t108 * t72 - t238 * t68;
t193 = t147 * t25 - t148 * t24;
t192 = t147 * t24 + t148 * t25;
t26 = t109 * t69 + t110 * t71 - t234 * t67;
t27 = t109 * t70 + t110 * t72 - t234 * t68;
t191 = t147 * t27 - t148 * t26;
t190 = t147 * t26 + t148 * t27;
t189 = t30 * t147 - t29 * t148;
t188 = t147 * t29 + t148 * t30;
t187 = t147 * t74 - t148 * t73;
t182 = t100 * t152 + t155 * t98;
t181 = t101 * t152 + t155 * t99;
t32 = Icges(7,5) * t64 + Icges(7,6) * t63 + Icges(7,3) * t164;
t174 = -t155 * t32 + t220 * t68;
t165 = t147 * t220 - t209;
t33 = Icges(7,5) * t66 + Icges(7,6) * t65 + Icges(7,3) * t165;
t173 = -t155 * t33 + t220 * t67;
t172 = t182 * t148;
t171 = t181 * t147;
t169 = qJD(5) * (Icges(6,1) * t155 - t252);
t168 = qJD(5) * (-Icges(6,2) * t152 + t251);
t167 = qJD(5) * (Icges(6,5) * t155 - Icges(6,6) * t152);
t163 = t155 * t271 + t201;
t161 = rSges(4,3) * t148 + t147 * t273 - t269;
t160 = rSges(5,2) * t148 + t147 * t217 - t269;
t103 = -t147 * rSges(6,3) + t227;
t102 = t148 * rSges(6,3) + t147 * t196;
t94 = -rSges(4,2) * t148 + rSges(4,3) * t147 + t211;
t93 = t143 + t161;
t92 = rSges(5,2) * t147 + rSges(5,3) * t148 + t202;
t91 = t143 + t160;
t89 = t228 * t147;
t83 = t140 + (-t149 + t273 * t148 + (-rSges(4,3) - qJ(3)) * t147) * qJD(1);
t82 = qJD(1) * t161 + t226;
t76 = t147 * t167 + t245;
t75 = -qJD(1) * t96 + t148 * t167;
t58 = -qJD(4) * t147 + t140 + (-t149 + (-rSges(5,2) - qJ(3)) * t147 + t217 * t148) * qJD(1);
t57 = qJD(1) * t160 + t212;
t54 = t114 * t234 + t152 * t74;
t53 = -t114 * t238 - t152 * t73;
t52 = -t147 * t97 + t181 * t148;
t51 = -t147 * t96 + t172;
t50 = t148 * t97 + t171;
t49 = t147 * t182 + t148 * t96;
t48 = -pkin(7) * t147 - t234 * t271 + t135 + t202 + t229;
t47 = t147 * t163 + t143 + t195 + t199;
t44 = qJD(1) * t90 + t147 * t253;
t43 = t148 * t253 - t224 * t228;
t40 = t187 * t155;
t39 = rSges(7,3) * t165 + t200;
t37 = Icges(7,1) * t66 + Icges(7,4) * t65 + Icges(7,5) * t165;
t36 = Icges(7,1) * t64 + Icges(7,4) * t63 + Icges(7,5) * t164;
t35 = Icges(7,4) * t66 + Icges(7,2) * t65 + Icges(7,6) * t165;
t34 = Icges(7,4) * t64 + Icges(7,2) * t63 + Icges(7,6) * t164;
t31 = t147 * t255 + t148 * t254;
t28 = -t148 * t123 + (-t131 * t145 + t146 * t263) * qJD(5) + (t282 * rSges(6,3) - t148 * t102 + t147 * t103) * qJD(1);
t23 = (-qJD(4) + (-t152 * t271 - t267) * qJD(5)) * t147 + (t148 * t163 - t283) * qJD(1) - t200 + t225;
t22 = (t147 * t201 + t199) * qJD(1) + t212 + t284;
t20 = (t114 * t221 - t39) * t152 + (-qJD(5) * t73 - t114 * t223 - t147 * t88) * t155;
t19 = (-qJD(5) * t114 * t148 + t38) * t152 + (qJD(5) * t74 - t114 * t224 + t148 * t88) * t155;
t18 = t107 * t86 + t108 * t87 + t111 * t165 + t112 * t65 + t113 * t66 - t238 * t85;
t17 = t109 * t86 + t110 * t87 + t111 * t164 + t112 * t63 + t113 * t64 - t234 * t85;
t14 = -t187 * t220 + (t147 * t38 - t148 * t39 + (t147 * t73 + t148 * t74) * qJD(1)) * t155;
t13 = t42 * t152 - t155 * t190;
t12 = t41 * t152 - t155 * t192;
t11 = (qJD(1) * t255 - t284) * t148 + (-t39 - t134 * t221 + (t216 - t254) * qJD(1)) * t147;
t10 = (qJD(5) * t185 + t32) * t152 + (qJD(5) * t68 - t151 * t34 + t154 * t36 + (-t151 * t72 - t154 * t70) * qJD(6)) * t155;
t9 = (qJD(5) * t186 + t33) * t152 + (qJD(5) * t67 - t151 * t35 + t154 * t37 + (-t151 * t71 - t154 * t69) * qJD(6)) * t155;
t8 = t107 * t34 + t108 * t36 + t147 * t174 - t209 * t68 + t65 * t70 + t66 * t72;
t7 = t107 * t35 + t108 * t37 + t147 * t173 - t209 * t67 + t65 * t69 + t66 * t71;
t6 = t109 * t34 + t110 * t36 + t148 * t174 + t210 * t68 + t63 * t70 + t64 * t72;
t5 = t109 * t35 + t110 * t37 + t148 * t173 + t210 * t67 + t63 * t69 + t64 * t71;
t4 = -qJD(1) * t192 - t147 * t8 + t148 * t7;
t3 = -qJD(1) * t190 - t147 * t6 + t148 * t5;
t2 = (qJD(5) * t192 + t18) * t152 + (qJD(1) * t193 + qJD(5) * t41 - t147 * t7 - t148 * t8) * t155;
t1 = (qJD(5) * t190 + t17) * t152 + (qJD(1) * t191 + qJD(5) * t42 - t147 * t5 - t148 * t6) * t155;
t15 = [(t22 * t48 + t23 * t47) * t278 - t152 * t169 - t180 * t219 + t178 * t220 + (t45 * t62 + t46 * t61) * t279 + 0.2e1 * m(4) * (t82 * t94 + t83 * t93) + 0.2e1 * m(5) * (t57 * t92 + t58 * t91) + t162 + t285 * t218 + (-t258 - t168) * t155; 0; 0; m(7) * (t147 * t23 - t148 * t22 + (t147 * t48 + t148 * t47) * qJD(1)) + m(6) * (t281 * qJD(1) + t147 * t46 - t148 * t45) + m(4) * (t147 * t83 - t148 * t82 + (t147 * t94 + t148 * t93) * qJD(1)) + m(5) * (t147 * t58 - t148 * t57 + (t147 * t92 + t148 * t91) * qJD(1)); 0; 0; m(7) * (t147 * t22 + t148 * t23 + (-t147 * t47 + t148 * t48) * qJD(1)) + m(6) * ((-t147 * t61 + t148 * t62) * qJD(1) + t280) + m(5) * (t147 * t57 + t148 * t58 + (-t147 * t91 + t148 * t92) * qJD(1)); 0; 0; 0; m(7) * (t22 * t89 + t23 * t90 + t43 * t47 + t44 * t48) + m(6) * (-t281 * t120 + t280 * t131) - (t146 / 0.2e1 + t145 / 0.2e1) * t176 * qJD(5) + ((t62 * t270 - t230 / 0.2e1 + t256 / 0.2e1 - t214) * t148 + (-t61 * t270 - t231 / 0.2e1 + t257 / 0.2e1 - t215) * t147) * qJD(1) + (-qJD(5) * t181 - t152 * (-t98 * qJD(1) + t148 * t168) + t155 * (-t100 * qJD(1) + t148 * t169) + t10 + t17) * t277 + (-qJD(5) * t182 - t152 * (qJD(1) * t99 + t147 * t168) + t155 * (qJD(1) * t101 + t147 * t169) + t18 + t9) * t148 / 0.2e1; m(6) * t28 + m(7) * t11; m(7) * (t147 * t43 - t148 * t44 + (t147 * t89 + t148 * t90) * qJD(1)); m(7) * (t44 * t147 + t43 * t148 + (-t147 * t90 + t148 * t89) * qJD(1)) - m(6) * t205; -t147 * ((t147 * t75 + (-t51 + t171) * qJD(1)) * t147 + (-t52 * qJD(1) + (t100 * t219 - t220 * t98) * t148 + (-t76 + (-t230 + t256) * qJD(5) + (-t182 + t97) * qJD(1)) * t147) * t148) + t148 * ((t148 * t76 + (-t50 + t172) * qJD(1)) * t148 + (-t49 * qJD(1) + (-t101 * t219 + t220 * t99 + t245) * t147 + (-t75 + (t231 - t257) * qJD(5) - t181 * qJD(1)) * t148) * t147) + ((-t102 * t147 - t103 * t148) * t28 - t131 * t205) * t279 - t147 * t3 + t148 * t4 + (t11 * t31 + t43 * t90 + t44 * t89) * t278 + (t147 * t50 - t49 * t148 + t193) * t224 + (t147 * t52 - t148 * t51 + t191) * t223; m(7) * (t19 * t48 + t20 * t47 + t22 * t54 + t23 * t53) + (t147 * t215 + t148 * t214) * t220 + ((-t10 / 0.2e1 - t17 / 0.2e1) * t148 + (-t9 / 0.2e1 - t18 / 0.2e1) * t147 + (t147 * t214 - t148 * t215) * qJD(1)) * t155 + t264; m(7) * t14; m(7) * (t147 * t20 - t148 * t19 + (t147 * t54 + t148 * t53) * qJD(1)); m(7) * (t147 * t19 + t148 * t20 + (-t147 * t53 + t148 * t54) * qJD(1)); m(7) * (t11 * t40 + t14 * t31 + t19 * t89 + t20 * t90 + t43 * t53 + t44 * t54) + (t13 * t266 - t191 * t206 + t2 / 0.2e1 + (-qJD(1) * t30 + t9) * t274) * t148 + (-t1 / 0.2e1 + t12 * t266 - t193 * t206 + (-qJD(1) * t29 - t10) * t274) * t147 + (t3 * t276 + t4 * t277 - qJD(5) * t189 / 0.2e1 + (t191 * t277 - t193 * t276) * qJD(1)) * t155; (t14 * t40 + t19 * t54 + t20 * t53) * t278 + ((t147 * t12 + t148 * t13 + t152 * t188) * qJD(5) + t264) * t152 + (-t148 * t1 - t147 * t2 + t152 * (-t10 * t148 - t147 * t9) + (t56 * t152 - t155 * t188) * qJD(5) + (-t148 * t12 + t147 * t13 + t152 * t189) * qJD(1)) * t155;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
