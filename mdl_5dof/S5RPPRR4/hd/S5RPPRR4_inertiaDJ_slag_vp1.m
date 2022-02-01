% Calculate time derivative of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:28
% DurationCPUTime: 5.24s
% Computational Cost: add. (12099->497), mult. (12764->738), div. (0->0), fcn. (12077->10), ass. (0->225)
t197 = cos(qJ(1));
t191 = sin(pkin(9));
t192 = sin(pkin(8));
t193 = cos(pkin(9));
t194 = cos(pkin(8));
t273 = pkin(1) - (-rSges(4,3) - qJ(3)) * t192 + (rSges(4,1) * t193 - rSges(4,2) * t191 + pkin(2)) * t194;
t277 = t273 * t197;
t189 = pkin(9) + qJ(4);
t178 = sin(t189);
t268 = pkin(3) * t191;
t165 = pkin(4) * t178 + t268;
t196 = sin(qJ(1));
t276 = (qJ(2) + t165) * t196;
t180 = qJ(5) + t189;
t174 = sin(t180);
t175 = cos(t180);
t247 = t196 * t175;
t249 = t194 * t197;
t152 = -t174 * t249 + t247;
t248 = t196 * t174;
t153 = t175 * t249 + t248;
t250 = t192 * t197;
t100 = t153 * rSges(6,1) + t152 * rSges(6,2) + rSges(6,3) * t250;
t176 = t193 * pkin(3) + pkin(2);
t179 = cos(t189);
t164 = pkin(4) * t179 + t176;
t183 = t196 * qJ(2);
t239 = t197 * pkin(1) + t183;
t275 = -t164 * t249 - t196 * t165 - t100 - t239;
t262 = Icges(5,4) * t179;
t141 = -Icges(5,6) * t194 + (-Icges(5,2) * t178 + t262) * t192;
t263 = Icges(5,4) * t178;
t142 = -Icges(5,5) * t194 + (Icges(5,1) * t179 - t263) * t192;
t274 = -t141 * t179 - t142 * t178;
t272 = 2 * m(5);
t271 = 2 * m(6);
t270 = m(5) / 0.2e1;
t269 = m(6) / 0.2e1;
t195 = qJ(3) + pkin(6);
t237 = qJD(1) * t197;
t177 = qJ(2) * t237;
t187 = -pkin(7) - t195;
t254 = t178 * t197;
t266 = pkin(4) * qJD(4);
t223 = t254 * t266;
t238 = qJD(1) * t196;
t232 = t192 * t238;
t234 = t196 * t266;
t201 = t165 * t237 + t179 * t234 + t187 * t232 - t194 * t223 + t177;
t161 = t176 * t194 + t192 * t195 + pkin(1);
t228 = -t164 * t194 - pkin(1);
t219 = t161 + t228;
t173 = qJ(2) + t268;
t256 = t173 * t197;
t190 = qJD(4) + qJD(5);
t209 = t175 * t197 + t194 * t248;
t103 = qJD(1) * t209 - t153 * t190;
t151 = -t174 * t197 + t194 * t247;
t104 = -qJD(1) * t151 + t152 * t190;
t243 = t104 * rSges(6,1) + t103 * rSges(6,2);
t65 = -rSges(6,3) * t232 + t243;
t267 = -t65 - (t196 * t219 - t256) * qJD(1) - t201;
t252 = t190 * t192;
t261 = Icges(6,4) * t174;
t131 = (-Icges(6,2) * t175 - t261) * t252;
t138 = -Icges(6,5) * t194 + (Icges(6,1) * t175 - t261) * t192;
t130 = (-Icges(6,5) * t174 - Icges(6,6) * t175) * t252;
t260 = Icges(6,4) * t175;
t132 = (-Icges(6,1) * t174 - t260) * t252;
t225 = t192 * t175 * t132 - t194 * t130;
t137 = -Icges(6,6) * t194 + (-Icges(6,2) * t174 + t260) * t192;
t233 = t190 * t175 * t137;
t42 = (-t233 + (-t138 * t190 - t131) * t174) * t192 + t225;
t265 = t42 * t194;
t253 = t187 * t192;
t257 = t173 * t196;
t264 = t257 - (-t161 - t253) * t197 + t275;
t139 = -rSges(6,3) * t194 + (rSges(6,1) * t175 - rSges(6,2) * t174) * t192;
t251 = t192 * t196;
t211 = -t151 * rSges(6,1) + rSges(6,2) * t209;
t99 = rSges(6,3) * t251 - t211;
t79 = t139 * t251 + t194 * t99;
t235 = qJD(4) * t192;
t147 = (-Icges(5,2) * t179 - t263) * t235;
t255 = t178 * t147;
t246 = t196 * t178;
t245 = t196 * t179;
t159 = t179 * t249 + t246;
t208 = t179 * t197 + t194 * t246;
t117 = qJD(1) * t208 - qJD(4) * t159;
t157 = t194 * t245 - t254;
t158 = -t178 * t249 + t245;
t118 = -qJD(1) * t157 + qJD(4) * t158;
t242 = t118 * rSges(5,1) + t117 * rSges(5,2);
t236 = qJD(3) * t192;
t170 = t197 * t236;
t181 = qJD(2) * t196;
t241 = t170 + t181;
t240 = t177 + t181;
t133 = (-rSges(6,1) * t174 - rSges(6,2) * t175) * t252;
t231 = t192 * t237;
t105 = qJD(1) * t152 - t151 * t190;
t106 = qJD(1) * t153 - t190 * t209;
t212 = -t106 * rSges(6,1) - t105 * rSges(6,2);
t66 = rSges(6,3) * t231 - t212;
t39 = t133 * t251 + t139 * t231 + t194 * t66;
t114 = t159 * rSges(5,1) + t158 * rSges(5,2) + rSges(5,3) * t250;
t226 = -rSges(5,3) * t192 - t161;
t146 = (-Icges(5,5) * t178 - Icges(5,6) * t179) * t235;
t148 = (-Icges(5,1) * t178 - t262) * t235;
t224 = t192 * t179 * t148 - t194 * t146;
t182 = qJD(2) * t197;
t218 = -t196 * t236 + t182;
t217 = rSges(3,1) * t194 - rSges(3,2) * t192;
t215 = rSges(4,1) * t191 + rSges(4,2) * t193;
t119 = qJD(1) * t158 - qJD(4) * t157;
t204 = t208 * qJD(4);
t120 = qJD(1) * t159 - t204;
t214 = -t120 * rSges(5,1) - t119 * rSges(5,2);
t213 = -t157 * rSges(5,1) + rSges(5,2) * t208;
t210 = -pkin(1) - t217;
t59 = Icges(6,5) * t104 + Icges(6,6) * t103 - Icges(6,3) * t232;
t61 = Icges(6,4) * t104 + Icges(6,2) * t103 - Icges(6,6) * t232;
t63 = Icges(6,1) * t104 + Icges(6,4) * t103 - Icges(6,5) * t232;
t96 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t250;
t98 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t250;
t10 = -t194 * t59 + ((-t190 * t96 + t63) * t175 + (-t190 * t98 - t61) * t174) * t192;
t136 = -Icges(6,3) * t194 + (Icges(6,5) * t175 - Icges(6,6) * t174) * t192;
t18 = t103 * t137 + t104 * t138 + t152 * t131 + t153 * t132 + (t130 * t197 - t136 * t238) * t192;
t19 = t105 * t137 + t106 * t138 - t209 * t131 + t151 * t132 + (t130 * t196 + t136 * t237) * t192;
t93 = Icges(6,5) * t151 - Icges(6,6) * t209 + Icges(6,3) * t251;
t95 = Icges(6,4) * t151 - Icges(6,2) * t209 + Icges(6,6) * t251;
t97 = Icges(6,1) * t151 - Icges(6,4) * t209 + Icges(6,5) * t251;
t40 = -t194 * t93 + (-t174 * t95 + t175 * t97) * t192;
t94 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t250;
t41 = -t194 * t94 + (-t174 * t96 + t175 * t98) * t192;
t52 = t136 * t251 - t137 * t209 + t138 * t151;
t53 = t136 * t250 + t152 * t137 + t153 * t138;
t60 = Icges(6,5) * t106 + Icges(6,6) * t105 + Icges(6,3) * t231;
t62 = Icges(6,4) * t106 + Icges(6,2) * t105 + Icges(6,6) * t231;
t64 = Icges(6,1) * t106 + Icges(6,4) * t105 + Icges(6,5) * t231;
t9 = -t194 * t60 + ((-t190 * t95 + t64) * t175 + (-t190 * t97 - t62) * t174) * t192;
t207 = (t19 + t9) * t251 / 0.2e1 - (t41 + t53) * t232 / 0.2e1 + (t10 + t18 + (t40 + t52) * qJD(1)) * t250 / 0.2e1;
t206 = -t219 - t253;
t205 = (-rSges(6,3) + t187) * t192 + t228;
t203 = t196 * t226 + t256;
t26 = t151 * t97 - t209 * t95 + t251 * t93;
t27 = t151 * t98 - t209 * t96 + t251 * t94;
t28 = t152 * t95 + t153 * t97 + t250 * t93;
t29 = t152 * t96 + t153 * t98 + t250 * t94;
t202 = -t194 * (-t265 + (t10 * t197 + t196 * t9 + (-t196 * t41 + t197 * t40) * qJD(1)) * t192) - t232 * (-t53 * t194 + (t196 * t28 + t197 * t29) * t192) + (-t19 * t194 + ((t105 * t96 + t106 * t98 + t151 * t63 - t209 * t61) * t197 - t27 * t238 + (t105 * t95 + t106 * t97 + t151 * t64 - t209 * t62) * t196 + t26 * t237 + ((t196 * t59 + t237 * t94) * t197 + (t196 * t60 + t237 * t93) * t196) * t192) * t192) * t251 + (-t18 * t194 + ((t103 * t96 + t104 * t98 + t152 * t61 + t153 * t63) * t197 - t29 * t238 + (t103 * t95 + t104 * t97 + t152 * t62 + t153 * t64) * t196 + t28 * t237 + ((t197 * t59 - t238 * t94) * t197 + (t197 * t60 - t238 * t93) * t196) * t192) * t192) * t250 + (-t52 * t194 + (t196 * t26 + t197 * t27) * t192) * t231;
t200 = rSges(3,3) * t197 + t196 * t210;
t198 = -t196 * t273 + t197 * t215;
t188 = t192 ^ 2;
t185 = t197 * qJ(2);
t167 = t173 * t238;
t149 = (-rSges(5,1) * t178 - rSges(5,2) * t179) * t235;
t143 = -rSges(5,3) * t194 + (rSges(5,1) * t179 - rSges(5,2) * t178) * t192;
t140 = -Icges(5,3) * t194 + (Icges(5,5) * t179 - Icges(5,6) * t178) * t192;
t135 = t196 * rSges(3,3) + t197 * t217 + t239;
t134 = t185 + t200;
t126 = t139 * t232;
t123 = t182 + ((-rSges(3,3) - qJ(2)) * t196 + t210 * t197) * qJD(1);
t122 = qJD(1) * t200 + t240;
t121 = (t187 + t195) * t194 + (t164 - t176) * t192;
t113 = rSges(5,3) * t251 - t213;
t112 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t250;
t111 = Icges(5,1) * t157 - Icges(5,4) * t208 + Icges(5,5) * t251;
t110 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t250;
t109 = Icges(5,4) * t157 - Icges(5,2) * t208 + Icges(5,6) * t251;
t108 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t250;
t107 = Icges(5,5) * t157 - Icges(5,6) * t208 + Icges(5,3) * t251;
t91 = t99 * t250;
t90 = t185 + t198;
t89 = t215 * t196 + t183 + t277;
t88 = t161 * t197 + t114 + t257;
t87 = t203 + t213;
t85 = -t185 + (-t165 + t173) * t197 + t206 * t196;
t84 = qJD(1) * t198 + t170 + t240;
t83 = ((-qJ(2) - t215) * t196 - t277) * qJD(1) + t218;
t82 = -t194 * t114 - t143 * t250;
t81 = t113 * t194 + t143 * t251;
t80 = -t194 * t100 - t139 * t250;
t78 = rSges(5,3) * t231 - t214;
t77 = -rSges(5,3) * t232 + t242;
t76 = Icges(5,1) * t120 + Icges(5,4) * t119 + Icges(5,5) * t231;
t75 = Icges(5,1) * t118 + Icges(5,4) * t117 - Icges(5,5) * t232;
t74 = Icges(5,4) * t120 + Icges(5,2) * t119 + Icges(5,6) * t231;
t73 = Icges(5,4) * t118 + Icges(5,2) * t117 - Icges(5,6) * t232;
t72 = Icges(5,5) * t120 + Icges(5,6) * t119 + Icges(5,3) * t231;
t71 = Icges(5,5) * t118 + Icges(5,6) * t117 - Icges(5,3) * t232;
t70 = -t187 * t250 - t275;
t69 = t165 * t197 + t196 * t205 + t185 + t211;
t68 = -t167 - pkin(4) * t204 + (t197 * t206 + t276) * qJD(1);
t58 = t140 * t250 + t158 * t141 + t159 * t142;
t57 = t140 * t251 - t141 * t208 + t142 * t157;
t55 = -t100 * t251 + t91;
t54 = t66 * t250;
t49 = qJD(1) * t203 + t241 + t242;
t48 = t226 * t237 - t167 + t214 + t218;
t47 = (t274 * qJD(4) - t255) * t192 + t224;
t46 = t194 * t78 + (t143 * t237 + t149 * t196) * t192;
t45 = -t194 * t77 + (t143 * t238 - t149 * t197) * t192;
t44 = -t108 * t194 + (-t110 * t178 + t112 * t179) * t192;
t43 = -t107 * t194 + (-t109 * t178 + t111 * t179) * t192;
t38 = -t133 * t250 - t194 * t65 + t126;
t35 = t108 * t250 + t158 * t110 + t159 * t112;
t34 = t107 * t250 + t158 * t109 + t159 * t111;
t33 = t108 * t251 - t110 * t208 + t112 * t157;
t32 = t107 * t251 - t109 * t208 + t111 * t157;
t31 = t208 * t266 + (t197 * t205 - t276) * qJD(1) + t212 + t218;
t30 = (-rSges(6,3) * t192 + t228) * t238 + t201 + t241 + t243;
t25 = t264 * t194 + (-t121 - t139) * t250;
t24 = t121 * t251 + t194 * t85 + t79;
t23 = t119 * t141 + t120 * t142 - t208 * t147 + t157 * t148 + (t140 * t237 + t146 * t196) * t192;
t22 = t117 * t141 + t118 * t142 + t158 * t147 + t159 * t148 + (-t140 * t238 + t146 * t197) * t192;
t21 = t91 + (t264 * t196 + t197 * t85) * t192;
t20 = (-t196 * t77 + t197 * t78 + (-t113 * t196 - t114 * t197) * qJD(1)) * t192;
t15 = -t178 * t188 * t234 + t121 * t231 + t194 * t68 + t39;
t14 = t188 * t223 + t126 + t267 * t194 + (t121 * t238 - t133 * t197) * t192;
t13 = t54 + (-t196 * t65 + (-t100 * t197 - t196 * t99) * qJD(1)) * t192;
t12 = -t194 * t71 + (-t178 * t73 + t179 * t75 + (-t110 * t179 - t112 * t178) * qJD(4)) * t192;
t11 = -t194 * t72 + (-t178 * t74 + t179 * t76 + (-t109 * t179 - t111 * t178) * qJD(4)) * t192;
t4 = t54 + (t197 * t68 + t267 * t196 + (t264 * t197 + (-t85 - t99) * t196) * qJD(1)) * t192;
t1 = [(t30 * t70 + t31 * t69) * t271 - t174 * t138 * t252 + (t48 * t87 + t49 * t88) * t272 + 0.2e1 * m(3) * (t122 * t135 + t123 * t134) + 0.2e1 * m(4) * (t83 * t90 + t84 * t89) + t224 + t225 + t274 * t235 + (-t174 * t131 - t233 - t255) * t192; m(6) * (t196 * t31 - t197 * t30 + (t196 * t70 + t197 * t69) * qJD(1)) + m(5) * (t196 * t48 - t197 * t49 + (t196 * t88 + t197 * t87) * qJD(1)) + m(3) * (-t122 * t197 + t196 * t123 + (t134 * t197 + t135 * t196) * qJD(1)) + m(4) * (t196 * t83 - t197 * t84 + (t196 * t89 + t197 * t90) * qJD(1)); 0; 0.2e1 * ((t196 * t30 + t197 * t31 + t237 * t70 - t238 * t69) * t269 + (t196 * t49 + t197 * t48 + t237 * t88 - t238 * t87) * t270 + m(4) * (t196 * t84 + t197 * t83 + t237 * t89 - t238 * t90) / 0.2e1) * t192; 0; 0; (-t42 - t47) * t194 + m(6) * (t14 * t70 + t15 * t69 + t24 * t31 + t25 * t30) + m(5) * (t45 * t88 + t46 * t87 + t48 * t81 + t49 * t82) + ((t12 / 0.2e1 + t22 / 0.2e1) * t197 + (t11 / 0.2e1 + t23 / 0.2e1) * t196 + ((t43 / 0.2e1 + t57 / 0.2e1) * t197 + (-t44 / 0.2e1 - t58 / 0.2e1) * t196) * qJD(1)) * t192 + t207; m(5) * (t46 * t196 - t197 * t45 + (t196 * t82 + t197 * t81) * qJD(1)) + m(6) * (-t14 * t197 + t15 * t196 + (t196 * t25 + t197 * t24) * qJD(1)); 0.2e1 * (-m(5) * t20 / 0.2e1 - m(6) * t4 / 0.2e1) * t194 + 0.2e1 * ((t196 * t45 + t197 * t46 + t237 * t82 - t238 * t81) * t270 + (t14 * t196 + t15 * t197 + t237 * t25 - t238 * t24) * t269) * t192; (t82 * t45 + t81 * t46 + (t113 * t197 - t114 * t196) * t20 * t192) * t272 - (-t58 * t194 + (t196 * t34 + t197 * t35) * t192) * t232 + (-t22 * t194 + ((t117 * t110 + t118 * t112 + t158 * t73 + t159 * t75) * t197 - t35 * t238 + (t117 * t109 + t118 * t111 + t158 * t74 + t159 * t76) * t196 + t34 * t237 + ((-t108 * t238 + t197 * t71) * t197 + (-t107 * t238 + t197 * t72) * t196) * t192) * t192) * t250 + (-t57 * t194 + (t196 * t32 + t197 * t33) * t192) * t231 + (-t23 * t194 + ((t119 * t110 + t120 * t112 + t157 * t75 - t208 * t73) * t197 - t33 * t238 + (t119 * t109 + t120 * t111 + t157 * t76 - t208 * t74) * t196 + t32 * t237 + ((t108 * t237 + t196 * t71) * t197 + (t107 * t237 + t196 * t72) * t196) * t192) * t192) * t251 - t194 * (-t47 * t194 + (t11 * t196 + t12 * t197 + (-t196 * t44 + t197 * t43) * qJD(1)) * t192) + (t14 * t25 + t15 * t24 + t21 * t4) * t271 + t202; m(6) * (t30 * t80 + t31 * t79 + t38 * t70 + t39 * t69) - t265 + t207; m(6) * (t39 * t196 - t197 * t38 + (t196 * t80 + t197 * t79) * qJD(1)); m(6) * (-t13 * t194 + (t196 * t38 + t197 * t39 + (-t196 * t79 + t197 * t80) * qJD(1)) * t192); m(6) * (t13 * t21 + t14 * t80 + t15 * t79 + t24 * t39 + t25 * t38 + t4 * t55) + t202; (t13 * t55 + t38 * t80 + t39 * t79) * t271 + t202;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
