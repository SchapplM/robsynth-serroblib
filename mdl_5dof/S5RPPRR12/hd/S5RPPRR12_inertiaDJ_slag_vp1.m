% Calculate time derivative of joint inertia matrix for
% S5RPPRR12
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:07:05
% DurationCPUTime: 4.32s
% Computational Cost: add. (8763->460), mult. (13081->688), div. (0->0), fcn. (12355->8), ass. (0->233)
t152 = pkin(8) + qJ(4);
t143 = sin(t152);
t144 = cos(t152);
t159 = sin(qJ(1));
t225 = qJD(4) * t159;
t209 = t144 * t225;
t161 = cos(qJ(1));
t228 = qJD(1) * t161;
t291 = t143 * t228 + t209;
t273 = rSges(6,3) + pkin(7);
t213 = t273 * t144;
t270 = pkin(4) * t143;
t290 = -t213 + t270;
t196 = rSges(5,1) * t143 + rSges(5,2) * t144;
t169 = t161 * t196;
t251 = Icges(5,4) * t143;
t179 = Icges(5,2) * t144 + t251;
t100 = Icges(5,6) * t161 + t179 * t159;
t250 = Icges(5,4) * t144;
t182 = Icges(5,1) * t143 + t250;
t102 = Icges(5,5) * t161 + t182 * t159;
t174 = t100 * t144 + t102 * t143;
t289 = t174 * t161;
t155 = sin(pkin(8));
t271 = pkin(3) * t155;
t140 = t161 * t271;
t149 = t161 * qJ(2);
t157 = -pkin(6) - qJ(3);
t215 = t159 * t157 + t140 + t149;
t274 = -rSges(5,3) - pkin(1);
t86 = t274 * t159 + t169 + t215;
t150 = t161 * rSges(5,3);
t238 = t144 * t159;
t239 = t143 * t159;
t104 = rSges(5,1) * t239 + rSges(5,2) * t238 + t150;
t151 = t161 * pkin(1);
t230 = t159 * qJ(2) + t151;
t171 = -t157 * t161 + t159 * t271 + t230;
t87 = t171 + t104;
t288 = -t159 * t86 + t161 * t87;
t158 = sin(qJ(5));
t160 = cos(qJ(5));
t248 = Icges(6,4) * t160;
t178 = -Icges(6,2) * t158 + t248;
t95 = Icges(6,6) * t143 + t178 * t144;
t249 = Icges(6,4) * t158;
t181 = Icges(6,1) * t160 - t249;
t96 = Icges(6,5) * t143 + t181 * t144;
t287 = -t158 * t96 - t160 * t95;
t231 = qJ(2) * t228 + qJD(2) * t159;
t214 = qJD(3) * t161 + t231;
t229 = qJD(1) * t159;
t184 = qJD(1) * t140 + t157 * t229 + t214;
t212 = t144 * t228;
t217 = t291 * rSges(5,1) + rSges(5,2) * t212;
t227 = qJD(4) * t143;
t54 = (-rSges(5,2) * t227 + t274 * qJD(1)) * t159 + t184 + t217;
t266 = rSges(5,2) * t143;
t124 = rSges(5,1) * t144 - t266;
t147 = qJD(2) * t161;
t203 = -qJD(3) * t159 + t147;
t198 = t157 * t228 + t203;
t206 = -qJ(2) - t271;
t223 = qJD(4) * t161;
t55 = t124 * t223 + (t274 * t161 + (-t196 + t206) * t159) * qJD(1) + t198;
t286 = t159 * t55 - t161 * t54;
t201 = qJD(1) * t143 + qJD(5);
t285 = -t144 * t223 + t201 * t159;
t176 = Icges(5,5) * t143 + Icges(5,6) * t144;
t284 = -Icges(5,3) * t159 + t176 * t161;
t283 = -Icges(5,6) * t159 + t179 * t161;
t282 = -Icges(5,5) * t159 + t182 * t161;
t281 = 2 * m(5);
t280 = 2 * m(6);
t153 = t159 ^ 2;
t154 = t161 ^ 2;
t279 = t143 / 0.2e1;
t278 = -t159 / 0.2e1;
t276 = t161 / 0.2e1;
t275 = rSges(3,2) - pkin(1);
t272 = m(5) * t124;
t269 = pkin(4) * t144;
t268 = t159 * pkin(1);
t224 = qJD(4) * t160;
t226 = qJD(4) * t144;
t264 = t158 * t95;
t175 = Icges(6,5) * t160 - Icges(6,6) * t158;
t222 = qJD(5) * t144;
t60 = (-Icges(6,5) * t158 - Icges(6,6) * t160) * t222 + (Icges(6,3) * t144 - t175 * t143) * qJD(4);
t62 = (-Icges(6,1) * t158 - t248) * t222 + (Icges(6,5) * t144 - t181 * t143) * qJD(4);
t94 = Icges(6,3) * t143 + t175 * t144;
t163 = t144 * t160 * t62 + t94 * t226 + t227 * t264 + (-t224 * t96 + t60) * t143;
t61 = (-Icges(6,2) * t160 - t249) * t222 + (Icges(6,6) * t144 - t178 * t143) * qJD(4);
t265 = t158 * t61;
t45 = t143 * t94 + (t160 * t96 - t264) * t144;
t267 = ((qJD(5) * t287 - t265) * t144 + t163) * t143 + t45 * t226;
t262 = t159 * rSges(5,3);
t256 = rSges(4,3) + qJ(3);
t135 = pkin(4) * t239;
t233 = t160 * t161;
t235 = t159 * t158;
t112 = -t143 * t235 + t233;
t234 = t159 * t160;
t236 = t158 * t161;
t113 = t143 * t234 + t236;
t232 = t113 * rSges(6,1) + t112 * rSges(6,2);
t82 = -rSges(6,3) * t238 + t232;
t255 = pkin(7) * t238 - t135 - t82;
t199 = pkin(7) * t144 - t270;
t114 = t143 * t236 + t234;
t115 = -t143 * t233 + t235;
t195 = -t115 * rSges(6,1) - t114 * rSges(6,2);
t237 = t144 * t161;
t83 = rSges(6,3) * t237 - t195;
t254 = t199 * t161 + t83;
t194 = rSges(6,1) * t160 - rSges(6,2) * t158;
t67 = (-rSges(6,1) * t158 - rSges(6,2) * t160) * t222 + (rSges(6,3) * t144 - t194 * t143) * qJD(4);
t253 = t199 * qJD(4) + t67;
t125 = pkin(7) * t143 + t269;
t97 = rSges(6,3) * t143 + t194 * t144;
t252 = t125 + t97;
t98 = Icges(5,3) * t161 + t176 * t159;
t244 = qJD(1) * t98;
t243 = t100 * t143;
t242 = t143 * t283;
t241 = t102 * t144;
t240 = t144 * t282;
t221 = -pkin(1) - t256;
t207 = t143 * t225;
t202 = qJD(5) * t143 + qJD(1);
t65 = -t202 * t234 + (-t201 * t161 - t209) * t158;
t66 = t201 * t233 + (t144 * t224 - t202 * t158) * t159;
t220 = t66 * rSges(6,1) + t65 * rSges(6,2) + rSges(6,3) * t207;
t78 = Icges(6,4) * t113 + Icges(6,2) * t112 - Icges(6,6) * t238;
t80 = Icges(6,1) * t113 + Icges(6,4) * t112 - Icges(6,5) * t238;
t191 = t158 * t78 - t160 * t80;
t76 = Icges(6,5) * t113 + Icges(6,6) * t112 - Icges(6,3) * t238;
t28 = t143 * t76 - t191 * t144;
t39 = t112 * t95 + t113 * t96 - t94 * t238;
t219 = t28 / 0.2e1 + t39 / 0.2e1;
t79 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t237;
t81 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t237;
t190 = t158 * t79 - t160 * t81;
t77 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t237;
t29 = t143 * t77 - t190 * t144;
t40 = t114 * t95 + t115 * t96 + t94 * t237;
t218 = -t29 / 0.2e1 - t40 / 0.2e1;
t216 = -t291 * pkin(4) - pkin(7) * t207;
t211 = t143 * t223;
t205 = qJD(1) * t252;
t119 = t196 * qJD(4);
t204 = t119 * (t153 + t154);
t172 = t202 * t161;
t63 = -t158 * t285 + t160 * t172;
t64 = t158 * t172 + t160 * t285;
t200 = t64 * rSges(6,1) + t63 * rSges(6,2);
t197 = rSges(4,1) * t155 + rSges(4,2) * cos(pkin(8));
t24 = t112 * t78 + t113 * t80 - t76 * t238;
t25 = t112 * t79 + t113 * t81 - t77 * t238;
t17 = t25 * t159 + t161 * t24;
t189 = t159 * t24 - t161 * t25;
t26 = t114 * t78 + t115 * t80 + t76 * t237;
t27 = t114 * t79 + t115 * t81 + t77 * t237;
t18 = t27 * t159 + t161 * t26;
t188 = t159 * t26 - t161 * t27;
t187 = t29 * t159 + t161 * t28;
t186 = t159 * t28 - t161 * t29;
t185 = t159 * t83 + t161 * t82;
t183 = Icges(5,1) * t144 - t251;
t180 = -Icges(5,2) * t143 + t250;
t177 = Icges(5,5) * t144 - Icges(5,6) * t143;
t173 = -t143 * t282 - t144 * t283;
t170 = rSges(3,3) * t161 + t275 * t159;
t168 = t173 * t159;
t167 = qJD(4) * t183;
t166 = qJD(4) * t180;
t165 = t207 - t212;
t164 = -t144 * t229 - t211;
t162 = t221 * t159 + t197 * t161;
t111 = -rSges(3,2) * t161 + t159 * rSges(3,3) + t230;
t110 = t149 + t170;
t105 = t262 - t169;
t93 = t147 + (t275 * t161 + (-rSges(3,3) - qJ(2)) * t159) * qJD(1);
t92 = t170 * qJD(1) + t231;
t91 = t197 * t159 + t256 * t161 + t230;
t90 = t149 + t162;
t85 = t252 * t161;
t84 = t252 * t159;
t71 = qJD(1) * t284 + t177 * t225;
t70 = -t177 * t223 + t244;
t69 = (t221 * t161 + (-qJ(2) - t197) * t159) * qJD(1) + t203;
t68 = t162 * qJD(1) + t214;
t53 = -t159 * t213 + t135 + t171 + t232;
t52 = t290 * t161 + t195 + t215 - t268;
t51 = -t143 * t83 + t97 * t237;
t50 = t143 * t82 + t97 * t238;
t49 = -t159 * t284 - t173 * t161;
t48 = t159 * t98 - t289;
t47 = -t161 * t284 + t168;
t46 = t174 * t159 + t161 * t98;
t43 = t185 * t144;
t42 = t253 * t159 + t161 * t205;
t41 = t159 * t205 - t253 * t161;
t38 = -rSges(6,3) * t212 + t220;
t37 = t164 * rSges(6,3) + t200;
t36 = Icges(6,1) * t66 + Icges(6,4) * t65 + t165 * Icges(6,5);
t35 = Icges(6,1) * t64 + Icges(6,4) * t63 + t164 * Icges(6,5);
t34 = Icges(6,4) * t66 + Icges(6,2) * t65 + t165 * Icges(6,6);
t33 = Icges(6,4) * t64 + Icges(6,2) * t63 + t164 * Icges(6,6);
t32 = Icges(6,5) * t66 + Icges(6,6) * t65 + t165 * Icges(6,3);
t31 = Icges(6,5) * t64 + Icges(6,6) * t63 + t164 * Icges(6,3);
t30 = t255 * t159 + t254 * t161;
t23 = (t273 * t143 + t269) * t223 + (-t151 + (t206 - t290) * t159) * qJD(1) + t198 - t200;
t22 = (-t161 * t213 - t268) * qJD(1) + t184 - t216 + t220;
t21 = (-t97 * t225 + t38) * t143 + (qJD(4) * t82 + t159 * t67 + t97 * t228) * t144;
t20 = (-t97 * t223 - t37) * t143 + (-qJD(4) * t83 + t161 * t67 - t97 * t229) * t144;
t16 = t94 * t207 + t112 * t61 + t113 * t62 + t65 * t95 + t66 * t96 + (-t159 * t60 - t94 * t228) * t144;
t15 = -t94 * t211 + t114 * t61 + t115 * t62 + t63 * t95 + t64 * t96 + (t161 * t60 - t94 * t229) * t144;
t14 = t185 * t227 + (-t159 * t37 - t161 * t38 + (t159 * t82 - t161 * t83) * qJD(1)) * t144;
t13 = (-t254 * qJD(1) + t216 - t38) * t159 + (t37 - t125 * t223 + (t135 + t255) * qJD(1)) * t161;
t12 = t40 * t143 - t188 * t144;
t11 = t39 * t143 - t189 * t144;
t10 = (t190 * qJD(4) + t31) * t143 + (qJD(4) * t77 - t158 * t33 + t160 * t35 + (-t158 * t81 - t160 * t79) * qJD(5)) * t144;
t9 = (t191 * qJD(4) + t32) * t143 + (qJD(4) * t76 - t158 * t34 + t160 * t36 + (-t158 * t80 - t160 * t78) * qJD(5)) * t144;
t8 = t77 * t207 + t112 * t33 + t113 * t35 + t65 * t79 + t66 * t81 + (-t159 * t31 - t77 * t228) * t144;
t7 = t76 * t207 + t112 * t34 + t113 * t36 + t65 * t78 + t66 * t80 + (-t159 * t32 - t76 * t228) * t144;
t6 = -t77 * t211 + t114 * t33 + t115 * t35 + t63 * t79 + t64 * t81 + (t161 * t31 - t77 * t229) * t144;
t5 = -t76 * t211 + t114 * t34 + t115 * t36 + t63 * t78 + t64 * t80 + (t161 * t32 - t76 * t229) * t144;
t4 = -t189 * qJD(1) + t8 * t159 + t161 * t7;
t3 = -t188 * qJD(1) + t6 * t159 + t161 * t5;
t2 = (t189 * qJD(4) + t16) * t143 + (-t17 * qJD(1) + qJD(4) * t39 - t159 * t7 + t161 * t8) * t144;
t1 = (t188 * qJD(4) + t15) * t143 + (-qJD(1) * t18 + qJD(4) * t40 - t159 * t5 + t161 * t6) * t144;
t19 = [(t22 * t53 + t23 * t52) * t280 - t143 * t167 - t182 * t226 + t179 * t227 + (t54 * t87 + t55 * t86) * t281 + 0.2e1 * m(4) * (t68 * t91 + t69 * t90) + 0.2e1 * m(3) * (t110 * t93 + t111 * t92) + t163 + t287 * t222 + (-t265 - t166) * t144; m(6) * (t159 * t23 - t161 * t22 + (t159 * t53 + t161 * t52) * qJD(1)) + m(5) * ((t159 * t87 + t161 * t86) * qJD(1) + t286) + m(4) * (t159 * t69 - t161 * t68 + (t159 * t91 + t161 * t90) * qJD(1)) + m(3) * (t159 * t93 - t161 * t92 + (t110 * t161 + t111 * t159) * qJD(1)); 0; m(6) * (t159 * t22 + t161 * t23 + (-t159 * t52 + t161 * t53) * qJD(1)) + m(5) * (qJD(1) * t288 + t159 * t54 + t161 * t55) + m(4) * (t159 * t68 + t161 * t69 + (-t159 * t90 + t161 * t91) * qJD(1)); 0; 0; m(6) * (-t22 * t85 + t23 * t84 + t41 * t53 + t42 * t52) + m(5) * (t119 * t288 + t124 * t286) - (t154 / 0.2e1 + t153 / 0.2e1) * t176 * qJD(4) + ((t87 * t272 + t243 / 0.2e1 - t241 / 0.2e1 - t219) * t159 + (t86 * t272 + t242 / 0.2e1 - t240 / 0.2e1 - t218) * t161) * qJD(1) + (-qJD(4) * t173 - t143 * (t100 * qJD(1) - t180 * t223) + t144 * (t102 * qJD(1) - t183 * t223) + t10 + t15) * t159 / 0.2e1 + (-qJD(4) * t174 - t143 * (qJD(1) * t283 + t159 * t166) + t144 * (qJD(1) * t282 + t159 * t167) + t16 + t9) * t276; m(6) * (t42 * t159 - t41 * t161 + (-t159 * t85 + t161 * t84) * qJD(1)) - m(5) * t204; m(6) * (t41 * t159 + t161 * t42 + (-t159 * t84 - t161 * t85) * qJD(1)); ((-t159 * t104 + t105 * t161) * (-t159 * t217 + (-t124 * t154 + t153 * t266) * qJD(4) + ((-t104 + t150) * t161 + (-t105 + t169 + t262) * t159) * qJD(1)) - t124 * t204) * t281 + t161 * ((t161 * t71 + (t47 + t289) * qJD(1)) * t161 + (-t46 * qJD(1) + (-t226 * t282 + t227 * t283) * t159 + (t70 + (t241 - t243) * qJD(4) + (t173 - t98) * qJD(1)) * t161) * t159) + t159 * ((t159 * t70 + (-t48 + t168) * qJD(1)) * t159 + (t49 * qJD(1) + (t100 * t227 - t102 * t226 + t244) * t161 + (t71 + (t240 - t242) * qJD(4) + t174 * qJD(1)) * t159) * t161) + (t30 * t13 - t41 * t85 + t42 * t84) * t280 + t161 * t4 + t159 * t3 + (-t47 * t159 - t161 * t46 - t17) * t229 + (t49 * t159 + t161 * t48 + t18) * t228; m(6) * (t20 * t52 + t21 * t53 + t22 * t50 + t23 * t51) + (t219 * t159 + t218 * t161) * t227 + ((t10 / 0.2e1 + t15 / 0.2e1) * t161 + (-t9 / 0.2e1 - t16 / 0.2e1) * t159 + (t218 * t159 - t219 * t161) * qJD(1)) * t144 + t267; m(6) * (t20 * t159 - t161 * t21 + (t159 * t50 + t161 * t51) * qJD(1)); m(6) * (t21 * t159 + t161 * t20 + (-t159 * t51 + t161 * t50) * qJD(1)); m(6) * (-t43 * t13 + t14 * t30 + t20 * t84 - t21 * t85 + t50 * t41 + t51 * t42) + (t2 / 0.2e1 + qJD(1) * t12 / 0.2e1 - t18 * t227 / 0.2e1 + (qJD(1) * t29 + t9) * t279) * t161 + (-qJD(1) * t11 / 0.2e1 + t17 * t227 / 0.2e1 + t1 / 0.2e1 + (-qJD(1) * t28 + t10) * t279) * t159 + (t4 * t278 + t3 * t276 + qJD(4) * t187 / 0.2e1 + (-t161 * t17 / 0.2e1 + t18 * t278) * qJD(1)) * t144; (-t14 * t43 + t20 * t51 + t21 * t50) * t280 + ((t159 * t11 - t161 * t12 + t186 * t143) * qJD(4) + t267) * t143 + (-t159 * t2 + t161 * t1 + t143 * (t10 * t161 - t159 * t9) + (t45 * t143 - t186 * t144) * qJD(4) + (-t161 * t11 - t159 * t12 - t143 * t187) * qJD(1)) * t144;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
