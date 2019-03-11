% Calculate time derivative of joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:39
% EndTime: 2019-03-09 01:31:47
% DurationCPUTime: 4.35s
% Computational Cost: add. (13953->466), mult. (13515->690), div. (0->0), fcn. (12693->10), ass. (0->240)
t158 = pkin(10) + qJ(5);
t153 = sin(t158);
t159 = qJ(1) + pkin(9);
t154 = sin(t159);
t155 = cos(t158);
t237 = qJD(5) * t155;
t156 = cos(t159);
t241 = qJD(1) * t156;
t305 = t153 * t241 + t154 * t237;
t284 = rSges(7,3) + pkin(8);
t224 = t284 * t155;
t280 = pkin(5) * t153;
t304 = -t224 + t280;
t163 = sin(qJ(6));
t165 = cos(qJ(6));
t262 = Icges(7,4) * t165;
t185 = -Icges(7,2) * t163 + t262;
t106 = Icges(7,6) * t153 + t155 * t185;
t263 = Icges(7,4) * t163;
t188 = Icges(7,1) * t165 - t263;
t107 = Icges(7,5) * t153 + t155 * t188;
t303 = -t106 * t165 - t107 * t163;
t264 = Icges(6,4) * t155;
t189 = Icges(6,1) * t153 + t264;
t101 = Icges(6,5) * t156 + t154 * t189;
t265 = Icges(6,4) * t153;
t186 = Icges(6,2) * t155 + t265;
t99 = Icges(6,6) * t156 + t154 * t186;
t191 = t101 * t153 + t155 * t99;
t302 = t156 * t191;
t204 = rSges(6,1) * t153 + rSges(6,2) * t155;
t177 = t156 * t204;
t157 = cos(qJ(1)) * pkin(1);
t301 = -t156 * pkin(2) - t157;
t282 = sin(qJ(1)) * pkin(1);
t285 = -rSges(6,3) - pkin(2);
t178 = t154 * t285 - t282;
t160 = sin(pkin(10));
t281 = pkin(4) * t160;
t141 = t156 * t281;
t148 = t156 * qJ(3);
t162 = -pkin(7) - qJ(4);
t298 = t154 * t162 + t141 + t148;
t59 = t177 + t178 + t298;
t149 = t156 * rSges(6,3);
t251 = t154 * t155;
t252 = t153 * t154;
t103 = rSges(6,1) * t252 + rSges(6,2) * t251 + t149;
t225 = t154 * qJ(3) - t301;
t176 = t154 * t281 - t156 * t162 + t225;
t60 = t176 + t103;
t300 = -t154 * t59 + t156 * t60;
t243 = qJ(3) * t241 + qJD(3) * t154;
t226 = qJD(4) * t156 + t243;
t240 = qJD(1) * t162;
t192 = qJD(1) * t141 + t154 * t240 + t226;
t238 = qJD(5) * t154;
t221 = t153 * t238;
t222 = t155 * t241;
t228 = rSges(6,1) * t305 + rSges(6,2) * t222;
t51 = -rSges(6,2) * t221 + qJD(1) * t178 + t192 + t228;
t277 = rSges(6,2) * t153;
t125 = rSges(6,1) * t155 - t277;
t146 = qJD(3) * t156;
t213 = -qJD(4) * t154 + t146;
t207 = t156 * t240 + t213;
t216 = -qJ(3) - t281;
t236 = qJD(5) * t156;
t52 = t125 * t236 + (-t157 + t285 * t156 + (-t204 + t216) * t154) * qJD(1) + t207;
t299 = t154 * t52 - t156 * t51;
t183 = Icges(6,5) * t153 + Icges(6,6) * t155;
t297 = -Icges(6,3) * t154 + t156 * t183;
t296 = -Icges(6,6) * t154 + t156 * t186;
t295 = -Icges(6,5) * t154 + t156 * t189;
t212 = qJD(6) * t153 + qJD(1);
t234 = qJD(5) * t165;
t294 = -t155 * t234 + t163 * t212;
t235 = qJD(5) * t163;
t293 = t155 * t235 + t165 * t212;
t292 = 2 * m(6);
t291 = 2 * m(7);
t151 = t154 ^ 2;
t152 = t156 ^ 2;
t290 = t153 / 0.2e1;
t289 = -t154 / 0.2e1;
t287 = t156 / 0.2e1;
t286 = rSges(4,2) - pkin(2);
t283 = m(6) * t125;
t279 = pkin(5) * t155;
t182 = Icges(7,5) * t165 - Icges(7,6) * t163;
t105 = Icges(7,3) * t153 + t155 * t182;
t233 = qJD(6) * t155;
t83 = (-Icges(7,5) * t163 - Icges(7,6) * t165) * t233 + (Icges(7,3) * t155 - t153 * t182) * qJD(5);
t85 = (-Icges(7,1) * t163 - t262) * t233 + (Icges(7,5) * t155 - t153 * t188) * qJD(5);
t169 = t155 * t165 * t85 + t105 * t237 + (t106 * t235 - t107 * t234 + t83) * t153;
t84 = (-Icges(7,2) * t165 - t263) * t233 + (Icges(7,6) * t155 - t153 * t185) * qJD(5);
t270 = t163 * t84;
t56 = t105 * t153 + (-t106 * t163 + t107 * t165) * t155;
t278 = ((qJD(6) * t303 - t270) * t155 + t169) * t153 + t56 * t237;
t276 = t153 * t99;
t275 = t154 * rSges(6,3);
t269 = rSges(5,3) + qJ(4);
t136 = pkin(5) * t252;
t246 = t156 * t165;
t250 = t154 * t163;
t113 = -t153 * t250 + t246;
t247 = t156 * t163;
t249 = t154 * t165;
t114 = t153 * t249 + t247;
t244 = t114 * rSges(7,1) + t113 * rSges(7,2);
t80 = -rSges(7,3) * t251 + t244;
t268 = pkin(8) * t251 - t136 - t80;
t208 = pkin(8) * t155 - t280;
t115 = t153 * t247 + t249;
t116 = -t153 * t246 + t250;
t203 = -rSges(7,1) * t116 - rSges(7,2) * t115;
t248 = t155 * t156;
t81 = rSges(7,3) * t248 - t203;
t267 = t208 * t156 + t81;
t202 = rSges(7,1) * t165 - rSges(7,2) * t163;
t86 = (-rSges(7,1) * t163 - rSges(7,2) * t165) * t233 + (rSges(7,3) * t155 - t153 * t202) * qJD(5);
t266 = t208 * qJD(5) + t86;
t97 = Icges(6,3) * t156 + t154 * t183;
t258 = qJD(1) * t97;
t257 = t153 * t296;
t256 = t101 * t155;
t255 = t155 * t295;
t108 = rSges(7,3) * t153 + t155 * t202;
t126 = pkin(8) * t153 + t279;
t245 = t108 + t126;
t242 = qJD(1) * t154;
t239 = qJD(5) * t153;
t232 = -pkin(2) - t269;
t211 = qJD(1) * t153 + qJD(6);
t179 = t211 * t163;
t65 = -t293 * t154 - t156 * t179;
t180 = t165 * t211;
t66 = -t294 * t154 + t156 * t180;
t231 = t66 * rSges(7,1) + t65 * rSges(7,2) + rSges(7,3) * t221;
t76 = Icges(7,4) * t114 + Icges(7,2) * t113 - Icges(7,6) * t251;
t78 = Icges(7,1) * t114 + Icges(7,4) * t113 - Icges(7,5) * t251;
t194 = t163 * t76 - t165 * t78;
t74 = Icges(7,5) * t114 + Icges(7,6) * t113 - Icges(7,3) * t251;
t29 = t153 * t74 - t155 * t194;
t41 = -t105 * t251 + t106 * t113 + t107 * t114;
t230 = t41 / 0.2e1 + t29 / 0.2e1;
t77 = Icges(7,4) * t116 + Icges(7,2) * t115 + Icges(7,6) * t248;
t79 = Icges(7,1) * t116 + Icges(7,4) * t115 + Icges(7,5) * t248;
t193 = t163 * t77 - t165 * t79;
t75 = Icges(7,5) * t116 + Icges(7,6) * t115 + Icges(7,3) * t248;
t30 = t153 * t75 - t155 * t193;
t42 = t105 * t248 + t106 * t115 + t107 * t116;
t229 = -t42 / 0.2e1 - t30 / 0.2e1;
t227 = -pkin(5) * t305 - pkin(8) * t221;
t220 = t153 * t236;
t120 = t204 * qJD(5);
t215 = (t151 + t152) * t120;
t214 = qJD(1) * t245;
t63 = -t154 * t179 + t293 * t156;
t64 = t154 * t180 + t294 * t156;
t210 = t64 * rSges(7,1) + t63 * rSges(7,2);
t209 = -pkin(2) * t154 - t282;
t205 = rSges(5,1) * t160 + rSges(5,2) * cos(pkin(10));
t24 = t113 * t76 + t114 * t78 - t251 * t74;
t25 = t113 * t77 + t114 * t79 - t251 * t75;
t15 = t154 * t25 + t156 * t24;
t199 = t154 * t24 - t156 * t25;
t26 = t115 * t76 + t116 * t78 + t248 * t74;
t27 = t115 * t77 + t116 * t79 + t248 * t75;
t16 = t154 * t27 + t156 * t26;
t198 = t154 * t26 - t156 * t27;
t197 = t30 * t154 + t29 * t156;
t196 = t154 * t29 - t156 * t30;
t195 = t154 * t81 + t156 * t80;
t190 = Icges(6,1) * t155 - t265;
t187 = -Icges(6,2) * t153 + t264;
t184 = Icges(6,5) * t155 - Icges(6,6) * t153;
t181 = -t153 * t295 - t155 * t296;
t175 = t181 * t154;
t174 = qJD(5) * t190;
t173 = qJD(5) * t187;
t171 = t221 - t222;
t170 = -t155 * t242 - t220;
t168 = rSges(4,3) * t156 + t154 * t286 - t282;
t167 = t154 * t232 + t156 * t205 - t282;
t104 = t275 - t177;
t96 = -rSges(4,2) * t156 + rSges(4,3) * t154 + t225;
t95 = t148 + t168;
t92 = t146 + (-t157 + t286 * t156 + (-rSges(4,3) - qJ(3)) * t154) * qJD(1);
t91 = qJD(1) * t168 + t243;
t90 = t154 * t205 + t156 * t269 + t225;
t89 = t148 + t167;
t88 = t245 * t156;
t87 = t245 * t154;
t68 = t297 * qJD(1) + t184 * t238;
t67 = -t184 * t236 + t258;
t58 = (-t157 + t232 * t156 + (-qJ(3) - t205) * t154) * qJD(1) + t213;
t57 = qJD(1) * t167 + t226;
t55 = t108 * t248 - t153 * t81;
t54 = t108 * t251 + t153 * t80;
t50 = -t154 * t224 + t136 + t176 + t244;
t49 = t156 * t304 + t203 + t209 + t298;
t48 = -t154 * t297 - t156 * t181;
t47 = t154 * t97 - t302;
t46 = -t156 * t297 + t175;
t45 = t154 * t191 + t156 * t97;
t44 = t154 * t266 + t156 * t214;
t43 = t154 * t214 - t156 * t266;
t40 = t195 * t155;
t39 = -rSges(7,3) * t222 + t231;
t38 = rSges(7,3) * t170 + t210;
t37 = Icges(7,1) * t66 + Icges(7,4) * t65 + Icges(7,5) * t171;
t36 = Icges(7,1) * t64 + Icges(7,4) * t63 + Icges(7,5) * t170;
t35 = Icges(7,4) * t66 + Icges(7,2) * t65 + Icges(7,6) * t171;
t34 = Icges(7,4) * t64 + Icges(7,2) * t63 + Icges(7,6) * t170;
t33 = Icges(7,5) * t66 + Icges(7,6) * t65 + Icges(7,3) * t171;
t32 = Icges(7,5) * t64 + Icges(7,6) * t63 + Icges(7,3) * t170;
t31 = t154 * t268 + t156 * t267;
t28 = -t154 * t228 + (-t125 * t152 + t151 * t277) * qJD(5) + ((-t103 + t149) * t156 + (-t104 + t177 + t275) * t154) * qJD(1);
t23 = (t153 * t284 + t279) * t236 + ((t216 - t304) * t154 + t301) * qJD(1) + t207 - t210;
t22 = (-t156 * t224 + t209) * qJD(1) + t192 - t227 + t231;
t20 = (-t108 * t238 + t39) * t153 + (qJD(5) * t80 + t108 * t241 + t154 * t86) * t155;
t19 = (-t108 * t236 - t38) * t153 + (-qJD(5) * t81 - t108 * t242 + t156 * t86) * t155;
t18 = t105 * t171 + t106 * t65 + t107 * t66 + t113 * t84 + t114 * t85 - t251 * t83;
t17 = t105 * t170 + t106 * t63 + t107 * t64 + t115 * t84 + t116 * t85 + t248 * t83;
t14 = t195 * t239 + (-t154 * t38 - t156 * t39 + (t154 * t80 - t156 * t81) * qJD(1)) * t155;
t13 = (-qJD(1) * t267 + t227 - t39) * t154 + (t38 - t126 * t236 + (t136 + t268) * qJD(1)) * t156;
t12 = t42 * t153 - t155 * t198;
t11 = t41 * t153 - t155 * t199;
t10 = (qJD(5) * t193 + t32) * t153 + (qJD(5) * t75 - t163 * t34 + t165 * t36 + (-t163 * t79 - t165 * t77) * qJD(6)) * t155;
t9 = (qJD(5) * t194 + t33) * t153 + (qJD(5) * t74 - t163 * t35 + t165 * t37 + (-t163 * t78 - t165 * t76) * qJD(6)) * t155;
t8 = -t75 * t222 + t113 * t34 + t114 * t36 + t65 * t77 + t66 * t79 + (-t155 * t32 + t239 * t75) * t154;
t7 = -t74 * t222 + t113 * t35 + t114 * t37 + t65 * t76 + t66 * t78 + (-t155 * t33 + t239 * t74) * t154;
t6 = -t75 * t220 + t115 * t34 + t116 * t36 + t63 * t77 + t64 * t79 + (t156 * t32 - t242 * t75) * t155;
t5 = -t74 * t220 + t115 * t35 + t116 * t37 + t63 * t76 + t64 * t78 + (t156 * t33 - t242 * t74) * t155;
t4 = -qJD(1) * t199 + t154 * t8 + t156 * t7;
t3 = -qJD(1) * t198 + t154 * t6 + t156 * t5;
t2 = (qJD(5) * t199 + t18) * t153 + (-qJD(1) * t15 + qJD(5) * t41 - t154 * t7 + t156 * t8) * t155;
t1 = (qJD(5) * t198 + t17) * t153 + (-qJD(1) * t16 + qJD(5) * t42 - t154 * t5 + t156 * t6) * t155;
t21 = [(t22 * t50 + t23 * t49) * t291 + (t51 * t60 + t52 * t59) * t292 - t153 * t174 - t189 * t237 + t186 * t239 + 0.2e1 * m(4) * (t91 * t96 + t92 * t95) + 0.2e1 * m(5) * (t57 * t90 + t58 * t89) + t169 + t303 * t233 + (-t270 - t173) * t155; 0; 0; m(7) * (t154 * t23 - t156 * t22 + (t154 * t50 + t156 * t49) * qJD(1)) + m(6) * ((t154 * t60 + t156 * t59) * qJD(1) + t299) + m(4) * (t154 * t92 - t156 * t91 + (t154 * t96 + t156 * t95) * qJD(1)) + m(5) * (t154 * t58 - t156 * t57 + (t154 * t90 + t156 * t89) * qJD(1)); 0; 0; m(7) * (t154 * t22 + t156 * t23 + (-t154 * t49 + t156 * t50) * qJD(1)) + m(6) * (t300 * qJD(1) + t154 * t51 + t156 * t52) + m(5) * (t154 * t57 + t156 * t58 + (-t154 * t89 + t156 * t90) * qJD(1)); 0; 0; 0; m(7) * (-t22 * t88 + t23 * t87 + t43 * t50 + t44 * t49) + m(6) * (t300 * t120 + t299 * t125) - (t151 / 0.2e1 + t152 / 0.2e1) * t183 * qJD(5) + ((-t256 / 0.2e1 + t276 / 0.2e1 + t60 * t283 - t230) * t154 + (t257 / 0.2e1 - t255 / 0.2e1 + t59 * t283 - t229) * t156) * qJD(1) + (-qJD(5) * t181 - t153 * (qJD(1) * t99 - t187 * t236) + t155 * (qJD(1) * t101 - t190 * t236) + t10 + t17) * t154 / 0.2e1 + (-qJD(5) * t191 - t153 * (t296 * qJD(1) + t154 * t173) + t155 * (t295 * qJD(1) + t154 * t174) + t18 + t9) * t287; m(6) * t28 + m(7) * t13; m(7) * (t44 * t154 - t43 * t156 + (-t154 * t88 + t156 * t87) * qJD(1)) - m(6) * t215; m(7) * (t154 * t43 + t156 * t44 + (-t154 * t87 - t156 * t88) * qJD(1)); t156 * ((t156 * t68 + (t46 + t302) * qJD(1)) * t156 + (-t45 * qJD(1) + (-t237 * t295 + t239 * t296) * t154 + (t67 + (t256 - t276) * qJD(5) + (t181 - t97) * qJD(1)) * t156) * t154) + t154 * ((t154 * t67 + (-t47 + t175) * qJD(1)) * t154 + (t48 * qJD(1) + (-t101 * t237 + t239 * t99 + t258) * t156 + (t68 + (t255 - t257) * qJD(5) + t191 * qJD(1)) * t154) * t156) + ((-t103 * t154 + t104 * t156) * t28 - t125 * t215) * t292 + t156 * t4 + t154 * t3 + (t13 * t31 - t43 * t88 + t44 * t87) * t291 + (-t154 * t46 - t156 * t45 - t15) * t242 + (t154 * t48 + t156 * t47 + t16) * t241; m(7) * (t19 * t49 + t20 * t50 + t22 * t54 + t23 * t55) + (t154 * t230 + t156 * t229) * t239 + ((t17 / 0.2e1 + t10 / 0.2e1) * t156 + (-t18 / 0.2e1 - t9 / 0.2e1) * t154 + (t154 * t229 - t156 * t230) * qJD(1)) * t155 + t278; m(7) * t14; m(7) * (t154 * t19 - t156 * t20 + (t154 * t54 + t156 * t55) * qJD(1)); m(7) * (t154 * t20 + t156 * t19 + (-t154 * t55 + t156 * t54) * qJD(1)); m(7) * (-t13 * t40 + t14 * t31 + t19 * t87 - t20 * t88 + t43 * t54 + t44 * t55) + (t2 / 0.2e1 + qJD(1) * t12 / 0.2e1 - t16 * t239 / 0.2e1 + (qJD(1) * t30 + t9) * t290) * t156 + (-qJD(1) * t11 / 0.2e1 + t15 * t239 / 0.2e1 + t1 / 0.2e1 + (-qJD(1) * t29 + t10) * t290) * t154 + (t4 * t289 + t3 * t287 + qJD(5) * t197 / 0.2e1 + (-t156 * t15 / 0.2e1 + t16 * t289) * qJD(1)) * t155; (-t14 * t40 + t19 * t55 + t20 * t54) * t291 + ((t154 * t11 - t156 * t12 + t153 * t196) * qJD(5) + t278) * t153 + (-t154 * t2 + t156 * t1 + t153 * (t10 * t156 - t154 * t9) + (t56 * t153 - t155 * t196) * qJD(5) + (-t156 * t11 - t154 * t12 - t153 * t197) * qJD(1)) * t155;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
