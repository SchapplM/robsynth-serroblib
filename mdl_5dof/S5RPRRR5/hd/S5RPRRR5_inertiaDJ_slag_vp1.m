% Calculate time derivative of joint inertia matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:53:57
% DurationCPUTime: 4.46s
% Computational Cost: add. (9111->349), mult. (5758->499), div. (0->0), fcn. (4324->10), ass. (0->204)
t165 = qJ(1) + pkin(9);
t232 = pkin(2) * cos(t165) + cos(qJ(1)) * pkin(1);
t158 = qJ(3) + t165;
t152 = sin(t158);
t167 = sin(qJ(4));
t169 = cos(qJ(4));
t153 = cos(t158);
t262 = Icges(5,4) * t169;
t197 = -Icges(5,2) * t167 + t262;
t86 = -Icges(5,6) * t152 - t153 * t197;
t263 = Icges(5,4) * t167;
t199 = Icges(5,1) * t169 - t263;
t88 = -Icges(5,5) * t152 - t153 * t199;
t200 = t167 * t86 - t169 * t88;
t287 = t152 * t200;
t166 = qJ(4) + qJ(5);
t159 = sin(t166);
t160 = cos(t166);
t260 = Icges(6,4) * t160;
t196 = -Icges(6,2) * t159 + t260;
t76 = -Icges(6,6) * t152 - t153 * t196;
t261 = Icges(6,4) * t159;
t198 = Icges(6,1) * t160 - t261;
t78 = -Icges(6,5) * t152 - t153 * t198;
t206 = t159 * t76 - t160 * t78;
t286 = t152 * t206;
t102 = t152 * rSges(4,1) + t153 * rSges(4,2);
t233 = pkin(2) * sin(t165) + sin(qJ(1)) * pkin(1);
t164 = qJD(1) + qJD(3);
t114 = Icges(6,1) * t159 + t260;
t163 = qJD(4) + qJD(5);
t251 = t114 * t163;
t284 = -Icges(6,5) * t164 + t251;
t113 = Icges(6,2) * t160 + t261;
t252 = t113 * t163;
t283 = -Icges(6,6) * t164 + t252;
t112 = Icges(6,5) * t159 + Icges(6,6) * t160;
t282 = -Icges(6,3) * t164 + t112 * t163;
t134 = Icges(5,5) * t167 + Icges(5,6) * t169;
t281 = -Icges(5,3) * t164 + qJD(4) * t134;
t135 = Icges(5,2) * t169 + t263;
t280 = -Icges(5,6) * t164 + qJD(4) * t135;
t136 = Icges(5,1) * t167 + t262;
t279 = -Icges(5,5) * t164 + qJD(4) * t136;
t117 = t197 * qJD(4);
t118 = t199 * qJD(4);
t278 = qJD(4) * (t135 * t169 + t136 * t167) + t117 * t167 - t118 * t169 - t134 * t164;
t100 = t198 * t163;
t99 = t196 * t163;
t277 = t159 * (t99 + t251) + t160 * (-t100 + t252) - t112 * t164;
t276 = 2 * m(4);
t275 = 2 * m(5);
t274 = 2 * m(6);
t273 = -t152 / 0.2e1;
t272 = -t153 / 0.2e1;
t267 = rSges(5,2) * t167;
t269 = rSges(5,1) * t169;
t123 = (-t267 + t269) * qJD(4);
t271 = m(5) * t123;
t137 = rSges(5,1) * t167 + rSges(5,2) * t169;
t270 = m(5) * t137;
t268 = rSges(6,1) * t160;
t266 = rSges(6,2) * t159;
t194 = Icges(6,5) * t160 - Icges(6,6) * t159;
t73 = -Icges(6,3) * t153 + t152 * t194;
t265 = t164 * t73;
t74 = -Icges(6,3) * t152 - t153 * t194;
t264 = t164 * t74;
t101 = (-t266 + t268) * t163;
t250 = t152 * t101;
t249 = t152 * t164;
t248 = t152 * t169;
t171 = -pkin(8) - pkin(7);
t247 = t152 * t171;
t246 = t153 * t164;
t245 = t159 * t163;
t244 = t160 * t163;
t243 = t164 * t167;
t228 = t152 * t266;
t242 = rSges(6,3) * t246 + t164 * t228;
t122 = t153 * t268;
t241 = rSges(6,3) * t249 + t164 * t122;
t227 = t152 * t243;
t240 = rSges(5,2) * t227 + rSges(5,3) * t246;
t133 = t153 * t269;
t239 = rSges(5,3) * t249 + t164 * t133;
t154 = pkin(4) * t169 + pkin(3);
t238 = t152 * t154 + t153 * t171;
t237 = -pkin(3) * t246 - pkin(7) * t249;
t236 = t232 * qJD(1);
t235 = t153 * pkin(3) + t152 * pkin(7);
t234 = t152 ^ 2 + t153 ^ 2;
t231 = qJD(4) * t167;
t230 = qJD(4) * t169;
t226 = t159 * t246;
t79 = -rSges(6,3) * t153 + t152 * t268 - t228;
t80 = -t152 * rSges(6,3) + t153 * t266 - t122;
t229 = t152 * (-t152 * rSges(6,1) * t245 + (-t152 * t244 - t226) * rSges(6,2) + t241) + t80 * t249 + t79 * t246;
t115 = rSges(6,1) * t159 + rSges(6,2) * t160;
t225 = t115 * t249;
t75 = -Icges(6,6) * t153 + t152 * t196;
t77 = -Icges(6,5) * t153 + t152 * t198;
t207 = t159 * t75 - t160 * t77;
t17 = -t152 * t207 - t153 * t73;
t18 = -t153 * t74 - t286;
t188 = t207 * t153;
t19 = -t152 * t73 + t188;
t20 = -t152 * t74 + t153 * t206;
t183 = t198 * t164;
t47 = t152 * t183 + t284 * t153;
t215 = -t163 * t76 + t47;
t181 = t196 * t164;
t45 = t152 * t181 + t283 * t153;
t217 = t163 * t78 + t45;
t179 = t194 * t164;
t43 = t152 * t179 + t282 * t153;
t44 = -t282 * t152 + t153 * t179;
t46 = -t283 * t152 + t153 * t181;
t48 = -t284 * t152 + t153 * t183;
t224 = -t152 * ((t152 * t43 + (t19 + t286) * t164) * t152 + (-t20 * t164 + (-t159 * t46 + t160 * t48 - t244 * t75 - t245 * t77 + t265) * t153 + (t264 + t44 + (-t164 * t77 + t215) * t160 + (t164 * t75 - t217) * t159) * t152) * t153) + (-t152 * t18 - t153 * t17) * t249;
t214 = -t163 * t75 + t48;
t216 = t163 * t77 + t46;
t2 = (t153 * t44 + (-t18 + t188) * t164) * t153 + (t17 * t164 + (t159 * t45 - t160 * t47 + t244 * t76 + t245 * t78 - t264) * t152 + (-t265 + t43 + (-t164 * t78 - t214) * t160 + (t164 * t76 + t216) * t159) * t153) * t152;
t5 = -t152 * t20 - t153 * t19;
t223 = -t164 * t5 - t2;
t222 = t152 * t231;
t221 = t153 * t230;
t220 = t249 / 0.2e1;
t219 = -t246 / 0.2e1;
t218 = pkin(4) * t167 + t115;
t103 = t153 * rSges(4,1) - rSges(4,2) * t152;
t211 = -t153 * t154 + t247;
t96 = rSges(4,1) * t246 - rSges(4,2) * t249;
t210 = rSges(5,1) * t248 - t152 * t267;
t85 = -Icges(5,6) * t153 + t152 * t197;
t87 = -Icges(5,5) * t153 + t152 * t199;
t203 = t167 * t87 + t169 * t85;
t202 = t167 * t85 - t169 * t87;
t201 = t167 * t88 + t169 * t86;
t195 = Icges(5,5) * t169 - Icges(5,6) * t167;
t193 = t113 * t159 - t114 * t160;
t191 = t135 * t167 - t136 * t169;
t90 = -t152 * rSges(5,3) + t153 * t267 - t133;
t190 = -pkin(4) * t231 - t164 * t171;
t189 = t233 * qJD(1);
t95 = t102 * t164;
t187 = t202 * t153;
t176 = -t194 * t163 - t164 * t193;
t186 = (t176 * t152 + t277 * t153 + t159 * t215 + t160 * t217) * t273 + (-t277 * t152 + t176 * t153 + t159 * t214 + t160 * t216) * t272 + (-t153 * t112 - t152 * t193 + t159 * t77 + t160 * t75) * t220 + (-t152 * t112 + t153 * t193 + t159 * t78 + t160 * t76) * t219;
t185 = qJD(4) * t137;
t184 = t199 * t164;
t182 = t197 * t164;
t180 = t195 * t164;
t178 = -t152 * t230 - t153 * t243;
t65 = -t90 + t235;
t175 = -t195 * qJD(4) - t164 * t191;
t62 = t79 + t238;
t63 = -t211 - t80;
t146 = t152 * pkin(3);
t64 = t146 + (-rSges(5,3) - pkin(7)) * t153 + t210;
t174 = t159 * t100 - t113 * t245 + t114 * t244 + t169 * t117 + t167 * t118 - t135 * t231 + t136 * t230 + t160 * t99;
t173 = -t115 * t163 + t190;
t172 = t186 + (-qJD(4) * t200 + t175 * t152 + t278 * t153 + t167 * (t152 * t184 + t279 * t153) + t169 * (t152 * t182 + t280 * t153)) * t273 + (-qJD(4) * t202 - t278 * t152 + t175 * t153 + t167 * (-t279 * t152 + t153 * t184) + t169 * (-t280 * t152 + t153 * t182)) * t272 + (-t153 * t134 - t152 * t191 + t203) * t220 + (-t152 * t134 + t153 * t191 + t201) * t219;
t34 = -rSges(5,1) * t222 + rSges(5,2) * t178 - t237 + t239;
t130 = pkin(7) * t246;
t33 = -rSges(5,2) * t221 - pkin(3) * t249 + t130 + (-t153 * t231 - t164 * t248) * rSges(5,1) + t240;
t104 = t154 * t246;
t28 = -rSges(6,2) * t226 + t152 * t173 + t104 + t241;
t27 = (-t154 - t268) * t249 + t173 * t153 + t242;
t92 = t103 + t232;
t91 = t233 + t102;
t89 = -rSges(5,3) * t153 + t210;
t84 = -Icges(5,3) * t152 - t153 * t195;
t83 = -Icges(5,3) * t153 + t152 * t195;
t82 = t218 * t153;
t81 = t218 * t152;
t72 = t96 + t236;
t71 = -t95 - t189;
t70 = t211 + t235;
t69 = pkin(7) * t153 - t146 + t238;
t68 = t152 * t79;
t61 = t65 + t232;
t60 = t64 + t233;
t55 = -t281 * t152 + t153 * t180;
t54 = t152 * t180 + t281 * t153;
t53 = t63 + t232;
t52 = t62 + t233;
t49 = t153 * rSges(6,2) * t244 + (t153 * t245 + t160 * t249) * rSges(6,1) - t242;
t42 = pkin(4) * t178 - t115 * t246 - t250;
t41 = -t225 + t153 * t101 + (t221 - t227) * pkin(4);
t38 = -t153 * t80 + t68;
t30 = t34 + t236;
t29 = -t189 + t33;
t26 = -t152 * t84 + t200 * t153;
t25 = -t152 * t83 + t187;
t24 = -t153 * t84 - t287;
t23 = -t202 * t152 - t153 * t83;
t22 = t28 + t236;
t21 = -t189 + t27;
t16 = t152 * t69 + t68 + (-t70 - t80) * t153;
t11 = (-t153 * t185 + t164 * t89 + t240) * t153 + (-t152 * t185 + (t90 + (-t267 - t269) * t153) * t164 + t239) * t152;
t8 = -t153 * t49 + t229;
t3 = (t190 * t153 + t164 * t69 - t130 - t49) * t153 + (-pkin(4) * t222 + t104 + (t70 - t153 * (-pkin(3) + t154) - t247) * t164 + t237) * t152 + t229;
t1 = [(t71 * t92 + t72 * t91) * t276 + (t29 * t61 + t30 * t60) * t275 + (t21 * t53 + t22 * t52) * t274 + t174; 0; 0; m(4) * (t102 * t72 + t103 * t71 + t91 * t96 - t92 * t95) + m(5) * (t29 * t65 + t30 * t64 + t33 * t61 + t34 * t60) + m(6) * (t21 * t63 + t22 * t62 + t27 * t53 + t28 * t52) + t174; 0; (t27 * t63 + t28 * t62) * t274 + (t33 * t65 + t34 * t64) * t275 + (t102 * t96 - t103 * t95) * t276 + t174; t172 + m(6) * (-t21 * t81 + t22 * t82 + t41 * t52 + t42 * t53) + ((-t164 * t61 + t30) * t153 + (-t164 * t60 - t29) * t152) * t270 + (-t152 * t61 + t153 * t60) * t271; m(5) * t11 + m(6) * t3; t172 + m(6) * (-t27 * t81 + t28 * t82 + t41 * t62 + t42 * t63) + ((-t164 * t65 + t34) * t153 + (-t164 * t64 - t33) * t152) * t270 + (-t152 * t65 + t153 * t64) * t271; ((t152 * t89 - t153 * t90) * t11 + t234 * t137 * t123) * t275 + (-t152 * t24 - t153 * t23) * t249 - t153 * ((t153 * t55 + (-t24 + t187) * t164) * t153 + (t23 * t164 + (t230 * t86 + t231 * t88) * t152 + (t203 * qJD(4) + t164 * t200 + t54) * t153) * t152) - t152 * ((t152 * t54 + (t25 + t287) * t164) * t152 + (-t26 * t164 + (-t230 * t85 - t231 * t87) * t153 + (-t201 * qJD(4) + t202 * t164 + t55) * t152) * t153) + (t16 * t3 + t41 * t82 - t42 * t81) * t274 - t153 * t2 + t224 + (t152 * t26 + t153 * t25 - t5) * t246; m(6) * ((-t152 * t53 + t153 * t52) * t101 + ((-t164 * t53 + t22) * t153 + (-t164 * t52 - t21) * t152) * t115) + t186; m(6) * t8; m(6) * ((-t152 * t63 + t153 * t62) * t101 + ((-t164 * t63 + t28) * t153 + (-t164 * t62 - t27) * t152) * t115) + t186; m(6) * (-t115 * t152 * t42 + t16 * t8 - t82 * t225 + t81 * t250 + t3 * t38) + (m(6) * (t101 * t82 + (t164 * t81 + t41) * t115) + t223) * t153 + t224; (t101 * t115 * t234 + t38 * t8) * t274 + t223 * t153 + t224;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
