% Calculate time derivative of joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:50
% EndTime: 2019-12-05 18:31:56
% DurationCPUTime: 3.19s
% Computational Cost: add. (9289->356), mult. (5874->504), div. (0->0), fcn. (4392->10), ass. (0->226)
t174 = qJ(1) + qJ(2);
t166 = pkin(9) + t174;
t162 = sin(t166);
t173 = qJ(4) + qJ(5);
t169 = cos(t173);
t279 = rSges(6,1) * t169;
t163 = cos(t166);
t167 = sin(t173);
t277 = rSges(6,2) * t167;
t306 = t163 * rSges(6,3) + t162 * t277;
t77 = -t162 * t279 + t306;
t177 = cos(qJ(4));
t244 = qJD(4) * t177;
t172 = qJD(1) + qJD(2);
t175 = sin(qJ(4));
t251 = t172 * t175;
t310 = t162 * t244 + t163 * t251;
t171 = qJD(4) + qJD(5);
t253 = t169 * t171;
t255 = t167 * t171;
t309 = rSges(6,1) * t255 + rSges(6,2) * t253;
t270 = Icges(5,4) * t177;
t205 = -Icges(5,2) * t175 + t270;
t85 = Icges(5,6) * t163 - t205 * t162;
t271 = Icges(5,4) * t175;
t207 = Icges(5,1) * t177 - t271;
t87 = Icges(5,5) * t163 - t207 * t162;
t210 = t175 * t85 - t177 * t87;
t308 = t210 * t163;
t268 = Icges(6,4) * t169;
t204 = -Icges(6,2) * t167 + t268;
t73 = Icges(6,6) * t163 - t204 * t162;
t269 = Icges(6,4) * t167;
t206 = Icges(6,1) * t169 - t269;
t75 = Icges(6,5) * t163 - t206 * t162;
t215 = t167 * t73 - t169 * t75;
t307 = t215 * t163;
t122 = Icges(6,1) * t167 + t268;
t259 = t122 * t171;
t305 = -Icges(6,5) * t172 + t259;
t121 = Icges(6,2) * t169 + t269;
t260 = t121 * t171;
t304 = -Icges(6,6) * t172 + t260;
t120 = Icges(6,5) * t167 + Icges(6,6) * t169;
t303 = -Icges(6,3) * t172 + t120 * t171;
t302 = -t309 * t163 + t77 * t172;
t146 = Icges(5,5) * t175 + Icges(5,6) * t177;
t301 = -Icges(5,3) * t172 + t146 * qJD(4);
t147 = Icges(5,2) * t177 + t271;
t300 = -Icges(5,6) * t172 + t147 * qJD(4);
t148 = Icges(5,1) * t175 + t270;
t299 = -Icges(5,5) * t172 + t148 * qJD(4);
t129 = t205 * qJD(4);
t130 = t207 * qJD(4);
t298 = t129 * t175 - t130 * t177 - t146 * t172 + (t147 * t177 + t148 * t175) * qJD(4);
t101 = t204 * t171;
t102 = t206 * t171;
t297 = (t101 + t259) * t167 + (-t102 + t260) * t169 - t120 * t172;
t296 = 2 * m(3);
t295 = 2 * m(4);
t294 = 2 * m(5);
t293 = 2 * m(6);
t292 = t162 / 0.2e1;
t291 = t163 / 0.2e1;
t290 = -rSges(5,3) - pkin(7);
t278 = rSges(5,2) * t175;
t280 = rSges(5,1) * t177;
t133 = (-t278 + t280) * qJD(4);
t289 = m(5) * t133;
t153 = rSges(5,1) * t175 + rSges(5,2) * t177;
t288 = m(5) * t153;
t168 = sin(t174);
t287 = pkin(2) * t168;
t170 = cos(t174);
t286 = pkin(2) * t170;
t285 = t162 * pkin(7);
t158 = t163 * pkin(7);
t176 = sin(qJ(1));
t284 = t176 * pkin(1);
t178 = cos(qJ(1));
t283 = t178 * pkin(1);
t164 = pkin(4) * t177 + pkin(3);
t282 = pkin(3) - t164;
t225 = t162 * t282;
t179 = -pkin(8) - pkin(7);
t256 = t163 * t179;
t281 = t158 - t225 + t256 - t77;
t276 = pkin(1) * qJD(1);
t275 = t162 * rSges(5,3);
t274 = t162 * rSges(6,3);
t157 = t163 * rSges(5,3);
t258 = t162 * t172;
t257 = t163 * t172;
t254 = t168 * t172;
t252 = t170 * t172;
t250 = t172 * t179;
t245 = qJD(4) * t175;
t232 = t162 * t245;
t249 = -pkin(4) * t232 - t162 * t250;
t144 = t162 * t278;
t247 = t144 + t157;
t104 = rSges(3,1) * t254 + rSges(3,2) * t252;
t246 = t162 ^ 2 + t163 ^ 2;
t243 = t162 * t280;
t240 = t178 * t276;
t132 = t163 * t277;
t230 = t163 * t245;
t237 = -pkin(4) * t230 - t163 * t250 - t164 * t258;
t236 = t172 * t132 + t309 * t162;
t229 = t163 * t244;
t234 = -rSges(5,1) * t230 - rSges(5,2) * t229 - t172 * t243;
t233 = rSges(5,1) * t232 + t310 * rSges(5,2);
t152 = pkin(2) * t254;
t79 = rSges(4,1) * t258 + rSges(4,2) * t257 + t152;
t228 = -t258 / 0.2e1;
t227 = t257 / 0.2e1;
t226 = -pkin(3) - t280;
t123 = rSges(6,1) * t167 + rSges(6,2) * t169;
t224 = pkin(4) * t175 + t123;
t189 = t204 * t172;
t45 = -t162 * t189 - t304 * t163;
t76 = Icges(6,5) * t162 + t206 * t163;
t223 = t171 * t76 + t45;
t46 = t304 * t162 - t163 * t189;
t222 = t171 * t75 + t46;
t191 = t206 * t172;
t47 = -t162 * t191 - t305 * t163;
t74 = Icges(6,6) * t162 + t204 * t163;
t221 = -t171 * t74 + t47;
t48 = t305 * t162 - t163 * t191;
t220 = -t171 * t73 + t48;
t125 = -rSges(3,1) * t170 + t168 * rSges(3,2);
t219 = -t164 - t279;
t105 = -rSges(3,1) * t252 + rSges(3,2) * t254;
t216 = -rSges(4,1) * t163 - t286;
t124 = -rSges(3,1) * t168 - rSges(3,2) * t170;
t214 = t167 * t74 - t169 * t76;
t211 = t175 * t87 + t177 * t85;
t86 = Icges(5,6) * t162 + t205 * t163;
t88 = Icges(5,5) * t162 + t207 * t163;
t209 = t175 * t88 + t177 * t86;
t208 = t175 * t86 - t177 * t88;
t203 = Icges(5,5) * t177 - Icges(5,6) * t175;
t202 = Icges(6,5) * t169 - Icges(6,6) * t167;
t201 = t121 * t167 - t122 * t169;
t199 = t147 * t175 - t148 * t177;
t187 = t202 * t162;
t71 = Icges(6,3) * t163 - t187;
t17 = t215 * t162 + t163 * t71;
t195 = t214 * t162;
t72 = Icges(6,3) * t162 + t202 * t163;
t18 = t163 * t72 + t195;
t19 = t162 * t71 - t307;
t20 = t162 * t72 - t214 * t163;
t43 = -t303 * t163 - t172 * t187;
t44 = t303 * t162 - t202 * t257;
t198 = -(t162 * t18 + t163 * t17) * t258 + t162 * ((t162 * t43 + (-t19 + t195) * t172) * t162 + (t20 * t172 + (-t167 * t46 + t169 * t48 - t73 * t253 - t75 * t255) * t163 + (t44 + (-t172 * t75 + t221) * t169 + (t172 * t73 - t223) * t167) * t162) * t163) + t163 * ((t163 * t44 + (t18 + t307) * t172) * t163 + (-t17 * t172 + (t167 * t45 - t169 * t47 + t74 * t253 + t76 * t255) * t162 + (t43 + (-t172 * t76 - t220) * t169 + (t172 * t74 + t222) * t167) * t163) * t162) + (t162 * t20 + t163 * t19) * t257;
t98 = t162 * rSges(4,2) + t216;
t197 = -t163 * t279 - t274;
t194 = t208 * t162;
t185 = t202 * t171 + t201 * t172;
t193 = (t185 * t162 - t297 * t163 + t221 * t167 + t223 * t169) * t292 + (t297 * t162 + t185 * t163 + t220 * t167 + t222 * t169) * t291 + (t163 * t120 + t201 * t162 + t167 * t75 + t169 * t73) * t228 + (t162 * t120 - t201 * t163 + t167 * t76 + t169 * t74) * t227;
t192 = t207 * t172;
t190 = t205 * t172;
t188 = t203 * t172;
t97 = -rSges(4,1) * t162 - rSges(4,2) * t163 - t287;
t80 = rSges(4,2) * t258 + t216 * t172;
t184 = t203 * qJD(4) + t199 * t172;
t183 = t219 * t163 - t274 - t286;
t64 = t226 * t162 + t158 + t247 - t287;
t182 = t290 * t162 + t226 * t163 - t286;
t145 = t163 * t278;
t65 = t145 + t182;
t154 = t162 * t179;
t61 = t132 + t154 + t183;
t60 = t219 * t162 - t256 - t287 + t306;
t181 = t169 * t101 + t167 * t102 - t121 * t255 + t122 * t253 + t177 * t129 + t175 * t130 - t147 * t245 + t148 * t244;
t180 = t193 + (-t208 * qJD(4) + t184 * t162 - t298 * t163 + t175 * (-t162 * t192 - t299 * t163) + t177 * (-t162 * t190 - t300 * t163)) * t292 + (-t210 * qJD(4) + t298 * t162 + t184 * t163 + t175 * (t299 * t162 - t163 * t192) + t177 * (t300 * t162 - t163 * t190)) * t291 + (t163 * t146 + t199 * t162 + t211) * t228 + (t162 * t146 - t199 * t163 + t209) * t227;
t143 = pkin(3) * t258;
t31 = t143 + t152 + (t290 * t163 - t144) * t172 - t234;
t23 = t152 - t237 - t302;
t32 = t182 * t172 + t233;
t24 = t183 * t172 + t236 - t249;
t165 = t176 * t276;
t107 = t125 - t283;
t106 = t124 - t284;
t103 = (-t277 + t279) * t171;
t94 = t105 - t240;
t93 = t165 + t104;
t92 = t98 - t283;
t91 = t97 - t284;
t90 = t163 * t280 - t145 + t275;
t89 = -t243 + t247;
t84 = Icges(5,3) * t162 + t203 * t163;
t83 = Icges(5,3) * t163 - t203 * t162;
t82 = t224 * t163;
t81 = t224 * t162;
t78 = -t132 - t197;
t70 = -t282 * t163 - t154 - t285;
t68 = t80 - t240;
t67 = t165 + t79;
t66 = t163 * t78;
t63 = t65 - t283;
t62 = t64 - t284;
t55 = t301 * t162 - t163 * t188;
t54 = -t162 * t188 - t301 * t163;
t53 = t61 - t283;
t52 = t60 - t284;
t49 = t197 * t172 + t236;
t42 = t310 * pkin(4) + t162 * t103 + t123 * t257;
t41 = t123 * t258 - t163 * t103 + (t162 * t251 - t229) * pkin(4);
t38 = -t162 * t77 + t66;
t37 = t163 * t302;
t30 = t32 - t240;
t29 = t165 + t31;
t28 = t162 * t84 - t208 * t163;
t27 = t162 * t83 - t308;
t26 = t163 * t84 + t194;
t25 = t210 * t162 + t163 * t83;
t22 = t24 - t240;
t21 = t165 + t23;
t16 = t281 * t162 + t163 * t70 + t66;
t11 = t163 * t234 - t162 * t233 + ((-t89 + t157) * t163 + (t275 - t90 + (t278 + t280) * t163) * t162) * t172;
t8 = -t77 * t257 + t37 + (-t172 * t78 - t49) * t162;
t3 = t163 * (t143 + t237) + t37 + (-t49 + t249) * t162 + ((-t70 - t78 - t285) * t162 + (-t225 + t281 - t158) * t163) * t172;
t1 = [(t106 * t94 + t107 * t93) * t296 + (t67 * t92 + t68 * t91) * t295 + (t29 * t63 + t30 * t62) * t294 + (t21 * t53 + t22 * t52) * t293 + t181; m(3) * (t104 * t107 + t105 * t106 + t124 * t94 + t125 * t93) + m(4) * (t67 * t98 + t68 * t97 + t79 * t92 + t80 * t91) + m(5) * (t29 * t65 + t30 * t64 + t31 * t63 + t32 * t62) + m(6) * (t21 * t61 + t22 * t60 + t23 * t53 + t24 * t52) + t181; (t23 * t61 + t24 * t60) * t293 + (t31 * t65 + t32 * t64) * t294 + (t79 * t98 + t80 * t97) * t295 + (t104 * t125 + t105 * t124) * t296 + t181; 0; 0; 0; m(6) * (t21 * t81 - t22 * t82 + t41 * t52 + t42 * t53) + ((t172 * t63 - t30) * t163 + (t172 * t62 + t29) * t162) * t288 + (t162 * t63 - t163 * t62) * t289 + t180; m(6) * (t23 * t81 - t24 * t82 + t41 * t60 + t42 * t61) + ((t172 * t65 - t32) * t163 + (t172 * t64 + t31) * t162) * t288 + (t162 * t65 - t163 * t64) * t289 + t180; m(5) * t11 + m(6) * t3; ((-t162 * t89 + t163 * t90) * t11 + t246 * t153 * t133) * t294 - (t162 * t26 + t163 * t25) * t258 + t163 * ((t163 * t55 + (t26 + t308) * t172) * t163 + (-t25 * t172 + (t86 * t244 + t88 * t245) * t162 + (t211 * qJD(4) + t208 * t172 + t54) * t163) * t162) + (t162 * t28 + t163 * t27) * t257 + t162 * ((t162 * t54 + (-t27 + t194) * t172) * t162 + (t28 * t172 + (-t85 * t244 - t87 * t245) * t163 + (-t209 * qJD(4) + t210 * t172 + t55) * t162) * t163) + (t16 * t3 - t41 * t82 + t42 * t81) * t293 + t198; m(6) * ((t162 * t53 - t163 * t52) * t103 + ((t172 * t53 - t22) * t163 + (t172 * t52 + t21) * t162) * t123) + t193; m(6) * ((t162 * t61 - t163 * t60) * t103 + ((t172 * t61 - t24) * t163 + (t172 * t60 + t23) * t162) * t123) + t193; m(6) * t8; m(6) * (t16 * t8 + t3 * t38 + (t162 * t81 + t163 * t82) * t103 + ((t172 * t81 - t41) * t163 + (-t172 * t82 + t42) * t162) * t123) + t198; (t246 * t123 * t103 + t38 * t8) * t293 + t198;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
