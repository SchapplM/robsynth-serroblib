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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:01:48
% EndTime: 2020-01-03 12:01:55
% DurationCPUTime: 3.25s
% Computational Cost: add. (9289->367), mult. (5874->516), div. (0->0), fcn. (4392->10), ass. (0->222)
t179 = sin(qJ(4));
t274 = rSges(5,2) * t179;
t181 = cos(qJ(4));
t277 = rSges(5,1) * t181;
t296 = -t274 + t277;
t178 = qJ(1) + qJ(2);
t168 = pkin(9) + t178;
t162 = sin(t168);
t163 = cos(t168);
t268 = Icges(5,4) * t181;
t209 = -Icges(5,2) * t179 + t268;
t88 = -Icges(5,6) * t162 - t163 * t209;
t269 = Icges(5,4) * t179;
t211 = Icges(5,1) * t181 - t269;
t90 = -Icges(5,5) * t162 - t163 * t211;
t212 = t179 * t88 - t181 * t90;
t295 = t162 * t212;
t177 = qJ(4) + qJ(5);
t169 = sin(t177);
t171 = cos(t177);
t266 = Icges(6,4) * t171;
t208 = -Icges(6,2) * t169 + t266;
t76 = -Icges(6,6) * t162 - t163 * t208;
t267 = Icges(6,4) * t169;
t210 = Icges(6,1) * t171 - t267;
t78 = -Icges(6,5) * t162 - t163 * t210;
t218 = t169 * t76 - t171 * t78;
t294 = t162 * t218;
t170 = sin(t178);
t172 = cos(t178);
t122 = t170 * rSges(3,1) + t172 * rSges(3,2);
t176 = qJD(1) + qJD(2);
t120 = Icges(6,1) * t169 + t266;
t175 = qJD(4) + qJD(5);
t257 = t120 * t175;
t293 = -Icges(6,5) * t176 + t257;
t119 = Icges(6,2) * t171 + t267;
t258 = t119 * t175;
t292 = -Icges(6,6) * t176 + t258;
t164 = pkin(2) * t170;
t99 = t162 * rSges(4,1) + t163 * rSges(4,2) + t164;
t118 = Icges(6,5) * t169 + Icges(6,6) * t171;
t291 = -Icges(6,3) * t176 + t118 * t175;
t142 = Icges(5,5) * t179 + Icges(5,6) * t181;
t290 = -Icges(5,3) * t176 + qJD(4) * t142;
t143 = Icges(5,2) * t181 + t269;
t289 = -Icges(5,6) * t176 + qJD(4) * t143;
t144 = Icges(5,1) * t179 + t268;
t288 = -Icges(5,5) * t176 + qJD(4) * t144;
t125 = t209 * qJD(4);
t126 = t211 * qJD(4);
t287 = t125 * t179 - t126 * t181 - t142 * t176 + (t143 * t181 + t144 * t179) * qJD(4);
t103 = t208 * t175;
t104 = t210 * t175;
t286 = t169 * (t103 + t257) + t171 * (-t104 + t258) - t118 * t176;
t285 = 2 * m(3);
t284 = 2 * m(4);
t283 = 2 * m(5);
t282 = 2 * m(6);
t281 = -t162 / 0.2e1;
t280 = -t163 / 0.2e1;
t131 = t296 * qJD(4);
t279 = m(5) * t131;
t147 = rSges(5,1) * t179 + rSges(5,2) * t181;
t278 = m(5) * t147;
t276 = rSges(6,1) * t171;
t275 = rSges(3,2) * t170;
t273 = rSges(6,2) * t169;
t272 = pkin(1) * qJD(1);
t206 = Icges(6,5) * t171 - Icges(6,6) * t169;
t73 = -Icges(6,3) * t163 + t162 * t206;
t271 = t176 * t73;
t74 = -Icges(6,3) * t162 - t163 * t206;
t270 = t176 * t74;
t105 = (-t273 + t276) * t175;
t256 = t162 * t105;
t255 = t162 * t176;
t183 = -pkin(8) - pkin(7);
t254 = t162 * t183;
t253 = t163 * t176;
t252 = t169 * t175;
t251 = t171 * t175;
t250 = t172 * t176;
t249 = t176 * t179;
t238 = t162 * t273;
t248 = rSges(6,3) * t253 + t176 * t238;
t130 = t163 * t276;
t247 = rSges(6,3) * t255 + t176 * t130;
t141 = t163 * t277;
t246 = rSges(5,3) * t255 + t176 * t141;
t166 = pkin(4) * t181 + pkin(3);
t245 = t162 * t166 + t163 * t183;
t244 = -pkin(3) * t253 - pkin(7) * t255;
t243 = t163 * pkin(3) + t162 * pkin(7);
t242 = t162 ^ 2 + t163 ^ 2;
t241 = qJD(4) * t179;
t240 = qJD(4) * t181;
t235 = t169 * t253;
t79 = -rSges(6,3) * t163 + t162 * t276 - t238;
t80 = -t162 * rSges(6,3) + t163 * t273 - t130;
t239 = t162 * (-t162 * rSges(6,1) * t252 + (-t162 * t251 - t235) * rSges(6,2) + t247) + t80 * t255 + t79 * t253;
t180 = sin(qJ(1));
t237 = t180 * t272;
t236 = t162 * t249;
t121 = rSges(6,1) * t169 + rSges(6,2) * t171;
t234 = t121 * t255;
t75 = -Icges(6,6) * t163 + t162 * t208;
t77 = -Icges(6,5) * t163 + t162 * t210;
t219 = t169 * t75 - t171 * t77;
t17 = -t162 * t219 - t163 * t73;
t18 = -t163 * t74 - t294;
t201 = t219 * t163;
t19 = -t162 * t73 + t201;
t20 = -t162 * t74 + t163 * t218;
t196 = t210 * t176;
t47 = t162 * t196 + t293 * t163;
t225 = -t175 * t76 + t47;
t194 = t208 * t176;
t45 = t162 * t194 + t292 * t163;
t227 = t175 * t78 + t45;
t192 = t206 * t176;
t43 = t162 * t192 + t291 * t163;
t44 = -t162 * t291 + t163 * t192;
t46 = -t162 * t292 + t163 * t194;
t48 = -t162 * t293 + t163 * t196;
t233 = -t162 * ((t162 * t43 + (t19 + t294) * t176) * t162 + (-t20 * t176 + (-t169 * t46 + t171 * t48 - t251 * t75 - t252 * t77 + t271) * t163 + (t270 + t44 + (-t176 * t77 + t225) * t171 + (t176 * t75 - t227) * t169) * t162) * t163) + (-t162 * t18 - t163 * t17) * t255;
t224 = -t175 * t75 + t48;
t226 = t175 * t77 + t46;
t2 = (t163 * t44 + (-t18 + t201) * t176) * t163 + (t17 * t176 + (t169 * t45 - t171 * t47 + t251 * t76 + t252 * t78 - t270) * t162 + (-t271 + t43 + (-t176 * t78 - t224) * t171 + (t176 * t76 + t226) * t169) * t163) * t162;
t5 = -t162 * t20 - t163 * t19;
t232 = -t176 * t5 - t2;
t231 = t162 * t241;
t230 = t255 / 0.2e1;
t229 = -t253 / 0.2e1;
t228 = pkin(4) * t179 + t121;
t123 = t172 * rSges(3,1) - t275;
t221 = -t163 * t166 + t254;
t107 = rSges(3,1) * t250 - t176 * t275;
t220 = t296 * t162;
t165 = pkin(2) * t172;
t100 = t163 * rSges(4,1) - rSges(4,2) * t162 + t165;
t87 = -Icges(5,6) * t163 + t162 * t209;
t197 = t211 * t162;
t89 = -Icges(5,5) * t163 + t197;
t215 = t179 * t89 + t181 * t87;
t214 = t179 * t87 - t181 * t89;
t213 = t179 * t90 + t181 * t88;
t207 = Icges(5,5) * t181 - Icges(5,6) * t179;
t205 = t119 * t169 - t120 * t171;
t203 = t143 * t179 - t144 * t181;
t146 = pkin(2) * t250;
t82 = rSges(4,1) * t253 - rSges(4,2) * t255 + t146;
t92 = -t162 * rSges(5,3) + t163 * t274 - t141;
t202 = -pkin(4) * t241 - t176 * t183;
t106 = t122 * t176;
t200 = t214 * t163;
t189 = -t206 * t175 - t176 * t205;
t199 = (t189 * t162 + t286 * t163 + t169 * t225 + t171 * t227) * t281 + (-t162 * t286 + t189 * t163 + t169 * t224 + t171 * t226) * t280 + (-t163 * t118 - t162 * t205 + t169 * t77 + t171 * t75) * t230 + (-t162 * t118 + t163 * t205 + t169 * t78 + t171 * t76) * t229;
t198 = qJD(4) * t147;
t195 = t209 * t176;
t193 = t207 * t176;
t191 = -t162 * t240 - t163 * t249;
t188 = -t207 * qJD(4) - t176 * t203;
t65 = t165 - t92 + t243;
t81 = t99 * t176;
t60 = t164 + t79 + t245;
t61 = t165 - t221 - t80;
t187 = rSges(5,2) * t236 + rSges(5,3) * t253 - t163 * t198;
t155 = t162 * pkin(3);
t64 = t155 + t164 + (-rSges(5,3) - pkin(7)) * t163 + t220;
t186 = -t121 * t175 + t202;
t185 = t171 * t103 + t169 * t104 - t119 * t252 + t120 * t251 + t181 * t125 + t179 * t126 - t143 * t241 + t144 * t240;
t184 = t199 + (-qJD(4) * t212 + t188 * t162 + t287 * t163 + t179 * (t288 * t163 + t176 * t197) + t181 * (t162 * t195 + t289 * t163)) * t281 + (-qJD(4) * t214 - t162 * t287 + t188 * t163 + t179 * (-t162 * t288 + t211 * t253) + t181 * (-t162 * t289 + t163 * t195)) * t280 + (-t163 * t142 - t162 * t203 + t215) * t230 + (-t162 * t142 + t163 * t203 + t213) * t229;
t32 = -rSges(5,1) * t231 + rSges(5,2) * t191 + t146 - t244 + t246;
t138 = pkin(7) * t253;
t31 = t138 + (-t164 + (-pkin(3) - t277) * t162) * t176 + t187;
t110 = t166 * t253;
t24 = -rSges(6,2) * t235 + t162 * t186 + t110 + t146 + t247;
t23 = (-t164 + (-t166 - t276) * t162) * t176 + t186 * t163 + t248;
t182 = cos(qJ(1));
t174 = t182 * pkin(1);
t173 = t180 * pkin(1);
t167 = t182 * t272;
t109 = t123 + t174;
t108 = t173 + t122;
t96 = t107 + t167;
t95 = -t106 - t237;
t94 = t100 + t174;
t93 = t173 + t99;
t91 = -rSges(5,3) * t163 + t220;
t86 = -Icges(5,3) * t162 - t163 * t207;
t85 = -Icges(5,3) * t163 + t162 * t207;
t84 = t228 * t163;
t83 = t228 * t162;
t72 = t221 + t243;
t71 = pkin(7) * t163 - t155 + t245;
t70 = t167 + t82;
t69 = -t81 - t237;
t68 = t162 * t79;
t63 = t174 + t65;
t62 = t173 + t64;
t55 = -t162 * t290 + t163 * t193;
t54 = t162 * t193 + t290 * t163;
t53 = t174 + t61;
t52 = t173 + t60;
t49 = t163 * rSges(6,2) * t251 + (t163 * t252 + t171 * t255) * rSges(6,1) - t248;
t42 = pkin(4) * t191 - t121 * t253 - t256;
t41 = -t234 + t163 * t105 + (t163 * t240 - t236) * pkin(4);
t38 = -t163 * t80 + t68;
t30 = t167 + t32;
t29 = t31 - t237;
t28 = -t162 * t86 + t212 * t163;
t27 = -t162 * t85 + t200;
t26 = -t163 * t86 - t295;
t25 = -t214 * t162 - t163 * t85;
t22 = t167 + t24;
t21 = t23 - t237;
t16 = t162 * t71 + t68 + (-t72 - t80) * t163;
t11 = (t176 * t91 + t187) * t163 + (-t162 * t198 + (t92 + (-t274 - t277) * t163) * t176 + t246) * t162;
t8 = -t163 * t49 + t239;
t3 = (t202 * t163 + t176 * t71 - t138 - t49) * t163 + (-pkin(4) * t231 + t110 + (t72 - t163 * (-pkin(3) + t166) - t254) * t176 + t244) * t162 + t239;
t1 = [(t108 * t96 + t109 * t95) * t285 + (t69 * t94 + t70 * t93) * t284 + (t29 * t63 + t30 * t62) * t283 + (t21 * t53 + t22 * t52) * t282 + t185; m(3) * (-t106 * t109 + t107 * t108 + t122 * t96 + t123 * t95) + m(4) * (t100 * t69 + t70 * t99 - t81 * t94 + t82 * t93) + m(5) * (t29 * t65 + t30 * t64 + t31 * t63 + t32 * t62) + m(6) * (t21 * t61 + t22 * t60 + t23 * t53 + t24 * t52) + t185; (t23 * t61 + t24 * t60) * t282 + (t31 * t65 + t32 * t64) * t283 + (-t100 * t81 + t82 * t99) * t284 + (-t106 * t123 + t107 * t122) * t285 + t185; 0; 0; 0; ((-t176 * t63 + t30) * t163 + (-t176 * t62 - t29) * t162) * t278 + (-t162 * t63 + t163 * t62) * t279 + t184 + m(6) * (-t21 * t83 + t22 * t84 + t41 * t52 + t42 * t53); m(6) * (-t23 * t83 + t24 * t84 + t41 * t60 + t42 * t61) + ((-t176 * t65 + t32) * t163 + (-t176 * t64 - t31) * t162) * t278 + (-t162 * t65 + t163 * t64) * t279 + t184; m(5) * t11 + m(6) * t3; ((t162 * t91 - t163 * t92) * t11 + t242 * t147 * t131) * t283 + (-t162 * t26 - t163 * t25) * t255 - t163 * ((t163 * t55 + (-t26 + t200) * t176) * t163 + (t25 * t176 + (t240 * t88 + t241 * t90) * t162 + (t215 * qJD(4) + t176 * t212 + t54) * t163) * t162) - t162 * ((t162 * t54 + (t27 + t295) * t176) * t162 + (-t28 * t176 + (-t240 * t87 - t241 * t89) * t163 + (-t213 * qJD(4) + t214 * t176 + t55) * t162) * t163) + (t16 * t3 + t41 * t84 - t42 * t83) * t282 - t163 * t2 + t233 + (t28 * t162 + t163 * t27 - t5) * t253; m(6) * ((-t162 * t53 + t163 * t52) * t105 + ((-t176 * t53 + t22) * t163 + (-t176 * t52 - t21) * t162) * t121) + t199; m(6) * ((-t162 * t61 + t163 * t60) * t105 + ((-t176 * t61 + t24) * t163 + (-t176 * t60 - t23) * t162) * t121) + t199; m(6) * t8; m(6) * (-t121 * t162 * t42 + t16 * t8 - t84 * t234 + t83 * t256 + t3 * t38) + (m(6) * (t105 * t84 + (t176 * t83 + t41) * t121) + t232) * t163 + t233; (t105 * t121 * t242 + t38 * t8) * t282 + t232 * t163 + t233;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
