% Calculate time derivative of joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:35
% DurationCPUTime: 3.92s
% Computational Cost: add. (7461->372), mult. (6266->525), div. (0->0), fcn. (4654->8), ass. (0->220)
t170 = qJ(4) + qJ(5);
t163 = sin(t170);
t165 = cos(t170);
t216 = rSges(6,1) * t163 + rSges(6,2) * t165;
t171 = qJ(1) + qJ(2);
t164 = sin(t171);
t166 = cos(t171);
t78 = t166 * rSges(6,3) + t164 * t216;
t307 = rSges(4,2) - pkin(2);
t252 = t216 * t166;
t169 = qJD(1) + qJD(2);
t174 = cos(qJ(4));
t245 = qJD(4) * t174;
t172 = sin(qJ(4));
t257 = t166 * t172;
t305 = -t164 * t245 - t169 * t257;
t271 = Icges(6,4) * t165;
t111 = -Icges(6,2) * t163 + t271;
t274 = Icges(5,4) * t172;
t205 = Icges(5,2) * t174 + t274;
t117 = t205 * qJD(4);
t273 = Icges(5,4) * t174;
t207 = Icges(5,1) * t172 + t273;
t118 = t207 * qJD(4);
t168 = qJD(4) + qJD(5);
t272 = Icges(6,4) * t163;
t204 = Icges(6,2) * t165 + t272;
t112 = Icges(6,1) * t165 - t272;
t256 = t168 * t112;
t221 = -t204 * t168 + t256;
t259 = t165 * t168;
t206 = Icges(6,1) * t163 + t271;
t99 = t206 * t168;
t304 = -t111 * t259 + t172 * t117 - t174 * t118 - t163 * t221 - t165 * t99;
t176 = -pkin(8) - pkin(7);
t260 = t164 * t172;
t303 = pkin(4) * t260 - t166 * t176 + t78;
t89 = Icges(5,6) * t166 + t164 * t205;
t193 = t207 * t164;
t91 = Icges(5,5) * t166 + t193;
t211 = t172 * t91 + t174 * t89;
t302 = t166 * t211;
t74 = Icges(6,6) * t166 + t164 * t204;
t76 = Icges(6,5) * t166 + t164 * t206;
t215 = t163 * t76 + t165 * t74;
t301 = t166 * t215;
t300 = -Icges(6,5) * t169 + t256;
t263 = t111 * t168;
t299 = -Icges(6,6) * t169 + t263;
t110 = Icges(6,5) * t165 - Icges(6,6) * t163;
t298 = -Icges(6,3) * t169 + t110 * t168;
t134 = Icges(5,5) * t174 - Icges(5,6) * t172;
t297 = -Icges(5,3) * t169 + qJD(4) * t134;
t135 = -Icges(5,2) * t172 + t273;
t296 = -Icges(5,6) * t169 + qJD(4) * t135;
t136 = Icges(5,1) * t174 - t274;
t295 = -Icges(5,5) * t169 + qJD(4) * t136;
t294 = (t135 * t172 - t136 * t174) * qJD(4) + t117 * t174 + t118 * t172 + t134 * t169;
t293 = t163 * (t99 + t263) - t165 * t221 + t110 * t169;
t292 = 2 * m(3);
t291 = 2 * m(4);
t290 = 2 * m(5);
t289 = 2 * m(6);
t288 = t164 / 0.2e1;
t287 = t166 / 0.2e1;
t286 = -rSges(6,3) - pkin(2);
t280 = rSges(5,2) * t174;
t217 = rSges(5,1) * t172 + t280;
t121 = t217 * qJD(4);
t285 = m(5) * t121;
t173 = sin(qJ(1));
t284 = pkin(1) * t173;
t283 = pkin(7) * t164;
t160 = t166 * pkin(2);
t159 = t166 * pkin(7);
t282 = t159 - t303;
t278 = pkin(1) * qJD(1);
t277 = t164 * rSges(5,3);
t262 = t163 * t168;
t261 = t164 * t169;
t258 = t166 * t169;
t255 = t169 * t174;
t254 = t169 * t176;
t230 = t166 * t245;
t253 = -pkin(4) * t230 - t166 * t254;
t251 = qJ(3) * t258 + qJD(3) * t164;
t250 = rSges(5,1) * t257 + t166 * t280;
t147 = pkin(4) * t257;
t249 = t164 * t176 + t147;
t248 = t164 * qJ(3) + t160;
t247 = t164 ^ 2 + t166 ^ 2;
t246 = qJD(4) * t172;
t244 = -rSges(5,3) - pkin(2) - pkin(7);
t243 = rSges(6,1) * t259;
t242 = rSges(6,2) * t262;
t241 = t173 * t278;
t175 = cos(qJ(1));
t240 = t175 * t278;
t239 = rSges(5,2) * t246;
t237 = t166 * t255;
t236 = t164 * t243 + t252 * t169;
t235 = t305 * rSges(5,1) - rSges(5,2) * t237;
t234 = t305 * pkin(4) - t164 * t254;
t93 = rSges(5,1) * t260 + t166 * rSges(5,3) + t164 * t280;
t233 = t164 * t246;
t231 = t166 * t246;
t229 = -t261 / 0.2e1;
t228 = t258 / 0.2e1;
t114 = rSges(6,1) * t165 - rSges(6,2) * t163;
t227 = pkin(4) * t174 + t114;
t190 = t204 * t169;
t41 = t164 * t190 - t299 * t166;
t77 = Icges(6,5) * t164 - t166 * t206;
t226 = -t168 * t77 - t41;
t42 = t299 * t164 + t166 * t190;
t225 = -t168 * t76 - t42;
t192 = t206 * t169;
t43 = t164 * t192 - t300 * t166;
t75 = Icges(6,6) * t164 - t166 * t204;
t224 = -t168 * t75 + t43;
t44 = t300 * t164 + t166 * t192;
t223 = -t168 * t74 + t44;
t115 = t166 * rSges(3,1) - rSges(3,2) * t164;
t100 = t216 * t168;
t220 = t247 * t100;
t219 = t247 * t121;
t102 = -rSges(3,1) * t258 + rSges(3,2) * t261;
t113 = -rSges(3,1) * t164 - rSges(3,2) * t166;
t214 = t163 * t77 + t165 * t75;
t210 = -t172 * t89 + t174 * t91;
t90 = Icges(5,6) * t164 - t166 * t205;
t92 = Icges(5,5) * t164 - t166 * t207;
t209 = t172 * t92 + t174 * t90;
t208 = t172 * t90 - t174 * t92;
t203 = Icges(5,5) * t172 + Icges(5,6) * t174;
t202 = Icges(6,5) * t163 + Icges(6,6) * t165;
t201 = t165 * t111 + t163 * t112;
t200 = t174 * t135 + t172 * t136;
t152 = t166 * qJ(3);
t82 = t166 * rSges(4,3) + t307 * t164 + t152;
t198 = (t242 - t243) * t166;
t83 = -rSges(4,2) * t166 + t164 * rSges(4,3) + t248;
t195 = t214 * t164;
t72 = Icges(6,3) * t166 + t164 * t202;
t20 = t164 * t215 + t166 * t72;
t73 = Icges(6,3) * t164 - t166 * t202;
t21 = t166 * t73 + t195;
t22 = t164 * t72 - t301;
t23 = t164 * t73 - t166 * t214;
t188 = t202 * t169;
t39 = t164 * t188 - t298 * t166;
t40 = t298 * t164 + t166 * t188;
t197 = -t261 * (t164 * t21 + t166 * t20) + t164 * ((t164 * t39 + (-t22 + t195) * t169) * t164 + (t23 * t169 + (-t163 * t44 - t165 * t42 - t259 * t76 + t262 * t74) * t166 + (t40 + (t169 * t74 + t226) * t165 + (t169 * t76 - t224) * t163) * t164) * t166) + t166 * ((t166 * t40 + (t21 + t301) * t169) * t166 + (-t20 * t169 + (t163 * t43 + t165 * t41 + t259 * t77 - t262 * t75) * t164 + (t39 + (t169 * t75 - t225) * t165 + (t169 * t77 + t223) * t163) * t166) * t164) + (t164 * t23 + t166 * t22) * t258;
t68 = t159 + t93 + t248;
t185 = -t202 * t168 + t169 * t201;
t196 = (t163 * t226 + t185 * t164 + t165 * t224 + t293 * t166) * t288 + (t163 * t225 - t293 * t164 + t165 * t223 + t185 * t166) * t287 + (t110 * t166 - t163 * t74 + t164 * t201 + t165 * t76) * t229 + (t110 * t164 - t163 * t75 + t165 * t77 - t166 * t201) * t228;
t101 = t113 * t169;
t194 = t209 * t164;
t191 = t205 * t169;
t189 = t203 * t169;
t65 = rSges(4,3) * t258 + t307 * t261 + t251;
t67 = t164 * t244 + t152 + t250;
t184 = -t203 * qJD(4) + t169 * t200;
t54 = t248 + t303;
t53 = t164 * t286 + t152 + t249 + t252;
t150 = qJD(3) * t166;
t66 = rSges(4,2) * t258 + t150 + (-t160 + (-rSges(4,3) - qJ(3)) * t164) * t169;
t17 = (t169 * t286 - t242) * t164 - t234 + t236 + t251;
t15 = t17 - t241;
t18 = t150 + (t286 * t166 + (-pkin(4) * t172 - qJ(3) - t216) * t164) * t169 - t198 - t253;
t16 = t18 - t240;
t49 = t53 - t284;
t167 = t175 * pkin(1);
t50 = t167 + t54;
t183 = (t169 * t49 - t15) * t166 + (t169 * t50 + t16) * t164;
t182 = (t169 * t53 - t17) * t166 + (t169 * t54 + t18) * t164;
t37 = t114 * t261 + t100 * t166 + (t164 * t255 + t231) * pkin(4);
t38 = t114 * t258 - t100 * t164 + (-t233 + t237) * pkin(4);
t80 = t227 * t164;
t81 = t227 * t166;
t181 = (t169 * t80 - t37) * t166 + (-t169 * t81 + t38) * t164;
t30 = (t169 * t244 - t239) * t164 - t235 + t251;
t28 = t30 - t241;
t130 = rSges(5,1) * t230;
t31 = -rSges(5,2) * t231 + t130 + t150 + (t244 * t166 + (-qJ(3) - t217) * t164) * t169;
t29 = t31 - t240;
t63 = t67 - t284;
t64 = t167 + t68;
t180 = m(5) * ((t169 * t63 - t28) * t166 + (t169 * t64 + t29) * t164);
t179 = m(5) * ((t169 * t67 - t30) * t166 + (t169 * t68 + t31) * t164);
t178 = t196 + (t134 * t164 - t166 * t200) * t228 - t208 * t258 / 0.2e1 + (-qJD(4) * t209 + t184 * t164 + t294 * t166 - t172 * (t164 * t191 - t296 * t166) + t174 * (-t295 * t166 + t169 * t193)) * t288 + (-qJD(4) * t211 - t294 * t164 + t184 * t166 - t172 * (t296 * t164 + t166 * t191) + t174 * (t295 * t164 + t207 * t258)) * t287 + (t134 * t166 + t164 * t200 + t210) * t229;
t177 = -qJD(4) * t200 + t304;
t141 = rSges(5,1) * t174 - rSges(5,2) * t172;
t104 = t115 + t167;
t103 = t113 - t284;
t95 = -t249 - t283;
t94 = -t250 + t277;
t88 = Icges(5,3) * t164 - t166 * t203;
t87 = Icges(5,3) * t166 + t164 * t203;
t85 = t102 - t240;
t84 = t101 - t241;
t79 = rSges(6,3) * t164 - t252;
t71 = t167 + t83;
t70 = t82 - t284;
t69 = t166 * t79;
t62 = t66 - t240;
t61 = t65 - t241;
t56 = t297 * t164 + t166 * t189;
t55 = t164 * t189 - t297 * t166;
t48 = (-rSges(6,3) * t169 - t242) * t164 + t236;
t45 = -t164 * t78 + t69;
t36 = t166 * (t78 * t169 + t198);
t27 = t164 * t88 - t209 * t166;
t26 = t164 * t87 - t302;
t25 = t166 * t88 + t194;
t24 = t211 * t164 + t166 * t87;
t19 = t164 * t282 + t166 * t95 + t69;
t10 = -t78 * t258 + t36 + (-t169 * t79 - t48) * t164;
t3 = t166 * t253 + t36 + (-t48 + t234) * t164 + ((t282 - t159) * t166 + (t147 - t79 - t95 - t283) * t164) * t169;
t1 = [(t15 * t50 + t16 * t49) * t289 - t136 * t246 - t135 * t245 + (t28 * t64 + t29 * t63) * t290 + (t61 * t71 + t62 * t70) * t291 + (t103 * t85 + t104 * t84) * t292 + t304; m(6) * (t15 * t54 + t16 * t53 + t17 * t50 + t18 * t49) + m(5) * (t28 * t68 + t29 * t67 + t30 * t64 + t31 * t63) + m(4) * (t61 * t83 + t62 * t82 + t65 * t71 + t66 * t70) + m(3) * (t101 * t104 + t102 * t103 + t113 * t85 + t115 * t84) + t177; (t17 * t54 + t18 * t53) * t289 + (t30 * t68 + t31 * t67) * t290 + (t65 * t83 + t66 * t82) * t291 + (t101 * t115 + t102 * t113) * t292 + t177; m(6) * t183 + t180 + m(4) * ((t169 * t70 - t61) * t166 + (t169 * t71 + t62) * t164); m(6) * t182 + t179 + m(4) * ((t169 * t82 - t65) * t166 + (t169 * t83 + t66) * t164); 0; t178 + m(6) * (-t15 * t81 + t16 * t80 + t37 * t50 + t38 * t49) + t141 * t180 - (t164 * t63 - t166 * t64) * t285; t178 + m(6) * (-t17 * t81 + t18 * t80 + t37 * t54 + t38 * t53) + t141 * t179 - (t164 * t67 - t166 * t68) * t285; -m(5) * t219 + m(6) * t181; ((-t164 * t93 + t166 * t94) * ((-t169 * t93 - t130 + (rSges(5,3) * t169 + t239) * t166) * t166 + (rSges(5,2) * t233 + (t166 * t217 + t277 - t94) * t169 + t235) * t164) - t141 * t219) * t290 - (t164 * t25 + t24 * t166) * t261 + t166 * ((t166 * t56 + (t25 + t302) * t169) * t166 + (-t24 * t169 + (t245 * t92 - t246 * t90) * t164 + (t210 * qJD(4) + t169 * t209 + t55) * t166) * t164) + (t27 * t164 + t166 * t26) * t258 + t164 * ((t164 * t55 + (-t26 + t194) * t169) * t164 + (t27 * t169 + (-t245 * t91 + t246 * t89) * t166 + (t208 * qJD(4) + t211 * t169 + t56) * t164) * t166) + (t19 * t3 - t37 * t81 + t38 * t80) * t289 + t197; m(6) * (-(t164 * t49 - t166 * t50) * t100 + t183 * t114) + t196; m(6) * (-(t164 * t53 - t166 * t54) * t100 + t182 * t114) + t196; -m(6) * t220; m(6) * (t10 * t19 + t3 * t45 - (t164 * t80 + t166 * t81) * t100 + t181 * t114) + t197; (t10 * t45 - t114 * t220) * t289 + t197;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
