% Calculate time derivative of joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:43
% EndTime: 2019-12-31 19:05:52
% DurationCPUTime: 4.86s
% Computational Cost: add. (7775->374), mult. (14196->541), div. (0->0), fcn. (14312->8), ass. (0->199)
t313 = sin(qJ(3));
t314 = sin(qJ(1));
t315 = cos(qJ(3));
t316 = cos(qJ(1));
t162 = t316 * t313 - t314 * t315;
t161 = -t314 * t313 - t316 * t315;
t188 = qJ(4) + qJ(5);
t181 = sin(t188);
t182 = cos(t188);
t225 = -Icges(6,5) * t182 + Icges(6,6) * t181;
t89 = -Icges(6,3) * t161 + t162 * t225;
t355 = t162 * t89;
t90 = Icges(6,3) * t162 + t161 * t225;
t357 = t161 * t90;
t359 = t355 - t357;
t358 = t161 ^ 2;
t347 = t162 ^ 2;
t338 = t162 / 0.2e1;
t189 = sin(qJ(4));
t336 = -t189 / 0.2e1;
t190 = cos(qJ(4));
t335 = -t190 / 0.2e1;
t340 = qJD(1) - qJD(3);
t123 = t340 * t161;
t187 = qJD(4) + qJD(5);
t282 = t181 * t187;
t260 = t161 * t282;
t124 = t340 * t162;
t287 = t124 * t182;
t218 = t260 + t287;
t281 = t182 * t187;
t288 = t124 * t181;
t219 = t161 * t281 - t288;
t356 = t162 * (Icges(6,5) * t218 + Icges(6,6) * t219 + Icges(6,3) * t123);
t354 = t123 * t189;
t353 = t123 * t190;
t352 = t124 * t189;
t351 = t124 * t190;
t294 = Icges(5,4) * t190;
t228 = Icges(5,2) * t189 - t294;
t100 = Icges(5,6) * t162 + t161 * t228;
t295 = Icges(5,4) * t189;
t230 = -Icges(5,1) * t190 + t295;
t101 = -Icges(5,5) * t161 + t162 * t230;
t102 = Icges(5,5) * t162 + t161 * t230;
t226 = -Icges(5,5) * t190 + Icges(5,6) * t189;
t158 = t226 * qJD(4);
t166 = -Icges(5,5) * t189 - Icges(5,6) * t190;
t159 = t228 * qJD(4);
t160 = t230 * qJD(4);
t167 = -Icges(5,2) * t190 - t295;
t168 = -Icges(5,1) * t189 - t294;
t194 = t159 * t189 - t160 * t190 + (t167 * t190 + t168 * t189) * qJD(4);
t130 = t225 * t187;
t152 = -Icges(6,5) * t181 - Icges(6,6) * t182;
t293 = Icges(6,4) * t181;
t229 = -Icges(6,1) * t182 + t293;
t132 = t229 * t187;
t153 = -Icges(6,2) * t182 - t293;
t292 = Icges(6,4) * t182;
t154 = -Icges(6,1) * t181 - t292;
t227 = Icges(6,2) * t181 - t292;
t327 = (-t154 - t227) * t187;
t195 = (t153 * t187 - t132) * t182 - t327 * t181;
t259 = t162 * t282;
t290 = t123 * t182;
t220 = t259 - t290;
t291 = t123 * t181;
t221 = t162 * t281 + t291;
t223 = t153 * t181 - t154 * t182;
t339 = t161 / 0.2e1;
t344 = -t162 / 0.2e1;
t345 = -t124 / 0.2e1;
t91 = -Icges(6,6) * t161 + t162 * t227;
t92 = Icges(6,6) * t162 + t161 * t227;
t93 = -Icges(6,5) * t161 + t162 * t229;
t94 = Icges(6,5) * t162 + t161 * t229;
t210 = -(t152 * t162 + t161 * t223 - t181 * t94 - t182 * t92) * t123 / 0.2e1 + (-t152 * t161 + t162 * t223 - t181 * t93 - t182 * t91) * t345 + (t123 * t152 - t124 * t223 + t130 * t162 + t161 * t195 + t181 * (-Icges(6,1) * t218 - Icges(6,4) * t219 - Icges(6,5) * t123 + t187 * t92) - t182 * (Icges(6,4) * t218 + Icges(6,2) * t219 + Icges(6,6) * t123 + t187 * t94)) * t344 + (t123 * t223 + t124 * t152 - t130 * t161 + t162 * t195 + t181 * (-Icges(6,1) * t220 - Icges(6,4) * t221 - Icges(6,5) * t124 + t187 * t91) - t182 * (Icges(6,4) * t220 + Icges(6,2) * t221 + Icges(6,6) * t124 + t187 * t93)) * t339;
t269 = qJD(4) * t189;
t255 = t161 * t269;
t214 = t255 + t351;
t268 = qJD(4) * t190;
t215 = t161 * t268 - t352;
t254 = t162 * t269;
t216 = t254 - t353;
t217 = t162 * t268 + t354;
t222 = t167 * t189 - t168 * t190;
t224 = t100 * t189 - t102 * t190;
t99 = -Icges(5,6) * t161 + t162 * t228;
t233 = t101 * t190 - t189 * t99;
t333 = qJD(4) / 0.2e1;
t337 = -t166 / 0.2e1;
t350 = -t210 + (t100 * t335 + t102 * t336 + t166 * t338) * t123 + (t101 * t336 + t222 * t338 + t335 * t99) * t124 + (-t123 * t337 + t158 * t338 + t194 * t339 + t222 * t345 + t224 * t333 + (Icges(5,4) * t214 + Icges(5,2) * t215 + Icges(5,6) * t123) * t335 + (Icges(5,1) * t214 + Icges(5,4) * t215 + Icges(5,5) * t123) * t336) * t162 - (-0.2e1 * t124 * t337 - t158 * t339 - t194 * t344 - t233 * t333 + t335 * (Icges(5,4) * t216 + Icges(5,2) * t217 + Icges(5,6) * t124) + t336 * (Icges(5,1) * t216 + Icges(5,4) * t217 + Icges(5,5) * t124)) * t161;
t236 = rSges(6,1) * t181 + rSges(6,2) * t182;
t198 = pkin(4) * t269 + t187 * t236;
t277 = rSges(6,1) * t287 + t123 * rSges(6,3);
t178 = pkin(4) * t190 + pkin(3);
t191 = -pkin(8) - pkin(7);
t279 = -t123 * t191 + t124 * t178;
t27 = rSges(6,2) * t288 - t161 * t198 - t277 - t279;
t276 = t123 * t178 + t124 * t191;
t278 = -rSges(6,1) * t290 + t124 * rSges(6,3);
t26 = rSges(6,2) * t291 + t162 * t198 - t276 + t278;
t349 = -t181 * t132 + t153 * t282 - t189 * t160 + t167 * t269 + t327 * t182;
t270 = t316 * pkin(1) + t314 * qJ(2);
t256 = t316 * pkin(2) + t270;
t342 = -rSges(3,1) * t316 - rSges(3,3) * t314 - t270;
t341 = -rSges(6,1) * t182 + rSges(6,2) * t181;
t96 = t162 * rSges(6,3) + t341 * t161;
t76 = t161 * t178 + t162 * t191 - t96;
t332 = t224 * t162;
t331 = t233 * t161;
t328 = t347 + t358;
t326 = t162 * t123 - t161 * t124;
t121 = t123 * pkin(7);
t211 = t214 * rSges(5,1) + t123 * rSges(5,3);
t35 = -rSges(5,2) * t215 - t124 * pkin(3) - t121 - t211;
t212 = t216 * rSges(5,1) + t124 * rSges(5,3);
t252 = -t123 * pkin(3) + t124 * pkin(7);
t34 = t217 * rSges(5,2) + t212 + t252;
t150 = t162 * pkin(7);
t305 = rSges(5,2) * t189;
t253 = -pkin(3) + t305;
t307 = rSges(5,1) * t190;
t272 = t162 * rSges(5,3) - t161 * t307;
t81 = -t161 * t253 - t150 - t272;
t149 = t161 * pkin(7);
t274 = t161 * rSges(5,3) + t162 * t307;
t80 = t162 * t253 - t149 - t274;
t321 = -(qJD(4) * t168 + t159) * t190 + t349;
t319 = 2 * m(4);
t318 = 2 * m(5);
t317 = 2 * m(6);
t164 = (t305 - t307) * qJD(4);
t312 = m(5) * t164;
t169 = -rSges(5,1) * t189 - rSges(5,2) * t190;
t311 = m(5) * t169;
t310 = pkin(3) * t161;
t95 = -t161 * rSges(6,3) + t341 * t162;
t309 = t162 * (rSges(6,1) * t259 + rSges(6,2) * t221 + t278) + t123 * t95;
t308 = -t150 + t310 - t76;
t273 = t161 * t191 - t162 * t178;
t184 = t316 * qJ(2);
t271 = qJD(1) * t184 + qJD(2) * t314;
t267 = t314 * pkin(1);
t41 = Icges(6,5) * t220 + Icges(6,6) * t221 + Icges(6,3) * t124;
t261 = t347 * t356 + (-t347 * t41 + (-t161 * t41 + t356) * t161) * t161 + (0.3e1 * t358 * t89 + (-t357 + t359) * t162) * t124 - (-0.3e1 * t347 * t90 + (t355 + t359) * t161) * t123;
t251 = pkin(4) * t189 + t236;
t86 = -t124 * rSges(4,1) + t123 * rSges(4,2);
t85 = -t123 * rSges(4,1) - t124 * rSges(4,2);
t125 = -t162 * rSges(4,1) + t161 * rSges(4,2);
t126 = rSges(4,1) * t161 + rSges(4,2) * t162;
t209 = -pkin(2) * t314 - t267;
t205 = -t161 * t314 + t162 * t316;
t75 = t95 + t273;
t202 = t184 + t209;
t197 = -rSges(3,1) * t314 + rSges(3,3) * t316 - t267;
t196 = qJD(1) * t209 + t271;
t180 = qJD(2) * t316;
t193 = -t256 * qJD(1) + t180;
t192 = t314 * t124 + t316 * t123 + (-t161 * t316 - t162 * t314) * qJD(1);
t138 = t184 + t197;
t133 = t341 * t187;
t128 = t342 * qJD(1) + t180;
t127 = qJD(1) * t197 + t271;
t108 = -t126 + t256;
t107 = -t125 + t202;
t106 = t251 * t161;
t105 = t251 * t162;
t104 = t161 * t305 + t272;
t103 = t162 * t305 - t274;
t98 = Icges(5,3) * t162 + t161 * t226;
t97 = -Icges(5,3) * t161 + t162 * t226;
t87 = t162 * pkin(3) + t149 + t273;
t84 = t162 * t95;
t78 = t193 - t85;
t77 = t196 - t86;
t74 = t256 - t81;
t73 = t202 - t80;
t70 = -t76 + t256;
t69 = t202 - t75;
t60 = pkin(4) * t215 - t124 * t236 - t133 * t161;
t59 = pkin(4) * t217 + t123 * t236 - t133 * t162;
t58 = t161 * t96 + t84;
t53 = Icges(5,5) * t214 + Icges(5,6) * t215 + Icges(5,3) * t123;
t52 = Icges(5,5) * t216 + Icges(5,6) * t217 + Icges(5,3) * t124;
t47 = rSges(6,1) * t260 + rSges(6,2) * t219 + t277;
t39 = t161 * t224 + t98 * t162;
t38 = t162 * t97 - t331;
t37 = -t161 * t98 + t332;
t36 = -t97 * t161 - t162 * t233;
t29 = t193 - t34;
t28 = t196 - t35;
t25 = t193 - t26;
t24 = t196 - t27;
t21 = t161 * t308 + t162 * t87 + t84;
t4 = -t124 * t96 + t161 * t47 + t309;
t3 = t123 * t87 + t162 * (pkin(4) * t254 - t252 - t276) + (pkin(4) * t255 - t121 + t279 + t47) * t161 + (-t308 - t310) * t124 + t309;
t1 = [(t24 * t70 + t25 * t69) * t317 + (t28 * t74 + t29 * t73) * t318 - t168 * t268 - t190 * t159 + (t107 * t78 + t108 * t77) * t319 + 0.2e1 * m(3) * (-t127 * t342 + t128 * t138) + t349; m(6) * (t314 * t25 - t316 * t24 + (t314 * t70 + t316 * t69) * qJD(1)) + m(5) * (t314 * t29 - t316 * t28 + (t314 * t74 + t316 * t73) * qJD(1)) + m(4) * (t314 * t78 - t316 * t77 + (t107 * t316 + t108 * t314) * qJD(1)) + m(3) * (t314 * t128 - t316 * t127 + (t138 * t316 - t314 * t342) * qJD(1)); 0; m(6) * (t24 * t76 + t25 * t75 + t26 * t69 + t27 * t70) + m(5) * (t28 * t81 + t29 * t80 + t34 * t73 + t35 * t74) + m(4) * (t107 * t85 + t108 * t86 + t125 * t78 + t126 * t77) - t321; m(4) * (t85 * t314 - t86 * t316 + (t125 * t316 + t126 * t314) * qJD(1)) + m(5) * (t34 * t314 - t35 * t316 + (t314 * t81 + t316 * t80) * qJD(1)) + m(6) * (t26 * t314 - t27 * t316 + (t314 * t76 + t316 * t75) * qJD(1)); (t26 * t75 + t27 * t76) * t317 + (t34 * t80 + t35 * t81) * t318 + (t125 * t85 + t126 * t86) * t319 + t321; m(6) * (t105 * t24 + t106 * t25 + t59 * t70 + t60 * t69) + (-t161 * t73 - t162 * t74) * t312 + (-t123 * t74 + t124 * t73 - t161 * t29 - t162 * t28) * t311 + t350; m(6) * (t60 * t314 - t59 * t316 + (t105 * t314 + t106 * t316) * qJD(1)) + t205 * t312 + t192 * t311; m(6) * (t105 * t27 + t106 * t26 + t59 * t76 + t60 * t75) + (-t161 * t80 - t162 * t81) * t312 + (-t123 * t81 + t124 * t80 - t161 * t34 - t162 * t35) * t311 - t350; ((t103 * t162 + t104 * t161) * (t123 * t103 + t162 * t212 - t124 * t104 + t161 * t211 + (t161 * t215 + t162 * t217) * rSges(5,2)) + (t328 * t164 + t169 * t326) * t169) * t318 + t123 * (-t38 * t161 + t39 * t162) + t162 * ((t123 * t98 + t162 * t53) * t162 + t39 * t123 + (t38 - t332) * t124 + (-t101 * t351 - t123 * t97 - t162 * t52 + t352 * t99) * t161) + t124 * (-t36 * t161 + t37 * t162) - t161 * (-(t124 * t97 - t161 * t52) * t161 + t36 * t124 - (-t37 - t331) * t123 + (t100 * t354 - t102 * t353 + t124 * t98 - t161 * t53) * t162) + (t105 * t59 + t106 * t60 + t21 * t3) * t317 + t261; m(6) * ((-t161 * t69 - t162 * t70) * t133 - (-t123 * t70 + t124 * t69 - t161 * t25 - t162 * t24) * t236) - t210; m(6) * (t133 * t205 - t192 * t236); m(6) * ((-t161 * t75 - t162 * t76) * t133 - (-t123 * t76 + t124 * t75 - t161 * t26 - t162 * t27) * t236) + t210; m(6) * (t21 * t4 + t3 * t58 + (-t105 * t162 - t106 * t161) * t133 - (-t105 * t123 + t106 * t124 - t161 * t60 - t162 * t59) * t236) + t261; (t4 * t58 - (t133 * t328 - t236 * t326) * t236) * t317 + t261;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
