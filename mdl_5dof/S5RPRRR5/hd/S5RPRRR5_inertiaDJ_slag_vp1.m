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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:48:35
% EndTime: 2022-01-20 09:48:43
% DurationCPUTime: 3.02s
% Computational Cost: add. (9111->334), mult. (5758->482), div. (0->0), fcn. (4324->10), ass. (0->200)
t160 = qJ(1) + pkin(9);
t154 = qJ(3) + t160;
t149 = sin(t154);
t150 = cos(t154);
t164 = cos(qJ(4));
t228 = qJD(4) * t164;
t159 = qJD(1) + qJD(3);
t162 = sin(qJ(4));
t239 = t159 * t162;
t275 = -t149 * t228 - t150 * t239;
t255 = rSges(5,2) * t162;
t257 = rSges(5,1) * t164;
t274 = -t255 + t257;
t250 = Icges(5,4) * t164;
t192 = -Icges(5,2) * t162 + t250;
t85 = -Icges(5,6) * t150 + t192 * t149;
t251 = Icges(5,4) * t162;
t194 = Icges(5,1) * t164 - t251;
t87 = -Icges(5,5) * t150 + t194 * t149;
t198 = t162 * t85 - t164 * t87;
t273 = t150 * t198;
t161 = qJ(4) + qJ(5);
t155 = sin(t161);
t156 = cos(t161);
t248 = Icges(6,4) * t156;
t191 = -Icges(6,2) * t155 + t248;
t75 = -Icges(6,6) * t150 + t191 * t149;
t249 = Icges(6,4) * t155;
t193 = Icges(6,1) * t156 - t249;
t77 = -Icges(6,5) * t150 + t193 * t149;
t203 = t155 * t75 - t156 * t77;
t272 = t150 * t203;
t181 = t192 * t150;
t86 = Icges(5,6) * t149 + t181;
t183 = t194 * t150;
t88 = Icges(5,5) * t149 + t183;
t196 = t162 * t86 - t164 * t88;
t271 = t196 * t149;
t180 = t191 * t150;
t76 = Icges(6,6) * t149 + t180;
t182 = t193 * t150;
t78 = Icges(6,5) * t149 + t182;
t202 = t155 * t76 - t156 * t78;
t270 = t202 * t149;
t138 = t149 * rSges(6,3);
t256 = rSges(6,1) * t156;
t269 = t150 * t256 + t138;
t268 = -pkin(2) * cos(t160) - cos(qJ(1)) * pkin(1);
t158 = qJD(4) + qJD(5);
t114 = Icges(6,2) * t156 + t249;
t115 = Icges(6,1) * t155 + t248;
t188 = t114 * t155 - t115 * t156;
t189 = Icges(6,5) * t156 - Icges(6,6) * t155;
t267 = t189 * t158 + t159 * t188;
t135 = Icges(5,2) * t164 + t251;
t136 = Icges(5,1) * t162 + t250;
t187 = t135 * t162 - t136 * t164;
t190 = Icges(5,5) * t164 - Icges(5,6) * t162;
t266 = t190 * qJD(4) + t159 * t187;
t265 = 2 * m(4);
t264 = 2 * m(5);
t263 = 2 * m(6);
t262 = t149 / 0.2e1;
t261 = -t150 / 0.2e1;
t124 = t274 * qJD(4);
t260 = m(5) * t124;
t137 = rSges(5,1) * t162 + rSges(5,2) * t164;
t259 = m(5) * t137;
t143 = t149 * pkin(7);
t254 = rSges(6,2) * t155;
t122 = t149 * t254;
t234 = t150 * rSges(6,3) + t122;
t79 = t149 * t256 - t234;
t225 = t150 * t254;
t80 = -t225 + t269;
t39 = t149 * t79 + t150 * t80;
t139 = t149 * rSges(5,3);
t247 = t114 * t158;
t246 = t115 * t158;
t245 = t149 * t159;
t244 = t149 * t164;
t243 = t150 * t159;
t166 = -pkin(8) - pkin(7);
t242 = t150 * t166;
t241 = t155 * t158;
t240 = t156 * t158;
t238 = t159 * t166;
t237 = rSges(6,3) * t243 + t159 * t122;
t224 = t149 * t239;
t236 = rSges(5,2) * t224 + rSges(5,3) * t243;
t229 = qJD(4) * t162;
t220 = t149 * t229;
t235 = -pkin(4) * t220 - t149 * t238;
t233 = t150 * rSges(5,3) + t149 * t255;
t232 = -t150 * pkin(3) - t143;
t231 = t149 ^ 2 + t150 ^ 2;
t226 = rSges(6,2) * t240;
t222 = -t159 * t225 + (-t241 * rSges(6,1) - t226) * t149;
t227 = t149 * (t269 * t159 + t222) + t150 * (-t150 * t226 + (-t150 * t241 - t156 * t245) * rSges(6,1) + t237) + t79 * t243;
t221 = -rSges(5,1) * t220 + t275 * rSges(5,2);
t218 = t150 * t229;
t217 = t150 * t228;
t216 = t245 / 0.2e1;
t215 = t243 / 0.2e1;
t214 = -pkin(3) - t257;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t213 = -pkin(4) * t162 - t116;
t175 = Icges(6,6) * t159 - t247;
t46 = t150 * t175 - t191 * t245;
t212 = t158 * t78 + t46;
t47 = t175 * t149 + t159 * t180;
t211 = t158 * t77 + t47;
t176 = Icges(6,5) * t159 - t246;
t48 = t150 * t176 - t193 * t245;
t210 = -t158 * t76 + t48;
t49 = t176 * t149 + t159 * t182;
t209 = -t158 * t75 + t49;
t151 = pkin(4) * t164 + pkin(3);
t208 = -t151 - t256;
t103 = t150 * rSges(4,1) - rSges(4,2) * t149;
t207 = -t149 * t166 + t150 * t151;
t96 = -rSges(4,1) * t243 + rSges(4,2) * t245;
t73 = -Icges(6,3) * t150 + t189 * t149;
t17 = -t203 * t149 - t150 * t73;
t178 = t189 * t150;
t74 = Icges(6,3) * t149 + t178;
t18 = -t150 * t74 - t270;
t19 = t149 * t73 - t272;
t20 = t149 * t74 - t150 * t202;
t113 = Icges(6,5) * t155 + Icges(6,6) * t156;
t174 = Icges(6,3) * t159 - t113 * t158;
t44 = t150 * t174 - t189 * t245;
t45 = t174 * t149 + t159 * t178;
t206 = -t150 * ((t150 * t45 + (t18 + t272) * t159) * t150 + (t17 * t159 + (-t155 * t46 + t156 * t48 - t240 * t76 - t241 * t78) * t149 + (-t44 + (t159 * t78 - t209) * t156 + (-t159 * t76 + t211) * t155) * t150) * t149) + t149 * ((t149 * t44 + (t19 + t270) * t159) * t149 + (t20 * t159 + (t155 * t47 - t156 * t49 + t240 * t75 + t241 * t77) * t150 + (-t45 + (t159 * t77 + t210) * t156 + (-t159 * t75 - t212) * t155) * t149) * t150) + (t149 * t18 - t150 * t17) * t245 + (t149 * t20 - t150 * t19) * t243;
t205 = -pkin(2) * sin(t160) - sin(qJ(1)) * pkin(1);
t102 = -rSges(4,1) * t149 - rSges(4,2) * t150;
t199 = t162 * t87 + t164 * t85;
t197 = t162 * t88 + t164 * t86;
t195 = t208 * t149;
t134 = Icges(5,5) * t162 + Icges(5,6) * t164;
t90 = t274 * t150 + t139;
t186 = t205 * qJD(1);
t185 = t268 * qJD(1);
t95 = t102 * t159;
t100 = t193 * t158;
t99 = t191 * t158;
t168 = t113 * t159 + (t100 - t247) * t156 + (-t99 - t246) * t155;
t184 = (t149 * t267 + t168 * t150 + t155 * t210 + t156 * t212) * t262 + (t168 * t149 - t150 * t267 + t155 * t209 + t156 * t211) * t261 + (-t150 * t113 - t188 * t149 + t155 * t77 + t156 * t75) * t216 + (t149 * t113 - t150 * t188 + t155 * t78 + t156 * t76) * t215;
t179 = t190 * t150;
t65 = t90 - t232;
t144 = t150 * pkin(7);
t64 = t214 * t149 + t144 + t233;
t63 = t207 + t80;
t173 = Icges(5,5) * t159 - qJD(4) * t136;
t172 = Icges(5,6) * t159 - qJD(4) * t135;
t171 = Icges(5,3) * t159 - qJD(4) * t134;
t62 = t195 + t234 - t242;
t119 = t192 * qJD(4);
t120 = t194 * qJD(4);
t170 = t155 * t100 - t114 * t241 + t115 * t240 + t164 * t119 + t162 * t120 - t135 * t229 + t136 * t228 + t156 * t99;
t167 = -t119 * t162 + t120 * t164 + t134 * t159 + (-t135 * t164 - t136 * t162) * qJD(4);
t169 = t184 + (-qJD(4) * t196 + t149 * t266 + t167 * t150 + t162 * (t150 * t173 - t194 * t245) + t164 * (t150 * t172 - t192 * t245)) * t262 + (-qJD(4) * t198 + t167 * t149 - t150 * t266 + t162 * (t173 * t149 + t159 * t183) + t164 * (t172 * t149 + t159 * t181)) * t261 + (-t150 * t134 - t187 * t149 + t199) * t216 + (t149 * t134 - t150 * t187 + t197) * t215;
t34 = (t214 * t150 + (-rSges(5,3) - pkin(7)) * t149) * t159 - t221;
t28 = (t150 * t208 - t138) * t159 - t222 - t235;
t131 = pkin(7) * t243;
t33 = -rSges(5,2) * t217 - pkin(3) * t245 + t131 + (-t159 * t244 - t218) * rSges(5,1) + t236;
t27 = t159 * t195 + (-pkin(4) * t229 - t116 * t158 - t238) * t150 + t237;
t101 = (-t254 + t256) * t158;
t92 = t103 - t268;
t91 = t102 + t205;
t89 = rSges(5,1) * t244 - t233;
t84 = Icges(5,3) * t149 + t179;
t83 = -Icges(5,3) * t150 + t190 * t149;
t82 = t213 * t150;
t81 = t213 * t149;
t72 = t185 + t96;
t71 = t95 + t186;
t70 = t207 + t232;
t69 = t242 + t144 + t149 * (-pkin(3) + t151);
t61 = t65 - t268;
t60 = t205 + t64;
t55 = t171 * t149 + t159 * t179;
t54 = t150 * t171 - t190 * t245;
t53 = t63 - t268;
t52 = t205 + t62;
t43 = t275 * pkin(4) - t149 * t101 - t116 * t243;
t42 = t116 * t245 - t150 * t101 + (-t217 + t224) * pkin(4);
t30 = t185 + t34;
t29 = t186 + t33;
t26 = t149 * t84 - t150 * t196;
t25 = t149 * t83 - t273;
t24 = -t150 * t84 - t271;
t23 = -t198 * t149 - t150 * t83;
t22 = t185 + t28;
t21 = t186 + t27;
t16 = t149 * t69 + t150 * t70 + t39;
t11 = ((-t90 + t139) * t159 + t221) * t149 + (-qJD(4) * t137 * t150 + t159 * t89 + t236) * t150;
t8 = -t245 * t80 + t227;
t3 = t149 * t235 + t150 * (-pkin(4) * t218 - t131) + ((t69 - t242) * t150 + (-t70 - t80 - t143) * t149) * t159 + t227;
t1 = [(t21 * t53 + t22 * t52) * t263 + (t29 * t61 + t30 * t60) * t264 + (t71 * t92 + t72 * t91) * t265 + t170; 0; 0; m(6) * (t21 * t63 + t22 * t62 + t27 * t53 + t28 * t52) + m(5) * (t29 * t65 + t30 * t64 + t33 * t61 + t34 * t60) + m(4) * (t102 * t72 + t103 * t71 + t91 * t96 + t92 * t95) + t170; 0; (t27 * t63 + t28 * t62) * t263 + (t33 * t65 + t34 * t64) * t264 + (t102 * t96 + t103 * t95) * t265 + t170; ((-t159 * t61 - t30) * t150 + (t159 * t60 - t29) * t149) * t259 + (-t149 * t61 - t150 * t60) * t260 + m(6) * (t21 * t81 + t22 * t82 + t42 * t52 + t43 * t53) + t169; m(5) * t11 + m(6) * t3; m(6) * (t27 * t81 + t28 * t82 + t42 * t62 + t43 * t63) + ((-t159 * t65 - t34) * t150 + (t159 * t64 - t33) * t149) * t259 + (-t149 * t65 - t150 * t64) * t260 + t169; ((t149 * t89 + t150 * t90) * t11 + t231 * t137 * t124) * t264 + (t149 * t26 - t150 * t25) * t243 + t149 * ((t149 * t54 + (t25 + t271) * t159) * t149 + (t26 * t159 + (t228 * t85 + t229 * t87) * t150 + (-t197 * qJD(4) - t198 * t159 - t55) * t149) * t150) + (t149 * t24 - t150 * t23) * t245 - t150 * ((t150 * t55 + (t24 + t273) * t159) * t150 + (t23 * t159 + (-t228 * t86 - t229 * t88) * t149 + (t199 * qJD(4) - t196 * t159 - t54) * t150) * t149) + (t16 * t3 + t42 * t82 + t43 * t81) * t263 + t206; m(6) * ((-t149 * t53 - t150 * t52) * t101 + ((-t159 * t53 - t22) * t150 + (t159 * t52 - t21) * t149) * t116) + t184; m(6) * t8; m(6) * ((-t149 * t63 - t150 * t62) * t101 + ((-t159 * t63 - t28) * t150 + (t159 * t62 - t27) * t149) * t116) + t184; m(6) * (t16 * t8 + t3 * t39 + (-t149 * t81 - t150 * t82) * t101 + ((-t159 * t81 - t42) * t150 + (t159 * t82 - t43) * t149) * t116) + t206; (t101 * t116 * t231 + t39 * t8) * t263 + t206;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
