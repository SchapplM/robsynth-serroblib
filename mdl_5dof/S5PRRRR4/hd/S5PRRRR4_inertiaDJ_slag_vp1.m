% Calculate time derivative of joint inertia matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:47
% DurationCPUTime: 2.80s
% Computational Cost: add. (9061->332), mult. (5676->481), div. (0->0), fcn. (4272->8), ass. (0->203)
t157 = pkin(9) + qJ(2);
t154 = qJ(3) + t157;
t149 = sin(t154);
t150 = cos(t154);
t162 = cos(qJ(4));
t222 = qJD(4) * t162;
t159 = qJD(2) + qJD(3);
t161 = sin(qJ(4));
t232 = t159 * t161;
t268 = -t149 * t222 - t150 * t232;
t249 = rSges(5,2) * t161;
t251 = rSges(5,1) * t162;
t267 = -t249 + t251;
t243 = Icges(5,4) * t162;
t186 = -Icges(5,2) * t161 + t243;
t85 = -Icges(5,6) * t150 + t186 * t149;
t244 = Icges(5,4) * t161;
t188 = Icges(5,1) * t162 - t244;
t87 = -Icges(5,5) * t150 + t188 * t149;
t192 = t161 * t85 - t162 * t87;
t266 = t150 * t192;
t160 = qJ(4) + qJ(5);
t155 = sin(t160);
t156 = cos(t160);
t241 = Icges(6,4) * t156;
t185 = -Icges(6,2) * t155 + t241;
t73 = -Icges(6,6) * t150 + t185 * t149;
t242 = Icges(6,4) * t155;
t187 = Icges(6,1) * t156 - t242;
t75 = -Icges(6,5) * t150 + t187 * t149;
t197 = t155 * t73 - t156 * t75;
t265 = t150 * t197;
t177 = t186 * t150;
t86 = Icges(5,6) * t149 + t177;
t179 = t188 * t150;
t88 = Icges(5,5) * t149 + t179;
t190 = t161 * t86 - t162 * t88;
t264 = t190 * t149;
t176 = t185 * t150;
t74 = Icges(6,6) * t149 + t176;
t178 = t187 * t150;
t76 = Icges(6,5) * t149 + t178;
t196 = t155 * t74 - t156 * t76;
t263 = t196 * t149;
t138 = t149 * rSges(6,3);
t250 = rSges(6,1) * t156;
t262 = t150 * t250 + t138;
t158 = qJD(4) + qJD(5);
t114 = Icges(6,2) * t156 + t242;
t115 = Icges(6,1) * t155 + t241;
t182 = t114 * t155 - t115 * t156;
t183 = Icges(6,5) * t156 - Icges(6,6) * t155;
t261 = t183 * t158 + t159 * t182;
t135 = Icges(5,2) * t162 + t244;
t136 = Icges(5,1) * t161 + t243;
t181 = t135 * t161 - t136 * t162;
t184 = Icges(5,5) * t162 - Icges(5,6) * t161;
t260 = t184 * qJD(4) + t159 * t181;
t259 = 2 * m(4);
t258 = 2 * m(5);
t257 = 2 * m(6);
t256 = t149 / 0.2e1;
t255 = -t150 / 0.2e1;
t124 = t267 * qJD(4);
t254 = m(5) * t124;
t137 = rSges(5,1) * t161 + rSges(5,2) * t162;
t253 = m(5) * t137;
t152 = sin(t157);
t252 = pkin(2) * t152;
t143 = t149 * pkin(7);
t248 = rSges(6,2) * t155;
t122 = t149 * t248;
t227 = t150 * rSges(6,3) + t122;
t77 = t149 * t250 - t227;
t217 = t150 * t248;
t78 = -t217 + t262;
t39 = t149 * t77 + t150 * t78;
t247 = pkin(2) * qJD(2);
t139 = t149 * rSges(5,3);
t240 = t114 * t158;
t239 = t115 * t158;
t238 = t149 * t159;
t237 = t149 * t162;
t236 = t150 * t159;
t163 = -pkin(8) - pkin(7);
t235 = t150 * t163;
t234 = t155 * t158;
t233 = t156 * t158;
t231 = t159 * t163;
t230 = rSges(6,3) * t236 + t159 * t122;
t216 = t149 * t232;
t229 = rSges(5,2) * t216 + rSges(5,3) * t236;
t223 = qJD(4) * t161;
t212 = t149 * t223;
t228 = -pkin(4) * t212 - t149 * t231;
t226 = t150 * rSges(5,3) + t149 * t249;
t225 = -t150 * pkin(3) - t143;
t224 = t149 ^ 2 + t150 ^ 2;
t220 = rSges(6,2) * t233;
t214 = -t159 * t217 + (-t234 * rSges(6,1) - t220) * t149;
t221 = t149 * (t262 * t159 + t214) + t150 * (-t150 * t220 + (-t150 * t234 - t156 * t238) * rSges(6,1) + t230) + t77 * t236;
t219 = t152 * t247;
t153 = cos(t157);
t218 = t153 * t247;
t213 = -rSges(5,1) * t212 + t268 * rSges(5,2);
t210 = t150 * t223;
t209 = t150 * t222;
t208 = t238 / 0.2e1;
t207 = t236 / 0.2e1;
t206 = -pkin(3) - t251;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t205 = -pkin(4) * t161 - t116;
t172 = Icges(6,6) * t159 - t240;
t46 = t150 * t172 - t185 * t238;
t204 = t158 * t76 + t46;
t47 = t172 * t149 + t159 * t176;
t203 = t158 * t75 + t47;
t173 = Icges(6,5) * t159 - t239;
t48 = t150 * t173 - t187 * t238;
t202 = -t158 * t74 + t48;
t49 = t173 * t149 + t159 * t178;
t201 = -t158 * t73 + t49;
t151 = pkin(4) * t162 + pkin(3);
t200 = -t151 - t250;
t103 = t150 * rSges(4,1) - rSges(4,2) * t149;
t199 = -t149 * t163 + t150 * t151;
t94 = -rSges(4,1) * t236 + rSges(4,2) * t238;
t71 = -Icges(6,3) * t150 + t183 * t149;
t17 = -t197 * t149 - t150 * t71;
t174 = t183 * t150;
t72 = Icges(6,3) * t149 + t174;
t18 = -t150 * t72 - t263;
t19 = t149 * t71 - t265;
t20 = t149 * t72 - t150 * t196;
t113 = Icges(6,5) * t155 + Icges(6,6) * t156;
t171 = Icges(6,3) * t159 - t113 * t158;
t44 = t150 * t171 - t183 * t238;
t45 = t171 * t149 + t159 * t174;
t198 = -t150 * ((t150 * t45 + (t18 + t265) * t159) * t150 + (t17 * t159 + (-t155 * t46 + t156 * t48 - t233 * t74 - t234 * t76) * t149 + (-t44 + (t159 * t76 - t201) * t156 + (-t159 * t74 + t203) * t155) * t150) * t149) + t149 * ((t149 * t44 + (t19 + t263) * t159) * t149 + (t20 * t159 + (t155 * t47 - t156 * t49 + t233 * t73 + t234 * t75) * t150 + (-t45 + (t159 * t75 + t202) * t156 + (-t159 * t73 - t204) * t155) * t149) * t150) + (t149 * t18 - t150 * t17) * t238 + (t149 * t20 - t150 * t19) * t236;
t102 = -rSges(4,1) * t149 - rSges(4,2) * t150;
t193 = t161 * t87 + t162 * t85;
t191 = t161 * t88 + t162 * t86;
t189 = t200 * t149;
t134 = Icges(5,5) * t161 + Icges(5,6) * t162;
t90 = t267 * t150 + t139;
t93 = t102 * t159;
t100 = t187 * t158;
t99 = t185 * t158;
t165 = t113 * t159 + (t100 - t240) * t156 + (-t99 - t239) * t155;
t180 = (t261 * t149 + t165 * t150 + t155 * t202 + t156 * t204) * t256 + (t165 * t149 - t261 * t150 + t155 * t201 + t156 * t203) * t255 + (-t150 * t113 - t182 * t149 + t155 * t75 + t156 * t73) * t208 + (t149 * t113 - t150 * t182 + t155 * t76 + t156 * t74) * t207;
t175 = t184 * t150;
t65 = t90 - t225;
t144 = t150 * pkin(7);
t64 = t206 * t149 + t144 + t226;
t61 = t199 + t78;
t170 = Icges(5,5) * t159 - qJD(4) * t136;
t169 = Icges(5,6) * t159 - qJD(4) * t135;
t168 = Icges(5,3) * t159 - qJD(4) * t134;
t60 = t189 + t227 - t235;
t119 = t186 * qJD(4);
t120 = t188 * qJD(4);
t167 = t155 * t100 - t114 * t234 + t115 * t233 + t162 * t119 + t161 * t120 - t135 * t223 + t136 * t222 + t156 * t99;
t164 = -t119 * t161 + t120 * t162 + t134 * t159 + (-t135 * t162 - t136 * t161) * qJD(4);
t166 = t180 + (-qJD(4) * t190 + t260 * t149 + t164 * t150 + t161 * (t150 * t170 - t188 * t238) + t162 * (t150 * t169 - t186 * t238)) * t256 + (-qJD(4) * t192 + t164 * t149 - t260 * t150 + t161 * (t170 * t149 + t159 * t179) + t162 * (t169 * t149 + t159 * t177)) * t255 + (-t150 * t134 - t181 * t149 + t193) * t208 + (t149 * t134 - t150 * t181 + t191) * t207;
t34 = (t206 * t150 + (-rSges(5,3) - pkin(7)) * t149) * t159 - t213;
t28 = (t150 * t200 - t138) * t159 - t214 - t228;
t131 = pkin(7) * t236;
t33 = -rSges(5,2) * t209 - pkin(3) * t238 + t131 + (-t159 * t237 - t210) * rSges(5,1) + t229;
t27 = t159 * t189 + (-pkin(4) * t223 - t116 * t158 - t231) * t150 + t230;
t148 = pkin(2) * t153;
t101 = (-t248 + t250) * t158;
t96 = t103 + t148;
t95 = t102 - t252;
t89 = rSges(5,1) * t237 - t226;
t84 = Icges(5,3) * t149 + t175;
t83 = -Icges(5,3) * t150 + t184 * t149;
t82 = t205 * t150;
t81 = t205 * t149;
t80 = t94 - t218;
t79 = t93 - t219;
t70 = t199 + t225;
t69 = t235 + t144 + t149 * (-pkin(3) + t151);
t63 = t148 + t65;
t62 = t64 - t252;
t59 = t148 + t61;
t58 = t60 - t252;
t53 = t168 * t149 + t159 * t175;
t52 = t150 * t168 - t184 * t238;
t43 = t268 * pkin(4) - t149 * t101 - t116 * t236;
t42 = t116 * t238 - t150 * t101 + (-t209 + t216) * pkin(4);
t30 = t34 - t218;
t29 = t33 - t219;
t26 = t149 * t84 - t150 * t190;
t25 = t149 * t83 - t266;
t24 = -t150 * t84 - t264;
t23 = -t192 * t149 - t150 * t83;
t22 = t28 - t218;
t21 = t27 - t219;
t16 = t149 * t69 + t150 * t70 + t39;
t11 = ((-t90 + t139) * t159 + t213) * t149 + (-qJD(4) * t137 * t150 + t159 * t89 + t229) * t150;
t8 = -t238 * t78 + t221;
t3 = t149 * t228 + t150 * (-pkin(4) * t210 - t131) + ((t69 - t235) * t150 + (-t70 - t78 - t143) * t149) * t159 + t221;
t1 = [0; 0; (t21 * t59 + t22 * t58) * t257 + (t29 * t63 + t30 * t62) * t258 + (t79 * t96 + t80 * t95) * t259 + t167; 0; m(6) * (t21 * t61 + t22 * t60 + t27 * t59 + t28 * t58) + m(5) * (t29 * t65 + t30 * t64 + t33 * t63 + t34 * t62) + m(4) * (t102 * t80 + t103 * t79 + t93 * t96 + t94 * t95) + t167; (t27 * t61 + t28 * t60) * t257 + (t33 * t65 + t34 * t64) * t258 + (t102 * t94 + t103 * t93) * t259 + t167; m(5) * t11 + m(6) * t3; m(6) * (t21 * t81 + t22 * t82 + t42 * t58 + t43 * t59) + ((-t159 * t63 - t30) * t150 + (t159 * t62 - t29) * t149) * t253 + (-t149 * t63 - t150 * t62) * t254 + t166; ((-t159 * t65 - t34) * t150 + (t159 * t64 - t33) * t149) * t253 + (-t149 * t65 - t150 * t64) * t254 + m(6) * (t27 * t81 + t28 * t82 + t42 * t60 + t43 * t61) + t166; ((t149 * t89 + t150 * t90) * t11 + t224 * t137 * t124) * t258 + (t26 * t149 - t150 * t25) * t236 + t149 * ((t149 * t52 + (t25 + t264) * t159) * t149 + (t26 * t159 + (t222 * t85 + t223 * t87) * t150 + (-t191 * qJD(4) - t192 * t159 - t53) * t149) * t150) + (t149 * t24 - t150 * t23) * t238 - t150 * ((t150 * t53 + (t24 + t266) * t159) * t150 + (t23 * t159 + (-t222 * t86 - t223 * t88) * t149 + (t193 * qJD(4) - t190 * t159 - t52) * t150) * t149) + (t16 * t3 + t42 * t82 + t43 * t81) * t257 + t198; m(6) * t8; m(6) * ((-t149 * t59 - t150 * t58) * t101 + ((-t159 * t59 - t22) * t150 + (t159 * t58 - t21) * t149) * t116) + t180; m(6) * ((-t149 * t61 - t150 * t60) * t101 + ((-t159 * t61 - t28) * t150 + (t159 * t60 - t27) * t149) * t116) + t180; m(6) * (t16 * t8 + t3 * t39 + (-t149 * t81 - t150 * t82) * t101 + ((-t159 * t81 - t42) * t150 + (t159 * t82 - t43) * t149) * t116) + t198; (t101 * t116 * t224 + t39 * t8) * t257 + t198;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
