% Calculate vector of inverse dynamics joint torques for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:11
% DurationCPUTime: 3.16s
% Computational Cost: add. (10669->331), mult. (5800->422), div. (0->0), fcn. (4382->8), ass. (0->204)
t165 = pkin(9) + qJ(2);
t158 = sin(t165);
t254 = pkin(2) * qJD(2);
t221 = t158 * t254;
t161 = qJ(3) + t165;
t153 = sin(t161);
t166 = qJD(2) + qJD(3);
t233 = t153 * t166;
t224 = pkin(3) * t233;
t155 = qJ(4) + t161;
t149 = sin(t155);
t150 = cos(t155);
t100 = rSges(5,1) * t149 + rSges(5,2) * t150;
t160 = qJD(4) + t166;
t92 = t160 * t100;
t62 = -t221 - t224 - t92;
t144 = t150 * pkin(8);
t103 = pkin(4) * t149 - t144;
t236 = t150 * t160;
t115 = pkin(8) * t236;
t167 = sin(qJ(5));
t238 = t149 * t167;
t125 = rSges(6,2) * t238;
t168 = cos(qJ(5));
t255 = rSges(6,2) * t168;
t218 = qJD(5) * t255;
t195 = rSges(6,3) * t236 + t160 * t125 - t150 * t218;
t225 = qJD(5) * t167;
t215 = t150 * t225;
t231 = t150 * rSges(6,3) + t125;
t237 = t149 * t168;
t73 = rSges(6,1) * t237 - t231;
t66 = t160 * t73;
t282 = -rSges(6,1) * t215 + t160 * t103 + t115 + t195 + t66;
t154 = cos(t161);
t109 = rSges(4,1) * t153 + rSges(4,2) * t154;
t243 = t109 * t166;
t78 = -t221 - t243;
t256 = rSges(6,1) * t168;
t136 = -rSges(6,2) * t167 + t256;
t121 = t136 * qJD(5);
t135 = rSges(6,1) * t167 + t255;
t164 = qJDD(2) + qJDD(3);
t157 = qJDD(4) + t164;
t163 = t166 ^ 2;
t159 = cos(t165);
t169 = qJD(2) ^ 2;
t191 = (-qJDD(2) * t158 - t159 * t169) * pkin(2);
t176 = t191 + (-t153 * t164 - t154 * t163) * pkin(3);
t227 = qJD(5) * t150;
t248 = -t103 - t73;
t273 = -t150 * pkin(4) - t149 * pkin(8);
t235 = t150 * t167;
t222 = rSges(6,2) * t235;
t217 = t160 * t222 + (rSges(6,1) * t225 + t218) * t149;
t234 = t150 * t168;
t274 = rSges(6,1) * t234 + t149 * rSges(6,3);
t46 = t274 * t160 - t217;
t226 = qJD(5) * t160;
t85 = -qJDD(5) * t150 + t149 * t226;
t12 = -t121 * t227 + t85 * t135 + t248 * t157 + (t273 * t160 - t46) * t160 + t176;
t281 = -g(1) + t12;
t239 = t149 * t160;
t81 = rSges(5,1) * t236 - rSges(5,2) * t239;
t280 = -t100 * t157 - t160 * t81 - g(1) + t176;
t232 = t154 * t166;
t97 = rSges(4,1) * t232 - rSges(4,2) * t233;
t279 = -t109 * t164 - t166 * t97 - g(1) + t191;
t148 = pkin(3) * t154;
t152 = pkin(2) * t159;
t264 = pkin(2) * t158;
t207 = qJDD(2) * t152 - t169 * t264;
t263 = pkin(3) * t153;
t185 = t164 * t148 - t163 * t263 + t207;
t228 = qJD(5) * t149;
t45 = (-t160 * t237 - t215) * rSges(6,1) + t195;
t74 = -t222 + t274;
t58 = t74 - t273;
t84 = qJDD(5) * t149 + t150 * t226;
t13 = -t121 * t228 - t84 * t135 + (-pkin(4) * t239 + t115 + t45) * t160 + t58 * t157 + t185;
t278 = -g(2) + t13;
t101 = t150 * rSges(5,1) - rSges(5,2) * t149;
t277 = t101 * t157 - t160 * t92 - g(2) + t185;
t110 = t154 * rSges(4,1) - rSges(4,2) * t153;
t276 = t110 * t164 - t166 * t243 - g(2) + t207;
t162 = Icges(6,4) * t168;
t198 = -Icges(6,2) * t167 + t162;
t133 = Icges(6,1) * t167 + t162;
t272 = -t224 + t282;
t200 = t135 * t228 - t160 * t58;
t130 = Icges(6,5) * t168 - Icges(6,6) * t167;
t129 = Icges(6,5) * t167 + Icges(6,6) * t168;
t181 = Icges(6,3) * t160 - qJD(5) * t129;
t193 = t198 * t150;
t70 = Icges(6,6) * t149 + t193;
t250 = t167 * t70;
t246 = Icges(6,4) * t167;
t134 = Icges(6,1) * t168 - t246;
t194 = t134 * t150;
t72 = Icges(6,5) * t149 + t194;
t201 = -t168 * t72 + t250;
t271 = -t130 * t239 + t150 * t181 + t160 * t201;
t192 = t130 * t150;
t69 = Icges(6,4) * t237 - Icges(6,2) * t238 - Icges(6,6) * t150;
t251 = t167 * t69;
t124 = Icges(6,4) * t238;
t71 = Icges(6,1) * t237 - Icges(6,5) * t150 - t124;
t202 = -t168 * t71 + t251;
t270 = t149 * t181 + (t192 + t202) * t160;
t131 = Icges(6,2) * t168 + t246;
t196 = t131 * t167 - t133 * t168;
t269 = t130 * qJD(5) + t160 * t196;
t67 = Icges(6,5) * t237 - Icges(6,6) * t238 - Icges(6,3) * t150;
t24 = -t149 * t202 - t150 * t67;
t214 = -pkin(4) - t256;
t216 = t135 * t227;
t188 = -t216 - t224;
t179 = t188 - t221;
t28 = t160 * t248 + t179;
t253 = t150 * t28;
t220 = t159 * t254;
t223 = pkin(3) * t232;
t189 = t220 + t223;
t29 = t189 - t200;
t170 = (t214 * t253 + (t28 * (-rSges(6,3) - pkin(8)) + t29 * t214) * t149) * t160;
t268 = t170 + (t217 - t200) * t28;
t258 = -Icges(6,2) * t237 - t124 + t71;
t260 = t133 * t149 + t69;
t267 = -t167 * t258 - t168 * t260;
t266 = t84 / 0.2e1;
t265 = t85 / 0.2e1;
t262 = -t149 * t67 - t71 * t234;
t68 = Icges(6,3) * t149 + t192;
t261 = t149 * t68 + t72 * t234;
t259 = -t133 * t150 - t70;
t257 = -t131 * t150 + t72;
t241 = t129 * t150;
t53 = -t149 * t196 - t241;
t249 = t53 * t160;
t244 = t101 * t160;
t242 = t129 * t149;
t240 = t130 * t160;
t230 = -t131 + t134;
t229 = t133 + t198;
t219 = m(2) + m(3) + m(4) + m(5);
t213 = -t228 / 0.2e1;
t212 = t228 / 0.2e1;
t211 = -t227 / 0.2e1;
t210 = t227 / 0.2e1;
t59 = t72 * t237;
t209 = t150 * t68 - t59;
t208 = -t67 + t250;
t83 = t101 + t148;
t112 = rSges(3,1) * t159 - rSges(3,2) * t158;
t111 = rSges(3,1) * t158 + rSges(3,2) * t159;
t25 = -t238 * t70 - t209;
t206 = t149 * t25 - t24 * t150;
t26 = -t235 * t69 - t262;
t27 = -t235 * t70 + t261;
t205 = t27 * t149 - t150 * t26;
t204 = -t149 * t29 - t253;
t203 = t149 * t73 + t150 * t74;
t43 = t167 * t71 + t168 * t69;
t44 = t167 * t72 + t168 * t70;
t197 = t131 * t168 + t167 * t133;
t82 = -t100 - t263;
t187 = -t167 * t257 + t168 * t259;
t186 = -t81 - t223;
t57 = t149 * t214 + t144 + t231;
t56 = t148 + t58;
t184 = (-t167 * t229 + t168 * t230) * t160;
t183 = Icges(6,5) * t160 - qJD(5) * t133;
t182 = Icges(6,6) * t160 - qJD(5) * t131;
t55 = t57 - t263;
t54 = -t150 * t196 + t242;
t52 = t54 * t160;
t10 = qJD(5) * t205 + t52;
t119 = t198 * qJD(5);
t120 = t134 * qJD(5);
t40 = t149 * t182 + t160 * t193;
t42 = t149 * t183 + t160 * t194;
t16 = -qJD(5) * t202 + t167 * t42 + t168 * t40;
t39 = t150 * t182 - t198 * t239;
t41 = -t134 * t239 + t150 * t183;
t17 = -qJD(5) * t201 + t167 * t41 + t168 * t39;
t172 = -qJD(5) * t197 - t119 * t167 + t120 * t168 + t129 * t160;
t20 = t269 * t149 + t172 * t150;
t21 = t172 * t149 - t269 * t150;
t9 = qJD(5) * t206 + t249;
t177 = (t52 + ((t25 - t59 + (t68 + t251) * t150 + t262) * t150 + t261 * t149) * qJD(5)) * t210 + (-qJD(5) * t196 + t119 * t168 + t120 * t167) * t160 + (t44 + t54) * t266 + (t43 + t53) * t265 + (-t249 + ((t150 * t208 - t261 + t27) * t150 + (t149 * t208 + t209 + t26) * t149) * qJD(5) + t9) * t213 + (t17 + t20) * t212 + (Icges(5,3) + t197) * t157 + (t16 + t21 + t10) * t211;
t175 = Icges(4,3) * t164 + t177;
t174 = -qJD(5) * t43 + t160 * t67 - t167 * t40 + t168 * t42;
t173 = -qJD(5) * t44 + t160 * t68 - t167 * t39 + t168 * t41;
t94 = t135 * t150;
t93 = t135 * t149;
t79 = t110 * t166 + t220;
t63 = t189 + t244;
t34 = qJD(5) * t203 + qJD(1);
t11 = t73 * t84 - t74 * t85 + qJDD(1) + (t149 * t46 + t150 * t45) * qJD(5);
t6 = t173 * t149 - t271 * t150;
t5 = t174 * t149 - t270 * t150;
t4 = t271 * t149 + t173 * t150;
t3 = t270 * t149 + t174 * t150;
t1 = [m(6) * t11 + t219 * qJDD(1) + (-m(6) - t219) * g(3); Icges(3,3) * qJDD(2) + t175 + (t276 * (t110 + t152) + t279 * (-t109 - t264) + (-t97 - t220 + t79) * t78) * m(4) + (g(1) * t111 - g(2) * t112 + (t111 ^ 2 + t112 ^ 2) * qJDD(2)) * m(3) + (t28 * (-t189 + t217) + t170 + t278 * (t152 + t56) + t281 * (t55 - t264) + (-t221 - t179 + t28 + t272) * t29) * m(6) + (t277 * (t152 + t83) + t280 * (t82 - t264) + (t63 + t186 - t220) * t62) * m(5); t175 + (t278 * t56 + t281 * t55 + (-t188 + t272) * t29 + t268) * m(6) + (t277 * t83 + t280 * t82 + (t186 + t223 + t244) * t62) * m(5) + (-t78 * t97 - t79 * t243 + (t78 * t166 + t276) * t110 + (t79 * t166 - t279) * t109) * m(4); t177 + (t278 * t58 + t281 * t57 + (t216 + t282) * t29 + t268) * m(6) + (-t62 * t81 - t63 * t92 + (t160 * t62 + t277) * t101 + (t160 * t63 - t280) * t100) * m(5); t10 * t236 / 0.2e1 + t149 * (t54 * t157 + t20 * t160 + t26 * t85 + t27 * t84 + (t149 * t4 - t3 * t150) * qJD(5)) / 0.2e1 + t205 * t266 + ((t160 * t27 - t3) * t150 + (t160 * t26 + t4) * t149) * t212 + t9 * t239 / 0.2e1 - t150 * (t53 * t157 + t21 * t160 + t24 * t85 + t25 * t84 + (t6 * t149 - t150 * t5) * qJD(5)) / 0.2e1 + t206 * t265 + ((t160 * t25 - t5) * t150 + (t160 * t24 + t6) * t149) * t211 + t157 * (t44 * t149 - t43 * t150) / 0.2e1 + t160 * ((t160 * t44 - t16) * t150 + (t160 * t43 + t17) * t149) / 0.2e1 + ((-t228 * t241 + t240) * t149 + (t184 + (-t267 * t150 + (t242 + t187) * t149) * qJD(5)) * t150) * t213 + ((-t227 * t242 - t240) * t150 + (t184 + (t187 * t149 + (-t267 + t241) * t150) * qJD(5)) * t149) * t210 - t160 * ((t167 * t230 + t168 * t229) * t160 + ((t149 * t257 - t150 * t258) * t168 + (t149 * t259 + t150 * t260) * t167) * qJD(5)) / 0.2e1 + (t11 * t203 + t34 * ((t45 + t66) * t150 + (-t160 * t74 + t46) * t149) + t204 * t121 + ((-t160 * t29 - t12) * t150 + (t160 * t28 - t13) * t149) * t135 - (t28 * t93 - t29 * t94) * t160 - (t34 * (-t149 * t93 - t150 * t94) + t204 * t136) * qJD(5) + g(1) * t94 + g(2) * t93 - g(3) * t136) * m(6);];
tau = t1;
