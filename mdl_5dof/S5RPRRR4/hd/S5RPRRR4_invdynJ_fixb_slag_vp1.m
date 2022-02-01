% Calculate vector of inverse dynamics joint torques for
% S5RPRRR4
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:34
% DurationCPUTime: 4.05s
% Computational Cost: add. (10808->359), mult. (6100->454), div. (0->0), fcn. (4528->10), ass. (0->218)
t176 = sin(qJ(5));
t178 = cos(qJ(5));
t272 = rSges(6,1) * t178;
t141 = -rSges(6,2) * t176 + t272;
t123 = t141 * qJD(5);
t270 = rSges(6,2) * t178;
t139 = rSges(6,1) * t176 + t270;
t173 = qJDD(1) + qJDD(3);
t165 = qJDD(4) + t173;
t174 = qJD(1) + qJD(3);
t168 = qJD(4) + t174;
t175 = qJ(1) + pkin(9);
t169 = qJ(3) + t175;
t160 = sin(t169);
t161 = cos(t169);
t172 = t174 ^ 2;
t166 = sin(t175);
t167 = cos(t175);
t180 = qJD(1) ^ 2;
t177 = sin(qJ(1));
t179 = cos(qJ(1));
t202 = (-qJDD(1) * t177 - t179 * t180) * pkin(1);
t186 = (-qJDD(1) * t166 - t167 * t180) * pkin(2) + t202;
t183 = (-t160 * t173 - t161 * t172) * pkin(3) + t186;
t162 = qJ(4) + t169;
t156 = cos(t162);
t241 = qJD(5) * t156;
t150 = t156 * pkin(8);
t155 = sin(t162);
t105 = pkin(4) * t155 - t150;
t254 = t155 * t176;
t129 = rSges(6,2) * t254;
t246 = t156 * rSges(6,3) + t129;
t253 = t155 * t178;
t77 = rSges(6,1) * t253 - t246;
t263 = -t105 - t77;
t289 = -t156 * pkin(4) - t155 * pkin(8);
t234 = qJD(5) * t270;
t251 = t156 * t176;
t235 = rSges(6,2) * t251;
t239 = qJD(5) * t176;
t233 = t168 * t235 + (rSges(6,1) * t239 + t234) * t155;
t250 = t156 * t178;
t290 = rSges(6,1) * t250 + t155 * rSges(6,3);
t48 = t168 * t290 - t233;
t240 = qJD(5) * t168;
t87 = -qJDD(5) * t156 + t155 * t240;
t12 = -t123 * t241 + t87 * t139 + t263 * t165 + (t168 * t289 - t48) * t168 + t183;
t295 = -g(1) + t12;
t102 = rSges(5,1) * t155 + rSges(5,2) * t156;
t252 = t156 * t168;
t255 = t155 * t168;
t81 = rSges(5,1) * t252 - rSges(5,2) * t255;
t294 = -t102 * t165 - t168 * t81 - g(1) + t183;
t117 = pkin(8) * t252;
t154 = pkin(3) * t161;
t159 = pkin(2) * t167;
t171 = t179 * pkin(1);
t163 = qJDD(1) * t171;
t280 = pkin(1) * t177;
t223 = -pkin(2) * t166 - t280;
t200 = qJDD(1) * t159 + t223 * t180 + t163;
t279 = pkin(3) * t160;
t192 = t173 * t154 - t172 * t279 + t200;
t242 = qJD(5) * t155;
t206 = rSges(6,3) * t252 + t168 * t129 - t156 * t234;
t231 = t156 * t239;
t47 = (-t168 * t253 - t231) * rSges(6,1) + t206;
t78 = -t235 + t290;
t58 = t78 - t289;
t86 = qJDD(5) * t155 + t156 * t240;
t13 = -t123 * t242 - t86 * t139 + (-pkin(4) * t255 + t117 + t47) * t168 + t58 * t165 + t192;
t293 = -g(2) + t13;
t103 = t156 * rSges(5,1) - rSges(5,2) * t155;
t94 = t168 * t102;
t292 = t103 * t165 - t168 * t94 - g(2) + t192;
t153 = t161 * rSges(4,1);
t112 = -rSges(4,2) * t160 + t153;
t111 = rSges(4,1) * t160 + rSges(4,2) * t161;
t247 = t174 * t111;
t298 = t112 * t173 - t174 * t247 - g(2) + t200;
t249 = t160 * t174;
t132 = rSges(4,2) * t249;
t248 = t161 * t174;
t99 = rSges(4,1) * t248 - t132;
t297 = t111 * t173 + t174 * t99 + g(1) - t186;
t236 = pkin(3) * t249;
t288 = t159 + t171;
t207 = t288 * qJD(1);
t237 = pkin(3) * t248;
t194 = t207 + t237;
t259 = t103 * t168;
t60 = t194 + t259;
t296 = t60 * (-t94 - t236);
t66 = t168 * t77;
t291 = -t168 * t105 - t66;
t208 = t223 * qJD(1);
t195 = t208 - t236;
t59 = t195 - t94;
t114 = t167 * rSges(3,1) - rSges(3,2) * t166;
t108 = t114 + t171;
t170 = Icges(6,4) * t178;
t212 = -Icges(6,2) * t176 + t170;
t137 = Icges(6,1) * t176 + t170;
t214 = t139 * t242 - t58 * t168;
t134 = Icges(6,5) * t178 - Icges(6,6) * t176;
t133 = Icges(6,5) * t176 + Icges(6,6) * t178;
t196 = Icges(6,3) * t168 - t133 * qJD(5);
t204 = t212 * t156;
t74 = Icges(6,6) * t155 + t204;
t265 = t176 * t74;
t261 = Icges(6,4) * t176;
t138 = Icges(6,1) * t178 - t261;
t205 = t138 * t156;
t76 = Icges(6,5) * t155 + t205;
t215 = -t178 * t76 + t265;
t286 = -t134 * t255 + t196 * t156 + t215 * t168;
t203 = t134 * t156;
t73 = Icges(6,4) * t253 - Icges(6,2) * t254 - Icges(6,6) * t156;
t266 = t176 * t73;
t128 = Icges(6,4) * t254;
t75 = Icges(6,1) * t253 - Icges(6,5) * t156 - t128;
t216 = -t178 * t75 + t266;
t285 = t196 * t155 + (t203 + t216) * t168;
t135 = Icges(6,2) * t178 + t261;
t210 = t135 * t176 - t137 * t178;
t284 = t134 * qJD(5) + t210 * t168;
t71 = Icges(6,5) * t253 - Icges(6,6) * t254 - Icges(6,3) * t156;
t24 = -t216 * t155 - t156 * t71;
t274 = -Icges(6,2) * t253 - t128 + t75;
t276 = t137 * t155 + t73;
t283 = -t274 * t176 - t276 * t178;
t282 = t86 / 0.2e1;
t281 = t87 / 0.2e1;
t278 = -t155 * t71 - t75 * t250;
t72 = Icges(6,3) * t155 + t203;
t277 = t155 * t72 + t76 * t250;
t275 = -t137 * t156 - t74;
t273 = -t135 * t156 + t76;
t70 = t112 * t174 + t207;
t269 = t111 * t70;
t232 = t139 * t241;
t190 = t195 - t232;
t30 = t263 * t168 + t190;
t268 = t156 * t30;
t257 = t133 * t156;
t53 = -t210 * t155 - t257;
t264 = t53 * t168;
t258 = t133 * t155;
t256 = t134 * t168;
t245 = -t135 + t138;
t244 = t137 + t212;
t238 = m(3) + m(4) + m(5);
t230 = -pkin(4) - t272;
t229 = -t242 / 0.2e1;
t228 = t242 / 0.2e1;
t227 = -t241 / 0.2e1;
t226 = t241 / 0.2e1;
t61 = t76 * t253;
t225 = t156 * t72 - t61;
t224 = -t71 + t265;
t85 = t103 + t154;
t142 = rSges(2,1) * t179 - rSges(2,2) * t177;
t140 = rSges(2,1) * t177 + rSges(2,2) * t179;
t113 = rSges(3,1) * t166 + rSges(3,2) * t167;
t25 = -t74 * t254 - t225;
t220 = t155 * t25 - t156 * t24;
t26 = -t73 * t251 - t278;
t27 = -t74 * t251 + t277;
t219 = t27 * t155 - t26 * t156;
t31 = t194 - t214;
t218 = -t155 * t31 - t268;
t217 = t155 * t77 + t156 * t78;
t45 = t176 * t75 + t178 * t73;
t46 = t176 * t76 + t178 * t74;
t211 = t135 * t178 + t137 * t176;
t209 = -t232 + t291;
t84 = -t102 - t279;
t201 = -t273 * t176 + t275 * t178;
t57 = t230 * t155 + t150 + t246;
t56 = t154 + t58;
t69 = t208 - t247;
t199 = (-t244 * t176 + t245 * t178) * t168;
t198 = Icges(6,5) * t168 - qJD(5) * t137;
t197 = Icges(6,6) * t168 - t135 * qJD(5);
t55 = t57 - t279;
t193 = -rSges(6,1) * t231 + t117 + t206;
t54 = -t210 * t156 + t258;
t52 = t54 * t168;
t10 = t219 * qJD(5) + t52;
t121 = t212 * qJD(5);
t122 = t138 * qJD(5);
t42 = t197 * t155 + t168 * t204;
t44 = t198 * t155 + t168 * t205;
t16 = -t216 * qJD(5) + t176 * t44 + t178 * t42;
t41 = t197 * t156 - t212 * t255;
t43 = -t138 * t255 + t198 * t156;
t17 = -t215 * qJD(5) + t176 * t43 + t178 * t41;
t185 = -t211 * qJD(5) - t121 * t176 + t122 * t178 + t133 * t168;
t20 = t155 * t284 + t185 * t156;
t21 = t185 * t155 - t156 * t284;
t9 = t220 * qJD(5) + t264;
t191 = (t52 + ((t25 - t61 + (t72 + t266) * t156 + t278) * t156 + t277 * t155) * qJD(5)) * t226 + (-t210 * qJD(5) + t121 * t178 + t122 * t176) * t168 + (t46 + t54) * t282 + (t45 + t53) * t281 + (-t264 + ((t224 * t156 + t27 - t277) * t156 + (t224 * t155 + t225 + t26) * t155) * qJD(5) + t9) * t229 + (t17 + t20) * t228 + (Icges(5,3) + t211) * t165 + (t16 + t21 + t10) * t227;
t189 = Icges(4,3) * t173 + t191;
t188 = -t45 * qJD(5) + t168 * t71 - t176 * t42 + t178 * t44;
t187 = -t46 * qJD(5) + t168 * t72 - t176 * t41 + t178 * t43;
t184 = (t230 * t268 + (t30 * (-rSges(6,3) - pkin(8)) + t31 * t230) * t155) * t168;
t182 = t59 * (-t81 - t237) + t296;
t181 = t30 * (t233 - t237) + t31 * (t193 - t236) + t184;
t96 = t139 * t156;
t95 = t139 * t155;
t34 = t217 * qJD(5) + qJD(2);
t11 = t77 * t86 - t78 * t87 + qJDD(2) + (t155 * t48 + t156 * t47) * qJD(5);
t6 = t187 * t155 - t156 * t286;
t5 = t188 * t155 - t156 * t285;
t4 = t155 * t286 + t187 * t156;
t3 = t155 * t285 + t188 * t156;
t1 = [t189 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + ((t60 * t223 - t288 * t59) * qJD(1) + t182 + t292 * (t85 + t288) + t294 * (t223 + t84)) * m(5) + (t69 * t132 + (-t69 * t153 - t269) * t174 + (t70 * t223 - t288 * t69) * qJD(1) + t298 * (t112 + t288) - t297 * (-t111 + t223)) * m(4) + ((qJDD(1) * t114 - g(2) + t163) * t108 + (-g(1) - qJDD(1) * t113 + t202 + (-0.2e1 * t114 + 0.2e1 * t108 - t171) * t180) * (-t113 - t280)) * m(3) + ((t140 ^ 2 + t142 ^ 2) * qJDD(1) + g(1) * t140 - g(2) * t142) * m(2) + (-(-t30 + t190 + t291) * t31 + (t31 * t223 - t288 * t30) * qJD(1) + t181 + t293 * (t56 + t288) + t295 * (t223 + t55)) * m(6); m(6) * t11 + t238 * qJDD(2) + (-m(6) - t238) * g(3); t189 + (t181 - t30 * (t214 - t237) - t31 * (t209 - t236) + t293 * t56 + t295 * t55) * m(6) + (t182 - t59 * (-t237 - t259) - t296 + t292 * t85 + t294 * t84) * m(5) + (t269 * t174 - t69 * t99 - t70 * t247 + (t69 * t174 + t298) * t112 + t297 * t111) * m(4); t191 + (t184 + t293 * t58 + t295 * t57 + (t193 - t209) * t31 + (-t214 + t233) * t30) * m(6) + (-t59 * t81 - t60 * t94 + (t59 * t168 + t292) * t103 + (t60 * t168 - t294) * t102) * m(5); t10 * t252 / 0.2e1 + t155 * (t54 * t165 + t20 * t168 + t26 * t87 + t27 * t86 + (t155 * t4 - t156 * t3) * qJD(5)) / 0.2e1 + t219 * t282 + ((t168 * t27 - t3) * t156 + (t168 * t26 + t4) * t155) * t228 + t9 * t255 / 0.2e1 - t156 * (t53 * t165 + t21 * t168 + t24 * t87 + t25 * t86 + (t155 * t6 - t156 * t5) * qJD(5)) / 0.2e1 + t220 * t281 + ((t168 * t25 - t5) * t156 + (t168 * t24 + t6) * t155) * t227 + t165 * (t46 * t155 - t45 * t156) / 0.2e1 + t168 * ((t168 * t46 - t16) * t156 + (t168 * t45 + t17) * t155) / 0.2e1 + ((-t242 * t257 + t256) * t155 + (t199 + (-t283 * t156 + (t258 + t201) * t155) * qJD(5)) * t156) * t229 + ((-t241 * t258 - t256) * t156 + (t199 + (t201 * t155 + (-t283 + t257) * t156) * qJD(5)) * t155) * t226 - t168 * ((t245 * t176 + t244 * t178) * t168 + ((t273 * t155 - t274 * t156) * t178 + (t275 * t155 + t276 * t156) * t176) * qJD(5)) / 0.2e1 + (t11 * t217 + t34 * ((t47 + t66) * t156 + (-t168 * t78 + t48) * t155) + t218 * t123 + ((-t168 * t31 - t12) * t156 + (t168 * t30 - t13) * t155) * t139 - (t30 * t95 - t31 * t96) * t168 - (t34 * (-t155 * t95 - t156 * t96) + t218 * t141) * qJD(5) + g(1) * t96 + g(2) * t95 - g(3) * t141) * m(6);];
tau = t1;
