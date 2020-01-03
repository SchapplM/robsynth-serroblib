% Calculate vector of inverse dynamics joint torques for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:32
% DurationCPUTime: 2.82s
% Computational Cost: add. (4806->337), mult. (5043->433), div. (0->0), fcn. (3776->6), ass. (0->191)
t166 = qJD(1) + qJD(2);
t168 = sin(qJ(4));
t170 = cos(qJ(4));
t240 = Icges(5,4) * t170;
t130 = -Icges(5,2) * t168 + t240;
t196 = Icges(5,1) * t168 + t240;
t226 = t130 + t196;
t241 = Icges(5,4) * t168;
t132 = Icges(5,1) * t170 - t241;
t195 = Icges(5,2) * t170 + t241;
t227 = -t195 + t132;
t282 = (t168 * t226 - t170 * t227) * t166;
t169 = sin(qJ(1));
t243 = pkin(1) * qJD(1);
t218 = t169 * t243;
t167 = qJ(1) + qJ(2);
t161 = sin(t167);
t162 = cos(t167);
t109 = rSges(3,1) * t161 + rSges(3,2) * t162;
t236 = t109 * t166;
t85 = -t218 - t236;
t204 = rSges(5,1) * t168 + rSges(5,2) * t170;
t119 = t204 * qJD(4);
t246 = rSges(5,2) * t168;
t247 = rSges(5,1) * t170;
t139 = -t246 + t247;
t156 = t162 * pkin(6);
t164 = t166 ^ 2;
t165 = qJDD(1) + qJDD(2);
t171 = cos(qJ(1));
t172 = qJD(1) ^ 2;
t185 = (-qJDD(1) * t169 - t171 * t172) * pkin(1);
t223 = qJD(3) * t166;
t177 = qJDD(3) * t161 + t162 * t223 + t185;
t148 = t162 * qJ(3);
t107 = pkin(2) * t161 - t148;
t253 = pkin(6) * t161;
t84 = -t161 * rSges(5,3) + t162 * t204;
t209 = -t107 + t84 - t253;
t222 = qJD(4) * t161;
t154 = t162 * rSges(5,3);
t215 = qJD(4) * t246;
t216 = qJD(4) * t247;
t190 = (t215 - t216) * t162;
t49 = (t161 * t204 + t154) * t166 + t190;
t157 = t162 * pkin(2);
t110 = t161 * qJ(3) + t157;
t146 = qJD(3) * t162;
t74 = t110 * t166 - t146;
t220 = qJD(4) * t166;
t98 = qJDD(4) * t161 + t162 * t220;
t11 = -t164 * t156 - t119 * t222 + t139 * t98 + (-t49 - t74) * t166 + t209 * t165 + t177;
t281 = -g(1) + t11;
t230 = t162 * t166;
t136 = rSges(4,2) * t230;
t244 = rSges(4,3) * t162;
t108 = rSges(4,2) * t161 + t244;
t229 = -t107 + t108;
t233 = t161 * t166;
t23 = (-rSges(4,3) * t233 + t136 - t74) * t166 + t229 * t165 + t177;
t280 = -g(1) + t23;
t88 = rSges(3,1) * t230 - rSges(3,2) * t233;
t279 = -t109 * t165 - t166 * t88 - g(1) + t185;
t163 = t171 * pkin(1);
t254 = pkin(1) * t169;
t207 = qJDD(1) * t163 - t172 * t254;
t126 = qJ(3) * t230;
t145 = qJD(3) * t161;
t228 = t126 + t145;
t184 = t161 * t223 + t207 + t166 * (-pkin(2) * t233 + t228) + t165 * t110;
t214 = t161 * t216 + t204 * t230;
t50 = (-rSges(5,3) * t166 - t215) * t161 + t214;
t231 = t161 * t170;
t232 = t161 * t168;
t83 = rSges(5,1) * t232 + rSges(5,2) * t231 + t154;
t99 = qJDD(4) * t162 - t161 * t220;
t12 = -t164 * t253 - t139 * t99 + t165 * t83 + t166 * t50 + (pkin(6) * t165 + qJD(4) * t119 - qJDD(3)) * t162 + t184;
t278 = -g(2) + t12;
t111 = -rSges(4,2) * t162 + t161 * rSges(4,3);
t225 = rSges(4,2) * t233 + rSges(4,3) * t230;
t24 = -qJDD(3) * t162 + t165 * t111 + t166 * t225 + t184;
t277 = -g(2) + t24;
t112 = t162 * rSges(3,1) - rSges(3,2) * t161;
t276 = t112 * t165 - t166 * t236 - g(2) + t207;
t76 = t110 + t111;
t275 = t166 * t76;
t73 = t166 * t84;
t274 = -t161 * t215 + t214 + t228 - t73;
t194 = Icges(5,5) * t168 + Icges(5,6) * t170;
t273 = t194 * t166;
t80 = -Icges(5,6) * t161 + t162 * t195;
t82 = -Icges(5,5) * t161 + t162 * t196;
t197 = t168 * t82 + t170 * t80;
t272 = t197 * t162;
t271 = -t166 * t108 + t225;
t270 = t156 + t83 + t110;
t40 = t168 * t80 - t170 * t82;
t187 = t195 * t166;
t266 = -Icges(5,6) * t166 + qJD(4) * t130;
t43 = t161 * t187 - t266 * t162;
t188 = t196 * t166;
t265 = -Icges(5,5) * t166 + qJD(4) * t132;
t45 = t161 * t188 - t265 * t162;
t78 = -Icges(5,3) * t161 + t162 * t194;
t269 = qJD(4) * t40 + t166 * t78 + t168 * t45 + t170 * t43;
t221 = qJD(4) * t162;
t268 = -t139 * t221 + t166 * t270;
t128 = Icges(5,5) * t170 - Icges(5,6) * t168;
t267 = -Icges(5,3) * t166 + qJD(4) * t128;
t115 = t195 * qJD(4);
t116 = t196 * qJD(4);
t192 = t130 * t168 - t132 * t170;
t264 = qJD(4) * t192 + t115 * t170 + t116 * t168 + t128 * t166;
t79 = Icges(5,6) * t162 + t161 * t195;
t141 = Icges(5,4) * t231;
t81 = Icges(5,1) * t232 + Icges(5,5) * t162 + t141;
t198 = t168 * t79 - t170 * t81;
t44 = t266 * t161 + t162 * t187;
t46 = t265 * t161 + t162 * t188;
t77 = Icges(5,3) * t162 + t161 * t194;
t263 = qJD(4) * t198 + t166 * t77 - t168 * t46 - t170 * t44;
t249 = t130 * t162 + t82;
t251 = -t132 * t162 + t80;
t261 = t168 * t251 - t170 * t249;
t250 = -Icges(5,2) * t232 + t141 + t81;
t252 = -t132 * t161 + t79;
t260 = t168 * t252 - t170 * t250;
t259 = t98 / 0.2e1;
t258 = t99 / 0.2e1;
t257 = -pkin(2) - pkin(6);
t256 = t161 / 0.2e1;
t255 = -t162 / 0.2e1;
t193 = t170 * t130 + t168 * t132;
t89 = t161 * t128;
t55 = t162 * t193 - t89;
t242 = t55 * t166;
t234 = t128 * t162;
t102 = t166 * t107;
t224 = t145 - t102;
t219 = -rSges(5,3) + t257;
t25 = t162 * t77 + t79 * t231 + t81 * t232;
t26 = -t162 * t78 - t80 * t231 - t82 * t232;
t217 = t171 * t243;
t213 = -t222 / 0.2e1;
t212 = t222 / 0.2e1;
t211 = -t221 / 0.2e1;
t210 = t221 / 0.2e1;
t206 = t145 - t218;
t205 = -t146 + t217;
t140 = rSges(2,1) * t171 - rSges(2,2) * t169;
t138 = rSges(2,1) * t169 + rSges(2,2) * t171;
t203 = t161 * t26 + t162 * t25;
t199 = t168 * t81 + t170 * t79;
t68 = t161 * t77;
t27 = -t162 * t199 + t68;
t28 = -t161 * t78 + t272;
t202 = t161 * t28 + t162 * t27;
t104 = t139 * t222;
t191 = t104 + t206;
t31 = t166 * t209 + t191;
t32 = t205 + t268;
t201 = t161 * t31 - t162 * t32;
t200 = -t161 * t83 - t162 * t84;
t189 = t146 - t190;
t75 = t244 + t148 + (rSges(4,2) - pkin(2)) * t161;
t180 = t161 * t273 - t267 * t162 - t166 * t197;
t179 = t267 * t161 + t162 * t273 + t166 * t199;
t178 = -t194 * qJD(4) + t166 * t193;
t60 = t161 * t257 + t148 + t84;
t56 = t166 * t229 + t206;
t57 = t205 + t275;
t175 = (-t56 * t157 + (t56 * (-rSges(4,3) - qJ(3)) - t57 * pkin(2)) * t161) * t166;
t10 = qJD(4) * t202 - t242;
t15 = -qJD(4) * t199 - t168 * t44 + t170 * t46;
t16 = qJD(4) * t197 - t168 * t43 + t170 * t45;
t19 = t178 * t161 + t264 * t162;
t20 = -t264 * t161 + t178 * t162;
t54 = t161 * t193 + t234;
t51 = t54 * t166;
t9 = qJD(4) * t203 + t51;
t174 = (t51 + ((-t27 + t68 + t26) * t161 + (t28 - t272 + (-t199 + t78) * t161 + t25) * t162) * qJD(4)) * t213 + t40 * t259 - t98 * t55 / 0.2e1 + (-qJD(4) * t193 + t115 * t168 - t116 * t170) * t166 + (-t198 + t54) * t258 + (t242 + (t161 ^ 2 * t78 + (-t68 + t26 + (t199 + t78) * t162) * t162) * qJD(4) + t10) * t211 + (t15 + t20) * t210 + (t16 + t19 + t9) * t212 + (Icges(3,3) + Icges(4,1) - t192) * t165;
t173 = (t31 * t219 * t162 + (t31 * (-qJ(3) - t204) + t32 * t219) * t161) * t166;
t96 = t139 * t162;
t95 = t139 * t161;
t86 = t112 * t166 + t217;
t38 = t200 * qJD(4);
t6 = t269 * t161 + t180 * t162;
t5 = -t263 * t161 + t179 * t162;
t4 = t180 * t161 - t269 * t162;
t3 = t179 * t161 + t263 * t162;
t1 = [Icges(2,3) * qJDD(1) + t174 + (t276 * (t112 + t163) + t279 * (-t109 - t254) + (-t88 - t217 + t86) * t85) * m(3) + (g(1) * t138 - g(2) * t140 + (t138 ^ 2 + t140 ^ 2) * qJDD(1)) * m(2) + (t31 * (t189 - t217) + t173 + t278 * (t163 + t270) + t281 * (t60 - t254) + (pkin(6) * t233 + t102 - t191 - t218 + t274 + t31) * t32) * m(5) + (t56 * (t136 - t205) + t175 + t277 * (t163 + t76) + t280 * (t75 - t254) + (t126 + t102 + t56 + t271) * t57) * m(4); t174 + (t173 + t278 * t270 + t281 * t60 + (t253 * t166 - t104 - t224 + t274) * t32 + (-t146 + t189 + t268) * t31) * m(5) + (t175 + t277 * t76 + t280 * t75 + (t228 - t224 + t271) * t57 + (t136 + t275) * t56) * m(4) + (-t85 * t88 - t86 * t236 + (t85 * t166 + t276) * t112 + (t86 * t166 - t279) * t109) * m(3); (-m(4) - m(5)) * (g(1) * t161 - g(2) * t162) + 0.2e1 * (t11 * t256 + t12 * t255) * m(5) + 0.2e1 * (t23 * t256 + t24 * t255) * m(4); -t9 * t233 / 0.2e1 + t162 * (t165 * t54 + t166 * t20 + t25 * t99 + t26 * t98 + (t161 * t6 + t162 * t5) * qJD(4)) / 0.2e1 + t203 * t258 + ((t166 * t26 + t5) * t162 + (-t166 * t25 + t6) * t161) * t210 + t10 * t230 / 0.2e1 + (-t165 * t55 + t166 * t19 + t27 * t99 + t28 * t98 + (t161 * t4 + t162 * t3) * qJD(4)) * t256 + t202 * t259 + ((t166 * t28 + t3) * t162 + (-t166 * t27 + t4) * t161) * t212 + t165 * (t161 * t40 - t162 * t198) / 0.2e1 + t166 * ((t166 * t40 + t15) * t162 + (t166 * t198 + t16) * t161) / 0.2e1 + ((t89 * t221 - t273) * t162 + (-t282 + (t261 * t161 + (-t260 - t234) * t162) * qJD(4)) * t161) * t211 + ((-t222 * t234 - t273) * t161 + (t282 + (t260 * t162 + (-t261 + t89) * t161) * qJD(4)) * t162) * t213 - t166 * ((-t168 * t227 - t170 * t226) * t166 + ((t161 * t251 - t162 * t252) * t170 + (t161 * t249 - t162 * t250) * t168) * qJD(4)) / 0.2e1 + ((-t83 * t98 - t84 * t99 + (-t161 * t50 + t162 * t49) * qJD(4)) * t200 + t38 * ((-t166 * t83 + t49) * t162 + (-t50 + t73) * t161) - t201 * t119 + ((t166 * t31 - t12) * t162 + (t166 * t32 + t11) * t161) * t139 - (t31 * t96 + t32 * t95) * t166 - (t38 * (-t161 * t95 - t162 * t96) - t201 * t204) * qJD(4) - g(1) * t95 + g(2) * t96 + g(3) * t204) * m(5);];
tau = t1;
