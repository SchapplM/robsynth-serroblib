% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:58
% EndTime: 2019-03-09 03:55:09
% DurationCPUTime: 6.02s
% Computational Cost: add. (9487->458), mult. (21249->625), div. (0->0), fcn. (15095->8), ass. (0->228)
t187 = qJD(3) + qJD(5);
t344 = t187 / 0.2e1;
t190 = sin(pkin(10));
t196 = cos(qJ(3));
t271 = cos(pkin(10));
t238 = t271 * t196;
t193 = sin(qJ(3));
t259 = qJD(1) * t193;
t152 = -qJD(1) * t238 + t190 * t259;
t161 = -t190 * t196 - t193 * t271;
t153 = t161 * qJD(1);
t192 = sin(qJ(5));
t195 = cos(qJ(5));
t218 = -t152 * t195 + t153 * t192;
t343 = t218 / 0.2e1;
t166 = pkin(3) * t259 + qJD(1) * qJ(2) + qJD(4);
t115 = -pkin(4) * t153 + t166;
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t299 = pkin(8) * t152;
t197 = -pkin(1) - pkin(7);
t172 = qJD(1) * t197 + qJD(2);
t145 = -qJ(4) * t259 + t172 * t193;
t126 = t190 * t145;
t164 = t196 * t172;
t258 = qJD(1) * t196;
t146 = -qJ(4) * t258 + t164;
t130 = qJD(3) * pkin(3) + t146;
t89 = t130 * t271 - t126;
t72 = qJD(3) * pkin(4) + t299 + t89;
t298 = pkin(8) * t153;
t239 = t271 * t145;
t90 = t130 * t190 + t239;
t75 = t90 + t298;
t44 = t192 * t72 + t195 * t75;
t40 = pkin(9) * t187 + t44;
t236 = t152 * t192 + t153 * t195;
t47 = -pkin(5) * t236 - pkin(9) * t218 + t115;
t15 = -t191 * t40 + t194 * t47;
t16 = t191 * t47 + t194 * t40;
t98 = qJD(6) - t236;
t312 = -t98 / 0.2e1;
t85 = t187 * t191 + t194 * t218;
t314 = -t85 / 0.2e1;
t84 = t187 * t194 - t191 * t218;
t315 = -t84 / 0.2e1;
t280 = Ifges(6,2) * t236;
t333 = t280 / 0.2e1;
t201 = -t115 * mrSges(6,1) - t15 * mrSges(7,1) + t16 * mrSges(7,2) + 0.2e1 * Ifges(6,4) * t343 + 0.2e1 * Ifges(7,5) * t314 + 0.2e1 * Ifges(6,6) * t344 + 0.2e1 * Ifges(7,6) * t315 + 0.2e1 * Ifges(7,3) * t312 + t333;
t342 = t44 * mrSges(6,3) + t201;
t162 = -t190 * t193 + t238;
t109 = t161 * t192 + t162 * t195;
t257 = qJD(3) * t193;
t154 = -qJD(3) * t238 + t190 * t257;
t155 = t161 * qJD(3);
t203 = qJD(5) * t109 - t154 * t195 + t155 * t192;
t43 = -t192 * t75 + t195 * t72;
t217 = t161 * t195 - t162 * t192;
t69 = qJD(5) * t217 + t154 * t192 + t155 * t195;
t341 = t203 * t44 + t43 * t69;
t340 = Ifges(6,1) * t343;
t244 = t271 * pkin(3);
t180 = t244 + pkin(4);
t300 = pkin(3) * t190;
t149 = t180 * t192 + t195 * t300;
t95 = -t146 * t190 - t239;
t214 = t95 - t298;
t96 = t146 * t271 - t126;
t81 = t96 + t299;
t330 = -qJD(5) * t149 + t192 * t81 - t195 * t214;
t251 = qJD(1) * qJD(3);
t139 = t162 * t251;
t140 = qJD(1) * t155;
t62 = qJD(5) * t236 - t139 * t192 + t140 * t195;
t41 = qJD(6) * t84 + t194 * t62;
t63 = qJD(5) * t218 + t139 * t195 + t140 * t192;
t19 = mrSges(7,1) * t63 - mrSges(7,3) * t41;
t42 = -qJD(6) * t85 - t191 * t62;
t20 = -mrSges(7,2) * t63 + mrSges(7,3) * t42;
t221 = -t191 * t19 + t194 * t20;
t252 = qJD(6) * t194;
t253 = qJD(6) * t191;
t53 = -mrSges(7,2) * t98 + mrSges(7,3) * t84;
t54 = mrSges(7,1) * t98 - mrSges(7,3) * t85;
t339 = -t252 * t54 - t253 * t53 + t221;
t337 = Ifges(6,5) * t344 + Ifges(6,4) * t236 / 0.2e1;
t336 = t340 + t337;
t225 = Ifges(7,5) * t194 - Ifges(7,6) * t191;
t283 = Ifges(7,4) * t194;
t227 = -Ifges(7,2) * t191 + t283;
t284 = Ifges(7,4) * t191;
t229 = Ifges(7,1) * t194 - t284;
t305 = t194 / 0.2e1;
t306 = -t191 / 0.2e1;
t313 = t85 / 0.2e1;
t304 = Ifges(7,4) * t85;
t35 = Ifges(7,2) * t84 + Ifges(7,6) * t98 + t304;
t83 = Ifges(7,4) * t84;
t36 = Ifges(7,1) * t85 + Ifges(7,5) * t98 + t83;
t335 = t98 * t225 / 0.2e1 + t229 * t313 + t84 * t227 / 0.2e1 + t36 * t305 + t35 * t306;
t290 = -mrSges(6,1) * t187 - mrSges(7,1) * t84 + mrSges(7,2) * t85 + mrSges(6,3) * t218;
t220 = -t191 * t54 + t194 * t53;
t87 = -mrSges(6,2) * t187 + mrSges(6,3) * t236;
t216 = -t220 - t87;
t329 = qJ(2) * (m(3) + m(4));
t240 = t196 * t251;
t256 = qJD(3) * t196;
t242 = t172 * t256;
t255 = qJD(4) * t193;
t204 = t242 + (-qJ(4) * t256 - t255) * qJD(1);
t243 = t172 * t257;
t254 = qJD(4) * t196;
t205 = -t243 + (qJ(4) * t257 - t254) * qJD(1);
t79 = -t190 * t204 + t205 * t271;
t200 = -t140 * pkin(8) + t79;
t80 = t190 * t205 + t204 * t271;
t67 = -pkin(8) * t139 + t80;
t12 = qJD(5) * t43 + t192 * t200 + t195 * t67;
t165 = pkin(3) * t240 + qJD(1) * qJD(2);
t110 = pkin(4) * t139 + t165;
t25 = pkin(5) * t63 - pkin(9) * t62 + t110;
t2 = qJD(6) * t15 + t12 * t194 + t191 * t25;
t3 = -qJD(6) * t16 - t12 * t191 + t194 * t25;
t327 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t41 + Ifges(7,6) * t42;
t295 = t3 * t191;
t326 = -t15 * t252 - t16 * t253 - t295;
t66 = pkin(5) * t218 - pkin(9) * t236;
t223 = t15 * t194 + t16 * t191;
t230 = mrSges(7,1) * t191 + mrSges(7,2) * t194;
t39 = -pkin(5) * t187 - t43;
t213 = t39 * t230;
t323 = t115 * mrSges(6,2) + t213 + t336;
t324 = -t223 * mrSges(7,3) + t323 + t335 + t337;
t321 = qJD(1) ^ 2;
t319 = t41 / 0.2e1;
t318 = t42 / 0.2e1;
t316 = t63 / 0.2e1;
t310 = -t152 / 0.2e1;
t308 = t154 / 0.2e1;
t307 = t155 / 0.2e1;
t13 = qJD(5) * t44 + t192 * t67 - t195 * t200;
t260 = qJ(4) - t197;
t167 = t260 * t193;
t168 = t260 * t196;
t113 = t167 * t190 - t168 * t271;
t93 = -pkin(8) * t162 + t113;
t114 = -t167 * t271 - t168 * t190;
t94 = pkin(8) * t161 + t114;
t51 = t192 * t94 - t195 * t93;
t297 = t13 * t51;
t296 = t194 * t2;
t294 = t43 * mrSges(6,3);
t292 = t62 * mrSges(6,3);
t291 = t63 * mrSges(6,3);
t289 = mrSges(5,3) * t153;
t287 = Ifges(4,4) * t193;
t286 = Ifges(4,4) * t196;
t278 = t109 * t13;
t277 = t152 * mrSges(5,3);
t276 = t152 * Ifges(5,4);
t270 = Ifges(4,5) * qJD(3);
t269 = Ifges(4,6) * qJD(3);
t268 = t236 * t191;
t267 = t236 * t194;
t265 = t161 * t139;
t263 = t162 * t140;
t170 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t258;
t262 = t193 * t170;
t181 = pkin(3) * t193 + qJ(2);
t141 = t257 * t260 - t254;
t142 = -qJD(3) * t168 - t255;
t92 = t141 * t190 + t142 * t271;
t173 = pkin(3) * t256 + qJD(2);
t184 = pkin(3) * t258;
t241 = t63 * mrSges(6,1) + mrSges(6,2) * t62;
t119 = -pkin(4) * t152 + t184;
t237 = t139 * mrSges(5,1) + mrSges(5,2) * t140;
t91 = t141 * t271 - t142 * t190;
t131 = -pkin(4) * t161 + t181;
t116 = -pkin(4) * t154 + t173;
t233 = -t191 * t2 - t194 * t3;
t232 = mrSges(4,1) * t193 + mrSges(4,2) * t196;
t231 = mrSges(7,1) * t194 - mrSges(7,2) * t191;
t228 = Ifges(7,1) * t191 + t283;
t226 = Ifges(7,2) * t194 + t284;
t224 = Ifges(7,5) * t191 + Ifges(7,6) * t194;
t222 = t15 * t191 - t16 * t194;
t52 = t192 * t93 + t195 * t94;
t55 = -pkin(5) * t217 - pkin(9) * t109 + t131;
t27 = t191 * t55 + t194 * t52;
t26 = -t191 * t52 + t194 * t55;
t219 = -t191 * t53 - t194 * t54;
t148 = t180 * t195 - t192 * t300;
t215 = -pkin(8) * t155 + t91;
t212 = qJ(2) * (mrSges(4,1) * t196 - mrSges(4,2) * t193);
t208 = -qJD(6) * t223 - t295;
t207 = -t154 * t90 + t155 * t89 - t161 * t80 + t162 * t79;
t206 = t208 + t296;
t8 = Ifges(7,4) * t41 + Ifges(7,2) * t42 + Ifges(7,6) * t63;
t9 = Ifges(7,1) * t41 + Ifges(7,4) * t42 + t63 * Ifges(7,5);
t202 = -t12 * mrSges(6,2) + mrSges(7,3) * t296 + t228 * t319 + t226 * t318 + t224 * t316 + t191 * t9 / 0.2e1 - Ifges(6,6) * t63 + Ifges(6,5) * t62 + t8 * t305 + (-t231 - mrSges(6,1)) * t13 + (t213 + t335) * qJD(6);
t169 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t259;
t163 = t232 * qJD(1);
t157 = t270 + (Ifges(4,1) * t196 - t287) * qJD(1);
t156 = t269 + (-Ifges(4,2) * t193 + t286) * qJD(1);
t147 = Ifges(5,4) * t153;
t143 = -pkin(5) - t148;
t125 = qJD(3) * mrSges(5,1) + t277;
t124 = -qJD(3) * mrSges(5,2) + t289;
t105 = -mrSges(5,1) * t153 - mrSges(5,2) * t152;
t100 = -t152 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t147;
t99 = t153 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t276;
t78 = pkin(8) * t154 + t92;
t65 = -mrSges(6,1) * t236 + mrSges(6,2) * t218;
t59 = Ifges(7,3) * t63;
t48 = t119 + t66;
t46 = t192 * t214 + t195 * t81;
t28 = pkin(5) * t203 - pkin(9) * t69 + t116;
t24 = t191 * t66 + t194 * t43;
t23 = -t191 * t43 + t194 * t66;
t22 = qJD(5) * t52 + t192 * t78 - t195 * t215;
t21 = -qJD(5) * t51 + t192 * t215 + t195 * t78;
t18 = t191 * t48 + t194 * t46;
t17 = -t191 * t46 + t194 * t48;
t14 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t5 = -qJD(6) * t27 - t191 * t21 + t194 * t28;
t4 = qJD(6) * t26 + t191 * t28 + t194 * t21;
t1 = [m(5) * (t113 * t79 + t114 * t80 + t165 * t181 + t166 * t173 + t89 * t91 + t90 * t92) + (-t113 * t140 - t114 * t139 - t207) * mrSges(5,3) + (-t139 * t162 + t140 * t161 + t153 * t307 + t154 * t310) * Ifges(5,4) + (-t280 / 0.2e1 - t201) * t203 - (-t12 * mrSges(6,3) + t59 / 0.2e1 - Ifges(6,4) * t62 + t110 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t63 + t327) * t217 + (t153 * t308 - t265) * Ifges(5,2) + (t340 + t324) * t69 + (((2 * mrSges(3,3)) + t232 + 0.2e1 * t329) * qJD(2) + (0.2e1 * t212 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t193 * t196 + (0.3e1 / 0.2e1 * t193 ^ 2 - 0.3e1 / 0.2e1 * t196 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1) + t100 * t307 + t99 * t308 + t26 * t19 + t27 * t20 + t181 * t237 + (t51 * t62 - t52 * t63 - t341) * mrSges(6,3) + t131 * t241 + t51 * t14 + t4 * t53 + t5 * t54 + t21 * t87 + t116 * t65 + t92 * t124 + t91 * t125 + t290 * t22 + m(7) * (t15 * t5 + t16 * t4 + t2 * t27 + t22 * t39 + t26 * t3 + t297) + qJD(2) * t163 + m(6) * (t110 * t131 + t115 * t116 + t12 * t52 + t21 * t44 - t22 * t43 + t297) + t165 * (-mrSges(5,1) * t161 + mrSges(5,2) * t162) + t166 * (-mrSges(5,1) * t154 + mrSges(5,2) * t155) + t173 * t105 + (Ifges(5,5) * t307 + Ifges(5,6) * t308 + (t197 * t169 - t156 / 0.2e1 - t269 / 0.2e1) * t196 + (-t197 * t170 - t157 / 0.2e1 - t270 / 0.2e1) * t193) * qJD(3) + (t155 * t310 + t263) * Ifges(5,1) + (t8 * t306 + t9 * t305 + Ifges(6,1) * t62 - Ifges(6,4) * t63 + t110 * mrSges(6,2) + t229 * t319 + t227 * t318 + t225 * t316 + (mrSges(6,3) + t230) * t13 + t233 * mrSges(7,3) + (t36 * t306 - t194 * t35 / 0.2e1 + t39 * t231 + t226 * t315 + t228 * t314 + t224 * t312 + t222 * mrSges(7,3)) * qJD(6)) * t109; -t154 * t124 + t155 * t125 - t290 * t69 - (t14 + t292) * t109 + (t169 * t196 - t262) * qJD(3) + (-t263 + t265) * mrSges(5,3) - t216 * t203 + (-mrSges(3,3) - t329) * t321 - (qJD(6) * t219 + t221 - t291) * t217 + m(6) * (-t12 * t217 - t278 + t341) + m(5) * t207 + m(7) * (-t203 * t222 - t206 * t217 - t39 * t69 - t278) + (-m(5) * t166 - m(6) * t115 - m(7) * t223 - t105 - t163 + t219 - t65) * qJD(1); (-t139 * t300 - t140 * t244) * mrSges(5,3) + (-t115 * t119 + t12 * t149 + t330 * t43 - t44 * t46) * m(6) - (Ifges(5,2) * t152 + t100 + t147) * t153 / 0.2e1 + t202 + (-t196 * (-Ifges(4,1) * t193 - t286) / 0.2e1 + t193 * (-Ifges(4,2) * t196 - t287) / 0.2e1 - t212) * t321 + ((m(6) * t44 - m(7) * t222 - t216) * qJD(5) - t13 * m(6) - t292) * t148 + (t15 * t267 + t16 * t268 + t326) * mrSges(7,3) + ((t190 * t80 + t271 * t79) * pkin(3) - t166 * t184 - t89 * t95 - t90 * t96) * m(5) + (t333 + t342) * t218 + t89 * t289 + t172 * t262 - t251 * Ifges(4,5) * t193 / 0.2e1 - t330 * t290 + (m(7) * t206 + t339) * (pkin(9) + t149) + (t225 * t312 + t227 * t315 + t229 * t314 + t294 - t323 - t336) * t236 + (t13 * t143 - t15 * t17 - t16 * t18 - t330 * t39) * m(7) + t99 * t310 - Ifges(4,6) * t240 / 0.2e1 - mrSges(4,2) * t242 - mrSges(4,1) * t243 - t18 * t53 - t17 * t54 - t105 * t184 + t79 * mrSges(5,1) - t80 * mrSges(5,2) + t156 * t258 / 0.2e1 + t157 * t259 / 0.2e1 - t46 * t87 - t169 * t164 - t36 * t267 / 0.2e1 + t35 * t268 / 0.2e1 - t119 * t65 - t96 * t124 - t95 * t125 + t152 * (Ifges(5,1) * t153 + t276) / 0.2e1 - t90 * t277 - Ifges(5,6) * t139 + Ifges(5,5) * t140 + t143 * t14 - qJD(3) * (Ifges(5,5) * t153 + Ifges(5,6) * t152) / 0.2e1 - t149 * t291 - t166 * (-mrSges(5,1) * t152 + mrSges(5,2) * t153); t220 * qJD(6) - t290 * t218 + t216 * t236 - t153 * t124 - t152 * t125 + t194 * t19 + t191 * t20 + t237 + t241 + (-t218 * t39 - t222 * t98 - t233) * m(7) + (t218 * t43 - t236 * t44 + t110) * m(6) + (-t152 * t89 - t153 * t90 + t165) * m(5); t202 + t339 * pkin(9) + t342 * t218 + (t294 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t218 - t324) * t236 + t208 * mrSges(7,3) - t290 * t44 - pkin(5) * t14 - t24 * t53 - t23 * t54 - t43 * t87 + (-t15 * t23 - t16 * t24 - t39 * t44 + (t296 + t326) * pkin(9) - pkin(5) * t13) * m(7); t59 - t39 * (mrSges(7,1) * t85 + mrSges(7,2) * t84) + (Ifges(7,1) * t84 - t304) * t314 + t35 * t313 + (Ifges(7,5) * t84 - Ifges(7,6) * t85) * t312 - t15 * t53 + t16 * t54 + (t15 * t84 + t16 * t85) * mrSges(7,3) + (-Ifges(7,2) * t85 + t36 + t83) * t315 + t327;];
tauc  = t1(:);
