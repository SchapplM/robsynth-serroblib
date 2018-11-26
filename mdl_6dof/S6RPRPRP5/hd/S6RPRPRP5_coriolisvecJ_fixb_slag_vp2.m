% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:59:16
% EndTime: 2018-11-23 15:59:26
% DurationCPUTime: 9.79s
% Computational Cost: add. (8339->524), mult. (22322->674), div. (0->0), fcn. (17021->8), ass. (0->226)
t332 = -Ifges(6,4) + Ifges(7,5);
t341 = t332 + Ifges(7,5);
t207 = sin(pkin(9));
t209 = cos(pkin(9));
t211 = sin(qJ(3));
t290 = cos(qJ(3));
t191 = t207 * t290 + t211 * t209;
t181 = t191 * qJD(1);
t206 = sin(pkin(10));
t208 = cos(pkin(10));
t158 = qJD(3) * t206 + t181 * t208;
t159 = qJD(3) * t208 - t181 * t206;
t210 = sin(qJ(5));
t212 = cos(qJ(5));
t104 = t210 * t158 - t159 * t212;
t235 = t290 * t209;
t250 = t211 * t207;
t216 = t235 - t250;
t184 = t216 * qJD(3);
t171 = qJD(1) * t184;
t188 = t206 * t210 - t212 * t208;
t67 = -qJD(5) * t104 - t171 * t188;
t313 = t67 / 0.2e1;
t190 = t206 * t212 + t208 * t210;
t338 = t212 * t158 + t159 * t210;
t68 = qJD(5) * t338 + t171 * t190;
t311 = t68 / 0.2e1;
t185 = t191 * qJD(3);
t172 = qJD(1) * t185;
t300 = t172 / 0.2e1;
t333 = Ifges(6,1) + Ifges(7,1);
t331 = Ifges(7,4) + Ifges(6,5);
t330 = Ifges(7,2) + Ifges(6,3);
t329 = -Ifges(6,6) + Ifges(7,6);
t200 = qJD(1) * t235;
t180 = qJD(1) * t250 - t200;
t129 = t190 * t180;
t183 = t190 * qJD(5);
t324 = t129 + t183;
t130 = t188 * t180;
t182 = t188 * qJD(5);
t248 = t130 + t182;
t297 = -t180 / 0.2e1;
t282 = pkin(7) + qJ(2);
t198 = t282 * t209;
t193 = qJD(1) * t198;
t178 = t211 * t193;
t196 = t282 * t207;
t192 = qJD(1) * t196;
t236 = t290 * t192;
t148 = -t178 - t236;
t142 = -qJD(3) * pkin(3) + qJD(4) - t148;
t237 = -pkin(2) * t209 - pkin(1);
t194 = qJD(1) * t237 + qJD(2);
t265 = Ifges(5,2) * t206;
t269 = Ifges(5,4) * t208;
t227 = -t265 + t269;
t270 = Ifges(5,4) * t206;
t228 = Ifges(5,1) * t208 - t270;
t229 = mrSges(5,1) * t206 + mrSges(5,2) * t208;
t263 = Ifges(4,5) * qJD(3);
t291 = t208 / 0.2e1;
t292 = -t206 / 0.2e1;
t121 = pkin(3) * t180 - qJ(4) * t181 + t194;
t149 = -t211 * t192 + t290 * t193;
t144 = qJD(3) * qJ(4) + t149;
t77 = t208 * t121 - t144 * t206;
t78 = t206 * t121 + t208 * t144;
t339 = (Ifges(5,1) * t158 + Ifges(5,4) * t159 + Ifges(5,5) * t180) * t291 + (Ifges(5,4) * t158 + Ifges(5,2) * t159 + Ifges(5,6) * t180) * t292 + t142 * t229 + t194 * mrSges(4,2) + (-t206 * t78 - t208 * t77) * mrSges(5,3) + t158 * t228 / 0.2e1 + t159 * t227 / 0.2e1 + t263 / 0.2e1;
t44 = pkin(4) * t180 - pkin(8) * t158 + t77;
t58 = pkin(8) * t159 + t78;
t14 = t210 * t44 + t212 * t58;
t257 = t171 * t208;
t101 = pkin(3) * t172 - qJ(4) * t171 - qJD(4) * t181;
t242 = qJD(1) * qJD(2);
t232 = t207 * t242;
t246 = qJD(2) * t200 - qJD(3) * t236;
t108 = -t211 * t232 + (qJD(4) - t178) * qJD(3) + t246;
t52 = t208 * t101 - t108 * t206;
t31 = pkin(4) * t172 - pkin(8) * t257 + t52;
t258 = t171 * t206;
t53 = t206 * t101 + t208 * t108;
t33 = -pkin(8) * t258 + t53;
t4 = -qJD(5) * t14 - t210 * t33 + t212 * t31;
t2 = -pkin(5) * t172 - t4;
t312 = -t68 / 0.2e1;
t214 = t191 * qJD(2);
t113 = qJD(1) * t214 + qJD(3) * t149;
t86 = pkin(4) * t258 + t113;
t9 = t68 * pkin(5) - t67 * qJ(6) - qJD(6) * t338 + t86;
t337 = mrSges(6,2) * t86 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t9 + Ifges(6,4) * t312 + 0.2e1 * t331 * t300 + t341 * t311 + 0.2e1 * t313 * t333;
t176 = qJD(5) + t180;
t13 = -t210 * t58 + t212 * t44;
t325 = qJD(6) - t13;
t11 = -pkin(5) * t176 + t325;
t12 = qJ(6) * t176 + t14;
t215 = Ifges(5,6) * t159;
t267 = Ifges(5,5) * t158;
t271 = Ifges(4,4) * t181;
t336 = t78 * mrSges(5,2) + t14 * mrSges(6,2) + t11 * mrSges(7,1) + Ifges(4,2) * t297 + Ifges(4,6) * qJD(3) + t271 / 0.2e1 - t77 * mrSges(5,1) - t13 * mrSges(6,1) - t12 * mrSges(7,3) - t194 * mrSges(4,1) - t215 / 0.2e1 - t267 / 0.2e1;
t335 = -Ifges(6,6) / 0.2e1;
t103 = Ifges(6,4) * t104;
t266 = Ifges(7,5) * t104;
t327 = t176 * t331 + t333 * t338 - t103 + t266;
t255 = t180 * t206;
t114 = -pkin(4) * t255 + t149;
t326 = pkin(5) * t324 + qJ(6) * t248 - qJD(6) * t190 - t114;
t323 = -t290 * t196 - t211 * t198;
t322 = t172 * t330 + t329 * t68 + t331 * t67;
t321 = (m(3) * qJ(2) + mrSges(3,3)) * (t207 ^ 2 + t209 ^ 2);
t252 = t191 * t208;
t147 = -pkin(3) * t216 - qJ(4) * t191 + t237;
t156 = -t211 * t196 + t198 * t290;
t92 = t208 * t147 - t156 * t206;
t73 = -pkin(4) * t216 - pkin(8) * t252 + t92;
t253 = t191 * t206;
t93 = t206 * t147 + t208 * t156;
t79 = -pkin(8) * t253 + t93;
t278 = t210 * t73 + t212 * t79;
t289 = pkin(8) * t208;
t116 = pkin(3) * t185 - qJ(4) * t184 - qJD(4) * t191;
t124 = t216 * qJD(2) + qJD(3) * t323;
t71 = t208 * t116 - t124 * t206;
t36 = pkin(4) * t185 - t184 * t289 + t71;
t254 = t184 * t206;
t72 = t206 * t116 + t208 * t124;
t51 = -pkin(8) * t254 + t72;
t8 = -qJD(5) * t278 - t210 * t51 + t212 * t36;
t243 = qJD(5) * t212;
t244 = qJD(5) * t210;
t3 = t210 * t31 + t212 * t33 + t44 * t243 - t244 * t58;
t1 = qJ(6) * t172 + qJD(6) * t176 + t3;
t320 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t316 = mrSges(6,1) * t86 + mrSges(7,1) * t9 - mrSges(7,2) * t1 - mrSges(6,3) * t3 + 0.2e1 * Ifges(7,3) * t311 - t67 * Ifges(6,4) / 0.2e1 + t172 * t335 + Ifges(7,6) * t300 + t341 * t313 + (-t312 + t311) * Ifges(6,2);
t305 = -t104 / 0.2e1;
t304 = t104 / 0.2e1;
t303 = -t338 / 0.2e1;
t302 = t338 / 0.2e1;
t299 = -t176 / 0.2e1;
t298 = t176 / 0.2e1;
t296 = t180 / 0.2e1;
t295 = -t181 / 0.2e1;
t281 = pkin(8) + qJ(4);
t45 = mrSges(6,1) * t172 - mrSges(6,3) * t67;
t46 = -t172 * mrSges(7,1) + t67 * mrSges(7,2);
t280 = t46 - t45;
t47 = -mrSges(6,2) * t172 - mrSges(6,3) * t68;
t48 = -mrSges(7,2) * t68 + mrSges(7,3) * t172;
t279 = t47 + t48;
t145 = pkin(3) * t181 + qJ(4) * t180;
t87 = t208 * t145 - t148 * t206;
t59 = pkin(4) * t181 + t180 * t289 + t87;
t88 = t206 * t145 + t208 * t148;
t75 = pkin(8) * t255 + t88;
t24 = t210 * t59 + t212 * t75;
t82 = -mrSges(7,2) * t104 + mrSges(7,3) * t176;
t273 = mrSges(6,3) * t104;
t83 = -mrSges(6,2) * t176 - t273;
t277 = t82 + t83;
t272 = mrSges(6,3) * t338;
t84 = mrSges(6,1) * t176 - t272;
t85 = -mrSges(7,1) * t176 + mrSges(7,2) * t338;
t276 = t85 - t84;
t275 = mrSges(4,3) * t180;
t274 = mrSges(4,3) * t181;
t268 = Ifges(6,4) * t338;
t264 = t181 * Ifges(4,1);
t261 = t113 * t323;
t120 = mrSges(5,1) * t258 + mrSges(5,2) * t257;
t247 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t159 - t158 * mrSges(5,2) - t274;
t203 = -pkin(4) * t208 - pkin(3);
t22 = t68 * mrSges(6,1) + t67 * mrSges(6,2);
t21 = t68 * mrSges(7,1) - t67 * mrSges(7,3);
t231 = t172 * mrSges(4,1) + t171 * mrSges(4,2);
t128 = pkin(4) * t253 - t323;
t226 = Ifges(5,5) * t208 - Ifges(5,6) * t206;
t225 = t206 * t53 + t208 * t52;
t224 = -t206 * t52 + t208 * t53;
t222 = t206 * t77 - t208 * t78;
t23 = -t210 * t75 + t212 * t59;
t27 = -t210 * t79 + t212 * t73;
t126 = -t180 * mrSges(5,2) + mrSges(5,3) * t159;
t127 = mrSges(5,1) * t180 - mrSges(5,3) * t158;
t218 = t126 * t208 - t127 * t206;
t195 = t281 * t206;
t197 = t281 * t208;
t217 = -t212 * t195 - t197 * t210;
t155 = -t195 * t210 + t197 * t212;
t7 = t210 * t36 + t212 * t51 + t73 * t243 - t244 * t79;
t97 = -pkin(4) * t159 + t142;
t125 = qJD(3) * t156 + t214;
t96 = pkin(4) * t254 + t125;
t174 = Ifges(4,4) * t180;
t163 = -qJD(3) * mrSges(4,2) - t275;
t146 = pkin(5) * t188 - qJ(6) * t190 + t203;
t141 = -t174 + t263 + t264;
t137 = t188 * t191;
t136 = t190 * t191;
t132 = mrSges(5,1) * t172 - mrSges(5,3) * t257;
t131 = -mrSges(5,2) * t172 - mrSges(5,3) * t258;
t123 = qJD(4) * t190 + qJD(5) * t155;
t122 = -qJD(4) * t188 + qJD(5) * t217;
t112 = (-qJD(3) * t193 - t232) * t211 + t246;
t102 = Ifges(7,5) * t338;
t95 = t172 * Ifges(5,5) + t171 * t228;
t94 = t172 * Ifges(5,6) + t171 * t227;
t89 = Ifges(5,3) * t180 + t215 + t267;
t81 = t184 * t190 + t243 * t252 - t244 * t253;
t80 = -t183 * t191 - t184 * t188;
t56 = mrSges(6,1) * t104 + mrSges(6,2) * t338;
t55 = mrSges(7,1) * t104 - mrSges(7,3) * t338;
t54 = pkin(5) * t338 + qJ(6) * t104;
t50 = pkin(5) * t136 + qJ(6) * t137 + t128;
t41 = -t104 * Ifges(6,2) + t176 * Ifges(6,6) + t268;
t40 = Ifges(7,4) * t338 + t176 * Ifges(7,2) + t104 * Ifges(7,6);
t39 = Ifges(6,5) * t338 - t104 * Ifges(6,6) + t176 * Ifges(6,3);
t38 = t176 * Ifges(7,6) + t104 * Ifges(7,3) + t102;
t29 = t104 * pkin(5) - qJ(6) * t338 + t97;
t26 = pkin(5) * t216 - t27;
t25 = -qJ(6) * t216 + t278;
t20 = -pkin(5) * t181 - t23;
t19 = qJ(6) * t181 + t24;
t10 = t81 * pkin(5) - t80 * qJ(6) + t137 * qJD(6) + t96;
t6 = -pkin(5) * t185 - t8;
t5 = qJ(6) * t185 - qJD(6) * t216 + t7;
t15 = [(t329 * t300 + t316) * t136 + (t331 * t298 + t333 * t302 + t327 / 0.2e1 + mrSges(6,2) * t97 + mrSges(7,2) * t11 - mrSges(6,3) * t13 - mrSges(7,3) * t29 + Ifges(6,4) * t305 + Ifges(7,5) * t304) * t80 + ((Ifges(4,1) + Ifges(5,1) * t208 ^ 2 / 0.2e1 + (-t269 + t265 / 0.2e1) * t206) * t171 - t172 * Ifges(4,4) - mrSges(5,3) * t225 + t226 * t300 + t291 * t95 + t292 * t94 + (mrSges(4,3) + t229) * t113) * t191 + (t329 * t298 + t332 * t302 - t12 * mrSges(7,2) - t14 * mrSges(6,3) + t38 / 0.2e1 + t29 * mrSges(7,1) - t41 / 0.2e1 + t97 * mrSges(6,1) + Ifges(7,3) * t304 - Ifges(6,2) * t305) * t81 + ((Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t180 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t176 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t338 + (Ifges(7,6) / 0.2e1 + t335) * t104 + t89 / 0.2e1 + t40 / 0.2e1 + t39 / 0.2e1 - t336) * t185 + (-t148 * t184 - t149 * t185 - t156 * t172 - t171 * t323) * mrSges(4,3) - t323 * t120 + 0.2e1 * t321 * t242 + (-mrSges(5,1) * t52 + mrSges(5,2) * t53 + mrSges(4,3) * t112 - Ifges(6,6) * t312 - Ifges(7,6) * t311 - t331 * t313 - t330 * t300 - (-Ifges(4,4) + t226) * t171 - t320 - (Ifges(5,3) + Ifges(4,2)) * t172 - t322 / 0.2e1) * t216 + m(5) * (t125 * t142 + t52 * t92 + t53 * t93 + t71 * t77 + t72 * t78 - t261) + m(4) * (t112 * t156 + t124 * t149 - t125 * t148 - t261) + (t184 * t297 + t185 * t295) * Ifges(4,4) + m(7) * (t1 * t25 + t10 * t29 + t11 * t6 + t12 * t5 + t2 * t26 + t50 * t9) + m(6) * (t128 * t86 + t13 * t8 + t14 * t7 + t27 * t4 + t278 * t3 + t96 * t97) + t278 * t47 - t337 * t137 + (t264 / 0.2e1 + t141 / 0.2e1 + t226 * t296 + t339) * t184 + t237 * t231 - t247 * t125 + t27 * t45 + t26 * t46 + t25 * t48 + t50 * t21 + t10 * t55 + t5 * t82 + t7 * t83 + t8 * t84 + t6 * t85 + t96 * t56 + t72 * t126 + t71 * t127 + t128 * t22 + t93 * t131 + t92 * t132 + t124 * t163; t206 * t131 + t208 * t132 + t279 * t190 + t280 * t188 + (t163 + t218) * t180 + (-t55 - t56 + t247) * t181 - m(4) * (-t148 * t181 - t149 * t180) + t231 - t248 * t277 + t324 * t276 - t321 * qJD(1) ^ 2 + (t1 * t190 + t11 * t324 - t12 * t248 - t181 * t29 + t188 * t2) * m(7) + (-t13 * t324 - t14 * t248 - t181 * t97 - t188 * t4 + t190 * t3) * m(6) + (-t142 * t181 - t180 * t222 + t225) * m(5); (-t174 + t141) * t296 + (-t129 * t329 + t130 * t331) * t299 + (-t182 * t331 + t183 * t329) * t298 + (-t129 * t332 + t130 * t333) * t303 + (-t182 * t333 + t183 * t332) * t302 + t327 * (-t130 / 0.2e1 - t182 / 0.2e1) + t326 * t55 + (-Ifges(4,2) * t296 + Ifges(6,6) * t304 + Ifges(7,6) * t305 + Ifges(5,3) * t297 + t330 * t299 + t331 * t303 + t336) * t181 + t224 * mrSges(5,3) + (Ifges(5,5) * t206 + Ifges(5,6) * t208 + t188 * t329) * t300 + (-mrSges(5,1) * t208 + mrSges(5,2) * t206 - mrSges(4,1)) * t113 + (-t11 * t248 - t12 * t324) * mrSges(7,2) + (t13 * t248 - t14 * t324) * mrSges(6,3) + (mrSges(7,1) * t324 + mrSges(7,3) * t248) * t29 + (mrSges(6,1) * t324 - mrSges(6,2) * t248) * t97 + (t131 * t208 - t132 * t206) * qJ(4) + (t1 * t155 + t146 * t9 - t217 * t2 + t326 * t29 + (t122 - t19) * t12 + (t123 - t20) * t11) * m(7) + (-t114 * t97 + t217 * t4 + t155 * t3 + t203 * t86 + (t122 - t24) * t14 + (-t123 - t23) * t13) * m(6) - t280 * t217 + (t41 - t38) * (-t129 / 0.2e1 - t183 / 0.2e1) + (-pkin(3) * t113 + qJ(4) * t224 - qJD(4) * t222 - t142 * t149 - t77 * t87 - t78 * t88) * m(5) + t218 * qJD(4) + t316 * t188 + (-t271 + t89 + t40 + t39) * t295 + t337 * t190 + (-Ifges(4,1) * t295 - t226 * t297 + t339) * t180 + (Ifges(6,4) * t130 - Ifges(7,5) * t182 + Ifges(6,2) * t129 + Ifges(7,3) * t183) * t304 + (-Ifges(6,4) * t182 + Ifges(7,5) * t130 - Ifges(6,2) * t183 - Ifges(7,3) * t129) * t305 - t19 * t82 - t24 * t83 - t23 * t84 - t20 * t85 + (t247 + t274) * t149 + (-t163 - t275) * t148 + t276 * t123 + t277 * t122 + t279 * t155 - t112 * mrSges(4,2) - t114 * t56 - pkin(3) * t120 + ((Ifges(5,1) * t206 + t269) * t291 + (Ifges(5,2) * t208 + t270) * t292 + Ifges(4,5)) * t171 - t88 * t126 - t87 * t127 + t146 * t21 - Ifges(4,6) * t172 + t203 * t22 + t206 * t95 / 0.2e1 + t94 * t291; -t276 * t338 + t277 * t104 - t159 * t126 + t158 * t127 + t120 + t21 + t22 + (t104 * t12 - t11 * t338 + t9) * m(7) + (t104 * t14 + t13 * t338 + t86) * m(6) + (t158 * t77 - t159 * t78 + t113) * m(5); t320 + (-t104 * t331 + t329 * t338) * t299 + (-Ifges(6,2) * t338 - t103 + t327) * t304 + (-t104 * t333 + t102 - t268 + t38) * t303 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t14 + t12 * t325 - t29 * t54) * m(7) + (t104 * t11 + t12 * t338) * mrSges(7,2) - pkin(5) * t46 + qJ(6) * t48 - t54 * t55 + qJD(6) * t82 + (t272 - t276) * t14 + (-t273 - t277) * t13 - t29 * (mrSges(7,1) * t338 + mrSges(7,3) * t104) - t97 * (mrSges(6,1) * t338 - mrSges(6,2) * t104) + t41 * t302 + (Ifges(7,3) * t338 - t266) * t305 + t322; t338 * t55 - t176 * t82 + 0.2e1 * (t2 / 0.2e1 + t29 * t302 + t12 * t299) * m(7) + t46;];
tauc  = t15(:);
