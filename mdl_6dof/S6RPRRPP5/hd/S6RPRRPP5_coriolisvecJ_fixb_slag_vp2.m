% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2018-11-23 16:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:13:36
% EndTime: 2018-11-23 16:13:43
% DurationCPUTime: 6.60s
% Computational Cost: add. (5629->496), mult. (14959->599), div. (0->0), fcn. (10579->6), ass. (0->226)
t347 = Ifges(7,4) + Ifges(6,5);
t337 = Ifges(6,4) + Ifges(5,5);
t348 = -Ifges(7,5) + t337;
t345 = Ifges(7,2) + Ifges(6,3);
t344 = Ifges(5,6) - Ifges(6,6);
t346 = Ifges(6,6) - Ifges(7,6);
t339 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t174 = sin(pkin(9));
t175 = cos(pkin(9));
t177 = sin(qJ(3));
t179 = cos(qJ(3));
t160 = t174 * t179 + t175 * t177;
t152 = t160 * qJD(1);
t176 = sin(qJ(4));
t178 = cos(qJ(4));
t129 = qJD(3) * t176 + t152 * t178;
t285 = pkin(7) + qJ(2);
t165 = t285 * t174;
t161 = qJD(1) * t165;
t166 = t285 * t175;
t162 = qJD(1) * t166;
t119 = -t161 * t177 + t162 * t179;
t312 = t160 * qJD(2);
t71 = qJD(1) * t312 + qJD(3) * t119;
t331 = -t174 * t177 + t179 * t175;
t153 = t331 * qJD(3);
t140 = qJD(1) * t153;
t191 = t178 * qJD(3) - t152 * t176;
t84 = qJD(4) * t191 + t140 * t178;
t343 = -qJ(5) * t84 - qJD(5) * t129 + t71;
t342 = t347 * t178;
t341 = t347 * t176;
t118 = -t179 * t161 - t177 * t162;
t151 = t331 * qJD(1);
t236 = -pkin(2) * t175 - pkin(1);
t163 = qJD(1) * t236 + qJD(2);
t145 = qJD(4) - t151;
t301 = pkin(4) + pkin(5);
t114 = qJD(3) * pkin(8) + t119;
t96 = -pkin(3) * t151 - pkin(8) * t152 + t163;
t40 = -t176 * t114 + t178 * t96;
t28 = qJ(6) * t129 + t40;
t322 = qJD(5) - t28;
t16 = -t145 * t301 + t322;
t139 = t145 * qJ(5);
t41 = t178 * t114 + t176 * t96;
t29 = -qJ(6) * t191 + t41;
t18 = t139 + t29;
t195 = t41 * t176 + t40 * t178;
t321 = qJD(5) - t40;
t31 = -pkin(4) * t145 + t321;
t32 = t139 + t41;
t197 = t32 * t176 - t31 * t178;
t203 = Ifges(7,5) * t178 + Ifges(7,6) * t176;
t273 = Ifges(5,4) * t178;
t213 = -Ifges(5,2) * t176 + t273;
t220 = -mrSges(7,1) * t176 + mrSges(7,2) * t178;
t222 = mrSges(6,1) * t176 - mrSges(6,3) * t178;
t224 = mrSges(5,1) * t176 + mrSges(5,2) * t178;
t229 = qJD(3) * pkin(3) + t118;
t186 = qJ(5) * t129 + t229;
t27 = t191 * t301 + qJD(6) + t186;
t289 = t178 / 0.2e1;
t291 = t176 / 0.2e1;
t292 = -t176 / 0.2e1;
t293 = t145 / 0.2e1;
t294 = -t145 / 0.2e1;
t297 = t129 / 0.2e1;
t299 = -t191 / 0.2e1;
t300 = t191 / 0.2e1;
t127 = Ifges(5,4) * t191;
t330 = t347 * t191;
t311 = t339 * t129 + t348 * t145 + t127 - t330;
t334 = t347 * t129;
t327 = t145 * t346 - t345 * t191 + t334;
t274 = Ifges(5,4) * t176;
t328 = t178 * t339 - t274 + t341;
t332 = -t344 * t176 + t178 * t337;
t333 = t176 * t345 + t342;
t42 = -pkin(4) * t191 - t186;
t275 = Ifges(5,4) * t129;
t55 = Ifges(5,2) * t191 + t145 * Ifges(5,6) + t275;
t306 = t197 * mrSges(6,2) + t195 * mrSges(5,3) + (t16 * t178 - t18 * t176) * mrSges(7,3) + t229 * t224 - t203 * t294 - t213 * t300 - t220 * t27 - t222 * t42 - t292 * t55 - t333 * t299 - t332 * t293 - t327 * t291 - t328 * t297 - t311 * t289;
t338 = t152 * Ifges(4,1) / 0.2e1;
t340 = t163 * mrSges(4,2) - t118 * mrSges(4,3) + Ifges(4,4) * t151 + Ifges(4,5) * qJD(3) - t306 + t338;
t154 = t160 * qJD(3);
t141 = qJD(1) * t154;
t85 = qJD(4) * t129 + t140 * t176;
t336 = t141 * t346 + t345 * t85 + t347 * t84;
t335 = qJD(5) * t176 + t119;
t124 = -t165 * t177 + t166 * t179;
t98 = t124 * qJD(3) + t312;
t329 = (-Ifges(5,4) + t347) * t85 + t339 * t84 + t348 * t141;
t298 = -t129 / 0.2e1;
t111 = t176 * t118;
t115 = pkin(3) * t152 - pkin(8) * t151;
t284 = pkin(8) - qJ(6);
t168 = t284 * t178;
t326 = -t111 - (-qJ(6) * t151 - t115) * t178 + t301 * t152 + qJD(4) * t168 - qJD(6) * t176;
t247 = qJD(4) * t176;
t252 = qJ(6) * t176;
t48 = t176 * t115 + t178 * t118;
t37 = t152 * qJ(5) + t48;
t325 = -qJD(6) * t178 - t151 * t252 - t247 * t284 - t37;
t254 = qJ(5) * t178;
t190 = -t176 * t301 + t254;
t324 = t145 * t190 + t335;
t200 = pkin(4) * t176 - t254;
t323 = t145 * t200 - t335;
t320 = t154 * qJ(5) - qJD(5) * t331;
t319 = -t179 * t165 - t166 * t177;
t318 = t176 * t337 + t344 * t178;
t317 = -t178 * t345 + t341;
t100 = pkin(3) * t141 - pkin(8) * t140;
t246 = qJD(4) * t178;
t187 = t331 * qJD(2);
t70 = qJD(1) * t187 + t118 * qJD(3);
t10 = t100 * t178 - t114 * t246 - t176 * t70 - t96 * t247;
t9 = t176 * t100 - t114 * t247 + t178 * t70 + t96 * t246;
t316 = -t10 * t176 + t9 * t178;
t6 = t141 * qJ(5) + t145 * qJD(5) + t9;
t7 = -pkin(4) * t141 - t10;
t315 = t7 * t176 + t6 * t178;
t310 = t176 * t339 + t273 - t342;
t309 = (m(3) * qJ(2) + mrSges(3,3)) * (t174 ^ 2 + t175 ^ 2);
t255 = qJ(5) * t176;
t308 = -t178 * t301 - t255;
t239 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t240 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t241 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t263 = t151 * Ifges(4,2);
t305 = (-Ifges(7,6) / 0.2e1 - t239) * t191 - (-Ifges(7,5) / 0.2e1 + t241) * t129 - (Ifges(7,3) / 0.2e1 + t240) * t145 - t32 * mrSges(6,3) + t31 * mrSges(6,1) - t163 * mrSges(4,1) - t18 * mrSges(7,2) + t16 * mrSges(7,1) + Ifges(4,6) * qJD(3) + t119 * mrSges(4,3) + t41 * mrSges(5,2) - t40 * mrSges(5,1) + Ifges(6,6) * t300 + Ifges(7,5) * t297 + t152 * Ifges(4,4) + t263 / 0.2e1 + (Ifges(5,6) + Ifges(7,6)) * t299 + t337 * t298 + (Ifges(5,3) + Ifges(6,2) + Ifges(7,3)) * t294;
t304 = t84 / 0.2e1;
t303 = -t85 / 0.2e1;
t302 = t85 / 0.2e1;
t296 = -t141 / 0.2e1;
t295 = t141 / 0.2e1;
t290 = -t178 / 0.2e1;
t60 = mrSges(5,1) * t141 - mrSges(5,3) * t84;
t61 = -t141 * mrSges(6,1) + t84 * mrSges(6,2);
t283 = -t60 + t61;
t63 = -mrSges(5,2) * t141 - mrSges(5,3) * t85;
t64 = -mrSges(6,2) * t85 + mrSges(6,3) * t141;
t282 = t63 + t64;
t73 = -mrSges(6,1) * t191 - mrSges(6,3) * t129;
t74 = mrSges(7,1) * t191 + mrSges(7,2) * t129;
t281 = t73 - t74;
t90 = mrSges(6,2) * t191 + mrSges(6,3) * t145;
t91 = mrSges(7,2) * t145 - mrSges(7,3) * t191;
t280 = t90 + t91;
t277 = mrSges(5,3) * t191;
t92 = -mrSges(5,2) * t145 + t277;
t279 = t90 + t92;
t276 = mrSges(5,3) * t129;
t94 = mrSges(5,1) * t145 - t276;
t95 = -mrSges(6,1) * t145 + mrSges(6,2) * t129;
t278 = t94 - t95;
t266 = t319 * t71;
t259 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t191 + mrSges(5,2) * t129 + t152 * mrSges(4,3);
t256 = qJ(5) * t191;
t253 = qJ(6) * t160;
t117 = -pkin(3) * t331 - pkin(8) * t160 + t236;
t66 = t176 * t117 + t178 * t124;
t244 = qJD(5) * t178;
t97 = t319 * qJD(3) + t187;
t238 = t117 * t247 + t124 * t246 + t176 * t97;
t116 = pkin(3) * t154 - pkin(8) * t153;
t237 = t176 * t116 + t117 * t246 + t178 * t97;
t43 = -qJ(5) * t331 + t66;
t235 = qJD(4) * t253;
t35 = -t85 * mrSges(7,1) + t84 * mrSges(7,2);
t59 = -t141 * mrSges(7,1) - t84 * mrSges(7,3);
t234 = t141 * mrSges(4,1) + t140 * mrSges(4,2);
t47 = t115 * t178 - t111;
t121 = t176 * t124;
t65 = t117 * t178 - t121;
t1 = -qJ(6) * t84 - qJD(6) * t129 - t141 * t301 - t10;
t2 = qJ(6) * t85 - qJD(6) * t191 + t6;
t228 = -t1 * t178 + t176 * t2;
t227 = -t176 * t6 + t178 * t7;
t226 = -t10 * t178 - t176 * t9;
t225 = mrSges(5,1) * t178 - mrSges(5,2) * t176;
t223 = mrSges(6,1) * t178 + mrSges(6,3) * t176;
t221 = mrSges(7,1) * t178 + mrSges(7,2) * t176;
t212 = Ifges(5,2) * t178 + t274;
t202 = Ifges(7,5) * t176 - Ifges(7,6) * t178;
t201 = pkin(4) * t178 + t255;
t198 = t16 * t176 + t178 * t18;
t196 = t176 * t31 + t178 * t32;
t194 = t176 * t40 - t178 * t41;
t192 = -qJ(6) * t153 - qJD(6) * t160;
t15 = t116 * t178 - t238;
t14 = -t124 * t247 + t237;
t185 = t10 * mrSges(5,1) - t7 * mrSges(6,1) - t1 * mrSges(7,1) - t9 * mrSges(5,2) + t2 * mrSges(7,2) + t6 * mrSges(6,3);
t167 = t284 * t176;
t164 = -pkin(3) - t201;
t157 = pkin(3) - t308;
t137 = Ifges(6,2) * t141;
t136 = Ifges(5,3) * t141;
t131 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t151;
t93 = -mrSges(7,1) * t145 - mrSges(7,3) * t129;
t81 = Ifges(6,4) * t84;
t80 = Ifges(5,5) * t84;
t79 = Ifges(5,6) * t85;
t78 = Ifges(6,6) * t85;
t72 = pkin(4) * t129 - t256;
t67 = t160 * t200 - t319;
t62 = mrSges(7,2) * t141 + mrSges(7,3) * t85;
t46 = -t129 * t301 + t256;
t45 = t160 * t190 + t319;
t44 = pkin(4) * t331 - t65;
t39 = -pkin(4) * t152 - t47;
t36 = mrSges(5,1) * t85 + mrSges(5,2) * t84;
t34 = mrSges(6,1) * t85 - mrSges(6,3) * t84;
t33 = t160 * t252 + t43;
t26 = t121 + (-t117 - t253) * t178 + t301 * t331;
t22 = t84 * Ifges(5,4) - t85 * Ifges(5,2) + t141 * Ifges(5,6);
t17 = t200 * t153 + (qJD(4) * t201 - t244) * t160 + t98;
t13 = -pkin(4) * t154 - t15;
t12 = t190 * t153 + (t308 * qJD(4) + t244) * t160 - t98;
t11 = pkin(4) * t85 + t343;
t8 = t14 + t320;
t5 = -t301 * t85 - t343;
t4 = t178 * t235 + (-qJD(4) * t124 - t192) * t176 + t237 + t320;
t3 = t176 * t235 - t301 * t154 + (-t116 + t192) * t178 + t238;
t19 = [0.2e1 * t309 * qJD(2) * qJD(1) + (-t263 / 0.2e1 - t305) * t154 + (t213 * t303 + t203 * t296 + t5 * t220 + t11 * t222 + t22 * t292 + Ifges(4,1) * t140 - Ifges(4,4) * t141 + (mrSges(4,3) + t224) * t71 + t228 * mrSges(7,3) + t226 * mrSges(5,3) + t227 * mrSges(6,2) + (-t196 * mrSges(6,2) + t194 * mrSges(5,3) + t198 * mrSges(7,3) + t202 * t293 + t212 * t299 - t27 * t221 + t42 * t223 - t225 * t229 + t55 * t290 + t311 * t292 + t318 * t294 + t310 * t298 + t317 * t300) * qJD(4) + t333 * t302 + t332 * t295 + t336 * t291 + t328 * t304 + (t327 * qJD(4) + t329) * t289) * t160 + (-t124 * t141 - t140 * t319) * mrSges(4,3) - (t136 / 0.2e1 + t137 / 0.2e1 + t80 / 0.2e1 + t81 / 0.2e1 - t79 / 0.2e1 + t78 / 0.2e1 - t70 * mrSges(4,3) - Ifges(4,4) * t140 + (-Ifges(7,6) - t239) * t85 + (-Ifges(7,5) + t241) * t84 + (Ifges(4,2) + Ifges(7,3) + t240) * t141 + t185) * t331 - t319 * t36 + m(7) * (t1 * t26 + t12 * t27 + t16 * t3 + t18 * t4 + t2 * t33 + t45 * t5) + m(6) * (t11 * t67 + t13 * t31 + t17 * t42 + t32 * t8 + t43 * t6 + t44 * t7) + m(5) * (t10 * t65 + t14 * t41 + t15 * t40 - t229 * t98 + t66 * t9 - t266) + m(4) * (-t118 * t98 + t119 * t97 + t124 * t70 - t266) + (t338 + t340) * t153 + t236 * t234 + t259 * t98 + t45 * t35 + t26 * t59 + t44 * t61 + t33 * t62 + t43 * t64 + t65 * t60 + t66 * t63 + t67 * t34 + t17 * t73 + t12 * t74 + t8 * t90 + t4 * t91 + t14 * t92 + t3 * t93 + t15 * t94 + t13 * t95 + t97 * t131; -t151 * t131 + (-t259 - t281) * t152 + (-t283 - t59 + t145 * (t91 + t279)) * t178 + (t282 + t62 + t145 * (t93 - t278)) * t176 - m(4) * (-t118 * t152 + t119 * t151) + t234 - t309 * qJD(1) ^ 2 + (t145 * t198 + t152 * t27 + t228) * m(7) + (t145 * t196 - t152 * t42 - t227) * m(6) + (-t145 * t194 + t152 * t229 - t226) * m(5); t305 * t152 + (-t1 * t176 - t2 * t178) * mrSges(7,3) + t310 * t304 + (t11 * t164 - t31 * t39 - t32 * t37 + t323 * t42) * m(6) + t323 * t73 + t324 * t74 + t325 * t91 + (t1 * t167 + t157 * t5 + t326 * t16 + t168 * t2 + t325 * t18 + t324 * t27) * m(7) + t326 * t93 + t315 * mrSges(6,2) + ((-m(5) * t195 - m(6) * t197 - t176 * t279 - t178 * t278) * qJD(4) + t176 * t283 + t178 * t282 + m(5) * t316 + m(6) * t315) * pkin(8) + t316 * mrSges(5,3) + t317 * t302 + t318 * t295 + t5 * t221 - t11 * t223 + (-mrSges(4,1) - t225) * t71 + t336 * t290 + t329 * t291 - t306 * qJD(4) + (-pkin(3) * t71 + t119 * t229 - t40 * t47 - t41 * t48) * m(5) + ((Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t152 - t340) * t151 - t259 * t119 + t22 * t289 + t202 * t296 + t212 * t303 - pkin(3) * t36 - t70 * mrSges(4,2) - t37 * t90 - t48 * t92 - t47 * t94 - t39 * t95 - t118 * t131 + Ifges(4,5) * t140 - Ifges(4,6) * t141 + t157 * t35 + t164 * t34 + t167 * t59 + t168 * t62; (-pkin(4) * t7 + qJ(5) * t6 - t31 * t41 + t321 * t32 - t42 * t72) * m(6) + t185 + (t277 - t279) * t40 + t136 + t137 + t80 + t81 - t79 + t78 + (t62 + t64) * qJ(5) + (t191 * t339 - t275 + t327 + t334) * t298 + (-t129 * t344 + t337 * t191) * t294 + (t129 * t345 + t330) * t300 + (t129 * t32 - t191 * t31) * mrSges(6,2) + (-t129 * t18 + t16 * t191) * mrSges(7,3) + (Ifges(7,5) * t191 + Ifges(7,6) * t129) * t293 - t27 * (-mrSges(7,1) * t129 + mrSges(7,2) * t191) - t42 * (mrSges(6,1) * t129 - mrSges(6,3) * t191) + t229 * (mrSges(5,1) * t129 + mrSges(5,2) * t191) + (qJ(5) * t2 - t1 * t301 - t16 * t29 + t322 * t18 - t27 * t46) * m(7) - t301 * t59 + (-Ifges(5,2) * t129 + t127 + t311) * t299 + t55 * t297 + (t276 + t278) * t41 + t280 * qJD(5) - pkin(4) * t61 - t72 * t73 - t46 * t74 - Ifges(7,5) * t84 - Ifges(7,6) * t85 - t28 * t91 - t29 * t93 + Ifges(7,3) * t141; t281 * t129 - t280 * t145 + t59 + t61 + (-t129 * t27 - t145 * t18 + t1) * m(7) + (t129 * t42 - t145 * t32 + t7) * m(6); t191 * t91 + t129 * t93 + 0.2e1 * (t5 / 0.2e1 + t18 * t300 + t16 * t297) * m(7) + t35;];
tauc  = t19(:);
