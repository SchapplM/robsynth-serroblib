% Calculate vector of inverse dynamics joint torques for
% S5RPRRR9
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
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:19
% EndTime: 2019-12-31 19:07:38
% DurationCPUTime: 11.20s
% Computational Cost: add. (7965->510), mult. (19743->694), div. (0->0), fcn. (15125->14), ass. (0->240)
t258 = qJD(1) * qJD(2);
t164 = qJ(2) * qJDD(1) + t258;
t186 = sin(pkin(9));
t187 = cos(pkin(9));
t265 = t186 ^ 2 + t187 ^ 2;
t184 = pkin(9) + qJ(3);
t179 = qJ(4) + t184;
t170 = sin(t179);
t171 = cos(t179);
t189 = sin(qJ(5));
t299 = mrSges(6,2) * t189;
t364 = -t170 * t299 + t171 * (-m(6) * pkin(8) - mrSges(6,3));
t363 = m(5) + m(6);
t341 = t171 * pkin(4) + t170 * pkin(8);
t362 = m(6) * t341;
t193 = cos(qJ(5));
t185 = qJD(3) + qJD(4);
t190 = sin(qJ(4));
t194 = cos(qJ(4));
t303 = pkin(6) + qJ(2);
t160 = t303 * t186;
t150 = qJD(1) * t160;
t161 = t303 * t187;
t151 = qJD(1) * t161;
t191 = sin(qJ(3));
t195 = cos(qJ(3));
t114 = -t150 * t191 + t151 * t195;
t271 = t187 * t195;
t146 = -t186 * t191 + t271;
t137 = t146 * qJD(1);
t89 = pkin(7) * t137 + t114;
t278 = t194 * t89;
t274 = t151 * t191;
t113 = -t195 * t150 - t274;
t147 = t186 * t195 + t187 * t191;
t138 = t147 * qJD(1);
t88 = -pkin(7) * t138 + t113;
t86 = qJD(3) * pkin(3) + t88;
t60 = t190 * t86 + t278;
t57 = pkin(8) * t185 + t60;
t172 = pkin(2) * t187 + pkin(1);
t154 = -qJD(1) * t172 + qJD(2);
t115 = -pkin(3) * t137 + t154;
t214 = t137 * t190 + t194 * t138;
t233 = t194 * t137 - t138 * t190;
t61 = -pkin(4) * t233 - pkin(8) * t214 + t115;
t20 = -t189 * t57 + t193 * t61;
t361 = t20 * mrSges(6,1);
t21 = t189 * t61 + t193 * t57;
t360 = t21 * mrSges(6,2);
t90 = t185 * t193 - t189 * t214;
t91 = t185 * t189 + t193 * t214;
t345 = -mrSges(5,1) * t185 - mrSges(6,1) * t90 + mrSges(6,2) * t91 + mrSges(5,3) * t214;
t359 = -t171 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t170;
t112 = t146 * t190 + t147 * t194;
t260 = qJD(5) * t193;
t139 = t146 * qJD(3);
t140 = t147 * qJD(3);
t213 = t194 * t146 - t147 * t190;
t77 = qJD(4) * t213 + t139 * t194 - t140 * t190;
t211 = t112 * t260 + t189 * t77;
t301 = mrSges(6,1) * t193;
t358 = t299 - t301;
t291 = Ifges(5,6) * t185;
t295 = Ifges(5,4) * t214;
t305 = t60 * mrSges(5,3);
t101 = qJD(5) - t233;
t290 = Ifges(6,3) * t101;
t316 = Ifges(6,6) * t90;
t317 = Ifges(6,5) * t91;
t39 = t290 + t316 + t317;
t343 = Ifges(5,2) * t233;
t71 = t291 + t295 + t343;
t356 = -t115 * mrSges(5,1) - t361 + t360 + t305 - t39 / 0.2e1 + t71 / 0.2e1 + t295 / 0.2e1 + t291 / 0.2e1;
t100 = Ifges(5,4) * t233;
t223 = mrSges(6,1) * t189 + mrSges(6,2) * t193;
t282 = t190 * t89;
t59 = t194 * t86 - t282;
t56 = -pkin(4) * t185 - t59;
t209 = t56 * t223;
t87 = Ifges(6,4) * t90;
t41 = Ifges(6,1) * t91 + Ifges(6,5) * t101 + t87;
t280 = t193 * t41;
t292 = Ifges(5,5) * t185;
t306 = t59 * mrSges(5,3);
t320 = t189 / 0.2e1;
t318 = Ifges(6,4) * t91;
t40 = Ifges(6,2) * t90 + Ifges(6,6) * t101 + t318;
t344 = Ifges(5,1) * t214;
t72 = t100 + t292 + t344;
t355 = -t115 * mrSges(5,2) - t209 + t40 * t320 - t280 / 0.2e1 + t306 - t72 / 0.2e1 - t292 / 0.2e1 - t100 / 0.2e1;
t240 = m(3) * qJ(2) + mrSges(3,3);
t353 = -m(4) * t303 + mrSges(2,2) - mrSges(4,3) - mrSges(5,3) - t240;
t262 = qJD(4) * t194;
t263 = qJD(4) * t190;
t109 = qJD(1) * t139 + qJDD(1) * t147;
t234 = pkin(6) * qJDD(1) + t164;
t131 = t234 * t186;
t132 = t234 * t187;
t76 = -qJD(3) * t114 - t195 * t131 - t132 * t191;
t53 = qJDD(3) * pkin(3) - pkin(7) * t109 + t76;
t110 = -qJD(1) * t140 + qJDD(1) * t146;
t264 = qJD(3) * t195;
t75 = -qJD(3) * t274 - t191 * t131 + t195 * t132 - t150 * t264;
t58 = pkin(7) * t110 + t75;
t15 = t190 * t53 + t194 * t58 + t86 * t262 - t263 * t89;
t180 = qJDD(3) + qJDD(4);
t12 = pkin(8) * t180 + t15;
t54 = qJD(4) * t233 + t109 * t194 + t110 * t190;
t55 = -qJD(4) * t214 - t109 * t190 + t110 * t194;
t153 = -qJDD(1) * t172 + qJDD(2);
t92 = -pkin(3) * t110 + t153;
t19 = -pkin(4) * t55 - pkin(8) * t54 + t92;
t2 = qJD(5) * t20 + t12 * t193 + t189 * t19;
t3 = -qJD(5) * t21 - t12 * t189 + t19 * t193;
t352 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t177 = sin(t184);
t178 = cos(t184);
t226 = mrSges(4,1) * t178 - mrSges(4,2) * t177;
t227 = -mrSges(3,1) * t187 + mrSges(3,2) * t186;
t351 = m(3) * pkin(1) + m(4) * t172 + mrSges(2,1) + t226 - t227 - t359;
t66 = -mrSges(6,2) * t101 + mrSges(6,3) * t90;
t67 = mrSges(6,1) * t101 - mrSges(6,3) * t91;
t217 = -t189 * t67 + t193 * t66;
t93 = -mrSges(5,2) * t185 + mrSges(5,3) * t233;
t342 = -t217 - t93;
t117 = -t191 * t160 + t195 * t161;
t340 = t171 * t358 + t359;
t261 = qJD(5) * t189;
t339 = -t20 * t260 - t21 * t261;
t16 = -qJD(4) * t60 - t190 * t58 + t194 * t53;
t34 = qJD(5) * t90 + t180 * t189 + t193 * t54;
t51 = qJDD(5) - t55;
t17 = mrSges(6,1) * t51 - mrSges(6,3) * t34;
t35 = -qJD(5) * t91 + t180 * t193 - t189 * t54;
t18 = -mrSges(6,2) * t51 + mrSges(6,3) * t35;
t337 = -t189 * t17 + t193 * t18 - t67 * t260 - t66 * t261;
t74 = pkin(4) * t214 - pkin(8) * t233;
t334 = m(6) * pkin(4);
t333 = t34 / 0.2e1;
t332 = t35 / 0.2e1;
t329 = t51 / 0.2e1;
t326 = -t90 / 0.2e1;
t325 = -t91 / 0.2e1;
t324 = t91 / 0.2e1;
t323 = -t101 / 0.2e1;
t321 = t138 / 0.2e1;
t315 = pkin(3) * t138;
t314 = pkin(3) * t140;
t313 = pkin(3) * t177;
t168 = pkin(3) * t178;
t312 = pkin(3) * t190;
t311 = pkin(3) * t194;
t308 = t189 * t3;
t307 = t193 * t2;
t300 = mrSges(5,2) * t171;
t298 = mrSges(6,3) * t189;
t297 = mrSges(6,3) * t193;
t296 = Ifges(4,4) * t138;
t294 = Ifges(6,4) * t189;
t293 = Ifges(6,4) * t193;
t289 = t113 * mrSges(4,3);
t288 = t114 * mrSges(4,3);
t277 = t112 * t189;
t276 = t112 * t193;
t192 = sin(qJ(1));
t270 = t189 * t192;
t196 = cos(qJ(1));
t269 = t189 * t196;
t268 = t192 * t193;
t267 = t193 * t196;
t257 = qJDD(1) * t186;
t256 = qJDD(1) * t187;
t255 = Ifges(6,5) * t34 + Ifges(6,6) * t35 + Ifges(6,3) * t51;
t244 = t280 / 0.2e1;
t238 = -t55 * mrSges(5,1) + t54 * mrSges(5,2);
t237 = -t261 / 0.2e1;
t236 = -t110 * mrSges(4,1) + t109 * mrSges(4,2);
t116 = -t195 * t160 - t161 * t191;
t231 = t364 * t192;
t230 = t364 * t196;
t228 = -mrSges(3,1) * t256 + mrSges(3,2) * t257;
t222 = Ifges(6,1) * t193 - t294;
t221 = -Ifges(6,2) * t189 + t293;
t220 = Ifges(6,5) * t193 - Ifges(6,6) * t189;
t219 = t21 * t189 + t20 * t193;
t218 = -t189 * t20 + t193 * t21;
t98 = -pkin(7) * t147 + t116;
t99 = pkin(7) * t146 + t117;
t69 = t190 * t98 + t194 * t99;
t123 = -pkin(3) * t146 - t172;
t70 = -pkin(4) * t213 - pkin(8) * t112 + t123;
t31 = t189 * t70 + t193 * t69;
t30 = -t189 * t69 + t193 * t70;
t68 = t190 * t99 - t194 * t98;
t210 = t112 * t261 - t193 * t77;
t208 = t90 * t221;
t207 = t91 * t222;
t206 = t101 * t220;
t202 = -qJD(5) * t219 - t308;
t95 = -t160 * t264 + qJD(2) * t271 + (-qJD(2) * t186 - qJD(3) * t161) * t191;
t201 = m(6) * (-pkin(4) * t170 - t313) - t170 * t301;
t200 = t300 + (mrSges(5,1) + t301 + t334) * t170;
t96 = -t147 * qJD(2) - qJD(3) * t117;
t13 = -pkin(4) * t180 - t16;
t8 = t34 * Ifges(6,4) + t35 * Ifges(6,2) + t51 * Ifges(6,6);
t9 = t34 * Ifges(6,1) + t35 * Ifges(6,4) + t51 * Ifges(6,5);
t199 = -t15 * mrSges(5,2) + t193 * t8 / 0.2e1 + t2 * t297 + t9 * t320 + t13 * t358 + t16 * mrSges(5,1) + Ifges(5,3) * t180 + (Ifges(6,1) * t189 + t293) * t333 + (Ifges(6,2) * t193 + t294) * t332 + t40 * t237 + (Ifges(6,5) * t189 + Ifges(6,6) * t193) * t329 + Ifges(5,6) * t55 + Ifges(5,5) * t54 + (t244 + t209) * qJD(5) + (t208 + t207 + t206) * qJD(5) / 0.2e1;
t198 = -pkin(7) * t139 + t96;
t181 = -pkin(7) - t303;
t176 = -qJDD(1) * pkin(1) + qJDD(2);
t174 = -pkin(4) - t311;
t152 = t168 + t172;
t133 = Ifges(4,4) * t137;
t130 = t171 * t267 + t270;
t129 = -t171 * t269 + t268;
t128 = -t171 * t268 + t269;
t127 = t171 * t270 + t267;
t119 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t138;
t118 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t137;
t103 = t138 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t133;
t102 = t137 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t296;
t82 = -pkin(7) * t140 + t95;
t78 = qJD(4) * t112 + t139 * t190 + t194 * t140;
t73 = -mrSges(5,1) * t233 + mrSges(5,2) * t214;
t65 = t315 + t74;
t63 = t194 * t88 - t282;
t62 = t190 * t88 + t278;
t45 = -mrSges(5,2) * t180 + mrSges(5,3) * t55;
t44 = mrSges(5,1) * t180 - mrSges(5,3) * t54;
t36 = pkin(4) * t78 - pkin(8) * t77 + t314;
t27 = t189 * t74 + t193 * t59;
t26 = -t189 * t59 + t193 * t74;
t25 = t189 * t65 + t193 * t63;
t24 = -t189 * t63 + t193 * t65;
t22 = -qJD(4) * t68 + t190 * t198 + t194 * t82;
t11 = -mrSges(6,1) * t35 + mrSges(6,2) * t34;
t5 = -qJD(5) * t31 - t189 * t22 + t193 * t36;
t4 = qJD(5) * t30 + t189 * t36 + t193 * t22;
t1 = [t78 * t361 + (-m(5) * t16 + m(6) * t13 + t11 - t44) * t68 + qJD(3) * (Ifges(4,5) * t139 - Ifges(4,6) * t140) / 0.2e1 + t137 * (Ifges(4,4) * t139 - Ifges(4,2) * t140) / 0.2e1 + t154 * (mrSges(4,1) * t140 + mrSges(4,2) * t139) + (Ifges(4,1) * t139 - Ifges(4,4) * t140) * t321 + t233 * (Ifges(5,4) * t77 - Ifges(5,2) * t78) / 0.2e1 + t214 * (Ifges(5,1) * t77 - Ifges(5,4) * t78) / 0.2e1 + m(5) * (t115 * t314 + t123 * t92 + t15 * t69 + t22 * t60) + (mrSges(4,2) * t153 - mrSges(4,3) * t76 + Ifges(4,1) * t109 + Ifges(4,4) * t110 + Ifges(4,5) * qJDD(3)) * t147 + (-t128 * mrSges(6,1) - t127 * mrSges(6,2) + (t181 * t363 + t353) * t196 + (m(5) * t152 - m(6) * (-t152 - t341) + t351) * t192) * g(1) + (-t130 * mrSges(6,1) - t129 * mrSges(6,2) - t363 * (t196 * t152 - t181 * t192) + t353 * t192 + (-t351 - t362) * t196) * g(2) - t78 * t360 + (-m(5) * t59 + m(6) * t56 + t345) * (qJD(4) * t69 + t190 * t82 - t194 * t198) - t211 * t40 / 0.2e1 + (-Ifges(6,6) * t332 - Ifges(6,5) * t333 + Ifges(5,6) * t180 + Ifges(5,2) * t55 + Ifges(5,4) * t54 - t92 * mrSges(5,1) - Ifges(6,3) * t329 - t255 / 0.2e1 + t15 * mrSges(5,3) + t352) * t213 + t56 * (mrSges(6,1) * t211 - mrSges(6,2) * t210) + t90 * (-Ifges(6,4) * t210 - Ifges(6,2) * t211 + Ifges(6,6) * t78) / 0.2e1 + t101 * (-Ifges(6,5) * t210 - Ifges(6,6) * t211 + Ifges(6,3) * t78) / 0.2e1 - t139 * t289 + m(6) * (t2 * t31 + t20 * t5 + t21 * t4 + t3 * t30) + t123 * t238 - t172 * t236 + (t92 * mrSges(5,2) - t16 * mrSges(5,3) + Ifges(5,1) * t54 + Ifges(5,4) * t55 + Ifges(5,5) * t180 + t13 * t223 + t220 * t329 + t221 * t332 + t222 * t333 + t237 * t41) * t112 + t176 * t227 + m(4) * (t113 * t96 + t114 * t95 + t116 * t76 + t117 * t75 - t153 * t172) + t77 * t72 / 0.2e1 + t78 * t39 / 0.2e1 - t78 * t71 / 0.2e1 + t69 * t45 + (-t2 * t277 + t20 * t210 - t21 * t211 - t3 * t276) * mrSges(6,3) + t4 * t66 + t5 * t67 + t30 * t17 + t31 * t18 + 0.2e1 * t265 * t164 * mrSges(3,3) + (-mrSges(4,1) * t153 + mrSges(4,3) * t75 + Ifges(4,4) * t109 + Ifges(4,2) * t110 + Ifges(4,6) * qJDD(3)) * t146 + t22 * t93 - t78 * t305 - t77 * t306 - pkin(1) * t228 + t73 * t314 + t115 * (mrSges(5,1) * t78 + mrSges(5,2) * t77) + t116 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t109) + t117 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t110) + t95 * t118 + t96 * t119 + t139 * t103 / 0.2e1 - t140 * t102 / 0.2e1 + t77 * t244 + m(3) * (-pkin(1) * t176 + (t164 + t258) * qJ(2) * t265) + (Ifges(3,4) * t186 + Ifges(3,2) * t187) * t256 + (Ifges(3,1) * t186 + Ifges(3,4) * t187) * t257 - t140 * t288 + (-Ifges(6,1) * t210 - Ifges(6,4) * t211 + Ifges(6,5) * t78) * t324 + t185 * (Ifges(5,5) * t77 - Ifges(5,6) * t78) / 0.2e1 + t9 * t276 / 0.2e1 - t8 * t277 / 0.2e1 + Ifges(2,3) * qJDD(1); -t345 * t214 + m(3) * t176 + t238 + t228 + t342 * t233 + t217 * qJD(5) + t236 - t137 * t118 + t138 * t119 + t189 * t18 + t193 * t17 + (-g(1) * t192 + g(2) * t196) * (m(3) + m(4) + t363) - t240 * t265 * qJD(1) ^ 2 + (t101 * t218 + t2 * t189 + t3 * t193 - t214 * t56) * m(6) + (t214 * t59 - t233 * t60 + t92) * m(5) + (t113 * t138 - t114 * t137 + t153) * m(4); (m(5) * t313 + mrSges(4,1) * t177 + mrSges(5,1) * t170 + mrSges(4,2) * t178 + t300) * (g(1) * t196 + g(2) * t192) + (-t226 - m(6) * (t168 + t341) - m(5) * t168 + t340) * g(3) - (-Ifges(4,2) * t138 + t103 + t133) * t137 / 0.2e1 + (t13 * t174 - t20 * t24 - t21 * t25 - t56 * t62) * m(6) + (-t115 * t315 + t59 * t62 - t60 * t63) * m(5) + (t345 * t263 - t342 * t262 + (t190 * t56 + t194 * t218) * qJD(4) * m(6) + (t15 * t190 + t16 * t194 + (-t190 * t59 + t194 * t60) * qJD(4)) * m(5)) * pkin(3) + (t20 * t297 + t21 * t298 - t344 / 0.2e1 + t220 * t323 + t222 * t325 + t221 * t326 + t355) * t233 + (t343 / 0.2e1 + Ifges(6,3) * t323 + Ifges(6,5) * t325 + Ifges(6,6) * t326 + t356) * t214 - t345 * t62 + t199 - t75 * mrSges(4,2) + t76 * mrSges(4,1) - t25 * t66 - t24 * t67 + (m(6) * (t202 + t307) + t337) * (pkin(8) + t312) - t138 * (Ifges(4,1) * t137 - t296) / 0.2e1 - t3 * t298 - t63 * t93 + t339 * mrSges(6,3) - g(1) * (t196 * t201 - t230) - g(2) * (t192 * t201 - t231) + Ifges(4,5) * t109 - t73 * t315 + Ifges(4,6) * t110 - t113 * t118 + t114 * t119 - qJD(3) * (Ifges(4,5) * t137 - Ifges(4,6) * t138) / 0.2e1 - t154 * (mrSges(4,1) * t138 + mrSges(4,2) * t137) + t174 * t11 + t138 * t288 + t137 * t289 + t44 * t311 + t45 * t312 + t102 * t321 + Ifges(4,3) * qJDD(3); (t192 * t200 + t231) * g(2) + (-t316 / 0.2e1 - t317 / 0.2e1 - t290 / 0.2e1 + t356) * t214 - m(6) * (t20 * t26 + t21 * t27 + t56 * t60) + (t340 - t362) * g(3) - t345 * t60 + t199 - t27 * t66 - t26 * t67 - pkin(4) * t11 + (m(6) * (t307 - t308 + t339) + t337) * pkin(8) - t13 * t334 + (t196 * t200 + t230) * g(1) + (-t208 / 0.2e1 - t207 / 0.2e1 - t206 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t214 + t219 * mrSges(6,3) + t355) * t233 + t202 * mrSges(6,3) - t59 * t93; -t56 * (mrSges(6,1) * t91 + mrSges(6,2) * t90) + (Ifges(6,1) * t90 - t318) * t325 + t40 * t324 + (Ifges(6,5) * t90 - Ifges(6,6) * t91) * t323 - t20 * t66 + t21 * t67 - g(1) * (mrSges(6,1) * t129 - mrSges(6,2) * t130) - g(2) * (-mrSges(6,1) * t127 + mrSges(6,2) * t128) + g(3) * t223 * t170 + (t20 * t90 + t21 * t91) * mrSges(6,3) + t255 + (-Ifges(6,2) * t91 + t41 + t87) * t326 - t352;];
tau = t1;
