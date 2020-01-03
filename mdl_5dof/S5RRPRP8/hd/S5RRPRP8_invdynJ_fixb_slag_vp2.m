% Calculate vector of inverse dynamics joint torques for
% S5RRPRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:45
% DurationCPUTime: 14.18s
% Computational Cost: add. (2590->438), mult. (5656->537), div. (0->0), fcn. (3322->6), ass. (0->198)
t325 = Ifges(5,4) + Ifges(6,4);
t167 = sin(qJ(2));
t234 = qJD(1) * t167;
t146 = pkin(6) * t234;
t298 = qJD(3) + t146;
t166 = sin(qJ(4));
t169 = cos(qJ(4));
t170 = cos(qJ(2));
t233 = qJD(1) * t170;
t78 = -t166 * t234 - t169 * t233;
t333 = t78 / 0.2e1;
t327 = Ifges(5,1) + Ifges(6,1);
t326 = Ifges(4,4) + Ifges(3,5);
t324 = Ifges(5,5) + Ifges(6,5);
t323 = Ifges(5,2) + Ifges(6,2);
t322 = Ifges(4,6) - Ifges(3,6);
t321 = Ifges(5,6) + Ifges(6,6);
t332 = -pkin(7) * t234 + t298;
t79 = -t166 * t233 + t169 * t234;
t282 = t79 / 0.2e1;
t329 = -mrSges(6,1) - mrSges(5,1);
t328 = mrSges(6,2) + mrSges(5,2);
t331 = t325 * t78;
t330 = t325 * t79;
t297 = qJD(2) - qJD(4);
t320 = -t297 * t321 + t323 * t78 + t330;
t319 = -t297 * t324 + t327 * t79 + t331;
t147 = pkin(6) * t233;
t106 = -pkin(7) * t233 + t147;
t172 = -pkin(2) - pkin(3);
t111 = -qJ(3) * t166 + t169 * t172;
t305 = qJD(4) * t111 - t166 * t106 + t169 * t332;
t112 = t169 * qJ(3) + t166 * t172;
t303 = -qJD(4) * t112 - t169 * t106 - t166 * t332;
t144 = Ifges(3,4) * t233;
t252 = Ifges(4,5) * t170;
t191 = t167 * Ifges(4,1) - t252;
t318 = Ifges(3,1) * t234 + qJD(1) * t191 + qJD(2) * t326 + t144;
t220 = mrSges(4,2) * t234;
t317 = mrSges(3,3) * t234 + t220 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t316 = t167 * t322 + t170 * t326;
t192 = t170 * mrSges(4,1) + t167 * mrSges(4,3);
t194 = mrSges(3,1) * t170 - mrSges(3,2) * t167;
t315 = t192 + t194;
t227 = qJD(1) * qJD(2);
t110 = qJDD(1) * t167 + t170 * t227;
t109 = -t170 * qJDD(1) + t167 * t227;
t242 = t167 * t166;
t184 = t170 * t169 + t242;
t239 = t170 * t166;
t185 = -t167 * t169 + t239;
t214 = mrSges(6,1) * t184 - mrSges(6,2) * t185;
t216 = mrSges(5,1) * t184 - mrSges(5,2) * t185;
t314 = -t214 - t216;
t93 = t109 * pkin(6);
t94 = t110 * pkin(6);
t313 = t167 * t94 - t170 * t93;
t57 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t93;
t208 = qJDD(3) + t94;
t61 = -qJDD(2) * pkin(2) + t208;
t312 = t167 * t61 + t170 * t57;
t217 = t172 * qJD(2);
t58 = t217 + t332;
t162 = qJD(2) * qJ(3);
t80 = t106 + t162;
t25 = t166 * t58 + t169 * t80;
t251 = qJ(5) * t78;
t10 = t25 + t251;
t311 = t25 * mrSges(5,3) + t10 * mrSges(6,3);
t228 = m(4) + m(5) + m(6);
t310 = mrSges(2,1) + t315;
t158 = -qJDD(2) + qJDD(4);
t177 = t184 * qJD(4);
t22 = -qJD(1) * t177 + t109 * t166 + t110 * t169;
t40 = -pkin(7) * t110 + qJDD(2) * t172 + t208;
t41 = pkin(7) * t109 + t57;
t4 = -qJD(4) * t25 - t166 * t41 + t169 * t40;
t1 = pkin(4) * t158 - qJ(5) * t22 - qJD(5) * t79 + t4;
t178 = t185 * qJD(4);
t23 = qJD(1) * t178 + t109 * t169 - t110 * t166;
t229 = qJD(4) * t169;
t230 = qJD(4) * t166;
t3 = t166 * t40 + t169 * t41 + t58 * t229 - t230 * t80;
t2 = qJ(5) * t23 + qJD(5) * t78 + t3;
t24 = -t166 * t80 + t169 * t58;
t262 = t24 * mrSges(5,3);
t272 = mrSges(6,3) * t78;
t81 = -qJD(1) * pkin(1) - pkin(2) * t233 - qJ(3) * t234;
t56 = pkin(3) * t233 - t81;
t32 = -pkin(4) * t78 + qJD(5) + t56;
t250 = qJ(5) * t79;
t9 = t24 - t250;
t8 = -pkin(4) * t297 + t9;
t309 = -t32 * (mrSges(6,1) * t79 + mrSges(6,2) * t78) - t56 * (mrSges(5,1) * t79 + mrSges(5,2) * t78) + t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) + t262 * t78 + t272 * t8 + t311 * t79 + t321 * t23 + t324 * t22 + (Ifges(5,3) + Ifges(6,3)) * t158 - (-t323 * t79 + t319 + t331) * t78 / 0.2e1 + (-t327 * t78 + t320 + t330) * t282;
t308 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t15 = -mrSges(6,2) * t158 + mrSges(6,3) * t23;
t16 = -mrSges(5,2) * t158 + mrSges(5,3) * t23;
t307 = t15 + t16;
t306 = -t250 + t305;
t304 = -t251 + t303;
t302 = -t321 * t79 + t324 * t78;
t171 = cos(qJ(1));
t238 = t170 * t171;
t240 = t167 * t171;
t65 = t166 * t238 - t169 * t240;
t66 = t184 * t171;
t301 = -t328 * t66 + t329 * t65;
t168 = sin(qJ(1));
t63 = t185 * t168;
t64 = t184 * t168;
t300 = -t328 * t64 + t329 * t63;
t299 = g(1) * t171 + g(2) * t168;
t296 = qJ(3) * t228;
t221 = t172 * t167;
t247 = t170 * mrSges(4,3);
t140 = pkin(4) * t169 + pkin(3);
t258 = -pkin(2) - t140;
t293 = -m(5) * t221 - m(6) * (pkin(4) * t239 + t167 * t258) - t247 - (-m(4) * pkin(2) - mrSges(4,1)) * t167;
t290 = m(6) * pkin(4);
t279 = pkin(6) - pkin(7);
t278 = pkin(4) * t79;
t274 = -t297 / 0.2e1;
t273 = t167 / 0.2e1;
t267 = pkin(6) * t167;
t266 = pkin(6) * t170;
t155 = t170 * pkin(2);
t255 = Ifges(3,4) * t167;
t254 = Ifges(3,4) * t170;
t253 = Ifges(4,5) * t167;
t121 = t279 * t167;
t122 = t279 * t170;
t48 = t166 * t121 + t169 * t122;
t163 = qJDD(1) * pkin(1);
t152 = t167 * qJ(3);
t151 = t167 * qJD(3);
t231 = qJD(2) * t170;
t237 = qJ(3) * t231 + t151;
t236 = t155 + t152;
t235 = t171 * pkin(1) + t168 * pkin(6);
t232 = qJD(2) * t167;
t219 = mrSges(4,2) * t233;
t218 = t170 * pkin(3) + t236;
t215 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t207 = -pkin(1) - t152;
t206 = -t227 / 0.2e1;
t67 = -qJDD(2) * mrSges(4,1) + t110 * mrSges(4,2);
t47 = t169 * t121 - t122 * t166;
t88 = pkin(1) + t218;
t36 = t109 * pkin(2) - t110 * qJ(3) - qJD(1) * t151 - t163;
t193 = mrSges(3,1) * t167 + mrSges(3,2) * t170;
t190 = t170 * Ifges(3,2) + t255;
t113 = -qJD(2) * pkin(2) + t298;
t117 = t147 + t162;
t186 = t113 * t170 - t117 * t167;
t135 = qJ(3) * t233;
t62 = qJD(1) * t221 + t135;
t183 = pkin(4) * t242 + t140 * t170;
t18 = -pkin(3) * t109 - t36;
t182 = pkin(1) * t193;
t181 = t81 * (t167 * mrSges(4,1) - t247);
t107 = t279 * t232;
t108 = qJD(2) * t122;
t11 = -t169 * t107 + t166 * t108 + t121 * t229 - t122 * t230;
t180 = t167 * (Ifges(3,1) * t170 - t255);
t179 = t170 * (Ifges(4,3) * t167 + t252);
t54 = t167 * t217 + t237;
t176 = t185 * t290;
t12 = -qJD(4) * t48 + t107 * t166 + t169 * t108;
t143 = Ifges(4,5) * t234;
t119 = qJD(2) * mrSges(4,3) + t219;
t118 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t233;
t114 = -pkin(1) - t236;
t104 = -pkin(4) + t111;
t103 = pkin(2) * t234 - t135;
t102 = t192 * qJD(1);
t75 = Ifges(3,6) * qJD(2) + qJD(1) * t190;
t74 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t233 + t143;
t73 = pkin(2) * t232 - t237;
t68 = -mrSges(4,2) * t109 + qJDD(2) * mrSges(4,3);
t53 = -mrSges(5,1) * t297 - mrSges(5,3) * t79;
t52 = -mrSges(6,1) * t297 - mrSges(6,3) * t79;
t51 = mrSges(5,2) * t297 + mrSges(5,3) * t78;
t50 = mrSges(6,2) * t297 + t272;
t46 = pkin(4) * t184 + t88;
t45 = qJD(2) * t184 - t177;
t44 = -qJD(2) * t185 + t178;
t39 = t62 - t278;
t38 = -mrSges(5,1) * t78 + mrSges(5,2) * t79;
t37 = -mrSges(6,1) * t78 + mrSges(6,2) * t79;
t31 = -qJ(5) * t184 + t48;
t30 = qJ(5) * t185 + t47;
t17 = -pkin(4) * t44 + t54;
t14 = mrSges(5,1) * t158 - mrSges(5,3) * t22;
t13 = mrSges(6,1) * t158 - mrSges(6,3) * t22;
t7 = -pkin(4) * t23 + qJDD(5) + t18;
t6 = -qJ(5) * t45 + qJD(5) * t185 + t12;
t5 = qJ(5) * t44 - qJD(5) * t184 + t11;
t19 = [((-t119 - t118) * pkin(6) - t117 * mrSges(4,2) - t75 / 0.2e1 + t74 / 0.2e1) * t232 + (-t329 * t64 - t328 * t63 + (-m(4) * (t207 - t155) - m(6) * (-pkin(1) + t258 * t170 + (-pkin(4) * t166 - qJ(3)) * t167) + m(3) * pkin(1) - m(5) * (t170 * t172 + t207) + t310) * t168 + ((-m(3) - t228) * pkin(6) + t308) * t171) * g(1) + (-mrSges(5,3) * t3 - mrSges(6,3) * t2 - t158 * t321 - t22 * t325 - t23 * t323) * t184 + t170 * (Ifges(3,4) * t110 + Ifges(3,6) * qJDD(2)) / 0.2e1 + m(5) * (t11 * t25 + t12 * t24 + t18 * t88 + t3 * t48 + t4 * t47 + t54 * t56) + m(6) * (t1 * t30 + t10 * t5 + t17 * t32 + t2 * t31 + t46 * t7 + t6 * t8) + (t317 * pkin(6) + t113 * mrSges(4,2) + t318 / 0.2e1) * t231 + t18 * t216 - t73 * t102 + t88 * (-mrSges(5,1) * t23 + mrSges(5,2) * t22) + t5 * t50 + t11 * t51 + t6 * t52 + t12 * t53 + t54 * t38 + t47 * t14 + t48 * t16 + t17 * t37 + t30 * t13 + t31 * t15 + t7 * t214 + t46 * t215 + t316 * qJD(2) ^ 2 / 0.2e1 + (-qJDD(2) * mrSges(3,1) + t67) * t267 + (t167 * (Ifges(4,1) * t170 + t253) + t170 * (-Ifges(3,2) * t167 + t254) + t180) * t227 / 0.2e1 + (t167 * Ifges(3,1) + t191 + t254) * t110 / 0.2e1 + t194 * t163 + (t324 * t274 + t325 * t333 - t8 * mrSges(6,3) + t56 * mrSges(5,2) + t32 * mrSges(6,2) - t262 + t319 / 0.2e1 + t327 * t282) * t45 + (t320 / 0.2e1 + t321 * t274 + t323 * t333 - t56 * mrSges(5,1) - t32 * mrSges(6,1) + t311 + t325 * t282) * t44 + m(4) * (t114 * t36 + t73 * t81 + (qJD(2) * t186 + t312) * pkin(6)) + (t110 * t267 + t313) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t313 * pkin(6)) + (-pkin(1) * t110 - qJDD(2) * t266) * mrSges(3,2) - t182 * t227 + (mrSges(5,3) * t4 + mrSges(6,3) * t1 - t158 * t324 - t22 * t327 - t23 * t325) * t185 - t36 * t192 + (t167 * t326 - t170 * t322) * qJDD(2) / 0.2e1 + ((Ifges(3,1) + Ifges(4,1)) * t110 + t326 * qJDD(2)) * t273 + (-m(5) * pkin(3) * t238 - m(3) * t235 + t329 * t66 + t328 * t65 - t228 * (pkin(2) * t238 + qJ(3) * t240 + t235) + (-m(6) * t183 - t310) * t171 + t308 * t168) * g(2) - t170 * (Ifges(4,5) * t110 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (t253 / 0.2e1 - pkin(1) * mrSges(3,1) + t114 * mrSges(4,1) - t190 / 0.2e1 - mrSges(3,3) * t266 + (Ifges(4,5) - Ifges(3,4)) * t273 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t170) * t109 + qJD(2) * t181 + t179 * t206 + t312 * mrSges(4,2) - t114 * mrSges(4,3) * t110 + t68 * t266 + Ifges(2,3) * qJDD(1); t75 * t234 / 0.2e1 + t118 * t146 + t117 * t220 - (Ifges(4,1) * t233 + t143 + t74) * t234 / 0.2e1 + t322 * t109 + t111 * t14 + t103 * t102 + t104 * t13 - t94 * mrSges(3,1) + t93 * mrSges(3,2) - pkin(2) * t67 + qJ(3) * t68 - t61 * mrSges(4,1) - t62 * t38 + t57 * mrSges(4,3) - t39 * t37 + (-m(6) * (t183 + t236) - m(4) * t236 - m(5) * t218 + t314 - t315) * g(3) + t316 * t206 - t317 * t147 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) - (-Ifges(3,2) * t234 + t144 + t318) * t233 / 0.2e1 + (t300 + (-t170 * t296 + t293) * t168) * g(2) + (t171 * t293 - t238 * t296 + t301) * g(1) + t326 * t110 + (-pkin(2) * t61 + qJ(3) * t57 + qJD(3) * t117 - t103 * t81) * m(4) + (-pkin(6) * t186 * m(4) - t181 + (t179 / 0.2e1 - t180 / 0.2e1 + t182) * qJD(1)) * qJD(1) - t113 * t219 - t309 + t298 * t119 + t299 * t193 + t302 * t274 + t303 * t53 + t304 * t52 + t305 * t51 + (t111 * t4 + t112 * t3 + t24 * t303 + t25 * t305 - t56 * t62) * m(5) + t306 * t50 + (t1 * t104 + t10 * t306 + t112 * t2 + t304 * t8 - t32 * t39) * m(6) + t307 * t112; -qJD(2) * t119 + t228 * t170 * g(3) + (t13 + t14 - t297 * (t50 + t51)) * t169 + (t297 * (t52 + t53) + t307) * t166 + ((-t102 - t37 - t38) * qJD(1) - t299 * t228) * t167 + t67 + (t1 * t169 + t166 * t2 - t234 * t32 - t297 * (t10 * t169 - t166 * t8)) * m(6) + (t166 * t3 + t169 * t4 - t234 * t56 - t297 * (-t166 * t24 + t169 * t25)) * m(5) + (-qJD(2) * t117 + t234 * t81 + t61) * m(4); -t37 * t278 - t9 * t50 - t24 * t51 + t10 * t52 + t25 * t53 + pkin(4) * t13 + t1 * t290 - m(6) * (t32 * t278 + (-t8 + t9) * t10) + t302 * t297 / 0.2e1 + (t184 * t290 - t314) * g(3) + (t168 * t176 - t300) * g(2) + (t171 * t176 - t301) * g(1) + t309; -t78 * t50 + t79 * t52 + (g(1) * t168 - g(2) * t171 - t10 * t78 + t8 * t79 + t7) * m(6) + t215;];
tau = t19;
