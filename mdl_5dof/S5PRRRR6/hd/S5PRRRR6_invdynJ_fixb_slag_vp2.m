% Calculate vector of inverse dynamics joint torques for
% S5PRRRR6
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:34
% EndTime: 2019-12-05 17:09:45
% DurationCPUTime: 4.22s
% Computational Cost: add. (3600->365), mult. (5936->507), div. (0->0), fcn. (3959->14), ass. (0->200)
t184 = cos(qJ(4));
t264 = t184 * mrSges(5,1);
t180 = sin(qJ(4));
t281 = mrSges(5,2) * t180;
t313 = t264 - t281;
t321 = -mrSges(4,1) - t313;
t320 = -mrSges(5,3) - mrSges(6,3);
t176 = qJ(2) + qJ(3);
t164 = sin(t176);
t175 = qJ(4) + qJ(5);
t163 = sin(t175);
t280 = mrSges(6,2) * t163;
t319 = (-t280 - t281) * t164;
t318 = -m(5) * pkin(7) + t320;
t170 = qJDD(2) + qJDD(3);
t172 = qJD(2) + qJD(3);
t238 = qJD(4) * t180;
t107 = t170 * t184 - t172 * t238;
t316 = t107 / 0.2e1;
t186 = cos(qJ(2));
t182 = sin(qJ(2));
t236 = qJD(1) * qJD(2);
t222 = t182 * t236;
t133 = t186 * qJDD(1) - t222;
t243 = qJD(1) * t182;
t315 = qJDD(2) * pkin(2) - qJD(3) * t243 + t133;
t165 = cos(t175);
t273 = t165 * mrSges(6,1);
t314 = (t264 + t273) * t164;
t221 = t186 * t236;
t134 = qJDD(1) * t182 + t221;
t242 = qJD(1) * t186;
t148 = qJD(2) * pkin(2) + t242;
t181 = sin(qJ(3));
t185 = cos(qJ(3));
t240 = qJD(3) * t185;
t50 = t185 * t134 + t148 * t240 + t315 * t181;
t43 = pkin(7) * t170 + t50;
t100 = t148 * t181 + t185 * t243;
t89 = pkin(7) * t172 + t100;
t19 = t184 * t43 - t238 * t89;
t237 = qJD(4) * t184;
t20 = -t180 * t43 - t237 * t89;
t205 = -t180 * t20 + t184 * t19;
t187 = -pkin(8) - pkin(7);
t227 = qJD(4) * t187;
t126 = t180 * t227;
t127 = t184 * t227;
t179 = sin(qJ(5));
t183 = cos(qJ(5));
t203 = t179 * t180 - t183 * t184;
t142 = t187 * t180;
t167 = t184 * pkin(8);
t143 = pkin(7) * t184 + t167;
t84 = t142 * t183 - t143 * t179;
t150 = t181 * t243;
t99 = t148 * t185 - t150;
t312 = qJD(5) * t84 + t126 * t183 + t127 * t179 + t203 * t99;
t124 = t179 * t184 + t180 * t183;
t85 = t142 * t179 + t143 * t183;
t311 = -qJD(5) * t85 + t124 * t99 - t126 * t179 + t127 * t183;
t111 = t185 * t242 - t150;
t291 = pkin(2) * t181;
t157 = pkin(7) + t291;
t284 = -pkin(8) - t157;
t115 = t284 * t180;
t259 = t157 * t184;
t116 = t167 + t259;
t71 = t115 * t179 + t116 * t183;
t231 = pkin(2) * t240;
t213 = t184 * t231;
t218 = qJD(4) * t284;
t86 = t180 * t218 + t213;
t214 = t180 * t231;
t87 = t184 * t218 - t214;
t310 = -qJD(5) * t71 + t124 * t111 - t179 * t86 + t183 * t87;
t70 = t115 * t183 - t116 * t179;
t309 = qJD(5) * t70 + t203 * t111 + t179 * t87 + t183 * t86;
t294 = m(6) * pkin(4);
t308 = -t294 - mrSges(5,1);
t125 = t181 * t186 + t182 * t185;
t65 = t203 * t125;
t158 = pkin(4) * t184 + pkin(3);
t166 = cos(t176);
t204 = -t158 * t164 - t166 * t187;
t288 = pkin(3) * t164;
t290 = pkin(2) * t182;
t307 = -m(6) * (t204 - t290) - m(5) * (-t288 - t290) + t314;
t256 = t172 * t180;
t131 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t256;
t255 = t172 * t184;
t132 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t255;
t306 = t131 * t184 + t132 * t180;
t305 = t321 * t172;
t177 = sin(pkin(9));
t178 = cos(pkin(9));
t304 = g(1) * t178 + g(2) * t177;
t171 = qJD(4) + qJD(5);
t103 = t203 * t172;
t104 = t124 * t172;
t62 = mrSges(6,1) * t103 + mrSges(6,2) * t104;
t303 = t62 + t305;
t302 = (-t273 + t280 + t321) * t166 + (mrSges(4,2) + t320) * t164;
t269 = t172 * mrSges(4,2);
t301 = -t180 * t131 + t184 * t132 - t269;
t258 = t166 * t177;
t300 = t319 * t177 + t318 * t258;
t257 = t166 * t178;
t299 = t319 * t178 + t318 * t257;
t298 = m(5) * t288 - m(6) * t204 + t314;
t108 = t170 * t180 + t172 * t237;
t266 = t180 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t108);
t297 = m(5) * t205 - t131 * t237 - t132 * t238 - t266;
t217 = (t180 ^ 2 + t184 ^ 2) * t89;
t296 = -m(5) * t217 - t301;
t73 = -t158 * t172 - t99;
t88 = -pkin(3) * t172 - t99;
t295 = -m(5) * t88 - m(6) * t73 - t303;
t292 = t104 / 0.2e1;
t289 = pkin(2) * t185;
t285 = g(3) * t164;
t168 = t186 * pkin(2);
t283 = (-t163 * t258 - t165 * t178) * mrSges(6,1) + (t163 * t178 - t165 * t258) * mrSges(6,2);
t282 = (-t163 * t257 + t165 * t177) * mrSges(6,1) + (-t163 * t177 - t165 * t257) * mrSges(6,2);
t279 = mrSges(6,3) * t103;
t278 = Ifges(5,4) * t180;
t277 = Ifges(5,4) * t184;
t276 = Ifges(6,4) * t104;
t274 = t104 * mrSges(6,3);
t272 = t170 * mrSges(4,1);
t271 = t170 * mrSges(4,2);
t223 = pkin(8) * t172 + t89;
t68 = t223 * t184;
t268 = t179 * t68;
t265 = t183 * t68;
t94 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t107;
t262 = t184 * t94;
t254 = t177 * t180;
t253 = t177 * t184;
t252 = t178 * t180;
t251 = t178 * t184;
t244 = t166 * pkin(3) + t164 * pkin(7);
t241 = qJD(3) * t181;
t239 = qJD(4) * t172;
t235 = pkin(4) * t256;
t232 = pkin(2) * t241;
t230 = pkin(4) * t238;
t215 = t166 * t158 - t164 * t187;
t67 = t223 * t180;
t211 = mrSges(4,1) * t164 + mrSges(4,2) * t166;
t209 = mrSges(5,1) * t180 + mrSges(5,2) * t184;
t208 = -mrSges(6,1) * t163 - mrSges(6,2) * t165;
t207 = Ifges(5,2) * t184 + t278;
t206 = Ifges(5,5) * t184 - Ifges(5,6) * t180;
t63 = qJD(4) * pkin(4) - t67;
t25 = t183 * t63 - t268;
t26 = t179 * t63 + t265;
t123 = t181 * t182 - t185 * t186;
t51 = -t181 * t134 - t148 * t241 + t315 * t185;
t198 = t88 * t209;
t197 = t180 * (Ifges(5,1) * t184 - t278);
t196 = t124 * qJD(5);
t195 = t203 * qJD(5);
t44 = -pkin(3) * t170 - t51;
t169 = qJDD(4) + qJDD(5);
t10 = qJDD(4) * pkin(4) - pkin(8) * t108 + t20;
t12 = pkin(8) * t107 + t19;
t3 = qJD(5) * t25 + t10 * t179 + t12 * t183;
t4 = -qJD(5) * t26 + t10 * t183 - t12 * t179;
t45 = t107 * t179 + t108 * t183 - t172 * t195;
t46 = t107 * t183 - t108 * t179 - t172 * t196;
t55 = -Ifges(6,2) * t103 + Ifges(6,6) * t171 + t276;
t96 = Ifges(6,4) * t103;
t56 = Ifges(6,1) * t104 + Ifges(6,5) * t171 - t96;
t190 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t25 * t279 + t55 * t292 - t73 * (mrSges(6,1) * t104 - mrSges(6,2) * t103) + Ifges(6,3) * t169 - t104 * (-Ifges(6,1) * t103 - t276) / 0.2e1 + Ifges(6,6) * t46 + Ifges(6,5) * t45 - t171 * (-Ifges(6,5) * t103 - Ifges(6,6) * t104) / 0.2e1 + (-Ifges(6,2) * t104 + t56 - t96) * t103 / 0.2e1;
t75 = -qJD(4) * t124 - t196;
t101 = Ifges(5,6) * qJD(4) + t172 * t207;
t149 = Ifges(5,4) * t255;
t102 = Ifges(5,1) * t256 + Ifges(5,5) * qJD(4) + t149;
t31 = -pkin(4) * t107 + t44;
t74 = -qJD(4) * t203 - t195;
t189 = t108 * t277 / 0.2e1 - t101 * t238 / 0.2e1 + t197 * t239 / 0.2e1 + t207 * t316 - t44 * t313 + (Ifges(6,1) * t74 + Ifges(6,4) * t75) * t292 + t184 * (Ifges(5,4) * t108 + Ifges(5,2) * t107) / 0.2e1 + Ifges(4,3) * t170 + t171 * (Ifges(6,5) * t74 + Ifges(6,6) * t75) / 0.2e1 - t103 * (Ifges(6,4) * t74 + Ifges(6,2) * t75) / 0.2e1 + t74 * t56 / 0.2e1 + t75 * t55 / 0.2e1 + t73 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) - t50 * mrSges(4,2) + t51 * mrSges(4,1) + (Ifges(5,1) * t108 + Ifges(5,4) * t316) * t180 + (t102 + t172 * (-Ifges(5,2) * t180 + t277)) * t237 / 0.2e1 + (t198 + t206 * qJD(4) / 0.2e1) * qJD(4) + (-t25 * t74 + t26 * t75) * mrSges(6,3) + t205 * mrSges(5,3) + qJDD(4) * (Ifges(5,5) * t180 + Ifges(5,6) * t184) + (t31 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t45 + Ifges(6,4) * t46 + Ifges(6,5) * t169) * t124 + (t31 * mrSges(6,1) - t3 * mrSges(6,3) - Ifges(6,4) * t45 - Ifges(6,2) * t46 - Ifges(6,6) * t169) * t203;
t188 = qJD(2) ^ 2;
t159 = -pkin(3) - t289;
t137 = -t158 - t289;
t128 = t230 + t232;
t81 = mrSges(6,1) * t171 - t274;
t80 = -mrSges(6,2) * t171 - t279;
t77 = t172 * t125;
t76 = t172 * t123;
t69 = -mrSges(5,1) * t107 + mrSges(5,2) * t108;
t64 = t124 * t125;
t37 = -mrSges(6,2) * t169 + mrSges(6,3) * t46;
t36 = mrSges(6,1) * t169 - mrSges(6,3) * t45;
t30 = -t183 * t67 - t268;
t29 = t179 * t67 - t265;
t9 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t6 = t124 * t76 + t171 * t65;
t5 = t125 * t75 + t203 * t76;
t1 = [m(2) * qJDD(1) - t64 * t36 - t65 * t37 + t5 * t80 + t6 * t81 + (-qJDD(2) * t182 - t186 * t188) * mrSges(3,2) + (qJDD(2) * t186 - t182 * t188) * mrSges(3,1) + t303 * t77 - t301 * t76 + (t69 + t9 - t272) * t123 + (-t306 * qJD(4) + t262 - t266 - t271) * t125 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + m(5) * (t123 * t44 + t125 * t205 - t217 * t76 + t77 * t88) + m(3) * (t133 * t186 + t134 * t182) + m(6) * (t123 * t31 + t25 * t6 + t26 * t5 - t3 * t65 - t4 * t64 + t73 * t77) + m(4) * (-t100 * t76 - t123 * t51 + t125 * t50 - t77 * t99); (m(4) * t99 + t295) * t125 * qJD(1) + m(5) * (t159 * t44 + (t181 * t88 + t185 * t217) * qJD(3) * pkin(2)) - t271 * t291 - t231 * t269 + (t221 - t134) * mrSges(3,2) + (t222 + t133) * mrSges(3,1) + t309 * t80 + t310 * t81 + (t128 * t73 + t137 * t31 + t310 * t25 + t309 * t26 + t3 * t71 + t4 * t70) * m(6) + (t307 * t178 + t299) * g(1) + (t307 * t177 + t300) * g(2) + t305 * t232 + (-m(4) * t168 - m(5) * (t168 + t244) - mrSges(3,1) * t186 + mrSges(3,2) * t182 - m(6) * (t168 + t215) + t302) * g(3) + (-m(4) * t100 + t296) * t111 + t297 * t157 + t189 + (m(4) * t290 + mrSges(3,1) * t182 + mrSges(3,2) * t186 + t211) * t304 + t272 * t289 + t94 * t259 + t159 * t69 + t128 * t62 + t137 * t9 + t71 * t37 + t70 * t36 + m(4) * (t181 * t50 + t185 * t51 + (t100 * t185 - t181 * t99) * qJD(3)) * pkin(2) + t132 * t213 + Ifges(3,3) * qJDD(2) - t131 * t214; -t158 * t9 + t62 * t230 + t84 * t36 + t85 * t37 + t189 + t311 * t81 + t312 * t80 + t304 * t211 + (t298 * t177 + t300) * g(2) + (t298 * t178 + t299) * g(1) + (-m(5) * t44 - t69) * pkin(3) + (-t158 * t31 + t73 * t230 + t311 * t25 + t312 * t26 + t3 * t85 + t4 * t84) * m(6) + t296 * t99 + t295 * t100 + (-m(5) * t244 - m(6) * t215 + t302) * g(3) + (t262 + t297) * pkin(7); t101 * t256 / 0.2e1 - t206 * t239 / 0.2e1 + t190 + (t179 * t3 + t183 * t4 + (-t179 * t25 + t183 * t26) * qJD(5)) * t294 + t26 * t274 + Ifges(5,5) * t108 + Ifges(5,6) * t107 - t30 * t80 - t29 * t81 - t19 * mrSges(5,2) + t20 * mrSges(5,1) - t62 * t235 - m(6) * (t235 * t73 + t25 * t29 + t26 * t30) + Ifges(5,3) * qJDD(4) + t306 * t89 + (t180 * t294 - t208 + t209) * t285 - (-Ifges(5,2) * t256 + t102 + t149) * t255 / 0.2e1 + (-t198 - t197 * t172 / 0.2e1) * t172 + (-t283 - (-t166 * t253 + t252) * mrSges(5,2) + t308 * (-t166 * t254 - t251)) * g(2) + (-t282 - (-t166 * t251 - t254) * mrSges(5,2) + t308 * (-t166 * t252 + t253)) * g(1) + ((-t179 * t81 + t183 * t80) * qJD(5) + t179 * t37 + t183 * t36) * pkin(4); t190 + (t81 + t274) * t26 - t208 * t285 - t25 * t80 - g(2) * t283 - g(1) * t282;];
tau = t1;
