% Calculate vector of inverse dynamics joint torques for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:49
% EndTime: 2019-12-05 18:31:55
% DurationCPUTime: 2.63s
% Computational Cost: add. (3753->352), mult. (5994->474), div. (0->0), fcn. (3531->16), ass. (0->180)
t178 = sin(qJ(4));
t182 = cos(qJ(4));
t136 = -t182 * mrSges(5,1) + mrSges(5,2) * t178;
t232 = Ifges(5,4) * t178;
t257 = t232 / 0.2e1;
t231 = Ifges(5,4) * t182;
t256 = t182 * Ifges(5,2);
t173 = qJ(4) + qJ(5);
t166 = cos(t173);
t148 = t166 * mrSges(6,1);
t242 = pkin(4) * t182;
t154 = pkin(3) + t242;
t255 = m(5) * pkin(3) + m(6) * t154 + mrSges(4,1) - t136 + t148;
t254 = -m(5) * pkin(7) - mrSges(6,3) - mrSges(5,3) + mrSges(4,2) + m(6) * (-pkin(8) - pkin(7));
t175 = sin(pkin(9));
t179 = sin(qJ(2));
t230 = qJD(1) * pkin(1);
t209 = t179 * t230;
t139 = t175 * t209;
t176 = cos(pkin(9));
t183 = cos(qJ(2));
t208 = t183 * t230;
t105 = t176 * t208 - t139;
t244 = pkin(2) * t175;
t151 = pkin(7) + t244;
t236 = -pkin(8) - t151;
t202 = qJD(4) * t236;
t109 = t178 * t202;
t110 = t182 * t202;
t177 = sin(qJ(5));
t181 = cos(qJ(5));
t125 = t177 * t182 + t178 * t181;
t119 = t236 * t178;
t168 = t182 * pkin(8);
t120 = t151 * t182 + t168;
t72 = t119 * t177 + t120 * t181;
t251 = -qJD(5) * t72 + t125 * t105 - t109 * t177 + t110 * t181;
t124 = -t177 * t178 + t181 * t182;
t71 = t119 * t181 - t120 * t177;
t250 = qJD(5) * t71 - t124 * t105 + t109 * t181 + t110 * t177;
t172 = qJD(1) + qJD(2);
t99 = t125 * t172;
t248 = t99 / 0.2e1;
t98 = t124 * t172;
t247 = mrSges(6,3) * t98;
t246 = Ifges(6,4) * t99;
t245 = pkin(1) * t183;
t243 = pkin(2) * t176;
t174 = qJ(1) + qJ(2);
t162 = pkin(9) + t174;
t150 = cos(t162);
t241 = g(3) * t150;
t163 = t182 * qJD(3);
t132 = pkin(2) * t172 + t208;
t87 = t175 * t132 + t176 * t209;
t77 = pkin(7) * t172 + t87;
t67 = -t178 * t77 + t163;
t240 = t67 * mrSges(5,3);
t213 = qJD(3) * t178;
t68 = t182 * t77 + t213;
t239 = t68 * mrSges(5,3);
t238 = t99 * mrSges(6,3);
t155 = pkin(2) + t245;
t218 = t176 * t179;
t108 = pkin(1) * t218 + t175 * t155;
t102 = pkin(7) + t108;
t237 = -pkin(8) - t102;
t164 = sin(t173);
t234 = mrSges(6,2) * t164;
t233 = mrSges(6,2) * t166;
t205 = pkin(8) * t172 + t77;
t58 = t182 * t205 + t213;
t229 = t177 * t58;
t228 = t181 * t58;
t212 = qJD(4) * t178;
t170 = qJDD(1) + qJDD(2);
t121 = -qJD(2) * t209 + qJDD(1) * t245;
t100 = pkin(2) * t170 + t121;
t214 = qJD(2) * t183;
t122 = (qJD(1) * t214 + qJDD(1) * t179) * pkin(1);
t66 = t175 * t100 + t176 * t122;
t57 = pkin(7) * t170 + t66;
t20 = qJD(4) * t163 + t178 * qJDD(3) + t182 * t57 - t212 * t77;
t226 = t182 * t20;
t225 = Ifges(5,5) * qJD(4);
t224 = Ifges(5,6) * qJD(4);
t223 = qJD(4) * t68;
t149 = sin(t162);
t222 = t149 * t164;
t221 = t172 * t178;
t220 = t172 * t182;
t219 = t175 * t179;
t106 = (t176 * t183 - t219) * qJD(2) * pkin(1);
t217 = t178 * t106;
t216 = t182 * t106;
t215 = mrSges(6,1) * t222 + t149 * t233;
t211 = qJD(4) * t182;
t210 = m(4) + m(5) + m(6);
t207 = pkin(4) * t212;
t206 = t172 * t212;
t203 = qJD(4) * t237;
t201 = t148 - t234;
t65 = t100 * t176 - t175 * t122;
t86 = t132 * t176 - t139;
t107 = -pkin(1) * t219 + t155 * t176;
t101 = -pkin(3) - t107;
t200 = pkin(2) * t210 + mrSges(3,1);
t199 = mrSges(5,1) * t178 + mrSges(5,2) * t182;
t198 = -mrSges(6,1) * t164 - t233;
t197 = t232 + t256;
t56 = -t178 * t205 + t163;
t50 = qJD(4) * pkin(4) + t56;
t15 = t181 * t50 - t229;
t16 = t177 * t50 + t228;
t82 = t237 * t178;
t83 = t102 * t182 + t168;
t43 = -t177 * t83 + t181 * t82;
t44 = t177 * t82 + t181 * t83;
t196 = -t178 * t67 + t182 * t68;
t195 = mrSges(2,1) + (m(3) + t210) * pkin(1);
t55 = -pkin(3) * t170 - t65;
t194 = pkin(1) * (t175 * t183 + t218);
t193 = t125 * qJD(5);
t192 = t124 * qJD(5);
t104 = qJD(2) * t194;
t21 = t182 * qJDD(3) - t178 * t57 - t223;
t169 = qJDD(4) + qJDD(5);
t171 = qJD(4) + qJD(5);
t115 = t170 * t178 + t172 * t211;
t12 = qJDD(4) * pkin(4) - pkin(8) * t115 + t21;
t114 = t170 * t182 - t206;
t13 = pkin(8) * t114 + t20;
t3 = qJD(5) * t15 + t12 * t177 + t13 * t181;
t4 = -qJD(5) * t16 + t12 * t181 - t13 * t177;
t41 = t114 * t177 + t115 * t181 + t172 * t192;
t42 = t114 * t181 - t115 * t177 - t172 * t193;
t51 = Ifges(6,2) * t98 + Ifges(6,6) * t171 + t246;
t93 = Ifges(6,4) * t98;
t52 = Ifges(6,1) * t99 + Ifges(6,5) * t171 + t93;
t69 = -t154 * t172 - t86;
t190 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t15 * t247 + t51 * t248 - t69 * (mrSges(6,1) * t99 + mrSges(6,2) * t98) + Ifges(6,3) * t169 - t99 * (Ifges(6,1) * t98 - t246) / 0.2e1 + Ifges(6,6) * t42 + Ifges(6,5) * t41 - t171 * (Ifges(6,5) * t98 - Ifges(6,6) * t99) / 0.2e1 - (-Ifges(6,2) * t99 + t52 + t93) * t98 / 0.2e1;
t189 = m(5) * (-t178 * t21 + t226 + (-t178 * t68 - t182 * t67) * qJD(4));
t165 = sin(t174);
t167 = cos(t174);
t188 = -t165 * mrSges(3,2) + t200 * t167 + (-t234 + t255) * t150 - t254 * t149;
t187 = mrSges(3,2) * t167 - mrSges(6,2) * t222 + t255 * t149 + t254 * t150 + t200 * t165;
t38 = -pkin(4) * t114 + t55;
t74 = qJD(4) * t124 + t192;
t75 = -qJD(4) * t125 - t193;
t76 = -pkin(3) * t172 - t86;
t96 = t172 * t197 + t224;
t140 = Ifges(5,4) * t220;
t97 = Ifges(5,1) * t221 + t140 + t225;
t186 = t171 * (Ifges(6,5) * t74 + Ifges(6,6) * t75) / 0.2e1 + t55 * t136 + t121 * mrSges(3,1) - t122 * mrSges(3,2) + t98 * (Ifges(6,4) * t74 + Ifges(6,2) * t75) / 0.2e1 + t69 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + t74 * t52 / 0.2e1 + t75 * t51 / 0.2e1 + t65 * mrSges(4,1) + qJD(4) ^ 2 * (Ifges(5,5) * t182 - Ifges(5,6) * t178) / 0.2e1 + (Ifges(6,1) * t74 + Ifges(6,4) * t75) * t248 + mrSges(5,3) * t226 + (Ifges(5,1) * t182 - t232) * t206 / 0.2e1 + t76 * t199 * qJD(4) - t96 * t212 / 0.2e1 + (t256 / 0.2e1 + t257 + t197 / 0.2e1) * t114 + t115 * (t178 * Ifges(5,1) + t231) + (t97 + t172 * (-Ifges(5,2) * t178 + t231)) * t211 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t170 + (-t15 * t74 + t16 * t75) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t178 + Ifges(5,6) * t182) + (t38 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t41 + Ifges(6,4) * t42 + Ifges(6,5) * t169) * t125 + (-t38 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t41 + Ifges(6,2) * t42 + Ifges(6,6) * t169) * t124;
t184 = cos(qJ(1));
t180 = sin(qJ(1));
t152 = -pkin(3) - t243;
t135 = -t154 - t243;
t131 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t220;
t130 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t221;
t113 = t136 * t172;
t103 = qJD(1) * t194;
t92 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t115;
t91 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t114;
t90 = t101 - t242;
t88 = t104 + t207;
t81 = mrSges(6,1) * t171 - t238;
t80 = -mrSges(6,2) * t171 + t247;
t70 = -mrSges(5,1) * t114 + mrSges(5,2) * t115;
t62 = -mrSges(6,1) * t98 + mrSges(6,2) * t99;
t47 = t182 * t203 - t217;
t46 = t178 * t203 + t216;
t35 = -mrSges(6,2) * t169 + mrSges(6,3) * t42;
t34 = mrSges(6,1) * t169 - mrSges(6,3) * t41;
t19 = t181 * t56 - t229;
t18 = -t177 * t56 - t228;
t10 = -mrSges(6,1) * t42 + mrSges(6,2) * t41;
t6 = -qJD(5) * t44 - t177 * t46 + t181 * t47;
t5 = qJD(5) * t43 + t177 * t47 + t181 * t46;
t1 = [(mrSges(2,2) * t184 + t180 * t195 + t187) * g(3) + (-t106 * t130 + (-qJD(4) * t131 - t92) * t102 + (-t21 - t223) * mrSges(5,3)) * t178 + (t102 * t91 + t106 * t131 + (-t102 * t130 - t240) * qJD(4)) * t182 + t186 + t104 * t113 + t101 * t70 + t88 * t62 + t90 * t10 + t5 * t80 + t6 * t81 + t43 * t34 + t44 * t35 + m(5) * (t101 * t55 + t104 * t76 + t216 * t68 - t217 * t67) + t102 * t189 + (-mrSges(2,2) * t180 + t184 * t195 + t188) * g(2) + m(4) * (-t104 * t86 + t106 * t87 + t107 * t65 + t108 * t66) + m(6) * (t15 * t6 + t16 * t5 + t3 * t44 + t38 * t90 + t4 * t43 + t69 * t88) + Ifges(2,3) * qJDD(1) + (-t104 * t172 + t107 * t170) * mrSges(4,1) + (-t106 * t172 - t108 * t170 - t66) * mrSges(4,2) + ((-t170 * t179 - t172 * t214) * mrSges(3,2) + (-qJD(2) * t172 * t179 + t170 * t183) * mrSges(3,1) + m(3) * (t121 * t183 + t122 * t179)) * pkin(1); t251 * t81 + t250 * t80 + t152 * t70 + t186 + t135 * t10 + t71 * t34 + t72 * t35 + (-t21 * mrSges(5,3) + t105 * t130 - t151 * t92 + (pkin(4) * t62 - t131 * t151 - t239) * qJD(4)) * t178 + t170 * mrSges(4,1) * t243 + t151 * t189 + (-t170 * t244 - t66) * mrSges(4,2) + (-t113 - t62) * t103 + (-t105 * t131 + t151 * t91 + (-t130 * t151 - t240) * qJD(4)) * t182 + (t103 * mrSges(4,1) + t105 * mrSges(4,2) + (mrSges(3,1) * t179 + mrSges(3,2) * t183) * t230) * t172 + t187 * g(3) + t188 * g(2) + (t135 * t38 + t3 * t72 + t4 * t71 + (-t103 + t207) * t69 + t250 * t16 + t251 * t15) * m(6) + (-t103 * t76 - t105 * t196 + t152 * t55) * m(5) + (t103 * t86 - t105 * t87 + (t175 * t66 + t176 * t65) * pkin(2)) * m(4); m(4) * qJDD(3) + t124 * t34 + t125 * t35 + t178 * t91 + t182 * t92 + t74 * t80 + t75 * t81 + (-t178 * t130 + t182 * t131) * qJD(4) + m(5) * (qJD(4) * t196 + t178 * t20 + t182 * t21) + m(6) * (t124 * t4 + t125 * t3 + t15 * t75 + t16 * t74) - t210 * g(1); t190 + t68 * t130 - t67 * t131 + Ifges(5,6) * t114 + Ifges(5,5) * t115 - t19 * t80 - t18 * t81 - t20 * mrSges(5,2) + t21 * mrSges(5,1) + t16 * t238 + ((-t76 * mrSges(5,2) - t225 / 0.2e1 - t140 / 0.2e1 - t97 / 0.2e1 + t240) * t182 + (-t76 * mrSges(5,1) + t224 / 0.2e1 + t96 / 0.2e1 + t239 + (t257 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t182) * t172 + (-m(6) * t69 - t62) * pkin(4)) * t178) * t172 + (t177 * t35 + t181 * t34 + (-g(1) * t182 + t177 * t3 + t181 * t4 + (-g(2) * t149 + t241) * t178) * m(6) + (-t177 * t81 + t181 * t80 + (-t15 * t177 + t16 * t181) * m(6)) * qJD(5)) * pkin(4) - m(6) * (t15 * t18 + t16 * t19) + (-t198 + t199) * t241 + Ifges(5,3) * qJDD(4) + (t136 - t201) * g(1) + (-t149 * t199 - t215) * g(2); (t81 + t238) * t16 - t198 * t241 - g(1) * t201 + t190 - g(2) * t215 - t15 * t80;];
tau = t1;
