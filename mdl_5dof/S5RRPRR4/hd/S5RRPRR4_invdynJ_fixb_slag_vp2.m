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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:01:49
% EndTime: 2020-01-03 12:01:54
% DurationCPUTime: 2.51s
% Computational Cost: add. (3753->379), mult. (5994->500), div. (0->0), fcn. (3531->16), ass. (0->191)
t190 = sin(qJ(4));
t248 = Ifges(5,4) * t190;
t269 = t248 / 0.2e1;
t194 = cos(qJ(4));
t247 = Ifges(5,4) * t194;
t268 = t194 * Ifges(5,2);
t249 = mrSges(5,1) * t194;
t267 = -mrSges(4,1) - t249;
t266 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t187 = sin(pkin(9));
t191 = sin(qJ(2));
t246 = qJD(1) * pkin(1);
t223 = t191 * t246;
t141 = t187 * t223;
t188 = cos(pkin(9));
t195 = cos(qJ(2));
t222 = t195 * t246;
t105 = t188 * t222 - t141;
t258 = pkin(2) * t187;
t159 = pkin(7) + t258;
t250 = -pkin(8) - t159;
t214 = qJD(4) * t250;
t109 = t190 * t214;
t110 = t194 * t214;
t189 = sin(qJ(5));
t193 = cos(qJ(5));
t125 = t189 * t194 + t190 * t193;
t119 = t250 * t190;
t179 = t194 * pkin(8);
t120 = t159 * t194 + t179;
t72 = t119 * t189 + t120 * t193;
t265 = -qJD(5) * t72 + t125 * t105 - t109 * t189 + t110 * t193;
t124 = -t189 * t190 + t193 * t194;
t71 = t119 * t193 - t120 * t189;
t264 = qJD(5) * t71 - t124 * t105 + t109 * t193 + t110 * t189;
t184 = qJD(1) + qJD(2);
t99 = t125 * t184;
t262 = t99 / 0.2e1;
t98 = t124 * t184;
t261 = mrSges(6,3) * t98;
t260 = Ifges(6,4) * t99;
t259 = pkin(1) * t195;
t186 = qJ(1) + qJ(2);
t175 = sin(t186);
t162 = pkin(2) * t175;
t177 = cos(t186);
t163 = pkin(2) * t177;
t257 = pkin(2) * t188;
t256 = pkin(4) * t194;
t172 = pkin(9) + t186;
t157 = sin(t172);
t255 = g(2) * t157;
t173 = t194 * qJD(3);
t134 = pkin(2) * t184 + t222;
t87 = t187 * t134 + t188 * t223;
t77 = pkin(7) * t184 + t87;
t67 = -t190 * t77 + t173;
t254 = t67 * mrSges(5,3);
t226 = qJD(3) * t190;
t68 = t194 * t77 + t226;
t253 = t68 * mrSges(5,3);
t252 = t99 * mrSges(6,3);
t165 = pkin(2) + t259;
t232 = t188 * t191;
t108 = pkin(1) * t232 + t187 * t165;
t102 = pkin(7) + t108;
t251 = -pkin(8) - t102;
t185 = qJ(4) + qJ(5);
t174 = sin(t185);
t245 = t174 * mrSges(6,2);
t176 = cos(t185);
t154 = t176 * mrSges(6,1);
t217 = pkin(8) * t184 + t77;
t58 = t194 * t217 + t226;
t244 = t189 * t58;
t243 = t190 * mrSges(5,2);
t182 = qJDD(1) + qJDD(2);
t121 = -qJD(2) * t223 + qJDD(1) * t259;
t100 = pkin(2) * t182 + t121;
t227 = qJD(2) * t195;
t122 = (qJD(1) * t227 + qJDD(1) * t191) * pkin(1);
t66 = t187 * t100 + t188 * t122;
t57 = pkin(7) * t182 + t66;
t21 = -qJD(4) * t68 + t194 * qJDD(3) - t190 * t57;
t242 = t190 * t21;
t241 = t193 * t58;
t225 = qJD(4) * t190;
t20 = qJD(4) * t173 + t190 * qJDD(3) + t194 * t57 - t225 * t77;
t240 = t194 * t20;
t239 = Ifges(5,5) * qJD(4);
t238 = Ifges(5,6) * qJD(4);
t158 = cos(t172);
t237 = t158 * t174;
t236 = t158 * t176;
t235 = t184 * t190;
t234 = t184 * t194;
t233 = t187 * t191;
t106 = (t188 * t195 - t233) * qJD(2) * pkin(1);
t231 = t190 * t106;
t230 = t194 * t106;
t229 = mrSges(6,1) * t237 + mrSges(6,2) * t236;
t192 = sin(qJ(1));
t178 = t192 * pkin(1);
t228 = t162 + t178;
t224 = qJD(4) * t194;
t221 = pkin(4) * t225;
t164 = pkin(3) + t256;
t197 = -pkin(8) - pkin(7);
t220 = t157 * t164 + t158 * t197 + t162;
t219 = t158 * pkin(3) + t157 * pkin(7) + t163;
t218 = t184 * t225;
t215 = qJD(4) * t251;
t213 = t154 - t245;
t65 = t100 * t188 - t187 * t122;
t86 = t134 * t188 - t141;
t107 = -pkin(1) * t233 + t165 * t188;
t101 = -pkin(3) - t107;
t212 = -t157 * t197 + t158 * t164 + t163;
t138 = t243 - t249;
t211 = mrSges(5,1) * t190 + mrSges(5,2) * t194;
t210 = -mrSges(6,1) * t174 - mrSges(6,2) * t176;
t209 = t243 + t245;
t208 = t248 + t268;
t56 = -t190 * t217 + t173;
t50 = qJD(4) * pkin(4) + t56;
t15 = t193 * t50 - t244;
t16 = t189 * t50 + t241;
t82 = t251 * t190;
t83 = t102 * t194 + t179;
t43 = -t189 * t83 + t193 * t82;
t44 = t189 * t82 + t193 * t83;
t207 = -t190 * t67 + t194 * t68;
t132 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t235;
t133 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t234;
t206 = -t190 * t132 + t194 * t133;
t55 = -pkin(3) * t182 - t65;
t205 = pkin(1) * (t187 * t195 + t232);
t204 = t125 * qJD(5);
t203 = t124 * qJD(5);
t104 = qJD(2) * t205;
t202 = -t242 + (-t190 * t68 - t194 * t67) * qJD(4);
t201 = -t177 * mrSges(3,1) - mrSges(6,1) * t236 + mrSges(3,2) * t175 + t266 * t157 + t267 * t158;
t200 = -t175 * mrSges(3,1) - t177 * mrSges(3,2) + (m(5) * pkin(7) - t266) * t158 + (-t154 + t267) * t157;
t181 = qJDD(4) + qJDD(5);
t183 = qJD(4) + qJD(5);
t115 = t182 * t190 + t184 * t224;
t12 = qJDD(4) * pkin(4) - pkin(8) * t115 + t21;
t114 = t182 * t194 - t218;
t13 = pkin(8) * t114 + t20;
t3 = qJD(5) * t15 + t12 * t189 + t13 * t193;
t4 = -qJD(5) * t16 + t12 * t193 - t13 * t189;
t41 = t114 * t189 + t115 * t193 + t184 * t203;
t42 = t114 * t193 - t115 * t189 - t184 * t204;
t51 = Ifges(6,2) * t98 + Ifges(6,6) * t183 + t260;
t93 = Ifges(6,4) * t98;
t52 = Ifges(6,1) * t99 + Ifges(6,5) * t183 + t93;
t69 = -t164 * t184 - t86;
t199 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t15 * t261 + t51 * t262 - t69 * (mrSges(6,1) * t99 + mrSges(6,2) * t98) + Ifges(6,3) * t181 - t99 * (Ifges(6,1) * t98 - t260) / 0.2e1 + Ifges(6,6) * t42 + Ifges(6,5) * t41 - t183 * (Ifges(6,5) * t98 - Ifges(6,6) * t99) / 0.2e1 - (-Ifges(6,2) * t99 + t52 + t93) * t98 / 0.2e1;
t38 = -pkin(4) * t114 + t55;
t74 = qJD(4) * t124 + t203;
t75 = -qJD(4) * t125 - t204;
t76 = -pkin(3) * t184 - t86;
t96 = t184 * t208 + t238;
t142 = Ifges(5,4) * t234;
t97 = Ifges(5,1) * t235 + t142 + t239;
t198 = t76 * t211 * qJD(4) + t183 * (Ifges(6,5) * t74 + Ifges(6,6) * t75) / 0.2e1 + t55 * t138 + t121 * mrSges(3,1) - t122 * mrSges(3,2) + t98 * (Ifges(6,4) * t74 + Ifges(6,2) * t75) / 0.2e1 + t74 * t52 / 0.2e1 + t75 * t51 / 0.2e1 + t69 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + t65 * mrSges(4,1) + (Ifges(5,1) * t194 - t248) * t218 / 0.2e1 - t96 * t225 / 0.2e1 + mrSges(5,3) * t240 + (Ifges(6,1) * t74 + Ifges(6,4) * t75) * t262 + qJD(4) ^ 2 * (Ifges(5,5) * t194 - Ifges(5,6) * t190) / 0.2e1 + (t268 / 0.2e1 + t269 + t208 / 0.2e1) * t114 + t115 * (t190 * Ifges(5,1) + t247) + (t184 * (-Ifges(5,2) * t190 + t247) + t97) * t224 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t182 + (-t15 * t74 + t16 * t75) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t190 + Ifges(5,6) * t194) + (t38 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t41 + Ifges(6,4) * t42 + Ifges(6,5) * t181) * t125 + (-t38 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t41 + Ifges(6,2) * t42 + Ifges(6,6) * t181) * t124;
t196 = cos(qJ(1));
t180 = t196 * pkin(1);
t160 = -pkin(3) - t257;
t149 = t157 * pkin(3);
t137 = -t164 - t257;
t113 = t138 * t184;
t103 = qJD(1) * t205;
t92 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t115;
t91 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t114;
t90 = t101 - t256;
t88 = t104 + t221;
t81 = mrSges(6,1) * t183 - t252;
t80 = -mrSges(6,2) * t183 + t261;
t70 = -mrSges(5,1) * t114 + mrSges(5,2) * t115;
t62 = -mrSges(6,1) * t98 + mrSges(6,2) * t99;
t47 = t194 * t215 - t231;
t46 = t190 * t215 + t230;
t35 = -mrSges(6,2) * t181 + mrSges(6,3) * t42;
t34 = mrSges(6,1) * t181 - mrSges(6,3) * t41;
t19 = t193 * t56 - t244;
t18 = -t189 * t56 - t241;
t10 = -mrSges(6,1) * t42 + mrSges(6,2) * t41;
t6 = -qJD(5) * t44 - t189 * t46 + t193 * t47;
t5 = qJD(5) * t43 + t189 * t47 + t193 * t46;
t1 = [(-mrSges(2,1) * t192 - mrSges(2,2) * t196 - m(5) * (t149 + t228) - m(6) * (t178 + t220) - m(4) * t228 + t209 * t157 + t200) * g(3) + t202 * mrSges(5,3) + t206 * t106 + t198 + m(5) * (t101 * t55 + t104 * t76 + t68 * t230 - t67 * t231) + t101 * t70 + t104 * t113 + t90 * t10 + t5 * t80 + t6 * t81 + t88 * t62 + t43 * t34 + t44 * t35 + (-t106 * t184 - t108 * t182 - t66) * mrSges(4,2) + m(4) * (-t104 * t86 + t106 * t87 + t107 * t65 + t108 * t66) + m(6) * (t15 * t6 + t16 * t5 + t3 * t44 + t38 * t90 + t4 * t43 + t69 * t88) + (-mrSges(2,1) * t196 + mrSges(2,2) * t192 - m(6) * (t180 + t212) - m(5) * (t180 + t219) - m(4) * (t163 + t180) + t209 * t158 + t201) * g(2) + (m(5) * (-t224 * t67 - t225 * t68 + t240 - t242) + t194 * t91 - t190 * t92 - t132 * t224 - t133 * t225) * t102 + ((-t191 * t182 - t184 * t227) * mrSges(3,2) + (-qJD(2) * t191 * t184 + t195 * t182) * mrSges(3,1) + (-g(2) * t196 - g(3) * t192 + t121 * t195 + t122 * t191) * m(3)) * pkin(1) + Ifges(2,3) * qJDD(1) + (-t104 * t184 + t107 * t182) * mrSges(4,1); (-t62 - t113) * t103 + (-t182 * t258 - t66) * mrSges(4,2) + t198 + t160 * t70 + t137 * t10 + t71 * t34 + t72 * t35 + t265 * t81 + t264 * t80 + (-t21 * mrSges(5,3) + t105 * t132 - t159 * t92 + (g(2) * t158 + g(3) * t157) * mrSges(5,2) + (pkin(4) * t62 - t133 * t159 - t253) * qJD(4)) * t190 + (t157 * t245 + t200) * g(3) + (t103 * mrSges(4,1) + t105 * mrSges(4,2) + (mrSges(3,1) * t191 + mrSges(3,2) * t195) * t246) * t184 + (-t105 * t133 + t159 * t91 + (-t132 * t159 - t254) * qJD(4)) * t194 + (mrSges(6,2) * t237 + t201) * g(2) + t182 * mrSges(4,1) * t257 + (-t212 * g(2) - t220 * g(3) + t137 * t38 + t3 * t72 + t4 * t71 + (-t103 + t221) * t69 + t264 * t16 + t265 * t15) * m(6) + (t103 * t86 - t105 * t87 + (t187 * t66 + t188 * t65) * pkin(2) - t162 * g(3) - t163 * g(2)) * m(4) + (t160 * t55 + (-t149 - t162) * g(3) - t103 * t76 - t105 * t207 - t219 * g(2) + (t202 + t240) * t159) * m(5); m(4) * qJDD(3) + t124 * t34 + t125 * t35 + t190 * t91 + t194 * t92 + t74 * t80 + t75 * t81 + t206 * qJD(4) + m(5) * (qJD(4) * t207 + t190 * t20 + t194 * t21) + m(6) * (t124 * t4 + t125 * t3 + t15 * t75 + t16 * t74) + (-m(4) - m(5) - m(6)) * g(1); t199 + t68 * t132 - t67 * t133 + Ifges(5,6) * t114 + Ifges(5,5) * t115 - t19 * t80 - t18 * t81 - t20 * mrSges(5,2) + t21 * mrSges(5,1) + t16 * t252 + (-t210 + t211) * t255 + Ifges(5,3) * qJDD(4) + (t138 - t213) * g(1) - m(6) * (t15 * t18 + t16 * t19) + ((-t76 * mrSges(5,2) - t239 / 0.2e1 - t142 / 0.2e1 - t97 / 0.2e1 + t254) * t194 + (-t76 * mrSges(5,1) + t238 / 0.2e1 + t96 / 0.2e1 + t253 + (t269 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t194) * t184 + (-m(6) * t69 - t62) * pkin(4)) * t190) * t184 + (t189 * t35 + t193 * t34 + (-g(1) * t194 + t189 * t3 + t193 * t4 + (-g(3) * t158 + t255) * t190) * m(6) + (-t189 * t81 + t193 * t80 + (-t15 * t189 + t16 * t193) * m(6)) * qJD(5)) * pkin(4) + (-t158 * t211 - t229) * g(3); t199 + (t81 + t252) * t16 - t210 * t255 - g(1) * t213 - g(3) * t229 - t15 * t80;];
tau = t1;
