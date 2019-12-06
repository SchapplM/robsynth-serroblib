% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:07
% EndTime: 2019-12-05 18:16:14
% DurationCPUTime: 2.25s
% Computational Cost: add. (3525->338), mult. (6055->448), div. (0->0), fcn. (3569->16), ass. (0->172)
t169 = sin(qJ(4));
t227 = Ifges(5,4) * t169;
t249 = t227 / 0.2e1;
t173 = cos(qJ(4));
t226 = Ifges(5,4) * t173;
t248 = t173 * Ifges(5,2);
t167 = cos(pkin(9));
t146 = pkin(1) * t167 + pkin(2);
t166 = sin(pkin(9));
t236 = pkin(1) * t166;
t204 = qJD(1) * t236;
t247 = -qJD(3) * t204 + t146 * qJDD(1);
t127 = t146 * qJD(1);
t246 = qJD(3) * t127 + qJDD(1) * t236;
t176 = -pkin(8) - pkin(7);
t240 = m(5) * pkin(7);
t245 = m(6) * t176 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - t240;
t168 = sin(qJ(5));
t172 = cos(qJ(5));
t114 = t168 * t173 + t169 * t172;
t200 = qJD(4) * t176;
t119 = t169 * t200;
t120 = t173 * t200;
t131 = t176 * t169;
t159 = t173 * pkin(8);
t132 = pkin(7) * t173 + t159;
t86 = t131 * t168 + t132 * t172;
t170 = sin(qJ(3));
t174 = cos(qJ(3));
t88 = t127 * t174 - t170 * t204;
t243 = -t86 * qJD(5) + t114 * t88 - t119 * t168 + t120 * t172;
t113 = -t168 * t169 + t172 * t173;
t85 = t131 * t172 - t132 * t168;
t242 = t85 * qJD(5) - t113 * t88 + t119 * t172 + t120 * t168;
t103 = t146 * t174 - t170 * t236;
t104 = t170 * t146 + t174 * t236;
t241 = m(5) * pkin(3);
t163 = qJD(1) + qJD(3);
t100 = t114 * t163;
t238 = t100 / 0.2e1;
t99 = t113 * t163;
t237 = mrSges(6,3) * t99;
t235 = pkin(4) * t173;
t164 = qJ(1) + pkin(9);
t155 = qJ(3) + t164;
t145 = cos(t155);
t234 = g(3) * t145;
t156 = t173 * qJD(2);
t89 = t127 * t170 + t174 * t204;
t77 = pkin(7) * t163 + t89;
t68 = -t169 * t77 + t156;
t233 = t68 * mrSges(5,3);
t209 = qJD(2) * t169;
t69 = t173 * t77 + t209;
t232 = t69 * mrSges(5,3);
t102 = pkin(7) + t104;
t231 = -pkin(8) - t102;
t230 = mrSges(5,2) * t169;
t165 = qJ(4) + qJ(5);
t157 = sin(t165);
t229 = mrSges(6,2) * t157;
t158 = cos(t165);
t228 = mrSges(6,2) * t158;
t225 = Ifges(6,4) * t100;
t224 = t100 * mrSges(6,3);
t143 = t158 * mrSges(6,1);
t198 = pkin(8) * t163 + t77;
t63 = t198 * t173 + t209;
t223 = t168 * t63;
t95 = t103 * qJD(3);
t222 = t169 * t95;
t207 = qJD(4) * t169;
t161 = qJDD(1) + qJDD(3);
t57 = t247 * t170 + t246 * t174;
t51 = pkin(7) * t161 + t57;
t17 = qJD(4) * t156 + t169 * qJDD(2) + t173 * t51 - t77 * t207;
t221 = t17 * t173;
t220 = t172 * t63;
t219 = t173 * mrSges(5,1);
t218 = t173 * t95;
t18 = -t69 * qJD(4) + t173 * qJDD(2) - t169 * t51;
t217 = t18 * t169;
t216 = Ifges(5,5) * qJD(4);
t215 = Ifges(5,6) * qJD(4);
t144 = sin(t155);
t214 = t144 * t157;
t212 = t163 * t169;
t211 = t163 * t173;
t210 = mrSges(6,1) * t214 + t144 * t228;
t206 = qJD(4) * t173;
t205 = m(4) + m(5) + m(6);
t203 = pkin(4) * t207;
t202 = m(3) + t205;
t148 = pkin(3) + t235;
t199 = t163 * t207;
t196 = qJD(4) * t231;
t195 = t143 - t229;
t101 = -pkin(3) - t103;
t193 = t205 * pkin(2) + mrSges(3,1);
t128 = -t219 + t230;
t192 = mrSges(5,1) * t169 + mrSges(5,2) * t173;
t191 = -mrSges(6,1) * t157 - t228;
t190 = t227 + t248;
t62 = -t198 * t169 + t156;
t59 = qJD(4) * pkin(4) + t62;
t19 = t172 * t59 - t223;
t20 = t168 * t59 + t220;
t80 = t231 * t169;
t81 = t102 * t173 + t159;
t42 = -t168 * t81 + t172 * t80;
t43 = t168 * t80 + t172 * t81;
t189 = -t169 * t68 + t173 * t69;
t121 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t212;
t122 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t211;
t188 = -t169 * t121 + t173 * t122;
t187 = t202 * pkin(1) + mrSges(2,1);
t58 = -t246 * t170 + t247 * t174;
t186 = t114 * qJD(5);
t185 = t113 * qJD(5);
t96 = t104 * qJD(3);
t184 = m(6) * t148 + mrSges(4,1) + t143 + t241;
t52 = -pkin(3) * t161 - t58;
t183 = -t217 + (-t169 * t69 - t173 * t68) * qJD(4);
t182 = t184 + t219;
t181 = (-t229 - t230) * t145 - t245 * t144;
t180 = -mrSges(6,2) * t214 - t144 * t230 + t245 * t145;
t160 = qJDD(4) + qJDD(5);
t162 = qJD(4) + qJD(5);
t107 = t161 * t169 + t163 * t206;
t11 = qJDD(4) * pkin(4) - pkin(8) * t107 + t18;
t106 = t161 * t173 - t199;
t13 = pkin(8) * t106 + t17;
t3 = t19 * qJD(5) + t11 * t168 + t13 * t172;
t4 = -t20 * qJD(5) + t11 * t172 - t13 * t168;
t40 = t106 * t168 + t107 * t172 + t163 * t185;
t41 = t106 * t172 - t107 * t168 - t163 * t186;
t55 = Ifges(6,2) * t99 + Ifges(6,6) * t162 + t225;
t92 = Ifges(6,4) * t99;
t56 = Ifges(6,1) * t100 + Ifges(6,5) * t162 + t92;
t71 = -t148 * t163 - t88;
t179 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t19 * t237 + t55 * t238 - t71 * (mrSges(6,1) * t100 + mrSges(6,2) * t99) + Ifges(6,3) * t160 - t100 * (Ifges(6,1) * t99 - t225) / 0.2e1 + Ifges(6,6) * t41 + Ifges(6,5) * t40 - t162 * (Ifges(6,5) * t99 - Ifges(6,6) * t100) / 0.2e1 - (-Ifges(6,2) * t100 + t56 + t92) * t99 / 0.2e1;
t31 = -pkin(4) * t106 + t52;
t72 = t113 * qJD(4) + t185;
t73 = -t114 * qJD(4) - t186;
t76 = -pkin(3) * t163 - t88;
t97 = t190 * t163 + t215;
t136 = Ifges(5,4) * t211;
t98 = Ifges(5,1) * t212 + t136 + t216;
t178 = t76 * t192 * qJD(4) + Ifges(4,3) * t161 + t162 * (Ifges(6,5) * t72 + Ifges(6,6) * t73) / 0.2e1 + t52 * t128 + t99 * (Ifges(6,4) * t72 + Ifges(6,2) * t73) / 0.2e1 + t73 * t55 / 0.2e1 + t71 * (-mrSges(6,1) * t73 + mrSges(6,2) * t72) + t72 * t56 / 0.2e1 + t58 * mrSges(4,1) + (Ifges(6,1) * t72 + Ifges(6,4) * t73) * t238 + mrSges(5,3) * t221 - t97 * t207 / 0.2e1 + (Ifges(5,1) * t173 - t227) * t199 / 0.2e1 + qJD(4) ^ 2 * (Ifges(5,5) * t173 - Ifges(5,6) * t169) / 0.2e1 + (t248 / 0.2e1 + t249 + t190 / 0.2e1) * t106 + t107 * (t169 * Ifges(5,1) + t226) + (t98 + t163 * (-Ifges(5,2) * t169 + t226)) * t206 / 0.2e1 + (-t19 * t72 + t20 * t73) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t169 + Ifges(5,6) * t173) + (t31 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t40 + Ifges(6,4) * t41 + Ifges(6,5) * t160) * t114 + (-t31 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t40 + Ifges(6,2) * t41 + Ifges(6,6) * t160) * t113;
t175 = cos(qJ(1));
t171 = sin(qJ(1));
t152 = cos(t164);
t151 = sin(t164);
t105 = t128 * t163;
t91 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t107;
t90 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t106;
t87 = t101 - t235;
t82 = t96 + t203;
t79 = mrSges(6,1) * t162 - t224;
t78 = -mrSges(6,2) * t162 + t237;
t67 = -mrSges(5,1) * t106 + mrSges(5,2) * t107;
t65 = -mrSges(6,1) * t99 + mrSges(6,2) * t100;
t50 = t173 * t196 - t222;
t49 = t169 * t196 + t218;
t34 = -mrSges(6,2) * t160 + mrSges(6,3) * t41;
t33 = mrSges(6,1) * t160 - mrSges(6,3) * t40;
t22 = t172 * t62 - t223;
t21 = -t168 * t62 - t220;
t10 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t6 = -t43 * qJD(5) - t168 * t49 + t172 * t50;
t5 = t42 * qJD(5) + t168 * t50 + t172 * t49;
t1 = [m(4) * (t103 * t58 + t104 * t57 - t88 * t96 + t89 * t95) + m(6) * (t19 * t6 + t20 * t5 + t3 * t43 + t31 * t87 + t4 * t42 + t71 * t82) + (-mrSges(2,2) * t171 - mrSges(3,2) * t151 + t182 * t145 + t193 * t152 + t187 * t175 + t181) * g(2) + t101 * t67 + t96 * t105 + t87 * t10 + t5 * t78 + t6 * t79 + t82 * t65 + t178 + t42 * t33 + t43 * t34 + (m(5) * (-t68 * t206 - t69 * t207 - t217 + t221) + t173 * t90 - t169 * t91 - t121 * t206 - t122 * t207) * t102 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t167 * mrSges(3,1) - 0.2e1 * t166 * mrSges(3,2) + m(3) * (t166 ^ 2 + t167 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-t104 * t161 - t95 * t163 - t57) * mrSges(4,2) + (t103 * t161 - t96 * t163) * mrSges(4,1) + m(5) * (t101 * t52 + t69 * t218 - t68 * t222 + t76 * t96) + t188 * t95 + t183 * mrSges(5,3) + (mrSges(2,2) * t175 + mrSges(3,2) * t152 + t182 * t144 + t193 * t151 + t187 * t171 + t180) * g(3); t113 * t33 + t114 * t34 + t169 * t90 + t173 * t91 + t72 * t78 + t73 * t79 + (m(3) + m(4)) * qJDD(2) + t188 * qJD(4) + m(5) * (t189 * qJD(4) + t169 * t17 + t173 * t18) + m(6) * (t113 * t4 + t114 * t3 + t19 * t73 + t20 * t72) - t202 * g(1); (t163 * mrSges(4,1) - t105 - t65) * t89 - t148 * t10 + t85 * t33 + t86 * t34 + t178 - pkin(3) * t67 + (t183 + t221) * t240 - t52 * t241 + (-t18 * mrSges(5,3) - pkin(7) * t91 + t88 * t121 + (pkin(4) * t65 - pkin(7) * t122 - t232) * qJD(4)) * t169 + (-qJD(4) * t233 - t88 * t122 + (-qJD(4) * t121 + t90) * pkin(7) + (g(2) * t145 + g(3) * t144) * mrSges(5,1)) * t173 + (t184 * t145 + t181) * g(2) + t243 * t79 + t242 * t78 + (t184 * t144 + t180) * g(3) - m(5) * (t189 * t88 + t76 * t89) + (t163 * t88 - t57) * mrSges(4,2) + (-t148 * t31 + t3 * t86 + t4 * t85 + (t203 - t89) * t71 + t242 * t20 + t243 * t19) * m(6); t20 * t224 + t69 * t121 - t68 * t122 + Ifges(5,6) * t106 + Ifges(5,5) * t107 - t22 * t78 - t21 * t79 + t179 - t17 * mrSges(5,2) + t18 * mrSges(5,1) - m(6) * (t19 * t21 + t20 * t22) + ((-t136 / 0.2e1 - t98 / 0.2e1 - t76 * mrSges(5,2) - t216 / 0.2e1 + t233) * t173 + (t97 / 0.2e1 - t76 * mrSges(5,1) + t215 / 0.2e1 + t232 + (t249 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t173) * t163 + (-m(6) * t71 - t65) * pkin(4)) * t169) * t163 + (t168 * t34 + t172 * t33 + (t168 * t3 + t172 * t4 - g(1) * t173 + (-g(2) * t144 + t234) * t169) * m(6) + (-t168 * t79 + t172 * t78 + (-t168 * t19 + t172 * t20) * m(6)) * qJD(5)) * pkin(4) + (-t191 + t192) * t234 + Ifges(5,3) * qJDD(4) + (t128 - t195) * g(1) + (-t192 * t144 - t210) * g(2); (t79 + t224) * t20 - t191 * t234 - g(1) * t195 - g(2) * t210 - t19 * t78 + t179;];
tau = t1;
