% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:48
% EndTime: 2019-12-05 15:59:51
% DurationCPUTime: 1.90s
% Computational Cost: add. (22293->215), mult. (41091->266), div. (0->0), fcn. (23850->8), ass. (0->94)
t192 = sin(pkin(8));
t218 = cos(pkin(8));
t174 = -t218 * g(1) - t192 * g(2);
t189 = -g(3) + qJDD(1);
t195 = sin(qJ(2));
t198 = cos(qJ(2));
t153 = t198 * t174 + t195 * t189;
t210 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t153;
t224 = m(3) + m(4);
t223 = -pkin(2) - pkin(6);
t199 = qJD(2) ^ 2;
t222 = pkin(4) * t199;
t221 = mrSges(3,1) - mrSges(4,2);
t220 = (Ifges(3,5) - Ifges(4,4));
t219 = -Ifges(3,6) + Ifges(4,5);
t152 = -t195 * t174 + t198 * t189;
t208 = -t199 * qJ(3) + qJDD(3) - t152;
t146 = t223 * qJDD(2) + t208;
t197 = cos(qJ(4));
t143 = t197 * t146;
t194 = sin(qJ(4));
t215 = qJD(2) * qJD(4);
t171 = t197 * qJDD(2) - t194 * t215;
t173 = t192 * g(1) - t218 * g(2);
t124 = (qJDD(4) * pkin(4)) - t171 * pkin(7) + t143 + (-pkin(7) * t215 - t197 * t222 + t173) * t194;
t134 = t194 * t146 - t197 * t173;
t170 = -t194 * qJDD(2) - t197 * t215;
t216 = qJD(2) * t197;
t177 = (qJD(4) * pkin(4)) - pkin(7) * t216;
t188 = t194 ^ 2;
t125 = t170 * pkin(7) - qJD(4) * t177 - t188 * t222 + t134;
t193 = sin(qJ(5));
t196 = cos(qJ(5));
t122 = t196 * t124 - t193 * t125;
t159 = (-t193 * t197 - t194 * t196) * qJD(2);
t136 = t159 * qJD(5) + t193 * t170 + t196 * t171;
t160 = (-t193 * t194 + t196 * t197) * qJD(2);
t141 = -t159 * mrSges(6,1) + t160 * mrSges(6,2);
t183 = qJD(4) + qJD(5);
t150 = -t183 * mrSges(6,2) + t159 * mrSges(6,3);
t182 = qJDD(4) + qJDD(5);
t119 = m(6) * t122 + t182 * mrSges(6,1) - t136 * mrSges(6,3) - t160 * t141 + t183 * t150;
t123 = t193 * t124 + t196 * t125;
t135 = -t160 * qJD(5) + t196 * t170 - t193 * t171;
t151 = t183 * mrSges(6,1) - t160 * mrSges(6,3);
t120 = m(6) * t123 - t182 * mrSges(6,2) + t135 * mrSges(6,3) + t159 * t141 - t183 * t151;
t109 = t196 * t119 + t193 * t120;
t133 = t194 * t173 + t143;
t169 = (mrSges(5,1) * t194 + mrSges(5,2) * t197) * qJD(2);
t217 = qJD(2) * t194;
t175 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t217;
t106 = m(5) * t133 + qJDD(4) * mrSges(5,1) - t171 * mrSges(5,3) + qJD(4) * t175 - t169 * t216 + t109;
t176 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t216;
t211 = -t193 * t119 + t196 * t120;
t107 = m(5) * t134 - qJDD(4) * mrSges(5,2) + t170 * mrSges(5,3) - qJD(4) * t176 - t169 * t217 + t211;
t103 = -t194 * t106 + t197 * t107;
t145 = t223 * t199 + t210;
t127 = t177 * t216 - t170 * pkin(4) + (-pkin(7) * t188 + t223) * t199 + t210;
t207 = m(6) * t127 - t135 * mrSges(6,1) + t136 * mrSges(6,2) - t159 * t150 + t160 * t151;
t115 = -m(5) * t145 + t170 * mrSges(5,1) - t171 * mrSges(5,2) - t175 * t217 - t176 * t216 - t207;
t147 = (t199 * pkin(2)) - t210;
t114 = -m(4) * t147 + (t199 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t115;
t113 = m(3) * t153 - (t199 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t114;
t102 = t197 * t106 + t194 * t107;
t149 = -qJDD(2) * pkin(2) + t208;
t206 = -m(4) * t149 + t199 * mrSges(4,3) - t102;
t97 = m(3) * t152 - t199 * mrSges(3,2) + t221 * qJDD(2) + t206;
t212 = t198 * t113 - t195 * t97;
t101 = -m(4) * t173 + t103;
t137 = Ifges(6,5) * t160 + Ifges(6,6) * t159 + Ifges(6,3) * t183;
t139 = Ifges(6,1) * t160 + Ifges(6,4) * t159 + Ifges(6,5) * t183;
t110 = -mrSges(6,1) * t127 + mrSges(6,3) * t123 + Ifges(6,4) * t136 + Ifges(6,2) * t135 + Ifges(6,6) * t182 - t160 * t137 + t183 * t139;
t138 = Ifges(6,4) * t160 + Ifges(6,2) * t159 + Ifges(6,6) * t183;
t111 = mrSges(6,2) * t127 - mrSges(6,3) * t122 + Ifges(6,1) * t136 + Ifges(6,4) * t135 + Ifges(6,5) * t182 + t159 * t137 - t183 * t138;
t156 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t197 - Ifges(5,6) * t194) * qJD(2);
t158 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t197 - Ifges(5,4) * t194) * qJD(2);
t92 = -mrSges(5,1) * t145 + mrSges(5,3) * t134 + Ifges(5,4) * t171 + Ifges(5,2) * t170 + Ifges(5,6) * qJDD(4) - pkin(4) * t207 + pkin(7) * t211 + qJD(4) * t158 + t196 * t110 + t193 * t111 - t156 * t216;
t157 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t197 - Ifges(5,2) * t194) * qJD(2);
t96 = mrSges(5,2) * t145 - mrSges(5,3) * t133 + Ifges(5,1) * t171 + Ifges(5,4) * t170 + Ifges(5,5) * qJDD(4) - pkin(7) * t109 - qJD(4) * t157 - t193 * t110 + t196 * t111 - t156 * t217;
t203 = -mrSges(4,1) * t147 - pkin(3) * t115 - pkin(6) * t103 - t194 * t96 - t197 * t92;
t88 = mrSges(3,3) * t153 - pkin(2) * t101 - t219 * qJDD(2) + t221 * t173 + (t220 * t199) + t203;
t204 = mrSges(6,1) * t122 - mrSges(6,2) * t123 + Ifges(6,5) * t136 + Ifges(6,6) * t135 + Ifges(6,3) * t182 + t160 * t138 - t159 * t139;
t202 = mrSges(5,1) * t133 - mrSges(5,2) * t134 + Ifges(5,5) * t171 + Ifges(5,6) * t170 + (Ifges(5,3) * qJDD(4)) + pkin(4) * t109 + t157 * t216 + t158 * t217 + t204;
t200 = mrSges(4,1) * t149 + pkin(3) * t102 + t202;
t90 = t200 + (-mrSges(3,2) + mrSges(4,3)) * t173 + t219 * t199 + t220 * qJDD(2) - mrSges(3,3) * t152 - qJ(3) * t101;
t209 = -mrSges(2,2) * t174 + pkin(5) * t212 + t195 * t90 + t198 * t88 + pkin(1) * (t224 * t173 - t103) + mrSges(2,1) * t173;
t205 = mrSges(4,2) * t149 - mrSges(4,3) * t147 + Ifges(4,1) * qJDD(2) - pkin(6) * t102 - t194 * t92 + t197 * t96;
t201 = mrSges(3,1) * t152 - mrSges(3,2) * t153 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) + t206) + qJ(3) * t114 + t205;
t99 = (m(2) + t224) * t173 - t103;
t94 = t195 * t113 + t198 * t97;
t91 = m(2) * t174 + t212;
t86 = -mrSges(2,1) * t189 + mrSges(2,3) * t174 - pkin(1) * t94 - t201;
t85 = mrSges(2,2) * t189 - mrSges(2,3) * t173 - pkin(5) * t94 - t195 * t88 + t198 * t90;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t218 * t85 - t192 * t86 - qJ(1) * (t192 * t91 + t218 * t99), t85, t90, t205, t96, t111; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t192 * t85 + t218 * t86 + qJ(1) * (-t192 * t99 + t218 * t91), t86, t88, -mrSges(4,3) * t173 + Ifges(4,4) * qJDD(2) - (t199 * Ifges(4,5)) - t200, t92, t110; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t201, mrSges(4,2) * t173 + t199 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t203, t202, t204;];
m_new = t1;
