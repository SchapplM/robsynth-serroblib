% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:58
% EndTime: 2019-12-31 18:26:00
% DurationCPUTime: 1.64s
% Computational Cost: add. (27359->182), mult. (35894->222), div. (0->0), fcn. (12906->8), ass. (0->79)
t169 = qJD(1) ^ 2;
t164 = sin(qJ(1));
t167 = cos(qJ(1));
t144 = -t167 * g(1) - t164 * g(2);
t179 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t144;
t190 = -pkin(1) - pkin(2);
t125 = t190 * t169 + t179;
t143 = t164 * g(1) - t167 * g(2);
t178 = -t169 * qJ(2) + qJDD(2) - t143;
t128 = t190 * qJDD(1) + t178;
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t120 = -t163 * t125 + t166 * t128;
t149 = -qJDD(1) + qJDD(3);
t117 = t149 * pkin(3) + t120;
t121 = t166 * t125 + t163 * t128;
t150 = -qJD(1) + qJD(3);
t148 = t150 ^ 2;
t118 = -t148 * pkin(3) + t121;
t160 = sin(pkin(8));
t161 = cos(pkin(8));
t114 = t160 * t117 + t161 * t118;
t110 = -(t148 * pkin(4)) + t149 * pkin(7) + t114;
t157 = g(3) + qJDD(4);
t162 = sin(qJ(5));
t165 = cos(qJ(5));
t107 = -t162 * t110 + t165 * t157;
t108 = t165 * t110 + t162 * t157;
t130 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t162 + Ifges(6,2) * t165) * t150;
t131 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t162 + Ifges(6,4) * t165) * t150;
t186 = qJD(5) * t150;
t138 = t162 * t149 + t165 * t186;
t139 = t165 * t149 - t162 * t186;
t191 = mrSges(6,1) * t107 - mrSges(6,2) * t108 + Ifges(6,5) * t138 + Ifges(6,6) * t139 + Ifges(6,3) * qJDD(5) + (t130 * t162 - t131 * t165) * t150;
t189 = mrSges(2,1) + mrSges(3,1);
t137 = (-mrSges(6,1) * t165 + mrSges(6,2) * t162) * t150;
t187 = t150 * t165;
t141 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t187;
t188 = t150 * t162;
t105 = m(6) * t107 + qJDD(5) * mrSges(6,1) - t138 * mrSges(6,3) + qJD(5) * t141 - t137 * t188;
t140 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t188;
t106 = m(6) * t108 - qJDD(5) * mrSges(6,2) + t139 * mrSges(6,3) - qJD(5) * t140 + t137 * t187;
t184 = -t162 * t105 + t165 * t106;
t88 = m(5) * t114 - (t148 * mrSges(5,1)) - t149 * mrSges(5,2) + t184;
t113 = t161 * t117 - t160 * t118;
t109 = -t149 * pkin(4) - t148 * pkin(7) - t113;
t175 = -m(6) * t109 + t139 * mrSges(6,1) - t138 * mrSges(6,2) - t140 * t188 + t141 * t187;
t99 = m(5) * t113 + t149 * mrSges(5,1) - t148 * mrSges(5,2) + t175;
t85 = t160 * t88 + t161 * t99;
t92 = t165 * t105 + t162 * t106;
t185 = -t160 * t99 + t161 * t88;
t82 = m(4) * t120 + t149 * mrSges(4,1) - (t148 * mrSges(4,2)) + t85;
t83 = m(4) * t121 - t148 * mrSges(4,1) - t149 * mrSges(4,2) + t185;
t79 = -t163 * t82 + t166 * t83;
t183 = m(5) * t157 + t92;
t78 = t163 * t83 + t166 * t82;
t132 = -t169 * pkin(1) + t179;
t181 = m(3) * t132 + qJDD(1) * mrSges(3,3) + t79;
t129 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t162 + Ifges(6,6) * t165) * t150;
t96 = -mrSges(6,1) * t109 + mrSges(6,3) * t108 + Ifges(6,4) * t138 + Ifges(6,2) * t139 + Ifges(6,6) * qJDD(5) + qJD(5) * t131 - t129 * t188;
t97 = mrSges(6,2) * t109 - mrSges(6,3) * t107 + Ifges(6,1) * t138 + Ifges(6,4) * t139 + Ifges(6,5) * qJDD(5) - qJD(5) * t130 + t129 * t187;
t180 = mrSges(5,1) * t113 - mrSges(5,2) * t114 + Ifges(5,3) * t149 + pkin(4) * t175 + pkin(7) * t184 + t162 * t97 + t165 * t96;
t136 = -qJDD(1) * pkin(1) + t178;
t177 = -m(3) * t136 + qJDD(1) * mrSges(3,1) + t169 * mrSges(3,3) - t78;
t80 = mrSges(5,2) * t157 - mrSges(5,3) * t113 + Ifges(5,5) * t149 - (t148 * Ifges(5,6)) - pkin(7) * t92 - t162 * t96 + t165 * t97;
t84 = -mrSges(5,1) * t157 + mrSges(5,3) * t114 + t148 * Ifges(5,5) + Ifges(5,6) * t149 - pkin(4) * t92 - t191;
t71 = -mrSges(4,1) * g(3) + mrSges(4,3) * t121 + (t148 * Ifges(4,5)) + Ifges(4,6) * t149 - pkin(3) * t183 + qJ(4) * t185 + t160 * t80 + t161 * t84;
t73 = mrSges(4,2) * g(3) - mrSges(4,3) * t120 + Ifges(4,5) * t149 - t148 * Ifges(4,6) - qJ(4) * t85 - t160 * t84 + t161 * t80;
t176 = mrSges(3,2) * t136 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t169 * Ifges(3,6) - pkin(6) * t78 - t163 * t71 + t166 * t73;
t174 = mrSges(3,2) * t132 - pkin(2) * (-m(4) * g(3) - t183) - pkin(6) * t79 - t163 * t73 - t166 * t71;
t172 = mrSges(4,1) * t120 - mrSges(4,2) * t121 + Ifges(4,3) * t149 + pkin(3) * t85 + t180;
t171 = -mrSges(3,1) * t136 + mrSges(3,3) * t132 + Ifges(3,2) * qJDD(1) - pkin(2) * t78 - t172;
t170 = -mrSges(2,2) * t144 + mrSges(2,1) * t143 + Ifges(2,3) * qJDD(1) + t171 + qJ(2) * (-t169 * mrSges(3,1) + t181) + pkin(1) * t177;
t89 = (-m(3) - m(4)) * g(3) - t183;
t75 = m(2) * t143 + qJDD(1) * mrSges(2,1) - t169 * mrSges(2,2) + t177;
t74 = m(2) * t144 - qJDD(1) * mrSges(2,2) - t189 * t169 + t181;
t70 = -mrSges(2,2) * g(3) - mrSges(2,3) * t143 + Ifges(2,5) * qJDD(1) - t169 * Ifges(2,6) - qJ(2) * t89 + t176;
t69 = mrSges(2,3) * t144 - pkin(1) * t89 + (Ifges(3,4) + Ifges(2,5)) * t169 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t189 * g(3) + t174;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t167 * t70 - t164 * t69 - pkin(5) * (t164 * t74 + t167 * t75), t70, t176, t73, t80, t97; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t164 * t70 + t167 * t69 + pkin(5) * (-t164 * t75 + t167 * t74), t69, t171, t71, t84, t96; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t170, t170, -mrSges(3,1) * g(3) - t169 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t174, t172, t180, t191;];
m_new = t1;
