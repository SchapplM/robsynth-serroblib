% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:07
% EndTime: 2019-12-31 18:01:09
% DurationCPUTime: 1.61s
% Computational Cost: add. (26037->182), mult. (35894->222), div. (0->0), fcn. (12906->8), ass. (0->79)
t170 = qJD(1) ^ 2;
t165 = sin(qJ(1));
t168 = cos(qJ(1));
t144 = -t168 * g(1) - t165 * g(2);
t180 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t144;
t190 = -pkin(1) - pkin(2);
t125 = t170 * t190 + t180;
t143 = t165 * g(1) - t168 * g(2);
t179 = -t170 * qJ(2) + qJDD(2) - t143;
t128 = qJDD(1) * t190 + t179;
t161 = sin(pkin(8));
t162 = cos(pkin(8));
t120 = -t125 * t161 + t162 * t128;
t117 = -qJDD(1) * pkin(3) + t120;
t121 = t162 * t125 + t161 * t128;
t118 = -pkin(3) * t170 + t121;
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t114 = t164 * t117 + t167 * t118;
t150 = -qJD(1) + qJD(4);
t148 = t150 ^ 2;
t149 = -qJDD(1) + qJDD(4);
t110 = -(pkin(4) * t148) + pkin(7) * t149 + t114;
t158 = g(3) + qJDD(3);
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t107 = -t110 * t163 + t158 * t166;
t108 = t110 * t166 + t158 * t163;
t130 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t163 + Ifges(6,2) * t166) * t150;
t131 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t163 + Ifges(6,4) * t166) * t150;
t186 = qJD(5) * t150;
t138 = t149 * t163 + t166 * t186;
t139 = t149 * t166 - t163 * t186;
t191 = mrSges(6,1) * t107 - mrSges(6,2) * t108 + Ifges(6,5) * t138 + Ifges(6,6) * t139 + Ifges(6,3) * qJDD(5) + (t130 * t163 - t131 * t166) * t150;
t189 = mrSges(2,1) + mrSges(3,1);
t137 = (-mrSges(6,1) * t166 + mrSges(6,2) * t163) * t150;
t187 = t150 * t166;
t141 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t187;
t188 = t150 * t163;
t105 = m(6) * t107 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t138 + qJD(5) * t141 - t137 * t188;
t140 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t188;
t106 = m(6) * t108 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t139 - qJD(5) * t140 + t137 * t187;
t184 = -t105 * t163 + t166 * t106;
t88 = m(5) * t114 - (mrSges(5,1) * t148) - mrSges(5,2) * t149 + t184;
t113 = t117 * t167 - t118 * t164;
t109 = -pkin(4) * t149 - pkin(7) * t148 - t113;
t176 = -m(6) * t109 + t139 * mrSges(6,1) - mrSges(6,2) * t138 - t140 * t188 + t141 * t187;
t99 = m(5) * t113 + mrSges(5,1) * t149 - mrSges(5,2) * t148 + t176;
t85 = t164 * t88 + t167 * t99;
t92 = t166 * t105 + t163 * t106;
t82 = m(4) * t120 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t170 + t85;
t185 = -t164 * t99 + t167 * t88;
t83 = m(4) * t121 - mrSges(4,1) * t170 + qJDD(1) * mrSges(4,2) + t185;
t79 = -t161 * t82 + t162 * t83;
t78 = t161 * t83 + t162 * t82;
t132 = -pkin(1) * t170 + t180;
t182 = m(3) * t132 + qJDD(1) * mrSges(3,3) + t79;
t90 = (-m(4) - m(5)) * t158 - t92;
t129 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t163 + Ifges(6,6) * t166) * t150;
t96 = -mrSges(6,1) * t109 + mrSges(6,3) * t108 + Ifges(6,4) * t138 + Ifges(6,2) * t139 + Ifges(6,6) * qJDD(5) + qJD(5) * t131 - t129 * t188;
t97 = mrSges(6,2) * t109 - mrSges(6,3) * t107 + Ifges(6,1) * t138 + Ifges(6,4) * t139 + Ifges(6,5) * qJDD(5) - qJD(5) * t130 + t129 * t187;
t181 = mrSges(5,1) * t113 - mrSges(5,2) * t114 + Ifges(5,3) * t149 + pkin(4) * t176 + pkin(7) * t184 + t163 * t97 + t166 * t96;
t136 = -qJDD(1) * pkin(1) + t179;
t178 = -m(3) * t136 + qJDD(1) * mrSges(3,1) + t170 * mrSges(3,3) - t78;
t80 = mrSges(5,2) * t158 - mrSges(5,3) * t113 + Ifges(5,5) * t149 - (Ifges(5,6) * t148) - pkin(7) * t92 - t163 * t96 + t166 * t97;
t84 = -mrSges(5,1) * t158 + mrSges(5,3) * t114 + t148 * Ifges(5,5) + Ifges(5,6) * t149 - pkin(4) * t92 - t191;
t71 = -Ifges(4,6) * qJDD(1) + t170 * Ifges(4,5) - mrSges(4,1) * t158 + mrSges(4,3) * t121 + t164 * t80 + t167 * t84 - pkin(3) * (m(5) * t158 + t92) + pkin(6) * t185;
t73 = mrSges(4,2) * t158 - mrSges(4,3) * t120 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t170 - pkin(6) * t85 - t164 * t84 + t167 * t80;
t177 = mrSges(3,2) * t136 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t170 * Ifges(3,6) - qJ(3) * t78 - t161 * t71 + t162 * t73;
t175 = mrSges(3,2) * t132 - pkin(2) * t90 - qJ(3) * t79 - t161 * t73 - t162 * t71;
t173 = mrSges(4,1) * t120 - mrSges(4,2) * t121 - Ifges(4,3) * qJDD(1) + pkin(3) * t85 + t181;
t172 = -mrSges(3,1) * t136 + mrSges(3,3) * t132 + Ifges(3,2) * qJDD(1) - pkin(2) * t78 - t173;
t171 = -mrSges(2,2) * t144 + mrSges(2,1) * t143 + Ifges(2,3) * qJDD(1) + t172 + qJ(2) * (-mrSges(3,1) * t170 + t182) + pkin(1) * t178;
t89 = -m(3) * g(3) + t90;
t75 = m(2) * t143 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t170 + t178;
t74 = m(2) * t144 - qJDD(1) * mrSges(2,2) - t170 * t189 + t182;
t70 = -mrSges(2,2) * g(3) - mrSges(2,3) * t143 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t170 - qJ(2) * t89 + t177;
t69 = mrSges(2,3) * t144 - pkin(1) * t89 + (Ifges(3,4) + Ifges(2,5)) * t170 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t189 * g(3) + t175;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t168 * t70 - t165 * t69 - pkin(5) * (t165 * t74 + t168 * t75), t70, t177, t73, t80, t97; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t165 * t70 + t168 * t69 + pkin(5) * (-t165 * t75 + t168 * t74), t69, t172, t71, t84, t96; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t171, t171, -mrSges(3,1) * g(3) - t170 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t175, t173, t181, t191;];
m_new = t1;
