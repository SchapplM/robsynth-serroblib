% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:41
% EndTime: 2019-12-05 15:12:44
% DurationCPUTime: 2.20s
% Computational Cost: add. (32410->159), mult. (42143->208), div. (0->0), fcn. (29126->10), ass. (0->78)
t152 = sin(pkin(8));
t154 = cos(pkin(8));
t144 = -t154 * g(1) - t152 * g(2);
t150 = -g(3) + qJDD(1);
t151 = sin(pkin(9));
t153 = cos(pkin(9));
t128 = -t151 * t144 + t153 * t150;
t129 = t153 * t144 + t151 * t150;
t157 = sin(qJ(3));
t160 = cos(qJ(3));
t123 = t160 * t128 - t157 * t129;
t120 = qJDD(3) * pkin(3) + t123;
t124 = t157 * t128 + t160 * t129;
t161 = qJD(3) ^ 2;
t121 = -t161 * pkin(3) + t124;
t156 = sin(qJ(4));
t159 = cos(qJ(4));
t117 = t156 * t120 + t159 * t121;
t148 = qJD(3) + qJD(4);
t146 = t148 ^ 2;
t147 = qJDD(3) + qJDD(4);
t114 = -t146 * pkin(4) + t147 * pkin(7) + t117;
t143 = t152 * g(1) - t154 * g(2);
t142 = qJDD(2) - t143;
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t111 = -t155 * t114 + t158 * t142;
t112 = t158 * t114 + t155 * t142;
t131 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t155 + Ifges(6,2) * t158) * t148;
t132 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t155 + Ifges(6,4) * t158) * t148;
t175 = qJD(5) * t148;
t136 = t155 * t147 + t158 * t175;
t137 = t158 * t147 - t155 * t175;
t178 = mrSges(6,1) * t111 - mrSges(6,2) * t112 + Ifges(6,5) * t136 + Ifges(6,6) * t137 + Ifges(6,3) * qJDD(5) + (t131 * t155 - t132 * t158) * t148;
t116 = t159 * t120 - t156 * t121;
t113 = -t147 * pkin(4) - t146 * pkin(7) - t116;
t177 = t148 * t155;
t139 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t177;
t176 = t148 * t158;
t140 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t176;
t165 = -m(6) * t113 + t137 * mrSges(6,1) - t136 * mrSges(6,2) - t139 * t177 + t140 * t176;
t103 = m(5) * t116 + t147 * mrSges(5,1) - t146 * mrSges(5,2) + t165;
t135 = (-mrSges(6,1) * t158 + mrSges(6,2) * t155) * t148;
t107 = m(6) * t111 + qJDD(5) * mrSges(6,1) - t136 * mrSges(6,3) + qJD(5) * t140 - t135 * t177;
t108 = m(6) * t112 - qJDD(5) * mrSges(6,2) + t137 * mrSges(6,3) - qJD(5) * t139 + t135 * t176;
t170 = -t155 * t107 + t158 * t108;
t92 = m(5) * t117 - t146 * mrSges(5,1) - t147 * mrSges(5,2) + t170;
t89 = t159 * t103 + t156 * t92;
t86 = m(4) * t123 + qJDD(3) * mrSges(4,1) - t161 * mrSges(4,2) + t89;
t171 = -t156 * t103 + t159 * t92;
t87 = m(4) * t124 - t161 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t171;
t80 = t157 * t87 + t160 * t86;
t96 = t158 * t107 + t155 * t108;
t174 = m(5) * t142 + t96;
t78 = m(3) * t128 + t80;
t172 = -t157 * t86 + t160 * t87;
t79 = m(3) * t129 + t172;
t173 = -t151 * t78 + t153 * t79;
t168 = (-m(3) - m(4)) * t142 - t174;
t130 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t155 + Ifges(6,6) * t158) * t148;
t100 = -mrSges(6,1) * t113 + mrSges(6,3) * t112 + Ifges(6,4) * t136 + Ifges(6,2) * t137 + Ifges(6,6) * qJDD(5) + qJD(5) * t132 - t130 * t177;
t101 = mrSges(6,2) * t113 - mrSges(6,3) * t111 + Ifges(6,1) * t136 + Ifges(6,4) * t137 + Ifges(6,5) * qJDD(5) - qJD(5) * t131 + t130 * t176;
t81 = mrSges(5,2) * t142 - mrSges(5,3) * t116 + Ifges(5,5) * t147 - t146 * Ifges(5,6) - pkin(7) * t96 - t155 * t100 + t158 * t101;
t82 = -mrSges(5,1) * t142 + mrSges(5,3) * t117 + t146 * Ifges(5,5) + Ifges(5,6) * t147 - pkin(4) * t96 - t178;
t75 = -mrSges(4,1) * t142 + mrSges(4,3) * t124 + t161 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t174 + pkin(6) * t171 + t156 * t81 + t159 * t82;
t76 = mrSges(4,2) * t142 - mrSges(4,3) * t123 + Ifges(4,5) * qJDD(3) - t161 * Ifges(4,6) - pkin(6) * t89 - t156 * t82 + t159 * t81;
t68 = -mrSges(3,1) * t142 + mrSges(3,3) * t129 + t157 * t76 + t160 * t75 - pkin(2) * (m(4) * t142 + t174) + pkin(5) * t172;
t70 = mrSges(3,2) * t142 - mrSges(3,3) * t128 - pkin(5) * t80 - t157 * t75 + t160 * t76;
t167 = mrSges(2,1) * t143 - mrSges(2,2) * t144 + pkin(1) * t168 + qJ(2) * t173 + t151 * t70 + t153 * t68;
t166 = mrSges(5,1) * t116 - mrSges(5,2) * t117 + Ifges(5,3) * t147 + pkin(4) * t165 + pkin(7) * t170 + t158 * t100 + t155 * t101;
t163 = mrSges(4,1) * t123 - mrSges(4,2) * t124 + Ifges(4,3) * qJDD(3) + pkin(3) * t89 + t166;
t162 = mrSges(3,1) * t128 - mrSges(3,2) * t129 + pkin(2) * t80 + t163;
t93 = m(2) * t143 + t168;
t74 = t151 * t79 + t153 * t78;
t72 = m(2) * t144 + t173;
t71 = -mrSges(2,1) * t150 + mrSges(2,3) * t144 - pkin(1) * t74 - t162;
t66 = mrSges(2,2) * t150 - mrSges(2,3) * t143 - qJ(2) * t74 - t151 * t68 + t153 * t70;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t154 * t66 - t152 * t71 - qJ(1) * (t152 * t72 + t154 * t93), t66, t70, t76, t81, t101; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t152 * t66 + t154 * t71 + qJ(1) * (-t152 * t93 + t154 * t72), t71, t68, t75, t82, t100; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t167, t167, t162, t163, t166, t178;];
m_new = t1;
