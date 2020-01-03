% Calculate vector of cutting torques with Newton-Euler for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:47
% EndTime: 2019-12-31 16:27:48
% DurationCPUTime: 0.83s
% Computational Cost: add. (5854->178), mult. (11119->227), div. (0->0), fcn. (5457->6), ass. (0->67)
t152 = sin(pkin(6));
t153 = cos(pkin(6));
t138 = t152 * g(1) - t153 * g(2);
t139 = -t153 * g(1) - t152 * g(2);
t155 = sin(qJ(2));
t157 = cos(qJ(2));
t108 = t155 * t138 + t157 * t139;
t159 = qJD(2) ^ 2;
t105 = -t159 * pkin(2) + qJDD(2) * pkin(5) + t108;
t154 = sin(qJ(3));
t151 = -g(3) + qJDD(1);
t156 = cos(qJ(3));
t176 = t156 * t151;
t101 = -t154 * t105 + t176;
t102 = t156 * t105 + t154 * t151;
t112 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t154 - Ifges(5,3) * t156) * qJD(2);
t115 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t154 + Ifges(4,2) * t156) * qJD(2);
t131 = (-mrSges(5,1) * t156 - mrSges(5,3) * t154) * qJD(2);
t171 = qJD(2) * qJD(3);
t133 = t154 * qJDD(2) + t156 * t171;
t134 = t156 * qJDD(2) - t154 * t171;
t130 = (-pkin(3) * t156 - qJ(4) * t154) * qJD(2);
t158 = qJD(3) ^ 2;
t100 = -qJDD(3) * pkin(3) - t158 * qJ(4) - t176 + qJDD(4) + (qJD(2) * t130 + t105) * t154;
t172 = qJD(2) * t156;
t98 = -t158 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t130 * t172 + t102;
t165 = -mrSges(5,1) * t100 + mrSges(5,3) * t98 + Ifges(5,4) * t133 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t134;
t143 = mrSges(5,2) * t172 + qJD(3) * mrSges(5,3);
t167 = -m(5) * t100 + qJDD(3) * mrSges(5,1) + qJD(3) * t143;
t173 = qJD(2) * t154;
t141 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t173;
t168 = m(5) * t98 + qJDD(3) * mrSges(5,3) + qJD(3) * t141 + t131 * t172;
t116 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t154 - Ifges(5,5) * t156) * qJD(2);
t174 = t116 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t154 + Ifges(4,4) * t156) * qJD(2);
t179 = -(t174 * t156 + (t112 - t115) * t154) * qJD(2) + mrSges(4,1) * t101 - mrSges(4,2) * t102 + Ifges(4,5) * t133 + Ifges(4,6) * t134 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t133 * mrSges(5,2) - t131 * t173 + t167) + qJ(4) * (t134 * mrSges(5,2) + t168) + t165;
t177 = mrSges(4,3) + mrSges(5,2);
t132 = (-mrSges(4,1) * t156 + mrSges(4,2) * t154) * qJD(2);
t140 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t173;
t90 = m(4) * t102 - qJDD(3) * mrSges(4,2) - qJD(3) * t140 + t132 * t172 + t177 * t134 + t168;
t142 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t172;
t91 = m(4) * t101 + qJDD(3) * mrSges(4,1) + qJD(3) * t142 - t177 * t133 + (-t131 - t132) * t173 + t167;
t170 = -t154 * t91 + t156 * t90;
t81 = m(3) * t108 - t159 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t170;
t107 = t157 * t138 - t155 * t139;
t104 = -qJDD(2) * pkin(2) - t159 * pkin(5) - t107;
t95 = -t134 * pkin(3) - t133 * qJ(4) + (-0.2e1 * qJD(4) * t154 + (pkin(3) * t154 - qJ(4) * t156) * qJD(3)) * qJD(2) + t104;
t92 = m(5) * t95 - t134 * mrSges(5,1) - t133 * mrSges(5,3) - t141 * t173 - t143 * t172;
t161 = -m(4) * t104 + t134 * mrSges(4,1) - t133 * mrSges(4,2) - t140 * t173 + t142 * t172 - t92;
t85 = m(3) * t107 + qJDD(2) * mrSges(3,1) - t159 * mrSges(3,2) + t161;
t74 = t155 * t81 + t157 * t85;
t83 = t154 * t90 + t156 * t91;
t169 = -t155 * t85 + t157 * t81;
t166 = -mrSges(5,1) * t95 + mrSges(5,2) * t98;
t113 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t154 + Ifges(4,6) * t156) * qJD(2);
t114 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t154 - Ifges(5,6) * t156) * qJD(2);
t77 = -mrSges(4,1) * t104 + mrSges(4,3) * t102 - pkin(3) * t92 + (Ifges(4,2) + Ifges(5,3)) * t134 + (Ifges(4,4) - Ifges(5,5)) * t133 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + t174 * qJD(3) + (-t113 - t114) * t173 + t166;
t163 = mrSges(5,2) * t100 - mrSges(5,3) * t95 + Ifges(5,1) * t133 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t134 + qJD(3) * t112 + t114 * t172;
t78 = mrSges(4,2) * t104 - mrSges(4,3) * t101 + Ifges(4,1) * t133 + Ifges(4,4) * t134 + Ifges(4,5) * qJDD(3) - qJ(4) * t92 - qJD(3) * t115 + t113 * t172 + t163;
t164 = mrSges(3,1) * t107 - mrSges(3,2) * t108 + Ifges(3,3) * qJDD(2) + pkin(2) * t161 + pkin(5) * t170 + t154 * t78 + t156 * t77;
t162 = mrSges(2,1) * t138 - mrSges(2,2) * t139 + pkin(1) * t74 + t164;
t72 = m(2) * t139 + t169;
t71 = m(2) * t138 + t74;
t70 = -mrSges(3,1) * t151 + mrSges(3,3) * t108 + t159 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t83 - t179;
t69 = mrSges(3,2) * t151 - mrSges(3,3) * t107 + Ifges(3,5) * qJDD(2) - t159 * Ifges(3,6) - pkin(5) * t83 - t154 * t77 + t156 * t78;
t68 = mrSges(2,2) * t151 - mrSges(2,3) * t138 - pkin(4) * t74 - t155 * t70 + t157 * t69;
t67 = -mrSges(2,1) * t151 + mrSges(2,3) * t139 + t155 * t69 + t157 * t70 - pkin(1) * (m(3) * t151 + t83) + pkin(4) * t169;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t153 * t68 - t152 * t67 - qJ(1) * (t152 * t72 + t153 * t71), t68, t69, t78, t163; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t152 * t68 + t153 * t67 + qJ(1) * (-t152 * t71 + t153 * t72), t67, t70, t77, (-t154 * t112 - t156 * t116) * qJD(2) + t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t162, t162, t164, t179, Ifges(5,5) * t133 + Ifges(5,6) * qJDD(3) - Ifges(5,3) * t134 - qJD(3) * t116 + t114 * t173 - t166;];
m_new = t1;
