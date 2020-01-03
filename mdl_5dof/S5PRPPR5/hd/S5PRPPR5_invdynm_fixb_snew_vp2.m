% Calculate vector of cutting torques with Newton-Euler for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:05
% EndTime: 2019-12-31 17:38:07
% DurationCPUTime: 1.13s
% Computational Cost: add. (13752->172), mult. (21381->210), div. (0->0), fcn. (9622->8), ass. (0->76)
t164 = qJD(2) ^ 2;
t158 = sin(pkin(7));
t182 = cos(pkin(7));
t142 = -t182 * g(1) - t158 * g(2);
t154 = -g(3) + qJDD(1);
t161 = sin(qJ(2));
t163 = cos(qJ(2));
t127 = t163 * t142 + t161 * t154;
t176 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t127;
t184 = -pkin(2) - pkin(3);
t119 = t184 * t164 + t176;
t126 = -t161 * t142 + t163 * t154;
t171 = -t164 * qJ(3) + qJDD(3) - t126;
t122 = t184 * qJDD(2) + t171;
t157 = sin(pkin(8));
t159 = cos(pkin(8));
t116 = t159 * t119 + t157 * t122;
t112 = -t164 * pkin(4) - qJDD(2) * pkin(6) + t116;
t141 = t158 * g(1) - t182 * g(2);
t140 = qJDD(4) + t141;
t160 = sin(qJ(5));
t162 = cos(qJ(5));
t109 = -t160 * t112 + t162 * t140;
t110 = t162 * t112 + t160 * t140;
t129 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t160 - Ifges(6,2) * t162) * qJD(2);
t130 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t160 - Ifges(6,4) * t162) * qJD(2);
t179 = qJD(2) * qJD(5);
t137 = -t160 * qJDD(2) - t162 * t179;
t138 = -t162 * qJDD(2) + t160 * t179;
t186 = mrSges(6,1) * t109 - mrSges(6,2) * t110 + Ifges(6,5) * t137 + Ifges(6,6) * t138 + Ifges(6,3) * qJDD(5) - (t129 * t160 - t130 * t162) * qJD(2);
t185 = m(3) + m(4);
t183 = mrSges(3,1) + mrSges(4,1);
t136 = (mrSges(6,1) * t162 - mrSges(6,2) * t160) * qJD(2);
t180 = qJD(2) * t162;
t144 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t180;
t181 = qJD(2) * t160;
t105 = m(6) * t109 + qJDD(5) * mrSges(6,1) - t137 * mrSges(6,3) + qJD(5) * t144 + t136 * t181;
t143 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t181;
t106 = m(6) * t110 - qJDD(5) * mrSges(6,2) + t138 * mrSges(6,3) - qJD(5) * t143 - t136 * t180;
t95 = t162 * t105 + t160 * t106;
t93 = -m(5) * t140 - t95;
t123 = -t164 * pkin(2) + t176;
t115 = -t157 * t119 + t159 * t122;
t111 = qJDD(2) * pkin(4) - t164 * pkin(6) - t115;
t169 = -m(6) * t111 + t138 * mrSges(6,1) - t137 * mrSges(6,2) + t143 * t181 - t144 * t180;
t101 = m(5) * t115 - qJDD(2) * mrSges(5,1) - t164 * mrSges(5,2) + t169;
t177 = -t160 * t105 + t162 * t106;
t89 = m(5) * t116 - t164 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + t177;
t87 = -t157 * t101 + t159 * t89;
t174 = m(4) * t123 + qJDD(2) * mrSges(4,3) + t87;
t81 = m(3) * t127 - qJDD(2) * mrSges(3,2) - t183 * t164 + t174;
t125 = -qJDD(2) * pkin(2) + t171;
t86 = t159 * t101 + t157 * t89;
t85 = -m(4) * t125 + qJDD(2) * mrSges(4,1) + t164 * mrSges(4,3) - t86;
t82 = m(3) * t126 + qJDD(2) * mrSges(3,1) - t164 * mrSges(3,2) + t85;
t178 = -t161 * t82 + t163 * t81;
t128 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t160 - Ifges(6,6) * t162) * qJD(2);
t100 = mrSges(6,2) * t111 - mrSges(6,3) * t109 + Ifges(6,1) * t137 + Ifges(6,4) * t138 + Ifges(6,5) * qJDD(5) - qJD(5) * t129 - t128 * t180;
t99 = -mrSges(6,1) * t111 + mrSges(6,3) * t110 + Ifges(6,4) * t137 + Ifges(6,2) * t138 + Ifges(6,6) * qJDD(5) + qJD(5) * t130 + t128 * t181;
t173 = mrSges(5,1) * t115 - mrSges(5,2) * t116 - Ifges(5,3) * qJDD(2) + pkin(4) * t169 + pkin(6) * t177 + t160 * t100 + t162 * t99;
t79 = mrSges(5,2) * t140 - mrSges(5,3) * t115 - Ifges(5,5) * qJDD(2) - t164 * Ifges(5,6) - pkin(6) * t95 + t162 * t100 - t160 * t99;
t83 = -mrSges(5,1) * t140 + mrSges(5,3) * t116 + t164 * Ifges(5,5) - Ifges(5,6) * qJDD(2) - pkin(4) * t95 - t186;
t168 = mrSges(4,2) * t123 - pkin(3) * t93 - qJ(4) * t87 - t157 * t79 - t159 * t83;
t92 = -m(4) * t141 + t93;
t72 = mrSges(3,3) * t127 - pkin(2) * t92 + (Ifges(4,4) + Ifges(3,5)) * t164 + t183 * t141 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t168;
t170 = mrSges(4,2) * t125 + Ifges(4,4) * qJDD(2) + t164 * Ifges(4,6) - qJ(4) * t86 - t157 * t83 + t159 * t79;
t74 = -mrSges(3,3) * t126 + Ifges(3,5) * qJDD(2) - t164 * Ifges(3,6) - qJ(3) * t92 + (-mrSges(3,2) + mrSges(4,3)) * t141 + t170;
t172 = -mrSges(2,2) * t142 + pkin(5) * t178 + t161 * t74 + t163 * t72 + mrSges(2,1) * t141 + pkin(1) * (t185 * t141 - t93);
t166 = -mrSges(4,1) * t125 + mrSges(4,3) * t123 + Ifges(4,2) * qJDD(2) - pkin(3) * t86 - t173;
t165 = mrSges(3,1) * t126 - mrSges(3,2) * t127 + Ifges(3,3) * qJDD(2) + pkin(2) * t85 + qJ(3) * (-t164 * mrSges(4,1) + t174) + t166;
t90 = (m(2) + t185) * t141 - t93;
t77 = t161 * t81 + t163 * t82;
t75 = m(2) * t142 + t178;
t70 = -mrSges(2,1) * t154 + mrSges(2,3) * t142 - pkin(1) * t77 - t165;
t69 = mrSges(2,2) * t154 - mrSges(2,3) * t141 - pkin(5) * t77 - t161 * t72 + t163 * t74;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t182 * t69 - t158 * t70 - qJ(1) * (t158 * t75 + t182 * t90), t69, t74, mrSges(4,3) * t141 + t170, t79, t100; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t158 * t69 + t182 * t70 + qJ(1) * (-t158 * t90 + t182 * t75), t70, t72, t166, t83, t99; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t172, t172, t165, -mrSges(4,1) * t141 - t164 * Ifges(4,4) + Ifges(4,6) * qJDD(2) - t168, t173, t186;];
m_new = t1;
