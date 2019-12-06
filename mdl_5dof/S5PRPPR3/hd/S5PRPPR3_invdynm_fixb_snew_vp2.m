% Calculate vector of cutting torques with Newton-Euler for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:22
% EndTime: 2019-12-05 15:26:25
% DurationCPUTime: 1.11s
% Computational Cost: add. (12658->169), mult. (19653->207), div. (0->0), fcn. (10502->8), ass. (0->77)
t184 = -pkin(3) - pkin(6);
t183 = mrSges(4,1) - mrSges(5,2);
t182 = -Ifges(5,4) + Ifges(4,5);
t181 = Ifges(5,5) - Ifges(4,6);
t155 = sin(pkin(8));
t157 = cos(pkin(8));
t156 = sin(pkin(7));
t158 = cos(pkin(7));
t142 = -t158 * g(1) - t156 * g(2);
t152 = -g(3) + qJDD(1);
t160 = sin(qJ(2));
t162 = cos(qJ(2));
t122 = -t160 * t142 + t162 * t152;
t120 = qJDD(2) * pkin(2) + t122;
t123 = t162 * t142 + t160 * t152;
t163 = qJD(2) ^ 2;
t121 = -t163 * pkin(2) + t123;
t115 = t157 * t120 - t155 * t121;
t173 = -t163 * qJ(4) + qJDD(4) - t115;
t113 = -qJDD(2) * pkin(3) + t173;
t110 = t184 * qJDD(2) + t173;
t141 = t156 * g(1) - t158 * g(2);
t140 = qJDD(3) - t141;
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t106 = t161 * t110 - t159 * t140;
t136 = (mrSges(6,1) * t159 + mrSges(6,2) * t161) * qJD(2);
t178 = qJD(2) * qJD(5);
t138 = t161 * qJDD(2) - t159 * t178;
t180 = qJD(2) * t159;
t143 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t180;
t179 = qJD(2) * t161;
t102 = m(6) * t106 + qJDD(5) * mrSges(6,1) - t138 * mrSges(6,3) + qJD(5) * t143 - t136 * t179;
t107 = t159 * t110 + t161 * t140;
t137 = -t159 * qJDD(2) - t161 * t178;
t144 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t179;
t103 = m(6) * t107 - qJDD(5) * mrSges(6,2) + t137 * mrSges(6,3) - qJD(5) * t144 - t136 * t180;
t90 = t161 * t102 + t159 * t103;
t171 = -m(5) * t113 + t163 * mrSges(5,3) - t90;
t84 = m(4) * t115 - t163 * mrSges(4,2) + t183 * qJDD(2) + t171;
t116 = t155 * t120 + t157 * t121;
t172 = qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) + t116;
t109 = t184 * t163 + t172;
t104 = -m(6) * t109 + t137 * mrSges(6,1) - t138 * mrSges(6,2) - t143 * t180 - t144 * t179;
t111 = t163 * pkin(3) - t172;
t168 = -m(5) * t111 + t163 * mrSges(5,2) + qJDD(2) * mrSges(5,3) - t104;
t94 = m(4) * t116 - t163 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t168;
t82 = t155 * t94 + t157 * t84;
t91 = -t159 * t102 + t161 * t103;
t89 = m(5) * t140 + t91;
t177 = -t155 * t84 + t157 * t94;
t80 = m(3) * t122 + qJDD(2) * mrSges(3,1) - t163 * mrSges(3,2) + t82;
t81 = m(3) * t123 - t163 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t177;
t176 = -t160 * t80 + t162 * t81;
t175 = m(4) * t140 + t89;
t126 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t161 - Ifges(6,6) * t159) * qJD(2);
t128 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t161 - Ifges(6,4) * t159) * qJD(2);
t96 = -mrSges(6,1) * t109 + mrSges(6,3) * t107 + Ifges(6,4) * t138 + Ifges(6,2) * t137 + Ifges(6,6) * qJDD(5) + qJD(5) * t128 - t126 * t179;
t127 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t161 - Ifges(6,2) * t159) * qJD(2);
t97 = mrSges(6,2) * t109 - mrSges(6,3) * t106 + Ifges(6,1) * t138 + Ifges(6,4) * t137 + Ifges(6,5) * qJDD(5) - qJD(5) * t127 - t126 * t180;
t167 = -mrSges(5,1) * t111 - pkin(4) * t104 - pkin(6) * t91 - t159 * t97 - t161 * t96;
t77 = mrSges(4,3) * t116 - pkin(3) * t89 - t181 * qJDD(2) - t183 * t140 + t182 * t163 + t167;
t170 = mrSges(6,1) * t106 - mrSges(6,2) * t107 + Ifges(6,5) * t138 + Ifges(6,6) * t137 + Ifges(6,3) * qJDD(5) + t127 * t179 + t128 * t180;
t166 = mrSges(5,1) * t113 + pkin(4) * t90 + t170;
t78 = -mrSges(4,3) * t115 - qJ(4) * t89 + t181 * t163 + (mrSges(4,2) - mrSges(5,3)) * t140 + t182 * qJDD(2) + t166;
t71 = mrSges(3,1) * t141 + mrSges(3,3) * t123 + t163 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t175 + qJ(3) * t177 + t155 * t78 + t157 * t77;
t73 = -mrSges(3,2) * t141 - mrSges(3,3) * t122 + Ifges(3,5) * qJDD(2) - t163 * Ifges(3,6) - qJ(3) * t82 - t155 * t77 + t157 * t78;
t174 = -mrSges(2,2) * t142 + pkin(5) * t176 + t160 * t73 + t162 * t71 + mrSges(2,1) * t141 + pkin(1) * (m(3) * t141 - t175);
t169 = mrSges(5,2) * t113 - mrSges(5,3) * t111 + Ifges(5,1) * qJDD(2) - pkin(6) * t90 - t159 * t96 + t161 * t97;
t165 = -mrSges(4,2) * t116 + mrSges(4,1) * t115 + Ifges(4,3) * qJDD(2) + t169 + pkin(3) * (-qJDD(2) * mrSges(5,2) + t171) + qJ(4) * t168;
t164 = mrSges(3,1) * t122 - mrSges(3,2) * t123 + Ifges(3,3) * qJDD(2) + pkin(2) * t82 + t165;
t87 = (m(2) + m(3)) * t141 - t175;
t76 = t160 * t81 + t162 * t80;
t74 = m(2) * t142 + t176;
t69 = -mrSges(2,1) * t152 + mrSges(2,3) * t142 - pkin(1) * t76 - t164;
t68 = mrSges(2,2) * t152 - mrSges(2,3) * t141 - pkin(5) * t76 - t160 * t71 + t162 * t73;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t158 * t68 - t156 * t69 - qJ(1) * (t156 * t74 + t158 * t87), t68, t73, t78, t169, t97; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t156 * t68 + t158 * t69 + qJ(1) * (-t156 * t87 + t158 * t74), t69, t71, t77, mrSges(5,3) * t140 + Ifges(5,4) * qJDD(2) - t163 * Ifges(5,5) - t166, t96; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t174, t174, t164, t165, -mrSges(5,2) * t140 + t163 * Ifges(5,4) + Ifges(5,5) * qJDD(2) - t167, t170;];
m_new = t1;
