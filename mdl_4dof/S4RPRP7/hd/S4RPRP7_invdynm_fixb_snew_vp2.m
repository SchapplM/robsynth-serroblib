% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:04
% EndTime: 2019-12-31 16:47:05
% DurationCPUTime: 0.62s
% Computational Cost: add. (3765->188), mult. (6953->228), div. (0->0), fcn. (2644->4), ass. (0->69)
t153 = sin(qJ(1));
t155 = cos(qJ(1));
t138 = -t155 * g(1) - t153 * g(2);
t184 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t138;
t183 = -pkin(1) - pkin(5);
t152 = sin(qJ(3));
t182 = t152 * g(3);
t181 = mrSges(2,1) - mrSges(3,2);
t180 = -mrSges(4,3) - mrSges(5,2);
t179 = Ifges(2,5) - Ifges(3,4);
t178 = (-Ifges(2,6) + Ifges(3,5));
t177 = qJD(1) * t152;
t154 = cos(qJ(3));
t176 = qJD(1) * t154;
t175 = qJD(1) * qJD(3);
t135 = -(qJD(3) * mrSges(5,1)) + mrSges(5,2) * t176;
t125 = (pkin(3) * t152 - qJ(4) * t154) * qJD(1);
t156 = qJD(3) ^ 2;
t137 = t153 * g(1) - t155 * g(2);
t157 = qJD(1) ^ 2;
t167 = -t157 * qJ(2) + qJDD(2) - t137;
t100 = t183 * qJDD(1) + t167;
t96 = -t154 * g(3) + t152 * t100;
t91 = -t156 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) - t125 * t177 + t96;
t173 = m(5) * t91 + qJDD(3) * mrSges(5,3) + qJD(3) * t135;
t128 = t152 * qJDD(1) + t154 * t175;
t134 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t176;
t126 = (mrSges(5,1) * t152 - mrSges(5,3) * t154) * qJD(1);
t171 = qJD(1) * (-t126 - (mrSges(4,1) * t152 + mrSges(4,2) * t154) * qJD(1));
t80 = m(4) * t96 - qJDD(3) * mrSges(4,2) - qJD(3) * t134 + t180 * t128 + t152 * t171 + t173;
t129 = t154 * qJDD(1) - t152 * t175;
t133 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t177;
t136 = -mrSges(5,2) * t177 + qJD(3) * mrSges(5,3);
t93 = -qJDD(3) * pkin(3) - t182 - t156 * qJ(4) + qJDD(4) + (qJD(1) * t125 - t100) * t154;
t170 = -m(5) * t93 + qJDD(3) * mrSges(5,1) + qJD(3) * t136;
t95 = t154 * t100 + t182;
t81 = m(4) * t95 + qJDD(3) * mrSges(4,1) + qJD(3) * t133 + t180 * t129 + t154 * t171 + t170;
t75 = -t152 * t81 + t154 * t80;
t110 = (Ifges(5,2) * qJD(3)) + (Ifges(5,4) * t154 + Ifges(5,6) * t152) * qJD(1);
t172 = qJD(1) * (-(Ifges(4,3) * qJD(3)) - (Ifges(4,5) * t154 - Ifges(4,6) * t152) * qJD(1) - t110);
t99 = t183 * t157 - t184;
t88 = t128 * pkin(3) - t129 * qJ(4) + (-0.2e1 * qJD(4) * t154 + (pkin(3) * t154 + qJ(4) * t152) * qJD(3)) * qJD(1) + t99;
t169 = -mrSges(5,1) * t88 + mrSges(5,2) * t91;
t74 = t152 * t80 + t154 * t81;
t108 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t154 + Ifges(5,3) * t152) * qJD(1);
t166 = mrSges(5,2) * t93 - mrSges(5,3) * t88 + Ifges(5,1) * t129 + Ifges(5,4) * qJDD(3) + Ifges(5,5) * t128 + qJD(3) * t108;
t107 = -qJDD(1) * pkin(1) + t167;
t165 = -m(3) * t107 + (t157 * mrSges(3,3)) - t74;
t82 = m(5) * t88 + t128 * mrSges(5,1) - t129 * mrSges(5,3) - t135 * t176 + t136 * t177;
t104 = t157 * pkin(1) + t184;
t112 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t154 + Ifges(5,5) * t152) * qJD(1);
t113 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t154 - Ifges(4,4) * t152) * qJD(1);
t69 = -mrSges(4,1) * t99 + mrSges(4,3) * t96 - pkin(3) * t82 + (Ifges(4,4) - Ifges(5,5)) * t129 + (-Ifges(4,2) - Ifges(5,3)) * t128 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + (t112 + t113) * qJD(3) + t154 * t172 + t169;
t111 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t154 - Ifges(4,2) * t152) * qJD(1);
t70 = mrSges(4,2) * t99 - mrSges(4,3) * t95 + Ifges(4,1) * t129 - Ifges(4,4) * t128 + Ifges(4,5) * qJDD(3) - qJ(4) * t82 - qJD(3) * t111 + t152 * t172 + t166;
t164 = mrSges(3,2) * t107 - mrSges(3,3) * t104 + Ifges(3,1) * qJDD(1) - pkin(5) * t74 - t152 * t69 + t154 * t70;
t163 = -mrSges(5,1) * t93 + mrSges(5,3) * t91 + Ifges(5,4) * t129 + Ifges(5,2) * qJDD(3) + Ifges(5,6) * t128 - t108 * t176 + t112 * t177;
t78 = -m(4) * t99 - t128 * mrSges(4,1) - t129 * mrSges(4,2) - t133 * t177 - t134 * t176 - t82;
t162 = -mrSges(3,1) * t104 - pkin(2) * t78 - pkin(5) * t75 - t152 * t70 - t154 * t69;
t160 = -m(3) * t104 + t157 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t78;
t161 = -mrSges(2,2) * t138 + mrSges(2,1) * t137 + Ifges(2,3) * qJDD(1) + t164 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t165) + qJ(2) * t160;
t159 = -mrSges(4,2) * t96 + t113 * t177 + t111 * t176 - Ifges(4,6) * t128 + Ifges(4,5) * t129 + Ifges(4,3) * qJDD(3) + t163 + qJ(4) * (-t128 * mrSges(5,2) - t126 * t177 + t173) + pkin(3) * (-t129 * mrSges(5,2) - t126 * t176 + t170) + mrSges(4,1) * t95;
t158 = mrSges(3,1) * t107 + pkin(2) * t74 + t159;
t76 = m(2) * t138 - t157 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t160;
t73 = -m(3) * g(3) + t75;
t71 = m(2) * t137 - t157 * mrSges(2,2) + t181 * qJDD(1) + t165;
t67 = t158 - mrSges(2,3) * t137 - qJ(2) * t73 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t179 * qJDD(1) + (t178 * t157);
t66 = mrSges(2,3) * t138 - pkin(1) * t73 + t181 * g(3) - t178 * qJDD(1) + t179 * t157 + t162;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t155 * t67 - t153 * t66 - pkin(4) * (t153 * t76 + t155 * t71), t67, t164, t70, -t110 * t177 + t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t153 * t67 + t155 * t66 + pkin(4) * (-t153 * t71 + t155 * t76), t66, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t157 * Ifges(3,5)) - t158, t69, t163; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t161, t161, mrSges(3,2) * g(3) + t157 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t162, t159, Ifges(5,5) * t129 + Ifges(5,6) * qJDD(3) + Ifges(5,3) * t128 - qJD(3) * t112 + t110 * t176 - t169;];
m_new = t1;
