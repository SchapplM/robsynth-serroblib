% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:40
% EndTime: 2019-12-31 17:35:42
% DurationCPUTime: 1.13s
% Computational Cost: add. (16859->157), mult. (21445->196), div. (0->0), fcn. (12906->8), ass. (0->72)
t140 = sin(pkin(8));
t141 = cos(pkin(8));
t133 = t140 * g(1) - t141 * g(2);
t131 = qJDD(2) - t133;
t134 = -t141 * g(1) - t140 * g(2);
t144 = sin(qJ(3));
t147 = cos(qJ(3));
t113 = t147 * t131 - t144 * t134;
t110 = qJDD(3) * pkin(3) + t113;
t114 = t144 * t131 + t147 * t134;
t148 = qJD(3) ^ 2;
t111 = -t148 * pkin(3) + t114;
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t107 = t143 * t110 + t146 * t111;
t138 = qJD(3) + qJD(4);
t136 = t138 ^ 2;
t137 = qJDD(3) + qJDD(4);
t103 = -t136 * pkin(4) + t137 * pkin(7) + t107;
t139 = g(3) - qJDD(1);
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t100 = -t142 * t103 + t145 * t139;
t101 = t145 * t103 + t142 * t139;
t116 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t142 + Ifges(6,2) * t145) * t138;
t117 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t142 + Ifges(6,4) * t145) * t138;
t162 = qJD(5) * t138;
t121 = t142 * t137 + t145 * t162;
t122 = t145 * t137 - t142 * t162;
t166 = mrSges(6,1) * t100 - mrSges(6,2) * t101 + Ifges(6,5) * t121 + Ifges(6,6) * t122 + Ifges(6,3) * qJDD(5) + (t116 * t142 - t117 * t145) * t138;
t165 = -m(4) - m(5);
t120 = (-mrSges(6,1) * t145 + mrSges(6,2) * t142) * t138;
t163 = t138 * t145;
t128 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t163;
t164 = t138 * t142;
t98 = m(6) * t100 + qJDD(5) * mrSges(6,1) - t121 * mrSges(6,3) + qJD(5) * t128 - t120 * t164;
t127 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t164;
t99 = m(6) * t101 - qJDD(5) * mrSges(6,2) + t122 * mrSges(6,3) - qJD(5) * t127 + t120 * t163;
t161 = -t142 * t98 + t145 * t99;
t81 = m(5) * t107 - t136 * mrSges(5,1) - t137 * mrSges(5,2) + t161;
t106 = t146 * t110 - t143 * t111;
t102 = -t137 * pkin(4) - t136 * pkin(7) - t106;
t153 = -m(6) * t102 + t122 * mrSges(6,1) - t121 * mrSges(6,2) - t127 * t164 + t128 * t163;
t92 = m(5) * t106 + t137 * mrSges(5,1) - t136 * mrSges(5,2) + t153;
t78 = t143 * t81 + t146 * t92;
t85 = t142 * t99 + t145 * t98;
t160 = -t143 * t92 + t146 * t81;
t76 = m(4) * t113 + qJDD(3) * mrSges(4,1) - t148 * mrSges(4,2) + t78;
t77 = m(4) * t114 - t148 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t160;
t72 = -t144 * t76 + t147 * t77;
t159 = m(3) * t134 + t72;
t71 = t144 * t77 + t147 * t76;
t115 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t142 + Ifges(6,6) * t145) * t138;
t89 = -mrSges(6,1) * t102 + mrSges(6,3) * t101 + Ifges(6,4) * t121 + Ifges(6,2) * t122 + Ifges(6,6) * qJDD(5) + qJD(5) * t117 - t115 * t164;
t90 = mrSges(6,2) * t102 - mrSges(6,3) * t100 + Ifges(6,1) * t121 + Ifges(6,4) * t122 + Ifges(6,5) * qJDD(5) - qJD(5) * t116 + t115 * t163;
t157 = mrSges(5,1) * t106 - mrSges(5,2) * t107 + Ifges(5,3) * t137 + pkin(4) * t153 + pkin(7) * t161 + t142 * t90 + t145 * t89;
t73 = mrSges(5,2) * t139 - mrSges(5,3) * t106 + Ifges(5,5) * t137 - t136 * Ifges(5,6) - pkin(7) * t85 - t142 * t89 + t145 * t90;
t75 = -mrSges(5,1) * t139 + mrSges(5,3) * t107 + t136 * Ifges(5,5) + Ifges(5,6) * t137 - pkin(4) * t85 - t166;
t64 = Ifges(4,6) * qJDD(3) + t148 * Ifges(4,5) - mrSges(4,1) * t139 + mrSges(4,3) * t114 + t143 * t73 + t146 * t75 - pkin(3) * (m(5) * t139 + t85) + pkin(6) * t160;
t66 = mrSges(4,2) * t139 - mrSges(4,3) * t113 + Ifges(4,5) * qJDD(3) - t148 * Ifges(4,6) - pkin(6) * t78 - t143 * t75 + t146 * t73;
t156 = mrSges(3,2) * t131 - pkin(5) * t71 - t144 * t64 + t147 * t66;
t155 = -m(3) * t131 - t71;
t154 = -pkin(2) * (t165 * t139 - t85) - pkin(5) * t72 - t144 * t66 - t147 * t64;
t151 = mrSges(4,1) * t113 - mrSges(4,2) * t114 + Ifges(4,3) * qJDD(3) + pkin(3) * t78 + t157;
t150 = -mrSges(3,1) * t131 + mrSges(3,3) * t134 - pkin(2) * t71 - t151;
t149 = mrSges(2,1) * t133 - mrSges(2,2) * t134 + pkin(1) * t155 + qJ(2) * t159 + t150;
t82 = (-m(3) + t165) * t139 - t85;
t68 = m(2) * t134 + t159;
t67 = m(2) * t133 + t155;
t63 = -mrSges(2,3) * t133 - qJ(2) * t82 + (-mrSges(2,2) + mrSges(3,3)) * t139 + t156;
t62 = -pkin(1) * t82 + (mrSges(2,1) + mrSges(3,1)) * t139 + (mrSges(3,2) + mrSges(2,3)) * t134 + t154;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t141 * t63 - t140 * t62 - qJ(1) * (t140 * t68 + t141 * t67), t63, mrSges(3,3) * t139 + t156, t66, t73, t90; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t140 * t63 + t141 * t62 + qJ(1) * (-t140 * t67 + t141 * t68), t62, t150, t64, t75, t89; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t149, t149, -mrSges(3,1) * t139 - mrSges(3,2) * t134 - t154, t151, t157, t166;];
m_new = t1;
