% Calculate vector of cutting torques with Newton-Euler for
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-05-04 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:24:06
% EndTime: 2019-05-04 19:24:06
% DurationCPUTime: 0.48s
% Computational Cost: add. (7149->119), mult. (8693->137), div. (0->0), fcn. (4338->6), ass. (0->54)
t128 = -mrSges(4,2) - mrSges(5,2);
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t109 = qJD(1) + qJD(2);
t107 = t109 ^ 2;
t108 = qJDD(1) + qJDD(2);
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t103 = qJD(3) + t109;
t101 = t103 ^ 2;
t102 = qJDD(3) + t108;
t114 = sin(qJ(1));
t117 = cos(qJ(1));
t92 = t114 * g(1) - t117 * g(2);
t89 = qJDD(1) * pkin(1) + t92;
t118 = qJD(1) ^ 2;
t93 = -t117 * g(1) - t114 * g(2);
t90 = -t118 * pkin(1) + t93;
t84 = -t113 * t90 + t116 * t89;
t81 = t108 * pkin(2) + t84;
t85 = t113 * t89 + t116 * t90;
t82 = -t107 * pkin(2) + t85;
t76 = -t112 * t82 + t115 * t81;
t73 = t102 * pkin(3) + t76;
t127 = m(5) * t73 + t102 * mrSges(5,1);
t65 = m(4) * t76 + t102 * mrSges(4,1) + t128 * t101 + t127;
t77 = t112 * t81 + t115 * t82;
t74 = -t101 * pkin(3) + t77;
t72 = m(5) * t74;
t66 = m(4) * t77 + t72 + t128 * t102 + (-mrSges(4,1) - mrSges(5,1)) * t101;
t59 = t112 * t66 + t115 * t65;
t56 = m(3) * t84 + t108 * mrSges(3,1) - t107 * mrSges(3,2) + t59;
t125 = -t112 * t65 + t115 * t66;
t57 = m(3) * t85 - t107 * mrSges(3,1) - t108 * mrSges(3,2) + t125;
t52 = t113 * t57 + t116 * t56;
t126 = mrSges(5,3) * t74 + t101 * Ifges(5,5) + Ifges(5,6) * t102;
t124 = -t113 * t56 + t116 * t57;
t123 = mrSges(5,1) * t73 - mrSges(5,2) * t74 + Ifges(5,3) * t102;
t111 = -g(3) + qJDD(4);
t122 = mrSges(5,2) * t111 - mrSges(5,3) * t73 + Ifges(5,5) * t102;
t68 = -t101 * mrSges(5,2) + t127;
t121 = mrSges(4,1) * t76 - mrSges(4,2) * t77 + Ifges(4,3) * t102 + pkin(3) * t68 + t123;
t120 = mrSges(3,1) * t84 - mrSges(3,2) * t85 + Ifges(3,3) * t108 + pkin(2) * t59 + t121;
t119 = mrSges(2,1) * t92 - mrSges(2,2) * t93 + Ifges(2,3) * qJDD(1) + pkin(1) * t52 + t120;
t105 = m(5) * t111;
t61 = -mrSges(4,2) * g(3) - mrSges(4,3) * t76 + Ifges(4,5) * t102 - qJ(4) * t68 + (-Ifges(4,6) - Ifges(5,6)) * t101 + t122;
t60 = Ifges(4,6) * t102 + t101 * Ifges(4,5) + mrSges(4,1) * g(3) + mrSges(4,3) * t77 + qJ(4) * (-t101 * mrSges(5,1) - t102 * mrSges(5,2) + t72) + (-pkin(3) * m(5) - mrSges(5,1)) * t111 + t126;
t50 = m(2) * t93 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t124;
t49 = m(2) * t92 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) + t52;
t48 = -mrSges(3,2) * g(3) - mrSges(3,3) * t84 + Ifges(3,5) * t108 - t107 * Ifges(3,6) - pkin(6) * t59 - t112 * t60 + t115 * t61;
t47 = Ifges(3,6) * t108 + t107 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t85 + t112 * t61 + t115 * t60 - pkin(2) * (-m(4) * g(3) + t105) + pkin(6) * t125;
t46 = -mrSges(2,2) * g(3) - mrSges(2,3) * t92 + Ifges(2,5) * qJDD(1) - t118 * Ifges(2,6) - pkin(5) * t52 - t113 * t47 + t116 * t48;
t45 = Ifges(2,6) * qJDD(1) + t118 * Ifges(2,5) + mrSges(2,3) * t93 + t113 * t48 + t116 * t47 - pkin(1) * t105 + pkin(5) * t124 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t117 * t46 - t114 * t45 - pkin(4) * (t114 * t50 + t117 * t49), t46, t48, t61, -t101 * Ifges(5,6) + t122; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t114 * t46 + t117 * t45 + pkin(4) * (-t114 * t49 + t117 * t50), t45, t47, t60, -mrSges(5,1) * t111 + t126; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t119, t119, t120, t121, t123;];
m_new  = t1;
