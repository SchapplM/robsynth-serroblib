% Calculate vector of cutting torques with Newton-Euler for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-05-04 19:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:03:23
% EndTime: 2019-05-04 19:03:23
% DurationCPUTime: 0.24s
% Computational Cost: add. (2829->103), mult. (3130->113), div. (0->0), fcn. (1344->4), ass. (0->46)
t118 = m(3) + m(4);
t117 = -mrSges(1,2) - mrSges(2,3);
t116 = -mrSges(4,2) - mrSges(5,2);
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t97 = -g(2) + qJDD(1);
t99 = sin(qJ(2));
t78 = t99 * g(1) + t101 * t97;
t74 = qJDD(2) * pkin(2) + t78;
t104 = qJD(2) ^ 2;
t79 = -t101 * g(1) + t99 * t97;
t75 = -t104 * pkin(2) + t79;
t98 = sin(qJ(3));
t69 = t100 * t74 - t98 * t75;
t93 = qJDD(2) + qJDD(3);
t66 = t93 * pkin(3) + t69;
t115 = m(5) * t66 + t93 * mrSges(5,1);
t94 = qJD(2) + qJD(3);
t92 = t94 ^ 2;
t58 = m(4) * t69 + t93 * mrSges(4,1) + t116 * t92 + t115;
t70 = t100 * t75 + t98 * t74;
t67 = -t92 * pkin(3) + t70;
t65 = m(5) * t67;
t59 = m(4) * t70 + t65 + t116 * t93 + (-mrSges(4,1) - mrSges(5,1)) * t92;
t52 = t100 * t58 + t98 * t59;
t114 = mrSges(5,3) * t67 + t92 * Ifges(5,5) + Ifges(5,6) * t93;
t49 = m(3) * t78 + qJDD(2) * mrSges(3,1) - t104 * mrSges(3,2) + t52;
t112 = t100 * t59 - t98 * t58;
t50 = m(3) * t79 - t104 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t112;
t113 = t101 * t50 - t99 * t49;
t111 = mrSges(5,1) * t66 - mrSges(5,2) * t67 + Ifges(5,3) * t93;
t96 = -g(3) + qJDD(4);
t110 = mrSges(5,2) * t96 - mrSges(5,3) * t66 + Ifges(5,5) * t93;
t53 = Ifges(4,6) * t93 + t92 * Ifges(4,5) + mrSges(4,1) * g(3) + mrSges(4,3) * t70 + qJ(4) * (-t92 * mrSges(5,1) - t93 * mrSges(5,2) + t65) + (-pkin(3) * m(5) - mrSges(5,1)) * t96 + t114;
t61 = -t92 * mrSges(5,2) + t115;
t54 = -mrSges(4,2) * g(3) - mrSges(4,3) * t69 + Ifges(4,5) * t93 - qJ(4) * t61 + (-Ifges(4,6) - Ifges(5,6)) * t92 + t110;
t90 = m(5) * t96;
t42 = Ifges(3,6) * qJDD(2) + t104 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t79 + t98 * t54 + t100 * t53 - pkin(2) * (-m(4) * g(3) + t90) + pkin(5) * t112;
t45 = -mrSges(3,2) * g(3) - mrSges(3,3) * t78 + Ifges(3,5) * qJDD(2) - t104 * Ifges(3,6) - pkin(5) * t52 + t100 * t54 - t98 * t53;
t109 = pkin(4) * t113 + t101 * t42 + t99 * t45 + mrSges(2,2) * g(1) + mrSges(2,1) * g(3) + pkin(1) * (t118 * g(3) - t90);
t47 = t101 * t49 + t99 * t50;
t108 = mrSges(2,2) * t97 - pkin(4) * t47 + t101 * t45 - t99 * t42;
t107 = mrSges(4,1) * t69 - mrSges(4,2) * t70 + Ifges(4,3) * t93 + pkin(3) * t61 + t111;
t106 = mrSges(3,1) * t78 - mrSges(3,2) * t79 + Ifges(3,3) * qJDD(2) + pkin(2) * t52 + t107;
t105 = mrSges(2,1) * t97 + pkin(1) * t47 + t106;
t1 = [mrSges(1,3) * g(2) + qJ(1) * t90 + (qJ(1) * (-m(2) - t118) + t117) * g(3) + t108, -mrSges(2,3) * g(3) + t108, t45, t54, -t92 * Ifges(5,6) + t110; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t109, -mrSges(2,3) * g(1) - t105, t42, t53, -mrSges(5,1) * t96 + t114; -mrSges(1,1) * g(2) + t105 - qJ(1) * t113 + (qJ(1) * m(2) - t117) * g(1), t109, t106, t107, t111;];
m_new  = t1;
