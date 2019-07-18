% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:23
% DurationCPUTime: 0.23s
% Computational Cost: add. (2673->99), mult. (3386->113), div. (0->0), fcn. (1718->6), ass. (0->45)
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t98 = sin(qJ(2));
t80 = t101 * g(1) + t98 * g(3);
t77 = qJDD(2) * pkin(1) + t80;
t102 = qJD(2) ^ 2;
t81 = t98 * g(1) - t101 * g(3);
t78 = -t102 * pkin(1) + t81;
t97 = sin(qJ(3));
t72 = t100 * t77 - t97 * t78;
t92 = qJDD(2) + qJDD(3);
t69 = t92 * pkin(2) + t72;
t73 = t100 * t78 + t97 * t77;
t93 = qJD(2) + qJD(3);
t91 = t93 ^ 2;
t70 = -t91 * pkin(2) + t73;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t67 = t99 * t69 - t96 * t70;
t86 = qJD(4) + t93;
t84 = t86 ^ 2;
t85 = qJDD(4) + t92;
t64 = m(5) * t67 + t85 * mrSges(5,1) - t84 * mrSges(5,2);
t68 = t96 * t69 + t99 * t70;
t65 = m(5) * t68 - t84 * mrSges(5,1) - t85 * mrSges(5,2);
t109 = t99 * t64 + t96 * t65;
t56 = m(4) * t72 + t92 * mrSges(4,1) - t91 * mrSges(4,2) + t109;
t57 = m(4) * t73 - t91 * mrSges(4,1) - t92 * mrSges(4,2) - t96 * t64 + t99 * t65;
t110 = t100 * t56 + t97 * t57;
t108 = qJ(1) * m(2) - mrSges(1,2) + mrSges(2,3);
t107 = -mrSges(5,1) * t67 + mrSges(5,2) * t68 - Ifges(5,3) * t85;
t95 = g(2) + qJDD(1);
t60 = -mrSges(5,1) * t95 + mrSges(5,3) * t68 + t84 * Ifges(5,5) + Ifges(5,6) * t85;
t61 = mrSges(5,2) * t95 - mrSges(5,3) * t67 + Ifges(5,5) * t85 - t84 * Ifges(5,6);
t52 = mrSges(4,3) * t73 + t91 * Ifges(4,5) + Ifges(4,6) * t92 + t99 * t60 + t96 * t61 + (-m(5) * pkin(2) - mrSges(4,1)) * t95;
t53 = mrSges(4,2) * t95 - mrSges(4,3) * t72 + Ifges(4,5) * t92 - t91 * Ifges(4,6) - t96 * t60 + t99 * t61;
t45 = mrSges(3,3) * t81 + t102 * Ifges(3,5) + Ifges(3,6) * qJDD(2) + t100 * t52 + t97 * t53 + (-mrSges(3,1) - pkin(1) * (m(4) + m(5))) * t95;
t47 = mrSges(3,2) * t95 - mrSges(3,3) * t80 + Ifges(3,5) * qJDD(2) - t102 * Ifges(3,6) + t100 * t53 - t97 * t52;
t106 = mrSges(2,2) * t95 + t101 * t47 - t98 * t45;
t105 = mrSges(2,1) * t95 - t101 * t45 - t98 * t47;
t104 = -mrSges(4,1) * t72 + mrSges(4,2) * t73 - Ifges(4,3) * t92 - pkin(2) * t109 + t107;
t103 = -mrSges(3,1) * t80 + mrSges(3,2) * t81 - Ifges(3,3) * qJDD(2) - pkin(1) * t110 + t104;
t49 = m(3) * t81 - t102 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t100 * t57 - t97 * t56;
t48 = m(3) * t80 + qJDD(2) * mrSges(3,1) - t102 * mrSges(3,2) + t110;
t1 = [mrSges(1,3) * g(2) - qJ(1) * (t101 * t49 - t98 * t48) + t108 * g(3) + t105, -mrSges(2,3) * g(1) + t106, t47, t53, t61; t103 + (mrSges(1,1) - mrSges(2,2)) * g(3) + (-mrSges(2,1) - mrSges(1,3)) * g(1), -mrSges(2,3) * g(3) - t105, t45, t52, t60; -mrSges(1,1) * g(2) + qJ(1) * (-t101 * t48 - t98 * t49) - t108 * g(1) + t106, mrSges(2,1) * g(1) + mrSges(2,2) * g(3) - t103, -t103, -t104, -t107;];
m_new  = t1;
