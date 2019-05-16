% Calculate vector of cutting torques with Newton-Euler for
% S4PPRP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2019-05-04 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:46:05
% EndTime: 2019-05-04 18:46:05
% DurationCPUTime: 0.21s
% Computational Cost: add. (963->89), mult. (1035->89), div. (0->0), fcn. (304->2), ass. (0->36)
t75 = -g(2) + qJDD(1);
t76 = -g(1) + qJDD(2);
t78 = sin(qJ(3));
t79 = cos(qJ(3));
t58 = -t78 * t75 + t79 * t76;
t82 = qJD(3) ^ 2;
t54 = qJDD(3) * pkin(3) + t58;
t95 = m(5) * t54 + qJDD(3) * mrSges(5,1);
t96 = -mrSges(4,2) - mrSges(5,2);
t45 = m(4) * t58 + qJDD(3) * mrSges(4,1) + t82 * t96 + t95;
t59 = t79 * t75 + t78 * t76;
t55 = -t82 * pkin(3) + t59;
t53 = m(5) * t55;
t46 = m(4) * t59 + t53 + (-mrSges(4,1) - mrSges(5,1)) * t82 + t96 * qJDD(3);
t92 = -t78 * t45 + t79 * t46;
t33 = m(3) * t75 + t92;
t74 = g(3) + qJDD(4);
t93 = mrSges(5,3) * t55 + t82 * Ifges(5,5) + Ifges(5,6) * qJDD(3);
t38 = Ifges(4,6) * qJDD(3) + t82 * Ifges(4,5) - mrSges(4,1) * g(3) + mrSges(4,3) * t59 + qJ(4) * (-t82 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t53) + (-m(5) * pkin(3) - mrSges(5,1)) * t74 + t93;
t48 = -t82 * mrSges(5,2) + t95;
t89 = mrSges(5,2) * t74 - mrSges(5,3) * t54 + Ifges(5,5) * qJDD(3);
t41 = mrSges(4,2) * g(3) - mrSges(4,3) * t58 + Ifges(4,5) * qJDD(3) - qJ(4) * t48 + (-Ifges(4,6) - Ifges(5,6)) * t82 + t89;
t99 = m(5) * t74;
t91 = pkin(4) * t92 + t79 * t38 + t78 * t41 + pkin(2) * (-m(4) * g(3) - t99);
t101 = -pkin(1) * t33 - (-mrSges(3,2) + mrSges(2,1)) * t75 - t91;
t100 = m(3) + m(4);
t98 = -mrSges(1,2) - mrSges(2,3);
t35 = t79 * t45 + t78 * t46;
t94 = mrSges(5,1) * t54 - mrSges(5,2) * t55 + Ifges(5,3) * qJDD(3);
t90 = m(3) * t76 + t35;
t87 = mrSges(4,1) * t58 - mrSges(4,2) * t59 + Ifges(4,3) * qJDD(3) + pkin(3) * t48 + t94;
t86 = mrSges(3,2) * t76 + mrSges(3,3) * g(3) - pkin(4) * t35 - t78 * t38 + t79 * t41;
t85 = -pkin(1) * t90 + qJ(2) * (g(3) * t100 + t99) + mrSges(2,1) * g(1) + t86;
t84 = -mrSges(3,1) * t76 + mrSges(3,3) * t75 - pkin(2) * t35 - t87;
t83 = -mrSges(2,2) * t75 + qJ(2) * t33 + t84;
t1 = [-qJ(1) * t99 + mrSges(1,3) * g(2) + (qJ(1) * (-m(2) - t100) - mrSges(3,1) + t98) * g(3) - t101, -mrSges(2,3) * g(1) - t83, t86, t41, -t82 * Ifges(5,6) + t89; -mrSges(1,3) * g(1) + (mrSges(1,1) - mrSges(2,2)) * g(3) + t85, (mrSges(3,1) + mrSges(2,3)) * g(3) + t101, t84, t38, -mrSges(5,1) * t74 + t93; t83 - qJ(1) * t90 + (qJ(1) * m(2) - t98) * g(1) - mrSges(1,1) * g(2), -mrSges(2,2) * g(3) + t85, -mrSges(3,1) * g(3) - mrSges(3,2) * t75 + t91, t87, t94;];
m_new  = t1;
