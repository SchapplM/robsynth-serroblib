% Calculate vector of cutting torques with Newton-Euler for
% S3RRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x4]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3RRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:32:32
% EndTime: 2019-05-04 18:32:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (3133->89), mult. (4111->110), div. (0->0), fcn. (2126->6), ass. (0->43)
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t70 = t84 * g(1) - t87 * g(2);
t67 = qJDD(1) * pkin(1) + t70;
t71 = -t87 * g(1) - t84 * g(2);
t88 = qJD(1) ^ 2;
t68 = -t88 * pkin(1) + t71;
t83 = sin(qJ(2));
t86 = cos(qJ(2));
t62 = t86 * t67 - t83 * t68;
t79 = qJDD(1) + qJDD(2);
t59 = t79 * pkin(2) + t62;
t63 = t83 * t67 + t86 * t68;
t80 = qJD(1) + qJD(2);
t78 = t80 ^ 2;
t60 = -t78 * pkin(2) + t63;
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t57 = t85 * t59 - t82 * t60;
t76 = qJD(3) + t80;
t74 = t76 ^ 2;
t75 = qJDD(3) + t79;
t54 = m(4) * t57 + t75 * mrSges(4,1) - t74 * mrSges(4,2);
t58 = t82 * t59 + t85 * t60;
t55 = m(4) * t58 - t74 * mrSges(4,1) - t75 * mrSges(4,2);
t48 = t85 * t54 + t82 * t55;
t45 = m(3) * t62 + t79 * mrSges(3,1) - t78 * mrSges(3,2) + t48;
t92 = -t82 * t54 + t85 * t55;
t46 = m(3) * t63 - t78 * mrSges(3,1) - t79 * mrSges(3,2) + t92;
t39 = t86 * t45 + t83 * t46;
t93 = -t83 * t45 + t86 * t46;
t91 = mrSges(4,1) * t57 - mrSges(4,2) * t58 + Ifges(4,3) * t75;
t90 = mrSges(3,1) * t62 - mrSges(3,2) * t63 + Ifges(3,3) * t79 + pkin(2) * t48 + t91;
t89 = mrSges(2,1) * t70 - mrSges(2,2) * t71 + Ifges(2,3) * qJDD(1) + pkin(1) * t39 + t90;
t53 = -mrSges(4,2) * g(3) - mrSges(4,3) * t57 + Ifges(4,5) * t75 - t74 * Ifges(4,6);
t52 = mrSges(4,1) * g(3) + mrSges(4,3) * t58 + t74 * Ifges(4,5) + Ifges(4,6) * t75;
t41 = -mrSges(3,2) * g(3) - mrSges(3,3) * t62 + Ifges(3,5) * t79 - t78 * Ifges(3,6) - pkin(5) * t48 - t82 * t52 + t85 * t53;
t40 = Ifges(3,6) * t79 + t78 * Ifges(3,5) + mrSges(3,3) * t63 + t82 * t53 + t85 * t52 + pkin(5) * t92 + (pkin(2) * m(4) + mrSges(3,1)) * g(3);
t37 = m(2) * t71 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t93;
t36 = m(2) * t70 + qJDD(1) * mrSges(2,1) - t88 * mrSges(2,2) + t39;
t35 = -mrSges(2,2) * g(3) - mrSges(2,3) * t70 + Ifges(2,5) * qJDD(1) - t88 * Ifges(2,6) - pkin(4) * t39 - t83 * t40 + t86 * t41;
t34 = Ifges(2,6) * qJDD(1) + t88 * Ifges(2,5) + mrSges(2,3) * t71 + t83 * t41 + t86 * t40 + pkin(4) * t93 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t87 * t35 - t84 * t34 - pkin(3) * (t87 * t36 + t84 * t37), t35, t41, t53; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t84 * t35 + t87 * t34 + pkin(3) * (-t84 * t36 + t87 * t37), t34, t40, t52; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t89, t89, t90, t91;];
m_new  = t1;
