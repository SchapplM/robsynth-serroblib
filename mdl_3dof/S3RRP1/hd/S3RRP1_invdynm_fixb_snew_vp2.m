% Calculate vector of cutting torques with Newton-Euler for
% S3RRP1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2019-05-04 18:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3RRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:31:14
% EndTime: 2019-05-04 18:31:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (1335->88), mult. (1668->101), div. (0->0), fcn. (660->4), ass. (0->35)
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t58 = t73 * g(1) - t75 * g(2);
t55 = qJDD(1) * pkin(1) + t58;
t59 = -t75 * g(1) - t73 * g(2);
t77 = qJD(1) ^ 2;
t56 = -t77 * pkin(1) + t59;
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t52 = t72 * t55 + t74 * t56;
t70 = (qJD(1) + qJD(2));
t68 = t70 ^ 2;
t69 = qJDD(1) + qJDD(2);
t47 = -t68 * pkin(2) + t69 * qJ(3) + (2 * qJD(3) * t70) + t52;
t86 = mrSges(4,2) * t47;
t85 = -mrSges(3,1) - mrSges(4,1);
t84 = m(4) * t47 + t69 * mrSges(4,3);
t40 = m(3) * t52 - t69 * mrSges(3,2) + t85 * t68 + t84;
t51 = t74 * t55 - t72 * t56;
t49 = -t69 * pkin(2) - t68 * qJ(3) + qJDD(3) - t51;
t81 = -m(4) * t49 + t69 * mrSges(4,1) + t68 * mrSges(4,3);
t42 = m(3) * t51 + t69 * mrSges(3,1) - t68 * mrSges(3,2) + t81;
t35 = t72 * t40 + t74 * t42;
t83 = t74 * t40 - t72 * t42;
t82 = mrSges(4,2) * t49 + mrSges(4,3) * g(3) + Ifges(4,4) * t69 + t68 * Ifges(4,6);
t80 = -mrSges(4,1) * t49 + mrSges(4,3) * t47 + Ifges(4,2) * t69;
t79 = -mrSges(3,2) * t52 + qJ(3) * (-t68 * mrSges(4,1) + t84) + pkin(2) * t81 + mrSges(3,1) * t51 + Ifges(3,3) * t69 + t80;
t78 = mrSges(2,1) * t58 - mrSges(2,2) * t59 + Ifges(2,3) * qJDD(1) + pkin(1) * t35 + t79;
t37 = -mrSges(3,3) * t51 + Ifges(3,5) * t69 - t68 * Ifges(3,6) + (m(4) * qJ(3) - mrSges(3,2)) * g(3) + t82;
t36 = t86 + mrSges(3,3) * t52 + (Ifges(3,6) - Ifges(4,6)) * t69 + (Ifges(4,4) + Ifges(3,5)) * t68 + (m(4) * pkin(2) - t85) * g(3);
t33 = m(2) * t59 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t83;
t32 = m(2) * t58 + qJDD(1) * mrSges(2,1) - t77 * mrSges(2,2) + t35;
t31 = -mrSges(2,2) * g(3) - mrSges(2,3) * t58 + Ifges(2,5) * qJDD(1) - t77 * Ifges(2,6) - pkin(4) * t35 - t72 * t36 + t74 * t37;
t30 = Ifges(2,6) * qJDD(1) + t77 * Ifges(2,5) + mrSges(2,3) * t59 + t72 * t37 + t74 * t36 + pkin(4) * t83 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t75 * t31 - t73 * t30 - pkin(3) * (t75 * t32 + t73 * t33), t31, t37, t82; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t73 * t31 + t75 * t30 + pkin(3) * (-t73 * t32 + t75 * t33), t30, t36, t80; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t78, t78, t79, -mrSges(4,1) * g(3) - t68 * Ifges(4,4) + Ifges(4,6) * t69 - t86;];
m_new  = t1;
