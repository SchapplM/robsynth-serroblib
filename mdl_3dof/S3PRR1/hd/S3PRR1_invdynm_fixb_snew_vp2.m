% Calculate vector of cutting torques with Newton-Euler for
% S3PRR1
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
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2019-05-04 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3PRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:24:49
% EndTime: 2019-05-04 18:24:50
% DurationCPUTime: 0.19s
% Computational Cost: add. (1180->72), mult. (1424->85), div. (0->0), fcn. (660->4), ass. (0->35)
t82 = m(3) + m(4);
t81 = -mrSges(1,2) - mrSges(2,3);
t66 = -g(2) + qJDD(1);
t68 = sin(qJ(2));
t70 = cos(qJ(2));
t54 = t68 * g(1) + t70 * t66;
t51 = qJDD(2) * pkin(2) + t54;
t55 = -t70 * g(1) + t68 * t66;
t73 = qJD(2) ^ 2;
t52 = -t73 * pkin(2) + t55;
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t49 = t69 * t51 - t67 * t52;
t64 = qJD(2) + qJD(3);
t62 = t64 ^ 2;
t63 = qJDD(2) + qJDD(3);
t46 = m(4) * t49 + t63 * mrSges(4,1) - t62 * mrSges(4,2);
t50 = t67 * t51 + t69 * t52;
t47 = m(4) * t50 - t62 * mrSges(4,1) - t63 * mrSges(4,2);
t40 = t69 * t46 + t67 * t47;
t37 = m(3) * t54 + qJDD(2) * mrSges(3,1) - t73 * mrSges(3,2) + t40;
t79 = -t67 * t46 + t69 * t47;
t38 = m(3) * t55 - t73 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t79;
t80 = -t68 * t37 + t70 * t38;
t78 = mrSges(4,1) * t49 - mrSges(4,2) * t50 + Ifges(4,3) * t63;
t44 = mrSges(4,1) * g(3) + mrSges(4,3) * t50 + t62 * Ifges(4,5) + Ifges(4,6) * t63;
t45 = -mrSges(4,2) * g(3) - mrSges(4,3) * t49 + Ifges(4,5) * t63 - t62 * Ifges(4,6);
t32 = Ifges(3,6) * qJDD(2) + t73 * Ifges(3,5) + mrSges(3,3) * t55 + t67 * t45 + t69 * t44 + pkin(4) * t79 + (pkin(2) * m(4) + mrSges(3,1)) * g(3);
t35 = -mrSges(3,2) * g(3) - mrSges(3,3) * t54 + Ifges(3,5) * qJDD(2) - t73 * Ifges(3,6) - pkin(4) * t40 - t67 * t44 + t69 * t45;
t77 = mrSges(2,2) * g(1) + pkin(3) * t80 + t70 * t32 + t68 * t35 + (pkin(1) * t82 + mrSges(2,1)) * g(3);
t30 = t70 * t37 + t68 * t38;
t76 = mrSges(2,2) * t66 - pkin(3) * t30 - t68 * t32 + t70 * t35;
t75 = mrSges(3,1) * t54 - mrSges(3,2) * t55 + Ifges(3,3) * qJDD(2) + pkin(2) * t40 + t78;
t74 = mrSges(2,1) * t66 + pkin(1) * t30 + t75;
t1 = [mrSges(1,3) * g(2) + (qJ(1) * (-m(2) - t82) + t81) * g(3) + t76, -mrSges(2,3) * g(3) + t76, t35, t45; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t77, -mrSges(2,3) * g(1) - t74, t32, t44; -qJ(1) * t80 - mrSges(1,1) * g(2) + (qJ(1) * m(2) - t81) * g(1) + t74, t77, t75, t78;];
m_new  = t1;
