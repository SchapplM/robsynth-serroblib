% Calculate vector of cutting torques with Newton-Euler for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2019-05-04 18:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:35:16
% EndTime: 2019-05-04 18:35:16
% DurationCPUTime: 0.17s
% Computational Cost: add. (1108->82), mult. (1397->90), div. (0->0), fcn. (828->4), ass. (0->39)
t70 = sin(pkin(5));
t71 = cos(pkin(5));
t64 = t70 * g(1) - t71 * g(2);
t61 = qJDD(2) - t64;
t65 = t71 * g(1) + t70 * g(2);
t62 = qJDD(3) - t65;
t72 = sin(qJ(4));
t73 = cos(qJ(4));
t53 = -t72 * t61 + t73 * t62;
t74 = qJD(4) ^ 2;
t50 = m(5) * t53 + qJDD(4) * mrSges(5,1) - t74 * mrSges(5,2);
t54 = t73 * t61 + t72 * t62;
t51 = m(5) * t54 - t74 * mrSges(5,1) - qJDD(4) * mrSges(5,2);
t85 = -t72 * t50 + t73 * t51;
t35 = m(4) * t61 + t85;
t69 = g(3) - qJDD(1);
t44 = mrSges(5,1) * t69 + mrSges(5,3) * t54 + t74 * Ifges(5,5) + Ifges(5,6) * qJDD(4);
t45 = -mrSges(5,2) * t69 - mrSges(5,3) * t53 + Ifges(5,5) * qJDD(4) - t74 * Ifges(5,6);
t83 = mrSges(4,1) * t69 + pkin(4) * t85 + t73 * t44 + t72 * t45;
t89 = pkin(2) * t35 + (mrSges(3,1) - mrSges(4,2)) * t61 + t83;
t88 = -m(3) - m(4);
t67 = m(5) * t69;
t86 = mrSges(3,2) - mrSges(4,3);
t39 = t73 * t50 + t72 * t51;
t36 = -m(4) * t62 - t39;
t84 = -m(5) * pkin(3) - mrSges(3,3);
t82 = mrSges(5,1) * t53 - mrSges(5,2) * t54 + Ifges(5,3) * qJDD(4);
t81 = t88 * t61 - t85;
t79 = -mrSges(4,2) * t62 + pkin(4) * t39 + t72 * t44 - t73 * t45;
t78 = -mrSges(4,1) * t62 + mrSges(4,3) * t61 - pkin(3) * t39 - t82;
t77 = pkin(2) * t36 + qJ(3) * (-m(4) * t69 - t67) - t79;
t76 = mrSges(3,2) * t61 - qJ(3) * t35 - t78;
t75 = pkin(1) * t81 + qJ(2) * (-m(3) * t65 - t36) + mrSges(2,1) * t64 + (mrSges(2,2) - mrSges(3,3)) * t65 + t76;
t55 = t88 * t69 - t67;
t32 = (-m(2) - m(3)) * t65 - t36;
t31 = m(2) * t64 + t81;
t30 = -mrSges(2,3) * t64 - qJ(2) * t55 + (-mrSges(2,2) - t84) * t69 + t89;
t29 = -pkin(1) * t55 + (-mrSges(3,1) - mrSges(2,3)) * t65 + (mrSges(2,1) - t86) * t69 - t77;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t71 * t30 - t70 * t29 - qJ(1) * (t71 * t31 + t70 * t32), t30, -mrSges(3,3) * t65 + t76, -mrSges(4,3) * t69 - t79, t45; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t70 * t30 + t71 * t29 + qJ(1) * (-t70 * t31 + t71 * t32), t29, t84 * t69 - t89, t78, t44; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t75, t75, mrSges(3,1) * t65 + t86 * t69 + t77, -mrSges(4,2) * t61 + pkin(3) * t67 + t83, t82;];
m_new  = t1;
