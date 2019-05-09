% Calculate vector of cutting torques with Newton-Euler for
% S4PPPR2
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
%   pkin=[a2,a3,a4,d4,theta2]';
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
% Datum: 2019-05-04 18:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:37:51
% EndTime: 2019-05-04 18:37:51
% DurationCPUTime: 0.16s
% Computational Cost: add. (1362->78), mult. (1534->89), div. (0->0), fcn. (906->4), ass. (0->39)
t70 = -g(2) + qJDD(1);
t71 = sin(pkin(5));
t72 = cos(pkin(5));
t63 = t71 * g(1) + t72 * t70;
t61 = qJDD(3) - t63;
t64 = -t72 * g(1) + t71 * t70;
t73 = sin(qJ(4));
t74 = cos(qJ(4));
t55 = t74 * t61 - t73 * t64;
t56 = t73 * t61 + t74 * t64;
t90 = mrSges(5,1) * t55 - mrSges(5,2) * t56 + Ifges(5,3) * qJDD(4);
t89 = -m(5) * pkin(3) - mrSges(4,1);
t69 = g(3) - qJDD(2);
t67 = m(4) * t69;
t65 = -m(5) * t69 - t67;
t77 = qJD(4) ^ 2;
t51 = m(5) * t55 + qJDD(4) * mrSges(5,1) - t77 * mrSges(5,2);
t52 = m(5) * t56 - t77 * mrSges(5,1) - qJDD(4) * mrSges(5,2);
t45 = t74 * t51 + t73 * t52;
t82 = -m(4) * t61 - t45;
t39 = m(3) * t63 + t82;
t46 = -t73 * t51 + t74 * t52;
t86 = m(4) * t64 + t46;
t40 = m(3) * t64 + t86;
t88 = -t71 * t39 + t72 * t40;
t87 = qJ(1) * m(2) + mrSges(1,2) + mrSges(2,3);
t48 = -mrSges(5,1) * t69 + mrSges(5,3) * t56 + t77 * Ifges(5,5) + Ifges(5,6) * qJDD(4);
t49 = mrSges(5,2) * t69 - mrSges(5,3) * t55 + Ifges(5,5) * qJDD(4) - t77 * Ifges(5,6);
t81 = pkin(4) * t46 + t74 * t48 + t73 * t49;
t32 = -pkin(2) * t65 + (mrSges(4,2) + mrSges(3,3)) * t64 + (mrSges(3,1) - t89) * t69 - t81;
t84 = mrSges(4,2) * t61 - pkin(4) * t45 - t73 * t48 + t74 * t49;
t37 = -mrSges(3,3) * t63 - qJ(3) * t65 + (-mrSges(3,2) + mrSges(4,3)) * t69 + t84;
t85 = qJ(2) * t88 + t72 * t32 + t71 * t37 + pkin(1) * (t67 + (m(3) + m(5)) * t69) + mrSges(2,2) * g(1) + mrSges(2,1) * g(3);
t36 = t72 * t39 + t71 * t40;
t83 = mrSges(2,2) * t70 - qJ(2) * t36 - t71 * t32 + t72 * t37;
t80 = -mrSges(4,1) * t61 + mrSges(4,3) * t64 - pkin(3) * t45 - t90;
t79 = mrSges(3,1) * t63 - mrSges(3,2) * t64 + pkin(2) * t82 + qJ(3) * t86 + t80;
t78 = mrSges(2,1) * t70 + pkin(1) * t36 + t79;
t1 = [mrSges(1,3) * g(2) + qJ(1) * (-m(3) * t69 + t65) - t87 * g(3) + t83, -mrSges(2,3) * g(3) + t83, t37, mrSges(4,3) * t69 + t84, t49; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t85, -mrSges(2,3) * g(1) - t78, t32, t80, t48; -mrSges(1,1) * g(2) + t87 * g(1) - qJ(1) * t88 + t78, t85, t79, -mrSges(4,2) * t64 + t89 * t69 + t81, t90;];
m_new  = t1;
