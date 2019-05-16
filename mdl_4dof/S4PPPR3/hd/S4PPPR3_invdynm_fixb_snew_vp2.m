% Calculate vector of cutting torques with Newton-Euler for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2019-05-04 18:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:40:36
% EndTime: 2019-05-04 18:40:37
% DurationCPUTime: 0.23s
% Computational Cost: add. (1783->79), mult. (1867->88), div. (0->0), fcn. (1120->4), ass. (0->37)
t66 = g(3) + qJDD(3);
t86 = (m(4) + m(5)) * t66;
t67 = -g(2) + qJDD(1);
t68 = -g(1) + qJDD(2);
t69 = sin(pkin(5));
t70 = cos(pkin(5));
t59 = -t69 * t67 + t70 * t68;
t60 = t70 * t67 + t69 * t68;
t71 = sin(qJ(4));
t72 = cos(qJ(4));
t54 = t72 * t59 - t71 * t60;
t75 = qJD(4) ^ 2;
t50 = m(5) * t54 + qJDD(4) * mrSges(5,1) - t75 * mrSges(5,2);
t55 = t71 * t59 + t72 * t60;
t51 = m(5) * t55 - t75 * mrSges(5,1) - qJDD(4) * mrSges(5,2);
t44 = t72 * t50 + t71 * t51;
t41 = m(4) * t59 + t44;
t84 = -t71 * t50 + t72 * t51;
t42 = m(4) * t60 + t84;
t85 = -t69 * t41 + t70 * t42;
t29 = m(3) * t67 + t85;
t45 = -mrSges(5,1) * t66 + mrSges(5,3) * t55 + t75 * Ifges(5,5) + Ifges(5,6) * qJDD(4);
t46 = mrSges(5,2) * t66 - mrSges(5,3) * t54 + Ifges(5,5) * qJDD(4) - t75 * Ifges(5,6);
t34 = mrSges(4,3) * t60 + t71 * t46 + t72 * t45 + pkin(4) * t84 + (-pkin(3) * m(5) - mrSges(4,1)) * t66;
t37 = mrSges(4,2) * t66 - mrSges(4,3) * t59 - pkin(4) * t44 - t71 * t45 + t72 * t46;
t83 = -pkin(2) * t86 + qJ(3) * t85 + t70 * t34 + t69 * t37;
t91 = -pkin(1) * t29 - (-mrSges(3,2) + mrSges(2,1)) * t67 - t83;
t89 = -mrSges(1,2) - mrSges(2,3);
t32 = t70 * t41 + t69 * t42;
t87 = mrSges(5,1) * t54 - mrSges(5,2) * t55 + Ifges(5,3) * qJDD(4);
t82 = m(3) * t68 + t32;
t81 = mrSges(4,1) * t59 - mrSges(4,2) * t60 + pkin(3) * t44 + t87;
t79 = mrSges(3,2) * t68 + mrSges(3,3) * g(3) - qJ(3) * t32 - t69 * t34 + t70 * t37;
t78 = -pkin(1) * t82 + qJ(2) * (m(3) * g(3) + t86) + mrSges(2,1) * g(1) + t79;
t77 = -mrSges(3,1) * t68 + mrSges(3,3) * t67 - pkin(2) * t32 - t81;
t76 = -mrSges(2,2) * t67 + qJ(2) * t29 + t77;
t1 = [mrSges(1,3) * g(2) - qJ(1) * t86 + (qJ(1) * (-m(2) - m(3)) - mrSges(3,1) + t89) * g(3) - t91, -mrSges(2,3) * g(1) - t76, t79, t37, t46; -mrSges(1,3) * g(1) + (mrSges(1,1) - mrSges(2,2)) * g(3) + t78, (mrSges(3,1) + mrSges(2,3)) * g(3) + t91, t77, t34, t45; -qJ(1) * t82 + (m(2) * qJ(1) - t89) * g(1) + t76 - mrSges(1,1) * g(2), -mrSges(2,2) * g(3) + t78, -mrSges(3,1) * g(3) - mrSges(3,2) * t67 + t83, t81, t87;];
m_new  = t1;
