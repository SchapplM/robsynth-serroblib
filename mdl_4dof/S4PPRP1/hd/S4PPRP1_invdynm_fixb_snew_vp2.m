% Calculate vector of cutting torques with Newton-Euler for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-05-04 18:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:43:08
% EndTime: 2019-05-04 18:43:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (1447->94), mult. (1984->103), div. (0->0), fcn. (974->4), ass. (0->39)
t94 = -m(4) - m(5);
t77 = sin(pkin(5));
t78 = cos(pkin(5));
t68 = t77 * g(1) - t78 * g(2);
t65 = qJDD(2) - t68;
t69 = -t78 * g(1) - t77 * g(2);
t79 = sin(qJ(3));
t80 = cos(qJ(3));
t58 = t79 * t65 + t80 * t69;
t81 = qJD(3) ^ 2;
t52 = -t81 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t58;
t93 = mrSges(5,2) * t52;
t92 = -mrSges(4,1) - mrSges(5,1);
t91 = m(5) * t52 + qJDD(3) * mrSges(5,3);
t57 = t80 * t65 - t79 * t69;
t55 = -qJDD(3) * pkin(3) - t81 * qJ(4) + qJDD(4) - t57;
t90 = -mrSges(5,1) * t55 + mrSges(5,3) * t52 + Ifges(5,2) * qJDD(3);
t89 = mrSges(5,2) * t55 + Ifges(5,4) * qJDD(3) + t81 * Ifges(5,6);
t46 = m(4) * t58 - qJDD(3) * mrSges(4,2) + t92 * t81 + t91;
t49 = -m(5) * t55 + qJDD(3) * mrSges(5,1) + t81 * mrSges(5,3);
t47 = m(4) * t57 + qJDD(3) * mrSges(4,1) - t81 * mrSges(4,2) + t49;
t41 = t80 * t46 - t79 * t47;
t88 = m(3) * t69 + t41;
t40 = t79 * t46 + t80 * t47;
t74 = g(3) - qJDD(1);
t43 = t93 + mrSges(4,3) * t58 + (Ifges(5,4) + Ifges(4,5)) * t81 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + (-m(5) * pkin(3) + t92) * t74;
t44 = -mrSges(4,3) * t57 + Ifges(4,5) * qJDD(3) - t81 * Ifges(4,6) + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t74 + t89;
t87 = mrSges(3,2) * t65 - pkin(4) * t40 - t79 * t43 + t80 * t44;
t86 = -m(3) * t65 - t40;
t85 = -pkin(2) * t94 * t74 - pkin(4) * t41 - t80 * t43 - t79 * t44;
t84 = mrSges(4,1) * t57 + Ifges(4,3) * qJDD(3) + pkin(3) * t49 + qJ(4) * (-t81 * mrSges(5,1) + t91) - mrSges(4,2) * t58 + t90;
t83 = -mrSges(3,1) * t65 + mrSges(3,3) * t69 - pkin(2) * t40 - t84;
t82 = mrSges(2,1) * t68 - mrSges(2,2) * t69 + pkin(1) * t86 + qJ(2) * t88 + t83;
t59 = (-m(3) + t94) * t74;
t37 = m(2) * t69 + t88;
t36 = m(2) * t68 + t86;
t35 = -mrSges(2,3) * t68 - qJ(2) * t59 + (-mrSges(2,2) + mrSges(3,3)) * t74 + t87;
t34 = -pkin(1) * t59 + (mrSges(2,1) + mrSges(3,1)) * t74 + (mrSges(3,2) + mrSges(2,3)) * t69 + t85;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t78 * t35 - t77 * t34 - qJ(1) * (t78 * t36 + t77 * t37), t35, mrSges(3,3) * t74 + t87, t44, -mrSges(5,3) * t74 + t89; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t77 * t35 + t78 * t34 + qJ(1) * (-t77 * t36 + t78 * t37), t34, t83, t43, t90; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t82, t82, -mrSges(3,1) * t74 - mrSges(3,2) * t69 - t85, t84, mrSges(5,1) * t74 - t81 * Ifges(5,4) + Ifges(5,6) * qJDD(3) - t93;];
m_new  = t1;
