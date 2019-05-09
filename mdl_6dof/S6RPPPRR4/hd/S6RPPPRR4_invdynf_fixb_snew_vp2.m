% Calculate vector of cutting forces with Newton-Euler
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:45:49
% EndTime: 2019-05-05 13:45:50
% DurationCPUTime: 0.62s
% Computational Cost: add. (5927->127), mult. (10069->154), div. (0->0), fcn. (4313->8), ass. (0->69)
t64 = qJD(1) ^ 2;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t84 = qJD(1) * qJD(5);
t79 = t61 * t84;
t41 = qJDD(1) * t58 + t79;
t80 = t58 * t84;
t42 = -qJDD(1) * t61 + t80;
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t77 = -g(1) * t62 - g(2) * t59;
t71 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t77;
t94 = -pkin(1) - pkin(2);
t31 = t94 * t64 + t71;
t81 = g(1) * t59 - t62 * g(2);
t68 = -qJ(2) * t64 + qJDD(2) - t81;
t33 = t94 * qJDD(1) + t68;
t55 = sin(pkin(9));
t56 = cos(pkin(9));
t88 = t56 * t31 + t55 * t33;
t96 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t88;
t67 = (-pkin(3) - pkin(7)) * t64 + t96;
t14 = (-t42 - t80) * pkin(8) + (-t41 - t79) * pkin(5) + t67;
t40 = (-pkin(5) * t58 + pkin(8) * t61) * qJD(1);
t63 = qJD(5) ^ 2;
t86 = qJD(1) * t58;
t78 = -t55 * t31 + t33 * t56;
t22 = qJDD(1) * pkin(3) - qJ(4) * t64 + qJDD(4) - t78;
t20 = qJDD(1) * pkin(7) + t22;
t53 = g(3) + qJDD(3);
t87 = t58 * t20 + t61 * t53;
t16 = -pkin(5) * t63 + qJDD(5) * pkin(8) + t40 * t86 + t87;
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t85 = qJD(1) * t61;
t37 = qJD(5) * t60 + t57 * t85;
t24 = qJD(6) * t37 + qJDD(5) * t57 + t42 * t60;
t38 = qJD(5) * t57 - t60 * t85;
t25 = -mrSges(7,1) * t37 + mrSges(7,2) * t38;
t45 = qJD(6) - t86;
t29 = -mrSges(7,2) * t45 + mrSges(7,3) * t37;
t36 = qJDD(6) - t41;
t12 = m(7) * (t14 * t60 - t16 * t57) - t24 * mrSges(7,3) + t36 * mrSges(7,1) - t38 * t25 + t45 * t29;
t23 = -qJD(6) * t38 + qJDD(5) * t60 - t42 * t57;
t30 = mrSges(7,1) * t45 - mrSges(7,3) * t38;
t13 = m(7) * (t14 * t57 + t16 * t60) + t23 * mrSges(7,3) - t36 * mrSges(7,2) + t37 * t25 - t45 * t30;
t43 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t86;
t44 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t85;
t65 = m(6) * t67 - t41 * mrSges(6,1) + t42 * mrSges(6,2) - (t58 * t43 + t61 * t44) * qJD(1) + t60 * t12 + t57 * t13;
t97 = m(5) * (pkin(3) * t64 - t96) - t65;
t95 = -m(2) - m(3);
t92 = t58 * t53;
t91 = mrSges(2,1) + mrSges(3,1);
t90 = -mrSges(4,1) + mrSges(5,2);
t89 = mrSges(4,2) - mrSges(5,3);
t39 = (-mrSges(6,1) * t58 - mrSges(6,2) * t61) * qJD(1);
t7 = m(6) * t87 - qJDD(5) * mrSges(6,2) + t41 * mrSges(6,3) - qJD(5) * t44 - t57 * t12 + t60 * t13 + t39 * t86;
t66 = m(7) * (-qJDD(5) * pkin(5) - pkin(8) * t63 + t92 + (-qJD(1) * t40 - t20) * t61) - t23 * mrSges(7,1) + t24 * mrSges(7,2) - t37 * t29 + t38 * t30;
t9 = m(6) * (t61 * t20 - t92) - t42 * mrSges(6,3) + qJDD(5) * mrSges(6,1) + t39 * t85 + qJD(5) * t43 - t66;
t82 = m(5) * t53 - t58 * t9 + t61 * t7;
t76 = m(4) * t53 + t82;
t69 = -m(5) * t22 - t58 * t7 - t61 * t9;
t3 = m(4) * t78 + t90 * qJDD(1) - t89 * t64 + t69;
t5 = m(4) * t88 + t89 * qJDD(1) + t90 * t64 - t97;
t73 = -t55 * t3 + t56 * t5 + m(3) * (-pkin(1) * t64 + t71) + qJDD(1) * mrSges(3,3);
t70 = -m(3) * (-qJDD(1) * pkin(1) + t68) - t56 * t3 - t55 * t5;
t2 = m(2) * t81 + (-mrSges(2,2) + mrSges(3,3)) * t64 + t91 * qJDD(1) + t70;
t1 = m(2) * t77 - qJDD(1) * mrSges(2,2) - t91 * t64 + t73;
t4 = [-m(1) * g(1) + t1 * t62 - t2 * t59, t1, -t64 * mrSges(3,1) + t73, t5, t82, t7, t13; -m(1) * g(2) + t1 * t59 + t2 * t62, t2, -m(3) * g(3) - t76, t3, -t64 * mrSges(5,2) + qJDD(1) * mrSges(5,3) + t97, t9, t12; (-m(1) + t95) * g(3) - t76, t95 * g(3) - t76, -qJDD(1) * mrSges(3,1) - t64 * mrSges(3,3) - t70, t76, -qJDD(1) * mrSges(5,2) - t64 * mrSges(5,3) - t69, t65, t66;];
f_new  = t4;
