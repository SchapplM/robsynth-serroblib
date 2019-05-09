% Calculate vector of cutting forces with Newton-Euler
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-05-05 13:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:31:57
% EndTime: 2019-05-05 13:31:58
% DurationCPUTime: 0.54s
% Computational Cost: add. (4800->124), mult. (8458->151), div. (0->0), fcn. (4074->8), ass. (0->68)
t69 = sin(qJ(1));
t72 = cos(qJ(1));
t90 = t69 * g(1) - t72 * g(2);
t43 = qJDD(1) * pkin(1) + t90;
t74 = qJD(1) ^ 2;
t85 = -t72 * g(1) - t69 * g(2);
t45 = -t74 * pkin(1) + t85;
t65 = sin(pkin(9));
t66 = cos(pkin(9));
t100 = t66 * t43 - t65 * t45;
t26 = -qJDD(1) * pkin(2) - t74 * qJ(3) + qJDD(3) - t100;
t105 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t26;
t99 = t65 * t43 + t66 * t45;
t104 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t99;
t63 = -g(3) + qJDD(2);
t68 = sin(qJ(5));
t103 = t68 * t63;
t102 = mrSges(4,2) - mrSges(5,3);
t101 = -mrSges(5,2) - mrSges(4,3);
t78 = qJDD(4) + (-pkin(2) - qJ(4)) * t74 + t104;
t20 = -qJDD(1) * pkin(7) + t78;
t71 = cos(qJ(5));
t98 = t68 * t20 + t71 * t63;
t97 = qJD(1) * t71;
t96 = t68 * qJD(1);
t95 = qJD(1) * qJD(5);
t92 = mrSges(3,1) - t102;
t89 = t68 * t95;
t88 = t71 * t95;
t47 = -t68 * qJDD(1) - t88;
t48 = t71 * qJDD(1) - t89;
t77 = -t74 * pkin(7) - t105;
t14 = (-t48 + t89) * pkin(8) + (-t47 + t88) * pkin(5) + t77;
t46 = (pkin(5) * t68 - pkin(8) * t71) * qJD(1);
t73 = qJD(5) ^ 2;
t16 = -t73 * pkin(5) + qJDD(5) * pkin(8) - t46 * t96 + t98;
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t41 = t70 * qJD(5) - t67 * t97;
t28 = t41 * qJD(6) + t67 * qJDD(5) + t70 * t48;
t42 = t67 * qJD(5) + t70 * t97;
t29 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t51 = qJD(6) + t96;
t30 = -t51 * mrSges(7,2) + t41 * mrSges(7,3);
t40 = qJDD(6) - t47;
t12 = m(7) * (t70 * t14 - t67 * t16) - t28 * mrSges(7,3) + t40 * mrSges(7,1) - t42 * t29 + t51 * t30;
t27 = -t42 * qJD(6) + t70 * qJDD(5) - t67 * t48;
t31 = t51 * mrSges(7,1) - t42 * mrSges(7,3);
t13 = m(7) * (t67 * t14 + t70 * t16) + t27 * mrSges(7,3) - t40 * mrSges(7,2) + t41 * t29 - t51 * t31;
t44 = (mrSges(6,1) * t68 + mrSges(6,2) * t71) * qJD(1);
t50 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t97;
t7 = m(6) * t98 - qJDD(5) * mrSges(6,2) + t47 * mrSges(6,3) - qJD(5) * t50 - t67 * t12 + t70 * t13 - t44 * t96;
t49 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t96;
t75 = m(7) * (-qJDD(5) * pkin(5) - t73 * pkin(8) + t103 + (qJD(1) * t46 - t20) * t71) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t41 * t30 + t42 * t31;
t9 = m(6) * (t71 * t20 - t103) - t48 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t44 * t97 + qJD(5) * t49 - t75;
t87 = m(5) * t78 + qJDD(1) * mrSges(5,2) + t68 * t7 + t71 * t9;
t86 = m(5) * t63 - t68 * t9 + t71 * t7;
t84 = m(4) * t63 + t86;
t83 = m(4) * (t74 * pkin(2) - t104) - t87;
t82 = m(3) * t63 + t84;
t80 = m(6) * t77 - t47 * mrSges(6,1) + t48 * mrSges(6,2) + t70 * t12 + t67 * t13 + t49 * t96 + t50 * t97;
t79 = m(5) * t105 - t80;
t76 = m(4) * t26 + t79;
t4 = m(3) * t100 + (-mrSges(3,2) - t101) * t74 + t92 * qJDD(1) - t76;
t3 = m(3) * t99 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(1) - t92 * t74 - t83;
t2 = m(2) * t85 - t74 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t66 * t3 - t65 * t4;
t1 = m(2) * t90 + qJDD(1) * mrSges(2,1) - t74 * mrSges(2,2) + t65 * t3 + t66 * t4;
t5 = [-m(1) * g(1) - t69 * t1 + t72 * t2, t2, t3, t84, t86, t7, t13; -m(1) * g(2) + t72 * t1 + t69 * t2, t1, t4, -qJDD(1) * mrSges(4,3) - t102 * t74 + t83, -t74 * mrSges(5,2) - qJDD(1) * mrSges(5,3) + t79, t9, t12; (-m(1) - m(2)) * g(3) + t82, -m(2) * g(3) + t82, t82, t102 * qJDD(1) + t101 * t74 + t76, -t74 * mrSges(5,3) + t87, t80, t75;];
f_new  = t5;
