% Calculate vector of cutting forces with Newton-Euler
% S6PRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-04 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:50:14
% EndTime: 2019-05-04 21:50:16
% DurationCPUTime: 0.99s
% Computational Cost: add. (12343->123), mult. (21762->163), div. (0->0), fcn. (14359->12), ass. (0->78)
t68 = sin(pkin(6));
t77 = cos(qJ(2));
t105 = t68 * t77;
t67 = sin(pkin(10));
t70 = cos(pkin(10));
t54 = t67 * g(1) - t70 * g(2);
t71 = cos(pkin(6));
t107 = t54 * t71;
t55 = -t70 * g(1) - t67 * g(2);
t65 = -g(3) + qJDD(1);
t74 = sin(qJ(2));
t87 = t65 * t105 + t77 * t107 - t74 * t55;
t30 = qJDD(2) * pkin(2) + t87;
t79 = qJD(2) ^ 2;
t106 = t68 * t74;
t94 = t65 * t106 + t74 * t107 + t77 * t55;
t31 = -t79 * pkin(2) + t94;
t66 = sin(pkin(11));
t69 = cos(pkin(11));
t101 = t66 * t30 + t69 * t31;
t109 = -qJDD(2) * qJ(4) - (2 * qJD(4) * qJD(2)) - t101;
t108 = -pkin(3) - pkin(8);
t89 = -t68 * t54 + t71 * t65;
t40 = qJDD(3) + t89;
t73 = sin(qJ(5));
t104 = t73 * t40;
t103 = mrSges(4,1) - mrSges(5,2);
t102 = -mrSges(4,2) + mrSges(5,3);
t90 = t69 * t30 - t66 * t31;
t83 = -t79 * qJ(4) + qJDD(4) - t90;
t24 = t108 * qJDD(2) + t83;
t76 = cos(qJ(5));
t100 = t73 * t24 + t76 * t40;
t99 = qJD(2) * t76;
t98 = t73 * qJD(2);
t97 = qJD(2) * qJD(5);
t51 = (pkin(5) * t73 - pkin(9) * t76) * qJD(2);
t78 = qJD(5) ^ 2;
t19 = -t78 * pkin(5) + qJDD(5) * pkin(9) - t51 * t98 + t100;
t91 = t76 * t97;
t52 = -t73 * qJDD(2) - t91;
t92 = t73 * t97;
t53 = t76 * qJDD(2) - t92;
t84 = t108 * t79 - t109;
t20 = (-t53 + t92) * pkin(9) + (-t52 + t91) * pkin(5) + t84;
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t48 = t75 * qJD(5) - t72 * t99;
t33 = t48 * qJD(6) + t72 * qJDD(5) + t75 * t53;
t49 = t72 * qJD(5) + t75 * t99;
t34 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t60 = qJD(6) + t98;
t38 = -t60 * mrSges(7,2) + t48 * mrSges(7,3);
t46 = qJDD(6) - t52;
t16 = m(7) * (-t72 * t19 + t75 * t20) - t33 * mrSges(7,3) + t46 * mrSges(7,1) - t49 * t34 + t60 * t38;
t32 = -t49 * qJD(6) + t75 * qJDD(5) - t72 * t53;
t39 = t60 * mrSges(7,1) - t49 * mrSges(7,3);
t17 = m(7) * (t75 * t19 + t72 * t20) + t32 * mrSges(7,3) - t46 * mrSges(7,2) + t48 * t34 - t60 * t39;
t56 = -(qJD(5) * mrSges(6,2)) - mrSges(6,3) * t98;
t57 = (qJD(5) * mrSges(6,1)) - mrSges(6,3) * t99;
t82 = m(6) * t84 - t52 * mrSges(6,1) + t53 * mrSges(6,2) + t75 * t16 + t72 * t17 + t56 * t98 + t57 * t99;
t81 = -m(5) * (t79 * pkin(3) + t109) + t82;
t10 = m(4) * t101 + t102 * qJDD(2) - t103 * t79 + t81;
t50 = (mrSges(6,1) * t73 + mrSges(6,2) * t76) * qJD(2);
t12 = m(6) * t100 - qJDD(5) * mrSges(6,2) + t52 * mrSges(6,3) - qJD(5) * t57 - t72 * t16 + t75 * t17 - t50 * t98;
t80 = m(7) * (-qJDD(5) * pkin(5) - t78 * pkin(9) + t104 + (qJD(2) * t51 - t24) * t76) - t32 * mrSges(7,1) + t33 * mrSges(7,2) - t48 * t38 + t49 * t39;
t13 = m(6) * (t76 * t24 - t104) - t53 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t50 * t99 + qJD(5) * t56 - t80;
t85 = -m(5) * (-qJDD(2) * pkin(3) + t83) - t73 * t12 - t76 * t13;
t7 = m(4) * t90 + t103 * qJDD(2) + t102 * t79 + t85;
t5 = m(3) * t87 + qJDD(2) * mrSges(3,1) - t79 * mrSges(3,2) + t66 * t10 + t69 * t7;
t6 = m(3) * t94 - t79 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t69 * t10 - t66 * t7;
t88 = m(5) * t40 + t76 * t12 - t73 * t13;
t86 = m(4) * t40 + t88;
t9 = m(3) * t89 + t86;
t93 = m(2) * t65 + t5 * t105 + t6 * t106 + t71 * t9;
t2 = m(2) * t55 - t74 * t5 + t77 * t6;
t1 = m(2) * t54 - t68 * t9 + (t5 * t77 + t6 * t74) * t71;
t3 = [-m(1) * g(1) - t67 * t1 + t70 * t2, t2, t6, t10, t88, t12, t17; -m(1) * g(2) + t70 * t1 + t67 * t2, t1, t5, t7, -t79 * mrSges(5,2) - qJDD(2) * mrSges(5,3) - t81, t13, t16; -m(1) * g(3) + t93, t93, t9, t86, qJDD(2) * mrSges(5,2) - t79 * mrSges(5,3) - t85, t82, t80;];
f_new  = t3;
