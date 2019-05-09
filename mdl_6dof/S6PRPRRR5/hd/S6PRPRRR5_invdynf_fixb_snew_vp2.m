% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:15:39
% EndTime: 2019-05-05 01:15:43
% DurationCPUTime: 1.67s
% Computational Cost: add. (21314->151), mult. (40300->202), div. (0->0), fcn. (27444->12), ass. (0->87)
t81 = sin(pkin(6));
t87 = sin(qJ(2));
t116 = t81 * t87;
t80 = sin(pkin(11));
t82 = cos(pkin(11));
t66 = t80 * g(1) - t82 * g(2);
t83 = cos(pkin(6));
t117 = t66 * t83;
t67 = -t82 * g(1) - t80 * g(2);
t79 = -g(3) + qJDD(1);
t91 = cos(qJ(2));
t107 = t79 * t116 + t87 * t117 + t91 * t67;
t103 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t107;
t85 = sin(qJ(5));
t86 = sin(qJ(4));
t89 = cos(qJ(5));
t90 = cos(qJ(4));
t55 = (t85 * t90 + t86 * t89) * qJD(2);
t121 = (t79 * t81 + t117) * t91 - t87 * t67;
t120 = -pkin(2) - pkin(8);
t92 = qJD(2) ^ 2;
t119 = pkin(4) * t92;
t114 = (-mrSges(3,2) + mrSges(4,3));
t115 = mrSges(3,1) - mrSges(4,2);
t110 = qJD(2) * t90;
t109 = qJD(2) * qJD(4);
t97 = -t92 * qJ(3) + qJDD(3) - t121;
t34 = t120 * qJDD(2) + t97;
t31 = t90 * t34;
t52 = -t81 * t66 + t83 * t79;
t65 = t90 * qJDD(2) - t86 * t109;
t23 = (qJDD(4) * pkin(4)) - t65 * pkin(9) + t31 + (-pkin(9) * t109 - t90 * t119 - t52) * t86;
t112 = t86 * t34 + t90 * t52;
t64 = -t86 * qJDD(2) - t90 * t109;
t71 = (qJD(4) * pkin(4)) - pkin(9) * t110;
t78 = t86 ^ 2;
t24 = t64 * pkin(9) - qJD(4) * t71 - t78 * t119 + t112;
t113 = t85 * t23 + t89 * t24;
t56 = (-t85 * t86 + t89 * t90) * qJD(2);
t45 = t55 * pkin(5) - t56 * pkin(10);
t76 = qJD(4) + qJD(5);
t74 = t76 ^ 2;
t75 = qJDD(4) + qJDD(5);
t19 = -t74 * pkin(5) + t75 * pkin(10) - t55 * t45 + t113;
t38 = -t56 * qJD(5) + t89 * t64 - t85 * t65;
t39 = -t55 * qJD(5) + t85 * t64 + t89 * t65;
t95 = -t64 * pkin(4) + t71 * t110 + (-pkin(9) * t78 + t120) * t92 + t103;
t20 = (t55 * t76 - t39) * pkin(10) + (t56 * t76 - t38) * pkin(5) + t95;
t84 = sin(qJ(6));
t88 = cos(qJ(6));
t46 = -t84 * t56 + t88 * t76;
t26 = t46 * qJD(6) + t88 * t39 + t84 * t75;
t47 = t88 * t56 + t84 * t76;
t29 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t36 = qJDD(6) - t38;
t53 = qJD(6) + t55;
t41 = -t53 * mrSges(7,2) + t46 * mrSges(7,3);
t16 = m(7) * (-t84 * t19 + t88 * t20) - t26 * mrSges(7,3) + t36 * mrSges(7,1) - t47 * t29 + t53 * t41;
t25 = -t47 * qJD(6) - t84 * t39 + t88 * t75;
t42 = t53 * mrSges(7,1) - t47 * mrSges(7,3);
t17 = m(7) * (t88 * t19 + t84 * t20) + t25 * mrSges(7,3) - t36 * mrSges(7,2) + t46 * t29 - t53 * t42;
t44 = t55 * mrSges(6,1) + t56 * mrSges(6,2);
t50 = t76 * mrSges(6,1) - t56 * mrSges(6,3);
t12 = m(6) * t113 - t75 * mrSges(6,2) + t38 * mrSges(6,3) - t84 * t16 + t88 * t17 - t55 * t44 - t76 * t50;
t102 = t89 * t23 - t85 * t24;
t49 = -t76 * mrSges(6,2) - t55 * mrSges(6,3);
t96 = m(7) * (-t75 * pkin(5) - t74 * pkin(10) + t56 * t45 - t102) - t25 * mrSges(7,1) + t26 * mrSges(7,2) - t46 * t41 + t47 * t42;
t13 = m(6) * t102 + t75 * mrSges(6,1) - t39 * mrSges(6,3) - t56 * t44 + t76 * t49 - t96;
t63 = (mrSges(5,1) * t86 + mrSges(5,2) * t90) * qJD(2);
t111 = qJD(2) * t86;
t68 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t111;
t8 = m(5) * (-t86 * t52 + t31) - t65 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t63 * t110 + qJD(4) * t68 + t85 * t12 + t89 * t13;
t69 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t110;
t9 = m(5) * t112 - qJDD(4) * mrSges(5,2) + t64 * mrSges(5,3) - qJD(4) * t69 - t63 * t111 + t89 * t12 - t85 * t13;
t99 = -m(4) * (-qJDD(2) * pkin(2) + t97) - t90 * t8 - t86 * t9;
t4 = m(3) * t121 + t115 * qJDD(2) + (t114 * t92) + t99;
t118 = t4 * t91;
t98 = m(6) * t95 - t38 * mrSges(6,1) + t39 * mrSges(6,2) + t88 * t16 + t84 * t17 + t55 * t49 + t56 * t50;
t94 = -t64 * mrSges(5,1) + m(5) * (t120 * t92 + t103) + t68 * t111 + t69 * t110 + t65 * mrSges(5,2) + t98;
t93 = -m(4) * (t92 * pkin(2) - t103) + t94;
t11 = m(3) * t107 + t114 * qJDD(2) - t115 * t92 + t93;
t104 = m(4) * t52 - t86 * t8 + t90 * t9;
t6 = m(3) * t52 + t104;
t106 = m(2) * t79 + t11 * t116 + t81 * t118 + t83 * t6;
t2 = m(2) * t67 + t91 * t11 - t87 * t4;
t1 = m(2) * t66 - t81 * t6 + (t11 * t87 + t118) * t83;
t3 = [-m(1) * g(1) - t80 * t1 + t82 * t2, t2, t11, t104, t9, t12, t17; -m(1) * g(2) + t82 * t1 + t80 * t2, t1, t4, -(t92 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t93, t8, t13, t16; -m(1) * g(3) + t106, t106, t6, qJDD(2) * mrSges(4,2) - t92 * mrSges(4,3) - t99, t94, t98, t96;];
f_new  = t3;
