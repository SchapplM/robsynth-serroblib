% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-05-05 06:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:33:27
% EndTime: 2019-05-05 06:33:33
% DurationCPUTime: 2.64s
% Computational Cost: add. (36085->172), mult. (71361->226), div. (0->0), fcn. (49950->12), ass. (0->92)
t86 = sin(pkin(10));
t88 = cos(pkin(10));
t76 = g(1) * t86 - g(2) * t88;
t89 = cos(pkin(6));
t123 = t76 * t89;
t77 = -g(1) * t88 - g(2) * t86;
t84 = -g(3) + qJDD(1);
t87 = sin(pkin(6));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t128 = (t84 * t87 + t123) * t95 - t92 * t77;
t94 = cos(qJ(3));
t115 = t94 * qJD(2);
t122 = t87 * t92;
t111 = t122 * t84 + t123 * t92 + t77 * t95;
t97 = qJD(2) ^ 2;
t44 = -pkin(2) * t97 + qJDD(2) * pkin(8) + t111;
t61 = -t76 * t87 + t84 * t89;
t91 = sin(qJ(3));
t118 = t44 * t94 + t61 * t91;
t73 = (-pkin(3) * t94 - pkin(9) * t91) * qJD(2);
t96 = qJD(3) ^ 2;
t27 = -pkin(3) * t96 + qJDD(3) * pkin(9) + t115 * t73 + t118;
t114 = qJD(2) * qJD(3);
t108 = t94 * t114;
t109 = t91 * t114;
t43 = -qJDD(2) * pkin(2) - t97 * pkin(8) - t128;
t74 = qJDD(2) * t91 + t108;
t75 = qJDD(2) * t94 - t109;
t30 = (-t74 - t108) * pkin(9) + (-t75 + t109) * pkin(3) + t43;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t107 = -t90 * t27 + t30 * t93;
t117 = cos(pkin(11));
t126 = -2 * qJD(5);
t116 = qJD(2) * t91;
t70 = qJD(3) * t93 - t116 * t90;
t53 = qJD(4) * t70 + qJDD(3) * t90 + t74 * t93;
t67 = qJDD(4) - t75;
t71 = qJD(3) * t90 + t116 * t93;
t82 = qJD(4) - t115;
t20 = (t70 * t82 - t53) * qJ(5) + (t70 * t71 + t67) * pkin(4) + t107;
t120 = t27 * t93 + t30 * t90;
t52 = -qJD(4) * t71 + qJDD(3) * t93 - t74 * t90;
t59 = pkin(4) * t82 - qJ(5) * t71;
t66 = t70 ^ 2;
t22 = -pkin(4) * t66 + qJ(5) * t52 - t59 * t82 + t120;
t85 = sin(pkin(11));
t54 = -t117 * t70 + t71 * t85;
t112 = t117 * t22 + t126 * t54 + t20 * t85;
t55 = t117 * t71 + t70 * t85;
t37 = pkin(5) * t54 - qJ(6) * t55;
t47 = -mrSges(7,1) * t82 + mrSges(7,2) * t55;
t81 = t82 ^ 2;
t113 = m(7) * (-pkin(5) * t81 + qJ(6) * t67 + 0.2e1 * qJD(6) * t82 - t37 * t54 + t112) + t82 * t47 + t67 * mrSges(7,3);
t38 = mrSges(7,1) * t54 - mrSges(7,3) * t55;
t119 = -mrSges(6,1) * t54 - mrSges(6,2) * t55 - t38;
t121 = -mrSges(6,3) - mrSges(7,2);
t33 = -t117 * t52 + t53 * t85;
t46 = mrSges(6,1) * t82 - mrSges(6,3) * t55;
t12 = m(6) * t112 - t67 * mrSges(6,2) + t119 * t54 + t121 * t33 - t82 * t46 + t113;
t103 = t117 * t20 - t85 * t22;
t125 = m(7) * (-t67 * pkin(5) - t81 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t37) * t55 - t103);
t34 = t117 * t53 + t52 * t85;
t45 = -mrSges(6,2) * t82 - mrSges(6,3) * t54;
t48 = -mrSges(7,2) * t54 + mrSges(7,3) * t82;
t13 = m(6) * t103 - t125 + (t45 + t48) * t82 + (mrSges(6,1) + mrSges(7,1)) * t67 + (m(6) * t126 + t119) * t55 + t121 * t34;
t56 = -mrSges(5,1) * t70 + mrSges(5,2) * t71;
t58 = -mrSges(5,2) * t82 + mrSges(5,3) * t70;
t10 = m(5) * t107 + mrSges(5,1) * t67 - mrSges(5,3) * t53 + t117 * t13 + t12 * t85 - t56 * t71 + t58 * t82;
t60 = mrSges(5,1) * t82 - mrSges(5,3) * t71;
t11 = m(5) * t120 - mrSges(5,2) * t67 + mrSges(5,3) * t52 + t117 * t12 - t13 * t85 + t56 * t70 - t60 * t82;
t78 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t116;
t79 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t115;
t127 = m(4) * t43 - t75 * mrSges(4,1) + t74 * mrSges(4,2) + t93 * t10 + t90 * t11 + (t78 * t91 - t79 * t94) * qJD(2);
t8 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t97 * mrSges(3,2) - t127;
t124 = t8 * t95;
t106 = -t44 * t91 + t94 * t61;
t72 = (-mrSges(4,1) * t94 + mrSges(4,2) * t91) * qJD(2);
t26 = -qJDD(3) * pkin(3) - t96 * pkin(9) + t116 * t73 - t106;
t99 = -t52 * pkin(4) - t66 * qJ(5) + t59 * t71 + qJDD(5) + t26;
t102 = t34 * mrSges(7,3) + t55 * t47 - m(7) * (-0.2e1 * qJD(6) * t55 + (t54 * t82 - t34) * qJ(6) + (t55 * t82 + t33) * pkin(5) + t99) - t33 * mrSges(7,1) - t54 * t48;
t100 = m(6) * t99 + mrSges(6,1) * t33 + mrSges(6,2) * t34 + t45 * t54 + t46 * t55 - t102;
t98 = m(5) * t26 - mrSges(5,1) * t52 + mrSges(5,2) * t53 - t58 * t70 + t60 * t71 + t100;
t14 = m(4) * t106 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t74 + qJD(3) * t79 - t116 * t72 - t98;
t9 = m(4) * t118 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t75 - qJD(3) * t78 - t10 * t90 + t11 * t93 + t115 * t72;
t4 = m(3) * t111 - mrSges(3,1) * t97 - qJDD(2) * mrSges(3,2) - t14 * t91 + t9 * t94;
t6 = m(3) * t61 + t14 * t94 + t9 * t91;
t110 = m(2) * t84 + t122 * t4 + t124 * t87 + t6 * t89;
t2 = m(2) * t77 + t4 * t95 - t8 * t92;
t1 = m(2) * t76 - t87 * t6 + (t4 * t92 + t124) * t89;
t3 = [-m(1) * g(1) - t1 * t86 + t2 * t88, t2, t4, t9, t11, t12, -t33 * mrSges(7,2) - t54 * t38 + t113; -m(1) * g(2) + t1 * t88 + t2 * t86, t1, t8, t14, t10, t13, -t102; -m(1) * g(3) + t110, t110, t6, t127, t98, t100, -t67 * mrSges(7,1) + t34 * mrSges(7,2) + t55 * t38 - t82 * t48 + t125;];
f_new  = t3;
