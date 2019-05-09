% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:56:50
% EndTime: 2019-05-06 02:57:00
% DurationCPUTime: 4.38s
% Computational Cost: add. (65940->179), mult. (129781->233), div. (0->0), fcn. (88868->12), ass. (0->97)
t102 = qJD(1) ^ 2;
t100 = cos(qJ(1));
t95 = sin(qJ(1));
t114 = g(1) * t95 - g(2) * t100;
t73 = qJDD(1) * pkin(1) + t114;
t109 = -g(1) * t100 - g(2) * t95;
t75 = -pkin(1) * t102 + t109;
t89 = sin(pkin(11));
t90 = cos(pkin(11));
t110 = t90 * t73 - t75 * t89;
t51 = -qJDD(1) * pkin(2) - t102 * pkin(7) - t110;
t117 = qJD(1) * qJD(3);
t99 = cos(qJ(3));
t115 = t99 * t117;
t94 = sin(qJ(3));
t77 = qJDD(1) * t94 + t115;
t85 = t94 * t117;
t78 = qJDD(1) * t99 - t85;
t38 = (-t77 - t115) * pkin(8) + (-t78 + t85) * pkin(3) + t51;
t101 = qJD(3) ^ 2;
t118 = t99 * qJD(1);
t120 = t73 * t89 + t75 * t90;
t52 = -pkin(2) * t102 + qJDD(1) * pkin(7) + t120;
t88 = -g(3) + qJDD(2);
t121 = t52 * t99 + t88 * t94;
t76 = (-pkin(3) * t99 - pkin(8) * t94) * qJD(1);
t44 = -pkin(3) * t101 + qJDD(3) * pkin(8) + t118 * t76 + t121;
t93 = sin(qJ(4));
t98 = cos(qJ(4));
t112 = t38 * t98 - t93 * t44;
t119 = qJD(1) * t94;
t71 = qJD(3) * t98 - t119 * t93;
t55 = qJD(4) * t71 + qJDD(3) * t93 + t77 * t98;
t70 = qJDD(4) - t78;
t72 = qJD(3) * t93 + t119 * t98;
t83 = qJD(4) - t118;
t24 = (t71 * t83 - t55) * pkin(9) + (t71 * t72 + t70) * pkin(4) + t112;
t122 = t38 * t93 + t44 * t98;
t54 = -qJD(4) * t72 + qJDD(3) * t98 - t77 * t93;
t62 = pkin(4) * t83 - pkin(9) * t72;
t69 = t71 ^ 2;
t26 = -pkin(4) * t69 + pkin(9) * t54 - t62 * t83 + t122;
t92 = sin(qJ(5));
t97 = cos(qJ(5));
t123 = t24 * t92 + t26 * t97;
t113 = t24 * t97 - t92 * t26;
t57 = t71 * t97 - t72 * t92;
t33 = qJD(5) * t57 + t54 * t92 + t55 * t97;
t58 = t71 * t92 + t72 * t97;
t68 = qJDD(5) + t70;
t82 = qJD(5) + t83;
t15 = (t57 * t82 - t33) * pkin(10) + (t57 * t58 + t68) * pkin(5) + t113;
t32 = -qJD(5) * t58 + t54 * t97 - t55 * t92;
t48 = pkin(5) * t82 - pkin(10) * t58;
t56 = t57 ^ 2;
t16 = -pkin(5) * t56 + pkin(10) * t32 - t48 * t82 + t123;
t91 = sin(qJ(6));
t96 = cos(qJ(6));
t41 = t57 * t96 - t58 * t91;
t21 = qJD(6) * t41 + t32 * t91 + t33 * t96;
t42 = t57 * t91 + t58 * t96;
t30 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t79 = qJD(6) + t82;
t34 = -mrSges(7,2) * t79 + mrSges(7,3) * t41;
t64 = qJDD(6) + t68;
t13 = m(7) * (t15 * t96 - t16 * t91) - t21 * mrSges(7,3) + t64 * mrSges(7,1) - t42 * t30 + t79 * t34;
t20 = -qJD(6) * t42 + t32 * t96 - t33 * t91;
t35 = mrSges(7,1) * t79 - mrSges(7,3) * t42;
t14 = m(7) * (t15 * t91 + t16 * t96) + t20 * mrSges(7,3) - t64 * mrSges(7,2) + t41 * t30 - t79 * t35;
t45 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t47 = mrSges(6,1) * t82 - mrSges(6,3) * t58;
t10 = m(6) * t123 - mrSges(6,2) * t68 + mrSges(6,3) * t32 - t13 * t91 + t14 * t96 + t45 * t57 - t47 * t82;
t59 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t60 = -mrSges(5,2) * t83 + mrSges(5,3) * t71;
t46 = -mrSges(6,2) * t82 + mrSges(6,3) * t57;
t9 = m(6) * t113 + mrSges(6,1) * t68 - mrSges(6,3) * t33 + t13 * t96 + t14 * t91 - t45 * t58 + t46 * t82;
t7 = m(5) * t112 + mrSges(5,1) * t70 - mrSges(5,3) * t55 + t10 * t92 - t59 * t72 + t60 * t83 + t9 * t97;
t61 = mrSges(5,1) * t83 - mrSges(5,3) * t72;
t8 = m(5) * t122 - mrSges(5,2) * t70 + mrSges(5,3) * t54 + t10 * t97 + t59 * t71 - t61 * t83 - t9 * t92;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t118;
t124 = m(4) * t51 - t78 * mrSges(4,1) + t77 * mrSges(4,2) + t98 * t7 + t93 * t8 + (t80 * t94 - t81 * t99) * qJD(1);
t111 = -t52 * t94 + t88 * t99;
t43 = -qJDD(3) * pkin(3) - pkin(8) * t101 + t119 * t76 - t111;
t105 = -pkin(4) * t54 - pkin(9) * t69 + t62 * t72 + t43;
t107 = t20 * mrSges(7,1) + t41 * t34 - m(7) * (-pkin(5) * t32 - pkin(10) * t56 + t48 * t58 + t105) - t21 * mrSges(7,2) - t42 * t35;
t104 = -m(6) * t105 + mrSges(6,1) * t32 - mrSges(6,2) * t33 + t46 * t57 - t47 * t58 + t107;
t103 = m(5) * t43 - mrSges(5,1) * t54 + mrSges(5,2) * t55 - t60 * t71 + t61 * t72 - t104;
t74 = (-mrSges(4,1) * t99 + mrSges(4,2) * t94) * qJD(1);
t12 = m(4) * t111 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t77 + qJD(3) * t81 - t119 * t74 - t103;
t6 = m(4) * t121 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t78 - qJD(3) * t80 + t118 * t74 - t7 * t93 + t8 * t98;
t116 = m(3) * t88 + t12 * t99 + t6 * t94;
t4 = m(3) * t110 + qJDD(1) * mrSges(3,1) - t102 * mrSges(3,2) - t124;
t3 = m(3) * t120 - t102 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t94 * t12 + t99 * t6;
t2 = m(2) * t109 - mrSges(2,1) * t102 - qJDD(1) * mrSges(2,2) + t3 * t90 - t4 * t89;
t1 = m(2) * t114 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t102 + t3 * t89 + t4 * t90;
t5 = [-m(1) * g(1) - t1 * t95 + t100 * t2, t2, t3, t6, t8, t10, t14; -m(1) * g(2) + t1 * t100 + t2 * t95, t1, t4, t12, t7, t9, t13; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t116, t124, t103, -t104, -t107;];
f_new  = t5;
