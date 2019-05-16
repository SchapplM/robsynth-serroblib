% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:32:57
% EndTime: 2019-05-05 18:33:06
% DurationCPUTime: 4.16s
% Computational Cost: add. (59206->178), mult. (123709->235), div. (0->0), fcn. (83663->12), ass. (0->95)
t102 = qJD(1) ^ 2;
t100 = cos(qJ(1));
t96 = sin(qJ(1));
t114 = g(1) * t96 - g(2) * t100;
t74 = qJDD(1) * pkin(1) + t114;
t109 = -g(1) * t100 - g(2) * t96;
t77 = -pkin(1) * t102 + t109;
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t111 = t92 * t74 - t77 * t90;
t51 = -qJDD(1) * pkin(2) - t102 * pkin(7) - t111;
t118 = qJD(1) * qJD(3);
t99 = cos(qJ(3));
t115 = t99 * t118;
t95 = sin(qJ(3));
t78 = qJDD(1) * t95 + t115;
t85 = t95 * t118;
t79 = qJDD(1) * t99 - t85;
t39 = (-t78 - t115) * qJ(4) + (-t79 + t85) * pkin(3) + t51;
t101 = qJD(3) ^ 2;
t119 = t99 * qJD(1);
t121 = t74 * t90 + t77 * t92;
t52 = -pkin(2) * t102 + qJDD(1) * pkin(7) + t121;
t88 = -g(3) + qJDD(2);
t122 = t52 * t99 + t88 * t95;
t75 = (-pkin(3) * t99 - qJ(4) * t95) * qJD(1);
t45 = -pkin(3) * t101 + qJDD(3) * qJ(4) + t119 * t75 + t122;
t120 = qJD(1) * t95;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t72 = qJD(3) * t89 + t120 * t91;
t110 = -0.2e1 * qJD(4) * t72 + t39 * t91 - t89 * t45;
t61 = qJDD(3) * t89 + t78 * t91;
t71 = qJD(3) * t91 - t120 * t89;
t22 = (-t119 * t71 - t61) * pkin(8) + (t71 * t72 - t79) * pkin(4) + t110;
t116 = 0.2e1 * qJD(4) * t71 + t39 * t89 + t45 * t91;
t60 = qJDD(3) * t91 - t78 * t89;
t62 = -pkin(4) * t119 - pkin(8) * t72;
t70 = t71 ^ 2;
t26 = -pkin(4) * t70 + pkin(8) * t60 + t119 * t62 + t116;
t94 = sin(qJ(5));
t98 = cos(qJ(5));
t113 = t22 * t98 - t94 * t26;
t55 = t71 * t98 - t72 * t94;
t35 = qJD(5) * t55 + t60 * t94 + t61 * t98;
t56 = t71 * t94 + t72 * t98;
t73 = qJDD(5) - t79;
t83 = qJD(5) - t119;
t15 = (t55 * t83 - t35) * pkin(9) + (t55 * t56 + t73) * pkin(5) + t113;
t123 = t22 * t94 + t26 * t98;
t34 = -qJD(5) * t56 + t60 * t98 - t61 * t94;
t48 = pkin(5) * t83 - pkin(9) * t56;
t54 = t55 ^ 2;
t16 = -pkin(5) * t54 + pkin(9) * t34 - t48 * t83 + t123;
t93 = sin(qJ(6));
t97 = cos(qJ(6));
t41 = t55 * t97 - t56 * t93;
t24 = qJD(6) * t41 + t34 * t93 + t35 * t97;
t42 = t55 * t93 + t56 * t97;
t28 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t82 = qJD(6) + t83;
t31 = -mrSges(7,2) * t82 + mrSges(7,3) * t41;
t69 = qJDD(6) + t73;
t11 = m(7) * (t15 * t97 - t16 * t93) - t24 * mrSges(7,3) + t69 * mrSges(7,1) - t42 * t28 + t82 * t31;
t23 = -qJD(6) * t42 + t34 * t97 - t35 * t93;
t32 = mrSges(7,1) * t82 - mrSges(7,3) * t42;
t12 = m(7) * (t15 * t93 + t16 * t97) + t23 * mrSges(7,3) - t69 * mrSges(7,2) + t41 * t28 - t82 * t32;
t44 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t47 = mrSges(6,1) * t83 - mrSges(6,3) * t56;
t10 = m(6) * t123 - mrSges(6,2) * t73 + mrSges(6,3) * t34 - t11 * t93 + t12 * t97 + t44 * t55 - t47 * t83;
t57 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t58 = mrSges(5,2) * t119 + mrSges(5,3) * t71;
t46 = -mrSges(6,2) * t83 + mrSges(6,3) * t55;
t9 = m(6) * t113 + mrSges(6,1) * t73 - mrSges(6,3) * t35 + t11 * t97 + t12 * t93 - t44 * t56 + t46 * t83;
t7 = m(5) * t110 - mrSges(5,1) * t79 - mrSges(5,3) * t61 + t10 * t94 - t119 * t58 - t57 * t72 + t9 * t98;
t59 = -mrSges(5,1) * t119 - mrSges(5,3) * t72;
t8 = m(5) * t116 + mrSges(5,2) * t79 + mrSges(5,3) * t60 + t10 * t98 + t119 * t59 + t57 * t71 - t9 * t94;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t120;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t119;
t124 = m(4) * t51 - t79 * mrSges(4,1) + t78 * mrSges(4,2) + t91 * t7 + t89 * t8 + (t80 * t95 - t81 * t99) * qJD(1);
t112 = -t52 * t95 + t88 * t99;
t43 = -qJDD(3) * pkin(3) - qJ(4) * t101 + t120 * t75 + qJDD(4) - t112;
t104 = -pkin(4) * t60 - pkin(8) * t70 + t62 * t72 + t43;
t107 = t23 * mrSges(7,1) + t41 * t31 - m(7) * (-pkin(5) * t34 - pkin(9) * t54 + t48 * t56 + t104) - t24 * mrSges(7,2) - t42 * t32;
t105 = -m(6) * t104 + mrSges(6,1) * t34 - mrSges(6,2) * t35 + t46 * t55 - t47 * t56 + t107;
t103 = m(5) * t43 - mrSges(5,1) * t60 + mrSges(5,2) * t61 - t58 * t71 + t59 * t72 - t105;
t76 = (-mrSges(4,1) * t99 + mrSges(4,2) * t95) * qJD(1);
t14 = m(4) * t112 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t78 + qJD(3) * t81 - t120 * t76 - t103;
t6 = m(4) * t122 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t79 - qJD(3) * t80 + t119 * t76 - t7 * t89 + t8 * t91;
t117 = m(3) * t88 + t14 * t99 + t6 * t95;
t4 = m(3) * t111 + qJDD(1) * mrSges(3,1) - t102 * mrSges(3,2) - t124;
t3 = m(3) * t121 - mrSges(3,1) * t102 - qJDD(1) * mrSges(3,2) - t14 * t95 + t6 * t99;
t2 = m(2) * t109 - mrSges(2,1) * t102 - qJDD(1) * mrSges(2,2) + t3 * t92 - t4 * t90;
t1 = m(2) * t114 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t102 + t3 * t90 + t4 * t92;
t5 = [-m(1) * g(1) - t1 * t96 + t100 * t2, t2, t3, t6, t8, t10, t12; -m(1) * g(2) + t1 * t100 + t2 * t96, t1, t4, t14, t7, t9, t11; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t117, t124, t103, -t105, -t107;];
f_new  = t5;
