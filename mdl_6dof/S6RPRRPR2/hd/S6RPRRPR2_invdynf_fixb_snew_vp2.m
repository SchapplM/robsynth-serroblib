% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:04:00
% EndTime: 2019-05-05 22:04:10
% DurationCPUTime: 4.29s
% Computational Cost: add. (62779->178), mult. (126343->235), div. (0->0), fcn. (85444->12), ass. (0->95)
t102 = qJD(1) ^ 2;
t100 = cos(qJ(1));
t96 = sin(qJ(1));
t114 = t96 * g(1) - t100 * g(2);
t74 = qJDD(1) * pkin(1) + t114;
t109 = -t100 * g(1) - t96 * g(2);
t76 = -t102 * pkin(1) + t109;
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t111 = t92 * t74 - t90 * t76;
t51 = -qJDD(1) * pkin(2) - t102 * pkin(7) - t111;
t118 = qJD(1) * qJD(3);
t99 = cos(qJ(3));
t115 = t99 * t118;
t95 = sin(qJ(3));
t78 = t95 * qJDD(1) + t115;
t85 = t95 * t118;
t79 = t99 * qJDD(1) - t85;
t36 = (-t78 - t115) * pkin(8) + (-t79 + t85) * pkin(3) + t51;
t101 = qJD(3) ^ 2;
t119 = t99 * qJD(1);
t121 = t90 * t74 + t92 * t76;
t52 = -t102 * pkin(2) + qJDD(1) * pkin(7) + t121;
t88 = -g(3) + qJDD(2);
t122 = t99 * t52 + t95 * t88;
t77 = (-pkin(3) * t99 - pkin(8) * t95) * qJD(1);
t44 = -t101 * pkin(3) + qJDD(3) * pkin(8) + t77 * t119 + t122;
t94 = sin(qJ(4));
t98 = cos(qJ(4));
t113 = t98 * t36 - t94 * t44;
t120 = qJD(1) * t95;
t72 = t98 * qJD(3) - t94 * t120;
t57 = t72 * qJD(4) + t94 * qJDD(3) + t98 * t78;
t71 = qJDD(4) - t79;
t73 = t94 * qJD(3) + t98 * t120;
t83 = qJD(4) - t119;
t21 = (t72 * t83 - t57) * qJ(5) + (t72 * t73 + t71) * pkin(4) + t113;
t123 = t94 * t36 + t98 * t44;
t56 = -t73 * qJD(4) + t98 * qJDD(3) - t94 * t78;
t63 = t83 * pkin(4) - t73 * qJ(5);
t70 = t72 ^ 2;
t26 = -t70 * pkin(4) + t56 * qJ(5) - t83 * t63 + t123;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t60 = t89 * t72 + t91 * t73;
t110 = -0.2e1 * qJD(5) * t60 + t91 * t21 - t89 * t26;
t38 = t89 * t56 + t91 * t57;
t59 = t91 * t72 - t89 * t73;
t15 = (t59 * t83 - t38) * pkin(9) + (t59 * t60 + t71) * pkin(5) + t110;
t116 = 0.2e1 * qJD(5) * t59 + t89 * t21 + t91 * t26;
t37 = t91 * t56 - t89 * t57;
t48 = t83 * pkin(5) - t60 * pkin(9);
t58 = t59 ^ 2;
t16 = -t58 * pkin(5) + t37 * pkin(9) - t83 * t48 + t116;
t93 = sin(qJ(6));
t97 = cos(qJ(6));
t41 = t97 * t59 - t93 * t60;
t25 = t41 * qJD(6) + t93 * t37 + t97 * t38;
t42 = t93 * t59 + t97 * t60;
t30 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t82 = qJD(6) + t83;
t31 = -t82 * mrSges(7,2) + t41 * mrSges(7,3);
t69 = qJDD(6) + t71;
t11 = m(7) * (t97 * t15 - t93 * t16) - t25 * mrSges(7,3) + t69 * mrSges(7,1) - t42 * t30 + t82 * t31;
t24 = -t42 * qJD(6) + t97 * t37 - t93 * t38;
t32 = t82 * mrSges(7,1) - t42 * mrSges(7,3);
t12 = m(7) * (t93 * t15 + t97 * t16) + t24 * mrSges(7,3) - t69 * mrSges(7,2) + t41 * t30 - t82 * t32;
t45 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t47 = t83 * mrSges(6,1) - t60 * mrSges(6,3);
t10 = m(6) * t116 - t71 * mrSges(6,2) + t37 * mrSges(6,3) - t93 * t11 + t97 * t12 + t59 * t45 - t83 * t47;
t61 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t62 = -t83 * mrSges(5,2) + t72 * mrSges(5,3);
t46 = -t83 * mrSges(6,2) + t59 * mrSges(6,3);
t9 = m(6) * t110 + t71 * mrSges(6,1) - t38 * mrSges(6,3) + t97 * t11 + t93 * t12 - t60 * t45 + t83 * t46;
t7 = m(5) * t113 + t71 * mrSges(5,1) - t57 * mrSges(5,3) + t89 * t10 - t73 * t61 + t83 * t62 + t91 * t9;
t64 = t83 * mrSges(5,1) - t73 * mrSges(5,3);
t8 = m(5) * t123 - t71 * mrSges(5,2) + t56 * mrSges(5,3) + t91 * t10 + t72 * t61 - t83 * t64 - t89 * t9;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t120;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t119;
t124 = m(4) * t51 - t79 * mrSges(4,1) + t78 * mrSges(4,2) + t98 * t7 + t94 * t8 + (t95 * t80 - t99 * t81) * qJD(1);
t112 = -t95 * t52 + t99 * t88;
t43 = -qJDD(3) * pkin(3) - t101 * pkin(8) + t77 * t120 - t112;
t104 = -t56 * pkin(4) - t70 * qJ(5) + t73 * t63 + qJDD(5) + t43;
t107 = t24 * mrSges(7,1) + t41 * t31 - m(7) * (-t37 * pkin(5) - t58 * pkin(9) + t60 * t48 + t104) - t25 * mrSges(7,2) - t42 * t32;
t105 = -m(6) * t104 + t37 * mrSges(6,1) - t38 * mrSges(6,2) + t59 * t46 - t60 * t47 + t107;
t103 = m(5) * t43 - t56 * mrSges(5,1) + t57 * mrSges(5,2) - t72 * t62 + t73 * t64 - t105;
t75 = (-mrSges(4,1) * t99 + mrSges(4,2) * t95) * qJD(1);
t14 = m(4) * t112 + qJDD(3) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(3) * t81 - t75 * t120 - t103;
t6 = m(4) * t122 - qJDD(3) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(3) * t80 + t75 * t119 - t94 * t7 + t98 * t8;
t117 = m(3) * t88 + t99 * t14 + t95 * t6;
t4 = m(3) * t111 + qJDD(1) * mrSges(3,1) - t102 * mrSges(3,2) - t124;
t3 = m(3) * t121 - t102 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t95 * t14 + t99 * t6;
t2 = m(2) * t109 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t92 * t3 - t90 * t4;
t1 = m(2) * t114 + qJDD(1) * mrSges(2,1) - t102 * mrSges(2,2) + t90 * t3 + t92 * t4;
t5 = [-m(1) * g(1) - t96 * t1 + t100 * t2, t2, t3, t6, t8, t10, t12; -m(1) * g(2) + t100 * t1 + t96 * t2, t1, t4, t14, t7, t9, t11; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t117, t124, t103, -t105, -t107;];
f_new  = t5;
