% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:23:36
% EndTime: 2019-05-05 10:23:46
% DurationCPUTime: 5.48s
% Computational Cost: add. (79003->179), mult. (164044->243), div. (0->0), fcn. (123247->14), ass. (0->101)
t101 = sin(qJ(2));
t106 = cos(qJ(2));
t93 = sin(pkin(12));
t95 = cos(pkin(12));
t79 = t93 * g(1) - t95 * g(2);
t96 = cos(pkin(6));
t131 = t79 * t96;
t80 = -t95 * g(1) - t93 * g(2);
t92 = -g(3) + qJDD(1);
t94 = sin(pkin(6));
t133 = -t101 * t80 + (t92 * t94 + t131) * t106;
t100 = sin(qJ(3));
t105 = cos(qJ(3));
t107 = qJD(2) ^ 2;
t113 = -qJDD(2) * pkin(2) - t133;
t125 = qJD(2) * t100;
t123 = qJD(2) * qJD(3);
t78 = t105 * qJDD(2) - t100 * t123;
t84 = qJD(3) * pkin(3) - pkin(9) * t125;
t91 = t105 ^ 2;
t110 = -t78 * pkin(3) + t84 * t125 + (-pkin(9) * t91 - pkin(8)) * t107 + t113;
t102 = cos(qJ(6));
t104 = cos(qJ(4));
t99 = sin(qJ(4));
t72 = (t100 * t104 + t105 * t99) * qJD(2);
t119 = t105 * t123;
t77 = t100 * qJDD(2) + t119;
t51 = -t72 * qJD(4) + t104 * t78 - t99 * t77;
t90 = qJD(3) + qJD(4);
t67 = t90 * pkin(4) - t72 * pkin(10);
t71 = (-t100 * t99 + t104 * t105) * qJD(2);
t70 = t71 ^ 2;
t108 = -t51 * pkin(4) - t70 * pkin(10) + t72 * t67 + t110;
t103 = cos(qJ(5));
t127 = t101 * t94;
t122 = t101 * t131 + t106 * t80 + t92 * t127;
t56 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t122;
t68 = -t94 * t79 + t96 * t92;
t118 = -t100 * t56 + t105 * t68;
t37 = (-t77 + t119) * pkin(9) + (t100 * t105 * t107 + qJDD(3)) * pkin(3) + t118;
t128 = t100 * t68 + t105 * t56;
t39 = -t91 * t107 * pkin(3) + t78 * pkin(9) - qJD(3) * t84 + t128;
t120 = t104 * t37 - t99 * t39;
t52 = t71 * qJD(4) + t104 * t77 + t99 * t78;
t89 = qJDD(3) + qJDD(4);
t23 = (t71 * t90 - t52) * pkin(10) + (t71 * t72 + t89) * pkin(4) + t120;
t129 = t104 * t39 + t99 * t37;
t25 = -t70 * pkin(4) + t51 * pkin(10) - t90 * t67 + t129;
t98 = sin(qJ(5));
t130 = t103 * t25 + t98 * t23;
t60 = t103 * t71 - t98 * t72;
t61 = t103 * t72 + t98 * t71;
t44 = -t60 * pkin(5) - t61 * pkin(11);
t87 = qJD(5) + t90;
t85 = t87 ^ 2;
t86 = qJDD(5) + t89;
t20 = -t85 * pkin(5) + t86 * pkin(11) + t60 * t44 + t130;
t32 = -t61 * qJD(5) + t103 * t51 - t98 * t52;
t33 = t60 * qJD(5) + t103 * t52 + t98 * t51;
t21 = (-t60 * t87 - t33) * pkin(11) + (t61 * t87 - t32) * pkin(5) + t108;
t97 = sin(qJ(6));
t47 = t102 * t87 - t97 * t61;
t27 = t47 * qJD(6) + t102 * t33 + t97 * t86;
t31 = qJDD(6) - t32;
t48 = t102 * t61 + t97 * t87;
t38 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t57 = qJD(6) - t60;
t40 = -t57 * mrSges(7,2) + t47 * mrSges(7,3);
t17 = m(7) * (t102 * t21 - t97 * t20) - t27 * mrSges(7,3) + t31 * mrSges(7,1) - t48 * t38 + t57 * t40;
t26 = -t48 * qJD(6) + t102 * t86 - t97 * t33;
t41 = t57 * mrSges(7,1) - t48 * mrSges(7,3);
t18 = m(7) * (t102 * t20 + t97 * t21) + t26 * mrSges(7,3) - t31 * mrSges(7,2) + t47 * t38 - t57 * t41;
t53 = -t87 * mrSges(6,2) + t60 * mrSges(6,3);
t54 = t87 * mrSges(6,1) - t61 * mrSges(6,3);
t114 = -m(6) * t108 + t32 * mrSges(6,1) - t33 * mrSges(6,2) - t102 * t17 - t97 * t18 + t60 * t53 - t61 * t54;
t65 = -t90 * mrSges(5,2) + t71 * mrSges(5,3);
t66 = t90 * mrSges(5,1) - t72 * mrSges(5,3);
t111 = -m(5) * t110 + t51 * mrSges(5,1) - t52 * mrSges(5,2) + t71 * t65 - t72 * t66 + t114;
t81 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t125;
t124 = qJD(2) * t105;
t82 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t124;
t132 = (t100 * t81 - t105 * t82) * qJD(2) + m(4) * (-t107 * pkin(8) + t113) - t78 * mrSges(4,1) + t77 * mrSges(4,2) - t111;
t12 = m(3) * t133 + qJDD(2) * mrSges(3,1) - t107 * mrSges(3,2) - t132;
t126 = t106 * t12;
t43 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t13 = m(6) * t130 - t86 * mrSges(6,2) + t32 * mrSges(6,3) + t102 * t18 - t97 * t17 + t60 * t43 - t87 * t54;
t116 = t103 * t23 - t98 * t25;
t112 = m(7) * (-t86 * pkin(5) - t85 * pkin(11) + t61 * t44 - t116) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t47 * t40 + t48 * t41;
t14 = m(6) * t116 + t86 * mrSges(6,1) - t33 * mrSges(6,3) - t61 * t43 + t87 * t53 - t112;
t62 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t10 = m(5) * t129 - t89 * mrSges(5,2) + t51 * mrSges(5,3) + t103 * t13 - t98 * t14 + t71 * t62 - t90 * t66;
t76 = (-mrSges(4,1) * t105 + mrSges(4,2) * t100) * qJD(2);
t9 = m(5) * t120 + t89 * mrSges(5,1) - t52 * mrSges(5,3) + t103 * t14 + t98 * t13 - t72 * t62 + t90 * t65;
t7 = m(4) * t118 + qJDD(3) * mrSges(4,1) - t77 * mrSges(4,3) + qJD(3) * t82 + t99 * t10 + t104 * t9 - t76 * t125;
t8 = m(4) * t128 - qJDD(3) * mrSges(4,2) + t78 * mrSges(4,3) - qJD(3) * t81 + t104 * t10 + t76 * t124 - t99 * t9;
t4 = m(3) * t122 - t107 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t100 * t7 + t105 * t8;
t6 = m(3) * t68 + t100 * t8 + t105 * t7;
t121 = m(2) * t92 + t94 * t126 + t4 * t127 + t96 * t6;
t2 = m(2) * t80 - t101 * t12 + t106 * t4;
t1 = m(2) * t79 - t94 * t6 + (t101 * t4 + t126) * t96;
t3 = [-m(1) * g(1) - t93 * t1 + t95 * t2, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t95 * t1 + t93 * t2, t1, t12, t7, t9, t14, t17; -m(1) * g(3) + t121, t121, t6, t132, -t111, -t114, t112;];
f_new  = t3;
