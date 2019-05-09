% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:55:50
% EndTime: 2019-05-05 11:56:04
% DurationCPUTime: 7.91s
% Computational Cost: add. (151929->184), mult. (312506->258), div. (0->0), fcn. (253000->16), ass. (0->109)
t102 = cos(pkin(7));
t114 = qJD(2) ^ 2;
t108 = sin(qJ(2));
t113 = cos(qJ(2));
t100 = sin(pkin(6));
t130 = t100 * t113;
t103 = cos(pkin(6));
t101 = cos(pkin(13));
t98 = sin(pkin(13));
t89 = g(1) * t98 - g(2) * t101;
t133 = t103 * t89;
t90 = -g(1) * t101 - g(2) * t98;
t97 = -g(3) + qJDD(1);
t121 = -t108 * t90 + t113 * t133 + t130 * t97;
t99 = sin(pkin(7));
t138 = pkin(9) * t99;
t57 = qJDD(2) * pkin(2) + t114 * t138 + t121;
t74 = -t100 * t89 + t103 * t97;
t140 = t102 * t57 + t74 * t99;
t107 = sin(qJ(3));
t112 = cos(qJ(3));
t129 = qJD(2) * qJD(3);
t82 = (-qJDD(2) * t112 + t107 * t129) * t99;
t131 = t100 * t108;
t127 = t108 * t133 + t113 * t90 + t131 * t97;
t58 = -pkin(2) * t114 + qJDD(2) * t138 + t127;
t139 = -t107 * t58 + t112 * t140;
t105 = sin(qJ(5));
t110 = cos(qJ(5));
t106 = sin(qJ(4));
t111 = cos(qJ(4));
t132 = qJD(2) * t99;
t125 = t112 * t132;
t128 = t107 * t140 + t112 * t58;
t80 = (-pkin(3) * t112 - pkin(10) * t107) * t132;
t95 = qJD(2) * t102 + qJD(3);
t93 = t95 ^ 2;
t94 = qJDD(2) * t102 + qJDD(3);
t35 = -pkin(3) * t93 + pkin(10) * t94 + t125 * t80 + t128;
t68 = t102 * t74;
t81 = (qJDD(2) * t107 + t112 * t129) * t99;
t38 = pkin(3) * t82 - pkin(10) * t81 + t68 + (-t57 + (pkin(3) * t107 - pkin(10) * t112) * t95 * qJD(2)) * t99;
t135 = t106 * t38 + t111 * t35;
t126 = t107 * t132;
t72 = -t106 * t126 + t111 * t95;
t73 = t106 * t95 + t111 * t126;
t60 = -pkin(4) * t72 - pkin(11) * t73;
t76 = qJDD(4) + t82;
t88 = qJD(4) - t125;
t86 = t88 ^ 2;
t24 = -pkin(4) * t86 + pkin(11) * t76 + t60 * t72 + t135;
t34 = -pkin(3) * t94 - pkin(10) * t93 + t126 * t80 - t139;
t52 = -qJD(4) * t73 - t106 * t81 + t111 * t94;
t53 = qJD(4) * t72 + t106 * t94 + t111 * t81;
t27 = (-t72 * t88 - t53) * pkin(11) + (t73 * t88 - t52) * pkin(4) + t34;
t136 = t105 * t27 + t110 * t24;
t104 = sin(qJ(6));
t109 = cos(qJ(6));
t123 = -t105 * t24 + t110 * t27;
t62 = -t105 * t73 + t110 * t88;
t40 = qJD(5) * t62 + t105 * t76 + t110 * t53;
t51 = qJDD(5) - t52;
t63 = t105 * t88 + t110 * t73;
t71 = qJD(5) - t72;
t18 = (t62 * t71 - t40) * pkin(12) + (t62 * t63 + t51) * pkin(5) + t123;
t39 = -qJD(5) * t63 - t105 * t53 + t110 * t76;
t49 = pkin(5) * t71 - pkin(12) * t63;
t61 = t62 ^ 2;
t19 = -pkin(5) * t61 + pkin(12) * t39 - t49 * t71 + t136;
t43 = -t104 * t63 + t109 * t62;
t30 = qJD(6) * t43 + t104 * t39 + t109 * t40;
t44 = t104 * t62 + t109 * t63;
t36 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t69 = qJD(6) + t71;
t41 = -mrSges(7,2) * t69 + mrSges(7,3) * t43;
t50 = qJDD(6) + t51;
t16 = m(7) * (-t104 * t19 + t109 * t18) - t30 * mrSges(7,3) + t50 * mrSges(7,1) - t44 * t36 + t69 * t41;
t29 = -qJD(6) * t44 - t104 * t40 + t109 * t39;
t42 = mrSges(7,1) * t69 - mrSges(7,3) * t44;
t17 = m(7) * (t104 * t18 + t109 * t19) + t29 * mrSges(7,3) - t50 * mrSges(7,2) + t43 * t36 - t69 * t42;
t45 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t47 = -mrSges(6,2) * t71 + mrSges(6,3) * t62;
t13 = m(6) * t123 + mrSges(6,1) * t51 - mrSges(6,3) * t40 + t104 * t17 + t109 * t16 - t45 * t63 + t47 * t71;
t48 = mrSges(6,1) * t71 - mrSges(6,3) * t63;
t14 = m(6) * t136 - mrSges(6,2) * t51 + mrSges(6,3) * t39 - t104 * t16 + t109 * t17 + t45 * t62 - t48 * t71;
t59 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t65 = mrSges(5,1) * t88 - mrSges(5,3) * t73;
t12 = m(5) * t135 - mrSges(5,2) * t76 + mrSges(5,3) * t52 - t105 * t13 + t110 * t14 + t59 * t72 - t65 * t88;
t122 = -t106 * t35 + t111 * t38;
t23 = -pkin(4) * t76 - pkin(11) * t86 + t60 * t73 - t122;
t117 = t29 * mrSges(7,1) + t43 * t41 - m(7) * (-pkin(5) * t39 - pkin(12) * t61 + t49 * t63 + t23) - t30 * mrSges(7,2) - t44 * t42;
t115 = m(6) * t23 - mrSges(6,1) * t39 + mrSges(6,2) * t40 - t47 * t62 + t48 * t63 - t117;
t64 = -mrSges(5,2) * t88 + mrSges(5,3) * t72;
t15 = m(5) * t122 + mrSges(5,1) * t76 - mrSges(5,3) * t53 - t59 * t73 + t64 * t88 - t115;
t77 = mrSges(4,1) * t95 - mrSges(4,3) * t126;
t78 = -mrSges(4,2) * t95 + mrSges(4,3) * t125;
t10 = m(4) * (-t57 * t99 + t68) + t81 * mrSges(4,2) + t82 * mrSges(4,1) + t106 * t12 + t111 * t15 + (t107 * t77 - t112 * t78) * t132;
t116 = m(5) * t34 - t52 * mrSges(5,1) + t53 * mrSges(5,2) + t105 * t14 + t110 * t13 - t72 * t64 + t73 * t65;
t79 = (-mrSges(4,1) * t112 + mrSges(4,2) * t107) * t132;
t11 = m(4) * t139 + t94 * mrSges(4,1) - t81 * mrSges(4,3) - t126 * t79 + t95 * t78 - t116;
t9 = m(4) * t128 - mrSges(4,2) * t94 - mrSges(4,3) * t82 - t106 * t15 + t111 * t12 + t125 * t79 - t77 * t95;
t119 = t107 * t9 + t11 * t112;
t4 = m(3) * t121 + qJDD(2) * mrSges(3,1) - t114 * mrSges(3,2) - t99 * t10 + t102 * t119;
t6 = m(3) * t74 + t10 * t102 + t119 * t99;
t8 = m(3) * t127 - mrSges(3,1) * t114 - qJDD(2) * mrSges(3,2) - t107 * t11 + t112 * t9;
t124 = m(2) * t97 + t103 * t6 + t130 * t4 + t131 * t8;
t2 = m(2) * t90 - t108 * t4 + t113 * t8;
t1 = m(2) * t89 - t100 * t6 + (t108 * t8 + t113 * t4) * t103;
t3 = [-m(1) * g(1) - t1 * t98 + t101 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t1 * t101 + t2 * t98, t1, t4, t11, t15, t13, t16; -m(1) * g(3) + t124, t124, t6, t10, t116, t115, -t117;];
f_new  = t3;
