% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 11:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:46:36
% EndTime: 2019-05-06 11:46:44
% DurationCPUTime: 3.39s
% Computational Cost: add. (41166->206), mult. (90670->259), div. (0->0), fcn. (57396->10), ass. (0->101)
t150 = -2 * qJD(3);
t106 = sin(qJ(2));
t110 = cos(qJ(2));
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t132 = t107 * g(1) - t111 * g(2);
t124 = -qJDD(1) * pkin(1) - t132;
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t135 = qJD(1) * qJD(2);
t129 = t110 * t135;
t130 = t106 * t135;
t83 = t106 * qJDD(1) + t129;
t96 = t106 * qJD(1);
t119 = pkin(2) * t130 + t96 * t150 + (-t83 - t129) * qJ(3) + t124;
t136 = qJD(1) * t110;
t113 = qJD(1) ^ 2;
t144 = t113 * pkin(7);
t105 = sin(qJ(5));
t109 = cos(qJ(5));
t101 = t110 ^ 2;
t84 = t110 * qJDD(1) - t130;
t88 = pkin(3) * t96 - qJD(2) * qJ(4);
t31 = -t88 * t96 + (-pkin(2) - qJ(4)) * t84 + (-pkin(3) * t101 - pkin(7)) * t113 + t119;
t112 = qJD(2) ^ 2;
t126 = -t111 * g(1) - t107 * g(2);
t72 = -t113 * pkin(1) + qJDD(1) * pkin(7) + t126;
t140 = -t110 * g(3) - t106 * t72;
t80 = (-pkin(2) * t110 - qJ(3) * t106) * qJD(1);
t49 = -qJDD(2) * pkin(2) - t112 * qJ(3) + t80 * t96 + qJDD(3) - t140;
t41 = (-t106 * t110 * t113 - qJDD(2)) * qJ(4) + (t83 - t129) * pkin(3) + t49;
t77 = t103 * qJD(2) - t102 * t136;
t127 = -0.2e1 * qJD(4) * t77 - t102 * t31 + t103 * t41;
t76 = -t102 * qJD(2) - t103 * t136;
t58 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t60 = -mrSges(5,2) * t96 + t76 * mrSges(5,3);
t63 = t103 * qJDD(2) - t102 * t84;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t20 = (t76 * t96 - t63) * pkin(8) + (t76 * t77 + t83) * pkin(4) + t127;
t133 = 0.2e1 * qJD(4) * t76 + t102 * t41 + t103 * t31;
t62 = -t102 * qJDD(2) - t103 * t84;
t64 = pkin(4) * t96 - t77 * pkin(8);
t75 = t76 ^ 2;
t22 = -t75 * pkin(4) + t62 * pkin(8) - t64 * t96 + t133;
t128 = -t105 * t22 + t109 * t20;
t56 = -t105 * t77 + t109 * t76;
t36 = t56 * qJD(5) + t105 * t62 + t109 * t63;
t57 = t105 * t76 + t109 * t77;
t79 = qJDD(5) + t83;
t93 = t96 + qJD(5);
t14 = (t56 * t93 - t36) * pkin(9) + (t56 * t57 + t79) * pkin(5) + t128;
t141 = t105 * t20 + t109 * t22;
t35 = -t57 * qJD(5) - t105 * t63 + t109 * t62;
t52 = t93 * pkin(5) - t57 * pkin(9);
t55 = t56 ^ 2;
t15 = -t55 * pkin(5) + t35 * pkin(9) - t93 * t52 + t141;
t43 = -t104 * t57 + t108 * t56;
t25 = t43 * qJD(6) + t104 * t35 + t108 * t36;
t44 = t104 * t56 + t108 * t57;
t29 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t91 = qJD(6) + t93;
t32 = -t91 * mrSges(7,2) + t43 * mrSges(7,3);
t73 = qJDD(6) + t79;
t10 = m(7) * (-t104 * t15 + t108 * t14) - t25 * mrSges(7,3) + t73 * mrSges(7,1) - t44 * t29 + t91 * t32;
t24 = -t44 * qJD(6) - t104 * t36 + t108 * t35;
t33 = t91 * mrSges(7,1) - t44 * mrSges(7,3);
t11 = m(7) * (t104 * t14 + t108 * t15) + t24 * mrSges(7,3) - t73 * mrSges(7,2) + t43 * t29 - t91 * t33;
t45 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t50 = -t93 * mrSges(6,2) + t56 * mrSges(6,3);
t8 = m(6) * t128 + t79 * mrSges(6,1) - t36 * mrSges(6,3) + t108 * t10 + t104 * t11 - t57 * t45 + t93 * t50;
t51 = t93 * mrSges(6,1) - t57 * mrSges(6,3);
t9 = m(6) * t141 - t79 * mrSges(6,2) + t35 * mrSges(6,3) - t104 * t10 + t108 * t11 + t56 * t45 - t93 * t51;
t6 = m(5) * t127 + t83 * mrSges(5,1) - t63 * mrSges(5,3) + t105 * t9 + t109 * t8 - t77 * t58 + t60 * t96;
t61 = mrSges(5,1) * t96 - t77 * mrSges(5,3);
t7 = m(5) * t133 - t83 * mrSges(5,2) + t62 * mrSges(5,3) - t105 * t8 + t109 * t9 + t76 * t58 - t61 * t96;
t89 = -mrSges(4,1) * t136 - qJD(2) * mrSges(4,3);
t125 = t102 * t6 - t103 * t7 - m(4) * (-t84 * pkin(2) + t119 - t144) - t89 * t136 + t83 * mrSges(4,3);
t90 = mrSges(4,1) * t96 + qJD(2) * mrSges(4,2);
t138 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t96 - t90;
t143 = mrSges(3,1) - mrSges(4,2);
t87 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t149 = (t138 * t106 - t110 * t87) * qJD(1) - t143 * t84 + m(3) * (t124 - t144) + t83 * mrSges(3,2) - t125;
t131 = -t106 * g(3) + t110 * t72;
t148 = t112 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t150 - t80 * t136 - t131;
t117 = -t101 * t113 * qJ(4) + t84 * pkin(3) + qJD(2) * t88 + qJDD(4) - t148;
t116 = -t62 * pkin(4) - t75 * pkin(8) + t77 * t64 + t117;
t122 = t24 * mrSges(7,1) + t43 * t32 - m(7) * (-t35 * pkin(5) - t55 * pkin(9) + t57 * t52 + t116) - t25 * mrSges(7,2) - t44 * t33;
t118 = -m(6) * t116 + t35 * mrSges(6,1) - t36 * mrSges(6,2) + t56 * t50 - t57 * t51 + t122;
t115 = -m(5) * t117 + t62 * mrSges(5,1) - t63 * mrSges(5,2) + t76 * t60 - t77 * t61 + t118;
t114 = m(4) * t148 + t115;
t81 = (mrSges(4,2) * t110 - mrSges(4,3) * t106) * qJD(1);
t139 = t81 + (-mrSges(3,1) * t110 + mrSges(3,2) * t106) * qJD(1);
t142 = mrSges(3,3) + mrSges(4,1);
t13 = m(3) * t131 - t114 + t139 * t136 - t138 * qJD(2) + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) + t142 * t84;
t123 = -m(4) * t49 - t102 * t7 - t103 * t6;
t4 = m(3) * t140 - t142 * t83 + t143 * qJDD(2) + (t87 - t89) * qJD(2) - t139 * t96 + t123;
t145 = t106 * t13 + t110 * t4;
t2 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t149;
t1 = m(2) * t126 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t110 * t13;
t3 = [-m(1) * g(1) + t111 * t1 - t107 * t2, t1, t13, t84 * mrSges(4,2) - t90 * t96 - t125, t7, t9, t11; -m(1) * g(2) + t107 * t1 + t111 * t2, t2, t4, -t84 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t90 - t81 * t136 + t114, t6, t8, t10; (-m(1) - m(2)) * g(3) + t145, -m(2) * g(3) + t145, t149, t83 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t89 + t81 * t96 - t123, -t115, -t118, -t122;];
f_new  = t3;
