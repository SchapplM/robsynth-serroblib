% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 00:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:17:10
% EndTime: 2019-05-07 00:17:20
% DurationCPUTime: 3.49s
% Computational Cost: add. (45646->207), mult. (94732->257), div. (0->0), fcn. (60873->10), ass. (0->103)
t150 = -2 * qJD(3);
t105 = sin(qJ(2));
t110 = cos(qJ(2));
t106 = sin(qJ(1));
t111 = cos(qJ(1));
t132 = t106 * g(1) - t111 * g(2);
t124 = -qJDD(1) * pkin(1) - t132;
t104 = sin(qJ(4));
t109 = cos(qJ(4));
t134 = qJD(1) * qJD(2);
t129 = t110 * t134;
t130 = t105 * t134;
t82 = t105 * qJDD(1) + t129;
t96 = t105 * qJD(1);
t119 = pkin(2) * t130 + t96 * t150 + (-t82 - t129) * qJ(3) + t124;
t135 = qJD(1) * t110;
t113 = qJD(1) ^ 2;
t144 = t113 * pkin(7);
t103 = sin(qJ(5));
t108 = cos(qJ(5));
t101 = t110 ^ 2;
t83 = t110 * qJDD(1) - t130;
t90 = pkin(3) * t96 - qJD(2) * pkin(8);
t34 = -t90 * t96 + (-pkin(2) - pkin(8)) * t83 + (-pkin(3) * t101 - pkin(7)) * t113 + t119;
t112 = qJD(2) ^ 2;
t126 = -t111 * g(1) - t106 * g(2);
t71 = -t113 * pkin(1) + qJDD(1) * pkin(7) + t126;
t139 = -t110 * g(3) - t105 * t71;
t79 = (-pkin(2) * t110 - qJ(3) * t105) * qJD(1);
t49 = -qJDD(2) * pkin(2) - t112 * qJ(3) + t79 * t96 + qJDD(3) - t139;
t41 = (-t105 * t110 * t113 - qJDD(2)) * pkin(8) + (t82 - t129) * pkin(3) + t49;
t127 = -t104 * t34 + t109 * t41;
t77 = -t104 * qJD(2) - t109 * t135;
t57 = t77 * qJD(4) + t109 * qJDD(2) - t104 * t83;
t78 = t109 * qJD(2) - t104 * t135;
t61 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t93 = t96 + qJD(4);
t62 = -t93 * mrSges(5,2) + t77 * mrSges(5,3);
t76 = qJDD(4) + t82;
t102 = sin(qJ(6));
t107 = cos(qJ(6));
t23 = (t77 * t93 - t57) * pkin(9) + (t77 * t78 + t76) * pkin(4) + t127;
t140 = t104 * t41 + t109 * t34;
t56 = -t78 * qJD(4) - t104 * qJDD(2) - t109 * t83;
t64 = t93 * pkin(4) - t78 * pkin(9);
t75 = t77 ^ 2;
t25 = -t75 * pkin(4) + t56 * pkin(9) - t93 * t64 + t140;
t128 = -t103 * t25 + t108 * t23;
t59 = -t103 * t78 + t108 * t77;
t32 = t59 * qJD(5) + t103 * t56 + t108 * t57;
t60 = t103 * t77 + t108 * t78;
t72 = qJDD(5) + t76;
t91 = qJD(5) + t93;
t14 = (t59 * t91 - t32) * pkin(10) + (t59 * t60 + t72) * pkin(5) + t128;
t141 = t103 * t23 + t108 * t25;
t31 = -t60 * qJD(5) - t103 * t57 + t108 * t56;
t52 = t91 * pkin(5) - t60 * pkin(10);
t58 = t59 ^ 2;
t15 = -t58 * pkin(5) + t31 * pkin(10) - t91 * t52 + t141;
t43 = -t102 * t60 + t107 * t59;
t20 = t43 * qJD(6) + t102 * t31 + t107 * t32;
t44 = t102 * t59 + t107 * t60;
t29 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t85 = qJD(6) + t91;
t35 = -t85 * mrSges(7,2) + t43 * mrSges(7,3);
t69 = qJDD(6) + t72;
t12 = m(7) * (-t102 * t15 + t107 * t14) - t20 * mrSges(7,3) + t69 * mrSges(7,1) - t44 * t29 + t85 * t35;
t19 = -t44 * qJD(6) - t102 * t32 + t107 * t31;
t36 = t85 * mrSges(7,1) - t44 * mrSges(7,3);
t13 = m(7) * (t102 * t14 + t107 * t15) + t19 * mrSges(7,3) - t69 * mrSges(7,2) + t43 * t29 - t85 * t36;
t45 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t50 = -t91 * mrSges(6,2) + t59 * mrSges(6,3);
t8 = m(6) * t128 + t72 * mrSges(6,1) - t32 * mrSges(6,3) + t102 * t13 + t107 * t12 - t60 * t45 + t91 * t50;
t51 = t91 * mrSges(6,1) - t60 * mrSges(6,3);
t9 = m(6) * t141 - t72 * mrSges(6,2) + t31 * mrSges(6,3) - t102 * t12 + t107 * t13 + t59 * t45 - t91 * t51;
t6 = m(5) * t127 + t76 * mrSges(5,1) - t57 * mrSges(5,3) + t103 * t9 + t108 * t8 - t78 * t61 + t93 * t62;
t63 = t93 * mrSges(5,1) - t78 * mrSges(5,3);
t7 = m(5) * t140 - t76 * mrSges(5,2) + t56 * mrSges(5,3) - t103 * t8 + t108 * t9 + t77 * t61 - t93 * t63;
t88 = -mrSges(4,1) * t135 - qJD(2) * mrSges(4,3);
t125 = t104 * t6 - t109 * t7 - m(4) * (-t83 * pkin(2) + t119 - t144) - t88 * t135 + t82 * mrSges(4,3);
t89 = mrSges(4,1) * t96 + qJD(2) * mrSges(4,2);
t137 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t96 - t89;
t143 = mrSges(3,1) - mrSges(4,2);
t87 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t149 = (t137 * t105 - t110 * t87) * qJD(1) - t143 * t83 + m(3) * (t124 - t144) + t82 * mrSges(3,2) - t125;
t131 = -t105 * g(3) + t110 * t71;
t148 = t112 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t150 - t79 * t135 - t131;
t117 = -t101 * t113 * pkin(8) + t83 * pkin(3) + qJD(2) * t90 - t148;
t116 = -t56 * pkin(4) - t75 * pkin(9) + t78 * t64 + t117;
t122 = t19 * mrSges(7,1) + t43 * t35 - m(7) * (-t31 * pkin(5) - t58 * pkin(10) + t60 * t52 + t116) - t20 * mrSges(7,2) - t44 * t36;
t118 = -m(6) * t116 + t31 * mrSges(6,1) - t32 * mrSges(6,2) + t59 * t50 - t60 * t51 + t122;
t115 = -m(5) * t117 + t56 * mrSges(5,1) - t57 * mrSges(5,2) + t77 * t62 - t78 * t63 + t118;
t114 = m(4) * t148 + t115;
t80 = (mrSges(4,2) * t110 - mrSges(4,3) * t105) * qJD(1);
t138 = t80 + (-mrSges(3,1) * t110 + mrSges(3,2) * t105) * qJD(1);
t142 = mrSges(3,3) + mrSges(4,1);
t11 = -t114 + t142 * t83 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t137 * qJD(2) + m(3) * t131 + t138 * t135;
t123 = -m(4) * t49 - t104 * t7 - t109 * t6;
t4 = m(3) * t139 - t142 * t82 + t143 * qJDD(2) + (t87 - t88) * qJD(2) - t138 * t96 + t123;
t145 = t105 * t11 + t110 * t4;
t2 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t149;
t1 = m(2) * t126 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t4 + t110 * t11;
t3 = [-m(1) * g(1) + t111 * t1 - t106 * t2, t1, t11, t83 * mrSges(4,2) - t89 * t96 - t125, t7, t9, t13; -m(1) * g(2) + t106 * t1 + t111 * t2, t2, t4, -t83 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t89 - t80 * t135 + t114, t6, t8, t12; (-m(1) - m(2)) * g(3) + t145, -m(2) * g(3) + t145, t149, t82 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t88 + t80 * t96 - t123, -t115, -t118, -t122;];
f_new  = t3;
