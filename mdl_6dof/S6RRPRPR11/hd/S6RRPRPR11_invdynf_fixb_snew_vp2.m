% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:59:36
% EndTime: 2019-05-06 15:59:45
% DurationCPUTime: 3.41s
% Computational Cost: add. (43533->206), mult. (92414->259), div. (0->0), fcn. (58569->10), ass. (0->101)
t150 = -2 * qJD(3);
t106 = sin(qJ(2));
t110 = cos(qJ(2));
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t132 = t107 * g(1) - t111 * g(2);
t124 = -qJDD(1) * pkin(1) - t132;
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t135 = qJD(1) * qJD(2);
t129 = t110 * t135;
t130 = t106 * t135;
t83 = t106 * qJDD(1) + t129;
t96 = t106 * qJD(1);
t119 = pkin(2) * t130 + t96 * t150 + (-t83 - t129) * qJ(3) + t124;
t136 = qJD(1) * t110;
t113 = qJD(1) ^ 2;
t144 = t113 * pkin(7);
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t101 = t110 ^ 2;
t84 = t110 * qJDD(1) - t130;
t90 = pkin(3) * t96 - qJD(2) * pkin(8);
t31 = -t90 * t96 + (-pkin(2) - pkin(8)) * t84 + (-pkin(3) * t101 - pkin(7)) * t113 + t119;
t112 = qJD(2) ^ 2;
t126 = -t111 * g(1) - t107 * g(2);
t72 = -t113 * pkin(1) + qJDD(1) * pkin(7) + t126;
t140 = -t110 * g(3) - t106 * t72;
t80 = (-pkin(2) * t110 - qJ(3) * t106) * qJD(1);
t49 = -qJDD(2) * pkin(2) - t112 * qJ(3) + t80 * t96 + qJDD(3) - t140;
t39 = (-t106 * t110 * t113 - qJDD(2)) * pkin(8) + (t83 - t129) * pkin(3) + t49;
t128 = -t105 * t31 + t109 * t39;
t78 = -t105 * qJD(2) - t109 * t136;
t59 = t78 * qJD(4) + t109 * qJDD(2) - t105 * t84;
t79 = t109 * qJD(2) - t105 * t136;
t63 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t93 = t96 + qJD(4);
t64 = -t93 * mrSges(5,2) + t78 * mrSges(5,3);
t77 = qJDD(4) + t83;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t20 = (t78 * t93 - t59) * qJ(5) + (t78 * t79 + t77) * pkin(4) + t128;
t141 = t105 * t39 + t109 * t31;
t58 = -t79 * qJD(4) - t105 * qJDD(2) - t109 * t84;
t65 = t93 * pkin(4) - t79 * qJ(5);
t76 = t78 ^ 2;
t22 = -t76 * pkin(4) + t58 * qJ(5) - t93 * t65 + t141;
t62 = t102 * t78 + t103 * t79;
t127 = -0.2e1 * qJD(5) * t62 - t102 * t22 + t103 * t20;
t41 = t102 * t58 + t103 * t59;
t61 = -t102 * t79 + t103 * t78;
t14 = (t61 * t93 - t41) * pkin(9) + (t61 * t62 + t77) * pkin(5) + t127;
t133 = 0.2e1 * qJD(5) * t61 + t102 * t20 + t103 * t22;
t40 = -t102 * t59 + t103 * t58;
t52 = t93 * pkin(5) - t62 * pkin(9);
t60 = t61 ^ 2;
t15 = -t60 * pkin(5) + t40 * pkin(9) - t93 * t52 + t133;
t43 = -t104 * t62 + t108 * t61;
t25 = t43 * qJD(6) + t104 * t40 + t108 * t41;
t44 = t104 * t61 + t108 * t62;
t29 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t91 = qJD(6) + t93;
t32 = -t91 * mrSges(7,2) + t43 * mrSges(7,3);
t73 = qJDD(6) + t77;
t10 = m(7) * (-t104 * t15 + t108 * t14) - t25 * mrSges(7,3) + t73 * mrSges(7,1) - t44 * t29 + t91 * t32;
t24 = -t44 * qJD(6) - t104 * t41 + t108 * t40;
t33 = t91 * mrSges(7,1) - t44 * mrSges(7,3);
t11 = m(7) * (t104 * t14 + t108 * t15) + t24 * mrSges(7,3) - t73 * mrSges(7,2) + t43 * t29 - t91 * t33;
t45 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t50 = -t93 * mrSges(6,2) + t61 * mrSges(6,3);
t8 = m(6) * t127 + t77 * mrSges(6,1) - t41 * mrSges(6,3) + t108 * t10 + t104 * t11 - t62 * t45 + t93 * t50;
t51 = t93 * mrSges(6,1) - t62 * mrSges(6,3);
t9 = m(6) * t133 - t77 * mrSges(6,2) + t40 * mrSges(6,3) - t104 * t10 + t108 * t11 + t61 * t45 - t93 * t51;
t6 = m(5) * t128 + t77 * mrSges(5,1) - t59 * mrSges(5,3) + t102 * t9 + t103 * t8 - t79 * t63 + t93 * t64;
t66 = t93 * mrSges(5,1) - t79 * mrSges(5,3);
t7 = m(5) * t141 - t77 * mrSges(5,2) + t58 * mrSges(5,3) - t102 * t8 + t103 * t9 + t78 * t63 - t93 * t66;
t88 = -mrSges(4,1) * t136 - qJD(2) * mrSges(4,3);
t125 = t105 * t6 - t109 * t7 - m(4) * (-t84 * pkin(2) + t119 - t144) - t88 * t136 + t83 * mrSges(4,3);
t89 = mrSges(4,1) * t96 + qJD(2) * mrSges(4,2);
t138 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t96 - t89;
t143 = mrSges(3,1) - mrSges(4,2);
t87 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t149 = qJD(1) * (t138 * t106 - t110 * t87) - t143 * t84 + m(3) * (t124 - t144) + t83 * mrSges(3,2) - t125;
t131 = -t106 * g(3) + t110 * t72;
t148 = t112 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t150 - t80 * t136 - t131;
t117 = -t101 * t113 * pkin(8) + t84 * pkin(3) + qJD(2) * t90 - t148;
t116 = -t58 * pkin(4) - t76 * qJ(5) + t79 * t65 + qJDD(5) + t117;
t122 = t24 * mrSges(7,1) + t43 * t32 - m(7) * (-t40 * pkin(5) - t60 * pkin(9) + t62 * t52 + t116) - t25 * mrSges(7,2) - t44 * t33;
t118 = -m(6) * t116 + t40 * mrSges(6,1) - t41 * mrSges(6,2) + t61 * t50 - t62 * t51 + t122;
t115 = -m(5) * t117 + t58 * mrSges(5,1) - t59 * mrSges(5,2) + t78 * t64 - t79 * t66 + t118;
t114 = m(4) * t148 + t115;
t81 = (mrSges(4,2) * t110 - mrSges(4,3) * t106) * qJD(1);
t139 = t81 + (-mrSges(3,1) * t110 + mrSges(3,2) * t106) * qJD(1);
t142 = mrSges(3,3) + mrSges(4,1);
t13 = -t114 + t142 * t84 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t138 * qJD(2) + m(3) * t131 + t139 * t136;
t123 = -m(4) * t49 - t105 * t7 - t109 * t6;
t4 = m(3) * t140 - t142 * t83 + t143 * qJDD(2) + (t87 - t88) * qJD(2) - t139 * t96 + t123;
t145 = t106 * t13 + t110 * t4;
t2 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t149;
t1 = m(2) * t126 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t110 * t13;
t3 = [-m(1) * g(1) + t111 * t1 - t107 * t2, t1, t13, t84 * mrSges(4,2) - t89 * t96 - t125, t7, t9, t11; -m(1) * g(2) + t107 * t1 + t111 * t2, t2, t4, -t84 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t89 - t81 * t136 + t114, t6, t8, t10; (-m(1) - m(2)) * g(3) + t145, -m(2) * g(3) + t145, t149, t83 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t88 + t81 * t96 - t123, -t115, -t118, -t122;];
f_new  = t3;
