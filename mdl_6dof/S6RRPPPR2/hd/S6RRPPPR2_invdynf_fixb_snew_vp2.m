% Calculate vector of cutting forces with Newton-Euler
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-05-06 08:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:22:09
% EndTime: 2019-05-06 08:22:14
% DurationCPUTime: 2.64s
% Computational Cost: add. (29894->205), mult. (71211->263), div. (0->0), fcn. (48316->10), ass. (0->100)
t151 = -2 * qJD(4);
t104 = sin(qJ(2));
t107 = cos(qJ(2));
t110 = qJD(1) ^ 2;
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t128 = t105 * g(1) - t108 * g(2);
t121 = -qJDD(1) * pkin(1) - t128;
t133 = qJD(1) * t104;
t131 = qJD(1) * qJD(2);
t91 = t107 * qJDD(1) - t104 * t131;
t92 = qJD(2) * pkin(2) - qJ(3) * t133;
t99 = t107 ^ 2;
t114 = -t91 * pkin(2) + qJDD(3) + t92 * t133 + (-qJ(3) * t99 - pkin(7)) * t110 + t121;
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t101 = sin(pkin(9));
t132 = qJD(1) * t107;
t136 = cos(pkin(9));
t83 = t101 * t133 - t136 * t132;
t135 = qJD(2) * t83;
t109 = qJD(2) ^ 2;
t125 = -t108 * g(1) - t105 * g(2);
t87 = -t110 * pkin(1) + qJDD(1) * pkin(7) + t125;
t137 = t104 * t87;
t145 = pkin(2) * t110;
t90 = t104 * qJDD(1) + t107 * t131;
t42 = qJDD(2) * pkin(2) - t90 * qJ(3) - t137 + (qJ(3) * t131 + t104 * t145 - g(3)) * t107;
t127 = -t104 * g(3) + t107 * t87;
t44 = t91 * qJ(3) - qJD(2) * t92 - t99 * t145 + t127;
t124 = -t101 * t44 + t136 * t42;
t84 = (t101 * t107 + t136 * t104) * qJD(1);
t57 = t83 * pkin(3) - t84 * qJ(4);
t26 = -qJDD(2) * pkin(3) - t109 * qJ(4) + qJDD(4) - t124 + ((2 * qJD(3)) + t57) * t84;
t64 = t101 * t91 + t136 * t90;
t20 = (t83 * t84 - qJDD(2)) * qJ(5) + (t64 + t135) * pkin(4) + t26;
t111 = (-t64 + t135) * qJ(4) + t114 + (qJD(2) * pkin(3) + t151) * t84;
t63 = t101 * t90 - t136 * t91;
t73 = t84 * pkin(4) - qJD(2) * qJ(5);
t82 = t83 ^ 2;
t24 = -t82 * pkin(4) - t84 * t73 + (pkin(3) + qJ(5)) * t63 + t111;
t69 = t102 * qJD(2) + t100 * t83;
t126 = -0.2e1 * qJD(5) * t69 - t100 * t24 + t102 * t20;
t54 = t102 * qJDD(2) + t100 * t63;
t68 = -t100 * qJD(2) + t102 * t83;
t14 = (t68 * t84 - t54) * pkin(8) + (t68 * t69 + t64) * pkin(5) + t126;
t129 = 0.2e1 * qJD(5) * t68 + t100 * t20 + t102 * t24;
t51 = t84 * pkin(5) - t69 * pkin(8);
t53 = -t100 * qJDD(2) + t102 * t63;
t67 = t68 ^ 2;
t15 = -t67 * pkin(5) + t53 * pkin(8) - t84 * t51 + t129;
t38 = -t103 * t69 + t106 * t68;
t31 = t38 * qJD(6) + t103 * t53 + t106 * t54;
t39 = t103 * t68 + t106 * t69;
t33 = -mrSges(7,1) * t38 + mrSges(7,2) * t39;
t81 = qJD(6) + t84;
t34 = -t81 * mrSges(7,2) + t38 * mrSges(7,3);
t62 = qJDD(6) + t64;
t12 = m(7) * (-t103 * t15 + t106 * t14) - t31 * mrSges(7,3) + t62 * mrSges(7,1) - t39 * t33 + t81 * t34;
t30 = -t39 * qJD(6) - t103 * t54 + t106 * t53;
t35 = t81 * mrSges(7,1) - t39 * mrSges(7,3);
t13 = m(7) * (t103 * t14 + t106 * t15) + t30 * mrSges(7,3) - t62 * mrSges(7,2) + t38 * t33 - t81 * t35;
t43 = -t68 * mrSges(6,1) + mrSges(6,2) * t69;
t50 = t84 * mrSges(6,1) - t69 * mrSges(6,3);
t10 = m(6) * t129 - t64 * mrSges(6,2) + t53 * mrSges(6,3) - t103 * t12 + t106 * t13 + t68 * t43 - t84 * t50;
t75 = t84 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t49 = -t84 * mrSges(6,2) + t68 * mrSges(6,3);
t9 = m(6) * t126 + t64 * mrSges(6,1) - t54 * mrSges(6,3) + t103 * t13 + t106 * t12 - t69 * t43 + t84 * t49;
t122 = -t102 * t10 + t100 * t9 - m(5) * (t63 * pkin(3) + t111) + t84 * t75 + t64 * mrSges(5,3);
t74 = t83 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t138 = -qJD(2) * mrSges(4,2) - t83 * mrSges(4,3) - t74;
t142 = mrSges(4,1) - mrSges(5,2);
t72 = qJD(2) * mrSges(4,1) - t84 * mrSges(4,3);
t113 = m(4) * t114 + t64 * mrSges(4,2) + t138 * t83 + t142 * t63 + t84 * t72 - t122;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t133;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t132;
t150 = t113 + (t104 * t93 - t107 * t94) * qJD(1) - t91 * mrSges(3,1) + t90 * mrSges(3,2) + m(3) * (-t110 * pkin(7) + t121);
t140 = t101 * t42 + t136 * t44;
t149 = t109 * pkin(3) - qJDD(2) * qJ(4) + qJD(2) * t151 + t83 * t57 - t140;
t134 = qJD(3) * t83;
t78 = -0.2e1 * t134;
t116 = -t63 * pkin(4) - t82 * qJ(5) + qJD(2) * t73 + qJDD(5) - t149 + t78;
t118 = -t30 * mrSges(7,1) - t38 * t34 + m(7) * (-t53 * pkin(5) - t67 * pkin(8) + t69 * t51 + t116) + t31 * mrSges(7,2) + t39 * t35;
t115 = m(6) * t116 - t53 * mrSges(6,1) + t54 * mrSges(6,2) - t68 * t49 + t69 * t50 + t118;
t112 = -m(5) * (0.2e1 * t134 + t149) + t115;
t59 = -t83 * mrSges(5,2) - t84 * mrSges(5,3);
t139 = -t83 * mrSges(4,1) - t84 * mrSges(4,2) - t59;
t141 = -mrSges(4,3) - mrSges(5,1);
t11 = t139 * t83 + t141 * t63 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(2) + (-t72 + t75) * qJD(2) + m(4) * (t78 + t140) + t112;
t119 = -m(5) * t26 - t100 * t10 - t102 * t9;
t7 = m(4) * t124 + (-0.2e1 * m(4) * qJD(3) + t139) * t84 + t141 * t64 + t142 * qJDD(2) + t138 * qJD(2) + t119;
t89 = (-mrSges(3,1) * t107 + mrSges(3,2) * t104) * qJD(1);
t4 = m(3) * (-t107 * g(3) - t137) - t90 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t89 * t133 + qJD(2) * t94 + t101 * t11 + t136 * t7;
t5 = m(3) * t127 - qJDD(2) * mrSges(3,2) + t91 * mrSges(3,3) - qJD(2) * t93 - t101 * t7 + t136 * t11 + t89 * t132;
t147 = t104 * t5 + t107 * t4;
t6 = m(2) * t128 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t150;
t1 = m(2) * t125 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t4 + t107 * t5;
t2 = [-m(1) * g(1) + t108 * t1 - t105 * t6, t1, t5, t11, -t63 * mrSges(5,2) - t83 * t74 - t122, t10, t13; -m(1) * g(2) + t105 * t1 + t108 * t6, t6, t4, t7, t63 * mrSges(5,1) - qJDD(2) * mrSges(5,3) - qJD(2) * t75 + t83 * t59 - t112, t9, t12; (-m(1) - m(2)) * g(3) + t147, -m(2) * g(3) + t147, t150, t113, t64 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + qJD(2) * t74 + t84 * t59 - t119, t115, t118;];
f_new  = t2;
