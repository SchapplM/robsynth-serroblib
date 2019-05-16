% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-05-06 09:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:02:18
% EndTime: 2019-05-06 09:02:26
% DurationCPUTime: 3.46s
% Computational Cost: add. (41063->200), mult. (96879->262), div. (0->0), fcn. (68028->10), ass. (0->96)
t142 = -2 * qJD(3);
t100 = sin(pkin(9));
t131 = cos(pkin(9));
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t127 = qJD(1) * qJD(2);
t108 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t120 = -g(1) * t106 - g(2) * t104;
t88 = -pkin(1) * t108 + qJDD(1) * pkin(7) + t120;
t132 = t103 * t88;
t136 = pkin(2) * t108;
t91 = qJDD(1) * t103 + t105 * t127;
t51 = qJDD(2) * pkin(2) - qJ(3) * t91 - t132 + (qJ(3) * t127 + t103 * t136 - g(3)) * t105;
t122 = -g(3) * t103 + t105 * t88;
t92 = qJDD(1) * t105 - t103 * t127;
t129 = qJD(1) * t103;
t93 = qJD(2) * pkin(2) - qJ(3) * t129;
t98 = t105 ^ 2;
t53 = qJ(3) * t92 - qJD(2) * t93 - t98 * t136 + t122;
t85 = (t100 * t105 + t131 * t103) * qJD(1);
t141 = -t100 * t53 + t131 * t51 + t85 * t142;
t102 = sin(qJ(5));
t137 = cos(qJ(5));
t101 = cos(pkin(10));
t107 = qJD(2) ^ 2;
t128 = qJD(1) * t105;
t84 = t100 * t129 - t131 * t128;
t124 = t100 * t51 + t131 * t53 + t84 * t142;
t64 = pkin(3) * t84 - qJ(4) * t85;
t27 = -pkin(3) * t107 + qJDD(2) * qJ(4) - t64 * t84 + t124;
t123 = g(1) * t104 - t106 * g(2);
t117 = -qJDD(1) * pkin(1) - t123;
t112 = -pkin(2) * t92 + qJDD(3) + t93 * t129 + (-qJ(3) * t98 - pkin(7)) * t108 + t117;
t69 = t100 * t91 - t131 * t92;
t70 = t100 * t92 + t131 * t91;
t30 = (qJD(2) * t84 - t70) * qJ(4) + (qJD(2) * t85 + t69) * pkin(3) + t112;
t99 = sin(pkin(10));
t76 = qJD(2) * t99 + t101 * t85;
t121 = -0.2e1 * qJD(4) * t76 + t101 * t30 - t27 * t99;
t62 = qJDD(2) * t99 + t101 * t70;
t75 = qJD(2) * t101 - t85 * t99;
t20 = (t75 * t84 - t62) * pkin(8) + (t75 * t76 + t69) * pkin(4) + t121;
t125 = 0.2e1 * qJD(4) * t75 + t101 * t27 + t99 * t30;
t59 = pkin(4) * t84 - pkin(8) * t76;
t61 = qJDD(2) * t101 - t70 * t99;
t74 = t75 ^ 2;
t22 = -pkin(4) * t74 + t61 * pkin(8) - t59 * t84 + t125;
t134 = t102 * t20 + t137 * t22;
t47 = t102 * t76 - t137 * t75;
t48 = t102 * t75 + t137 * t76;
t37 = pkin(5) * t47 - qJ(6) * t48;
t83 = qJD(5) + t84;
t44 = -mrSges(7,1) * t83 + t48 * mrSges(7,2);
t68 = qJDD(5) + t69;
t82 = t83 ^ 2;
t126 = m(7) * (-pkin(5) * t82 + qJ(6) * t68 + 0.2e1 * qJD(6) * t83 - t47 * t37 + t134) + t83 * t44 + t68 * mrSges(7,3);
t38 = mrSges(7,1) * t47 - mrSges(7,3) * t48;
t133 = -mrSges(6,1) * t47 - mrSges(6,2) * t48 - t38;
t135 = -mrSges(6,3) - mrSges(7,2);
t33 = t48 * qJD(5) + t102 * t62 - t137 * t61;
t43 = mrSges(6,1) * t83 - t48 * mrSges(6,3);
t12 = m(6) * t134 - t68 * mrSges(6,2) + t133 * t47 + t135 * t33 - t83 * t43 + t126;
t116 = -t102 * t22 + t137 * t20;
t138 = m(7) * (-t68 * pkin(5) - t82 * qJ(6) + t48 * t37 + qJDD(6) - t116);
t34 = -t47 * qJD(5) + t102 * t61 + t137 * t62;
t41 = -t47 * mrSges(7,2) + mrSges(7,3) * t83;
t42 = -mrSges(6,2) * t83 - t47 * mrSges(6,3);
t13 = m(6) * t116 - t138 + (t42 + t41) * t83 + (mrSges(6,1) + mrSges(7,1)) * t68 + t133 * t48 + t135 * t34;
t52 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t57 = -mrSges(5,2) * t84 + mrSges(5,3) * t75;
t10 = m(5) * t121 + t69 * mrSges(5,1) - t62 * mrSges(5,3) + t102 * t12 + t137 * t13 - t76 * t52 + t84 * t57;
t58 = mrSges(5,1) * t84 - mrSges(5,3) * t76;
t11 = m(5) * t125 - t69 * mrSges(5,2) + t61 * mrSges(5,3) - t102 * t13 + t137 * t12 + t75 * t52 - t84 * t58;
t77 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t84;
t78 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t85;
t114 = m(4) * t112 + t69 * mrSges(4,1) + t70 * mrSges(4,2) + t101 * t10 + t99 * t11 + t84 * t77 + t85 * t78;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t140 = (t103 * t94 - t105 * t95) * qJD(1) + m(3) * (-pkin(7) * t108 + t117) - t92 * mrSges(3,1) + t91 * mrSges(3,2) + t114;
t26 = -qJDD(2) * pkin(3) - t107 * qJ(4) + t85 * t64 + qJDD(4) - t141;
t110 = -t61 * pkin(4) - t74 * pkin(8) + t76 * t59 + t26;
t115 = t34 * mrSges(7,3) + t48 * t44 - m(7) * (-0.2e1 * qJD(6) * t48 + (t47 * t83 - t34) * qJ(6) + (t48 * t83 + t33) * pkin(5) + t110) - t33 * mrSges(7,1) - t47 * t41;
t113 = m(6) * t110 + t33 * mrSges(6,1) + t34 * mrSges(6,2) + t47 * t42 + t48 * t43 - t115;
t109 = m(5) * t26 - t61 * mrSges(5,1) + t62 * mrSges(5,2) - t75 * t57 + t76 * t58 + t113;
t65 = mrSges(4,1) * t84 + mrSges(4,2) * t85;
t14 = m(4) * t141 + qJDD(2) * mrSges(4,1) - t70 * mrSges(4,3) + qJD(2) * t77 - t85 * t65 - t109;
t7 = m(4) * t124 - qJDD(2) * mrSges(4,2) - t69 * mrSges(4,3) - qJD(2) * t78 - t99 * t10 + t101 * t11 - t84 * t65;
t90 = (-mrSges(3,1) * t105 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * (-g(3) * t105 - t132) - t91 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t90 * t129 + qJD(2) * t95 + t100 * t7 + t131 * t14;
t5 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t94 - t100 * t14 + t90 * t128 + t131 * t7;
t139 = t103 * t5 + t105 * t4;
t6 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t108 * mrSges(2,2) - t140;
t1 = m(2) * t120 - t108 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t1 * t106 - t104 * t6, t1, t5, t7, t11, t12, -t33 * mrSges(7,2) - t47 * t38 + t126; -m(1) * g(2) + t1 * t104 + t106 * t6, t6, t4, t14, t10, t13, -t115; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, t114, t109, t113, -t68 * mrSges(7,1) + t34 * mrSges(7,2) + t48 * t38 - t83 * t41 + t138;];
f_new  = t2;
