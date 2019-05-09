% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:03:45
% EndTime: 2019-05-06 04:03:52
% DurationCPUTime: 2.79s
% Computational Cost: add. (37533->179), mult. (78422->231), div. (0->0), fcn. (55425->10), ass. (0->94)
t92 = sin(qJ(1));
t97 = cos(qJ(1));
t111 = -t97 * g(1) - t92 * g(2);
t108 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t111;
t127 = -m(2) - m(3);
t126 = -pkin(1) - pkin(7);
t125 = (mrSges(2,1) - mrSges(3,2));
t124 = -mrSges(2,2) + mrSges(3,3);
t118 = qJD(1) * qJD(3);
t91 = sin(qJ(3));
t113 = t91 * t118;
t115 = t92 * g(1) - t97 * g(2);
t98 = qJD(1) ^ 2;
t106 = -t98 * qJ(2) + qJDD(2) - t115;
t62 = qJDD(1) * t126 + t106;
t96 = cos(qJ(3));
t121 = t91 * g(3) + t96 * t62;
t74 = t96 * qJDD(1) - t113;
t38 = (-t74 - t113) * pkin(8) + (-t91 * t96 * t98 + qJDD(3)) * pkin(3) + t121;
t114 = -t96 * g(3) + t91 * t62;
t73 = -t91 * qJDD(1) - t118 * t96;
t119 = qJD(1) * t96;
t77 = qJD(3) * pkin(3) - pkin(8) * t119;
t87 = t91 ^ 2;
t39 = -t87 * t98 * pkin(3) + t73 * pkin(8) - qJD(3) * t77 + t114;
t90 = sin(qJ(4));
t95 = cos(qJ(4));
t112 = t95 * t38 - t90 * t39;
t66 = (-t90 * t96 - t91 * t95) * qJD(1);
t46 = t66 * qJD(4) + t90 * t73 + t95 * t74;
t67 = (-t90 * t91 + t95 * t96) * qJD(1);
t84 = qJDD(3) + qJDD(4);
t85 = qJD(3) + qJD(4);
t18 = (t66 * t85 - t46) * pkin(9) + (t66 * t67 + t84) * pkin(4) + t112;
t122 = t90 * t38 + t95 * t39;
t45 = -t67 * qJD(4) + t95 * t73 - t90 * t74;
t60 = t85 * pkin(4) - t67 * pkin(9);
t65 = t66 ^ 2;
t20 = -t65 * pkin(4) + t45 * pkin(9) - t85 * t60 + t122;
t89 = sin(qJ(5));
t94 = cos(qJ(5));
t123 = t89 * t18 + t94 * t20;
t120 = qJD(1) * t91;
t54 = -t66 * mrSges(5,1) + t67 * mrSges(5,2);
t58 = -t85 * mrSges(5,2) + t66 * mrSges(5,3);
t52 = t94 * t66 - t89 * t67;
t53 = t89 * t66 + t94 * t67;
t34 = -t52 * pkin(5) - t53 * pkin(10);
t80 = qJD(5) + t85;
t78 = t80 ^ 2;
t79 = qJDD(5) + t84;
t15 = -t78 * pkin(5) + t79 * pkin(10) + t52 * t34 + t123;
t103 = -t73 * pkin(3) + t77 * t119 + (-pkin(8) * t87 + t126) * t98 + t108;
t101 = -t45 * pkin(4) - t65 * pkin(9) + t67 * t60 + t103;
t27 = -t53 * qJD(5) + t94 * t45 - t89 * t46;
t28 = t52 * qJD(5) + t89 * t45 + t94 * t46;
t16 = (-t52 * t80 - t28) * pkin(10) + (t53 * t80 - t27) * pkin(5) + t101;
t88 = sin(qJ(6));
t93 = cos(qJ(6));
t42 = -t88 * t53 + t93 * t80;
t22 = t42 * qJD(6) + t93 * t28 + t88 * t79;
t26 = qJDD(6) - t27;
t43 = t93 * t53 + t88 * t80;
t29 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t49 = qJD(6) - t52;
t30 = -t49 * mrSges(7,2) + t42 * mrSges(7,3);
t12 = m(7) * (-t88 * t15 + t93 * t16) - t22 * mrSges(7,3) + t26 * mrSges(7,1) - t43 * t29 + t49 * t30;
t21 = -t43 * qJD(6) - t88 * t28 + t93 * t79;
t31 = t49 * mrSges(7,1) - t43 * mrSges(7,3);
t13 = m(7) * (t93 * t15 + t88 * t16) + t21 * mrSges(7,3) - t26 * mrSges(7,2) + t42 * t29 - t49 * t31;
t33 = -t52 * mrSges(6,1) + t53 * mrSges(6,2);
t48 = t80 * mrSges(6,1) - t53 * mrSges(6,3);
t8 = m(6) * t123 - t79 * mrSges(6,2) + t27 * mrSges(6,3) - t88 * t12 + t93 * t13 + t52 * t33 - t80 * t48;
t110 = t94 * t18 - t89 * t20;
t104 = m(7) * (-t79 * pkin(5) - t78 * pkin(10) + t53 * t34 - t110) - t21 * mrSges(7,1) + t22 * mrSges(7,2) - t42 * t30 + t43 * t31;
t47 = -t80 * mrSges(6,2) + t52 * mrSges(6,3);
t9 = m(6) * t110 + t79 * mrSges(6,1) - t28 * mrSges(6,3) - t53 * t33 + t80 * t47 - t104;
t5 = m(5) * t112 + t84 * mrSges(5,1) - t46 * mrSges(5,3) - t67 * t54 + t85 * t58 + t89 * t8 + t94 * t9;
t59 = t85 * mrSges(5,1) - t67 * mrSges(5,3);
t6 = m(5) * t122 - t84 * mrSges(5,2) + t45 * mrSges(5,3) + t66 * t54 - t85 * t59 + t94 * t8 - t89 * t9;
t72 = (mrSges(4,1) * t91 + mrSges(4,2) * t96) * qJD(1);
t75 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t120;
t3 = m(4) * t121 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t75 - t119 * t72 + t95 * t5 + t90 * t6;
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t4 = m(4) * t114 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t76 - t120 * t72 - t90 * t5 + t95 * t6;
t116 = -t91 * t3 + t96 * t4;
t107 = -m(3) * (-qJDD(1) * pkin(1) + t106) - t96 * t3 - t91 * t4;
t105 = m(6) * t101 - t27 * mrSges(6,1) + t28 * mrSges(6,2) + t93 * t12 + t88 * t13 - t52 * t47 + t53 * t48;
t102 = m(5) * t103 - t45 * mrSges(5,1) + t46 * mrSges(5,2) - t66 * t58 + t67 * t59 + t105;
t100 = -t73 * mrSges(4,1) + t102 + m(4) * (t126 * t98 + t108) + t75 * t120 + t76 * t119 + t74 * mrSges(4,2);
t99 = -m(3) * (t98 * pkin(1) - t108) + t100;
t7 = m(2) * t111 + qJDD(1) * t124 - (t125 * t98) + t99;
t1 = m(2) * t115 + qJDD(1) * t125 + t124 * t98 + t107;
t2 = [-m(1) * g(1) - t92 * t1 + t97 * t7, t7, -m(3) * g(3) + t116, t4, t6, t8, t13; -m(1) * g(2) + t97 * t1 + t92 * t7, t1, -(t98 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t99, t3, t5, t9, t12; (-m(1) + t127) * g(3) + t116, g(3) * t127 + t116, qJDD(1) * mrSges(3,2) - t98 * mrSges(3,3) - t107, t100, t102, t105, t104;];
f_new  = t2;
