% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:45:05
% EndTime: 2019-05-05 04:45:16
% DurationCPUTime: 6.68s
% Computational Cost: add. (116647->181), mult. (283231->260), div. (0->0), fcn. (225775->16), ass. (0->108)
t103 = cos(pkin(7));
t113 = qJD(2) ^ 2;
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t100 = sin(pkin(6));
t133 = t100 * t112;
t104 = cos(pkin(6));
t102 = cos(pkin(12));
t98 = sin(pkin(12));
t86 = t98 * g(1) - t102 * g(2);
t138 = t104 * t86;
t87 = -t102 * g(1) - t98 * g(2);
t96 = -g(3) + qJDD(1);
t121 = -t108 * t87 + t112 * t138 + t96 * t133;
t99 = sin(pkin(7));
t142 = pkin(9) * t99;
t53 = qJDD(2) * pkin(2) + t113 * t142 + t121;
t76 = -t100 * t86 + t104 * t96;
t144 = t103 * t53 + t76 * t99;
t101 = cos(pkin(13));
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t136 = qJD(2) * t99;
t97 = sin(pkin(13));
t74 = (t101 * t111 - t107 * t97) * t136;
t134 = t100 * t108;
t129 = t108 * t138 + t112 * t87 + t96 * t134;
t54 = -t113 * pkin(2) + qJDD(2) * t142 + t129;
t122 = -t107 * t54 + t144 * t111;
t126 = t111 * t136;
t137 = t113 * t99 ^ 2;
t132 = qJD(2) * qJD(3);
t81 = (qJDD(2) * t107 + t111 * t132) * t99;
t91 = t103 * qJDD(2) + qJDD(3);
t92 = t103 * qJD(2) + qJD(3);
t33 = (t126 * t92 - t81) * qJ(4) + (t107 * t111 * t137 + t91) * pkin(3) + t122;
t128 = t111 ^ 2 * t137;
t130 = t144 * t107 + t111 * t54;
t127 = t107 * t136;
t77 = t92 * pkin(3) - qJ(4) * t127;
t82 = (qJDD(2) * t111 - t107 * t132) * t99;
t34 = -pkin(3) * t128 + t82 * qJ(4) - t92 * t77 + t130;
t75 = (t101 * t107 + t111 * t97) * t136;
t143 = -0.2e1 * qJD(4) * t75 + t101 * t33 - t97 * t34;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t131 = 0.2e1 * qJD(4) * t74 + t101 * t34 + t97 * t33;
t56 = -t74 * pkin(4) - t75 * pkin(10);
t90 = t92 ^ 2;
t25 = -t90 * pkin(4) + t91 * pkin(10) + t74 * t56 + t131;
t124 = t103 * t76 - t99 * t53;
t115 = -t82 * pkin(3) - qJ(4) * t128 + t77 * t127 + qJDD(4) + t124;
t59 = t101 * t82 - t97 * t81;
t60 = t101 * t81 + t97 * t82;
t27 = (-t74 * t92 - t60) * pkin(10) + (t75 * t92 - t59) * pkin(4) + t115;
t140 = t106 * t27 + t110 * t25;
t105 = sin(qJ(6));
t109 = cos(qJ(6));
t62 = -t106 * t75 + t110 * t92;
t63 = t106 * t92 + t110 * t75;
t43 = -t62 * pkin(5) - t63 * pkin(11);
t58 = qJDD(5) - t59;
t73 = qJD(5) - t74;
t72 = t73 ^ 2;
t21 = -t72 * pkin(5) + t58 * pkin(11) + t62 * t43 + t140;
t24 = -t91 * pkin(4) - t90 * pkin(10) + t75 * t56 - t143;
t40 = -t63 * qJD(5) - t106 * t60 + t110 * t91;
t41 = t62 * qJD(5) + t106 * t91 + t110 * t60;
t22 = (-t62 * t73 - t41) * pkin(11) + (t63 * t73 - t40) * pkin(5) + t24;
t45 = -t105 * t63 + t109 * t73;
t29 = t45 * qJD(6) + t105 * t58 + t109 * t41;
t46 = t105 * t73 + t109 * t63;
t35 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t61 = qJD(6) - t62;
t37 = -t61 * mrSges(7,2) + t45 * mrSges(7,3);
t39 = qJDD(6) - t40;
t18 = m(7) * (-t105 * t21 + t109 * t22) - t29 * mrSges(7,3) + t39 * mrSges(7,1) - t46 * t35 + t61 * t37;
t28 = -t46 * qJD(6) - t105 * t41 + t109 * t58;
t38 = t61 * mrSges(7,1) - t46 * mrSges(7,3);
t19 = m(7) * (t105 * t22 + t109 * t21) + t28 * mrSges(7,3) - t39 * mrSges(7,2) + t45 * t35 - t61 * t38;
t42 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t48 = t73 * mrSges(6,1) - t63 * mrSges(6,3);
t15 = m(6) * t140 - t58 * mrSges(6,2) + t40 * mrSges(6,3) - t105 * t18 + t109 * t19 + t62 * t42 - t73 * t48;
t119 = -t106 * t25 + t110 * t27;
t116 = m(7) * (-t58 * pkin(5) - t72 * pkin(11) + t63 * t43 - t119) - t28 * mrSges(7,1) + t29 * mrSges(7,2) - t45 * t37 + t46 * t38;
t47 = -t73 * mrSges(6,2) + t62 * mrSges(6,3);
t17 = m(6) * t119 + t58 * mrSges(6,1) - t41 * mrSges(6,3) - t63 * t42 + t73 * t47 - t116;
t64 = -t92 * mrSges(5,2) + t74 * mrSges(5,3);
t65 = t92 * mrSges(5,1) - t75 * mrSges(5,3);
t117 = m(5) * t115 - t59 * mrSges(5,1) + t60 * mrSges(5,2) + t106 * t15 + t110 * t17 - t74 * t64 + t75 * t65;
t78 = t92 * mrSges(4,1) - mrSges(4,3) * t127;
t79 = -t92 * mrSges(4,2) + mrSges(4,3) * t126;
t12 = m(4) * t124 + t81 * mrSges(4,2) - t82 * mrSges(4,1) + (t107 * t78 - t111 * t79) * t136 + t117;
t55 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t11 = m(5) * t131 - t91 * mrSges(5,2) + t59 * mrSges(5,3) - t106 * t17 + t110 * t15 + t74 * t55 - t92 * t65;
t114 = m(6) * t24 - t40 * mrSges(6,1) + t41 * mrSges(6,2) + t105 * t19 + t109 * t18 - t62 * t47 + t63 * t48;
t13 = m(5) * t143 + t91 * mrSges(5,1) - t60 * mrSges(5,3) - t75 * t55 + t92 * t64 - t114;
t80 = (-mrSges(4,1) * t111 + mrSges(4,2) * t107) * t136;
t10 = m(4) * t130 - t91 * mrSges(4,2) + t82 * mrSges(4,3) + t101 * t11 + t126 * t80 - t97 * t13 - t92 * t78;
t9 = m(4) * t122 + t91 * mrSges(4,1) - t81 * mrSges(4,3) + t101 * t13 + t97 * t11 - t127 * t80 + t92 * t79;
t120 = t107 * t10 + t111 * t9;
t4 = m(3) * t121 + qJDD(2) * mrSges(3,1) - t113 * mrSges(3,2) + t103 * t120 - t99 * t12;
t6 = m(3) * t76 + t103 * t12 + t120 * t99;
t8 = m(3) * t129 - t113 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t111 * t10 - t107 * t9;
t125 = m(2) * t96 + t104 * t6 + t4 * t133 + t8 * t134;
t2 = m(2) * t87 - t108 * t4 + t112 * t8;
t1 = m(2) * t86 - t100 * t6 + (t108 * t8 + t112 * t4) * t104;
t3 = [-m(1) * g(1) - t98 * t1 + t102 * t2, t2, t8, t10, t11, t15, t19; -m(1) * g(2) + t102 * t1 + t98 * t2, t1, t4, t9, t13, t17, t18; -m(1) * g(3) + t125, t125, t6, t12, t117, t114, t116;];
f_new  = t3;
