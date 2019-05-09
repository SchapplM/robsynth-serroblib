% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 09:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:05:25
% EndTime: 2019-05-05 09:05:32
% DurationCPUTime: 3.84s
% Computational Cost: add. (52386->182), mult. (108390->243), div. (0->0), fcn. (84270->14), ass. (0->103)
t106 = qJD(2) ^ 2;
t102 = sin(qJ(2));
t105 = cos(qJ(2));
t95 = sin(pkin(6));
t126 = t105 * t95;
t93 = sin(pkin(12));
t96 = cos(pkin(12));
t84 = g(1) * t93 - g(2) * t96;
t98 = cos(pkin(6));
t133 = t84 * t98;
t85 = -g(1) * t96 - g(2) * t93;
t92 = -g(3) + qJDD(1);
t118 = -t102 * t85 + t105 * t133 + t92 * t126;
t94 = sin(pkin(7));
t138 = pkin(9) * t94;
t47 = qJDD(2) * pkin(2) + t106 * t138 + t118;
t68 = -t84 * t95 + t92 * t98;
t97 = cos(pkin(7));
t142 = t47 * t97 + t68 * t94;
t101 = sin(qJ(3));
t104 = cos(qJ(3));
t124 = qJD(2) * qJD(3);
t75 = (-qJDD(2) * t104 + t101 * t124) * t94;
t127 = t102 * t95;
t122 = t102 * t133 + t105 * t85 + t92 * t127;
t48 = -pkin(2) * t106 + qJDD(2) * t138 + t122;
t141 = -t101 * t48 + t142 * t104;
t103 = cos(qJ(6));
t100 = sin(qJ(4));
t125 = qJD(2) * t94;
t121 = t101 * t125;
t137 = cos(qJ(4));
t90 = qJD(2) * t97 + qJD(3);
t66 = t100 * t121 - t137 * t90;
t120 = t104 * t125;
t82 = -qJD(4) + t120;
t135 = t66 * t82;
t139 = -2 * qJD(5);
t73 = (-pkin(3) * t104 - pkin(10) * t101) * t125;
t88 = t90 ^ 2;
t89 = qJDD(2) * t97 + qJDD(3);
t27 = -pkin(3) * t89 - pkin(10) * t88 + t73 * t121 - t141;
t74 = (qJDD(2) * t101 + t104 * t124) * t94;
t43 = -t66 * qJD(4) + t100 * t89 + t137 * t74;
t67 = t100 * t90 + t137 * t121;
t107 = (-t43 - t135) * qJ(5) + t27 + (-t82 * pkin(4) + t139) * t67;
t123 = t142 * t101 + t104 * t48;
t28 = -pkin(3) * t88 + pkin(10) * t89 + t73 * t120 + t123;
t62 = t97 * t68;
t30 = pkin(3) * t75 - pkin(10) * t74 + t62 + (-t47 + (pkin(3) * t101 - pkin(10) * t104) * t90 * qJD(2)) * t94;
t117 = -t100 * t28 + t137 * t30;
t49 = pkin(4) * t66 - qJ(5) * t67;
t69 = qJDD(4) + t75;
t81 = t82 ^ 2;
t22 = -t69 * pkin(4) - t81 * qJ(5) + t67 * t49 + qJDD(5) - t117;
t17 = (t66 * t67 - t69) * pkin(11) + (t43 - t135) * pkin(5) + t22;
t42 = t67 * qJD(4) + t100 * t74 - t137 * t89;
t58 = t67 * pkin(5) + pkin(11) * t82;
t65 = t66 ^ 2;
t20 = -t65 * pkin(5) - t67 * t58 + (pkin(4) + pkin(11)) * t42 + t107;
t99 = sin(qJ(6));
t52 = t103 * t66 + t82 * t99;
t33 = t52 * qJD(6) + t103 * t69 + t42 * t99;
t53 = -t103 * t82 + t66 * t99;
t35 = -mrSges(7,1) * t52 + mrSges(7,2) * t53;
t64 = qJD(6) + t67;
t37 = -mrSges(7,2) * t64 + mrSges(7,3) * t52;
t40 = qJDD(6) + t43;
t15 = m(7) * (t103 * t17 - t20 * t99) - t33 * mrSges(7,3) + t40 * mrSges(7,1) - t53 * t35 + t64 * t37;
t32 = -t53 * qJD(6) + t103 * t42 - t69 * t99;
t38 = mrSges(7,1) * t64 - mrSges(7,3) * t53;
t16 = m(7) * (t103 * t20 + t17 * t99) + t32 * mrSges(7,3) - t40 * mrSges(7,2) + t52 * t35 - t64 * t38;
t57 = t67 * mrSges(6,1) - mrSges(6,2) * t82;
t114 = -t103 * t16 + t99 * t15 - m(6) * (t42 * pkin(4) + t107) + t43 * mrSges(6,3) + t67 * t57;
t56 = t66 * mrSges(6,1) + mrSges(6,3) * t82;
t128 = mrSges(5,2) * t82 - t66 * mrSges(5,3) - t56;
t132 = mrSges(5,1) - mrSges(6,2);
t55 = -mrSges(5,1) * t82 - t67 * mrSges(5,3);
t140 = m(5) * t27 + t43 * mrSges(5,2) + t128 * t66 + t132 * t42 + t67 * t55 - t114;
t131 = -mrSges(5,3) - mrSges(6,1);
t130 = t100 * t30 + t137 * t28;
t51 = -mrSges(6,2) * t66 - mrSges(6,3) * t67;
t129 = -mrSges(5,1) * t66 - mrSges(5,2) * t67 - t51;
t112 = -m(6) * t22 - t103 * t15 - t99 * t16;
t12 = m(5) * t117 - t128 * t82 + t129 * t67 + t131 * t43 + t132 * t69 + t112;
t110 = -pkin(4) * t81 + t69 * qJ(5) - t66 * t49 + t130;
t111 = -t32 * mrSges(7,1) - t52 * t37 + m(7) * (-t42 * pkin(5) - t65 * pkin(11) + (t139 - t58) * t82 + t110) + t33 * mrSges(7,2) + t53 * t38;
t109 = -m(6) * (0.2e1 * qJD(5) * t82 - t110) + t111;
t13 = m(5) * t130 + (t55 - t57) * t82 + (-mrSges(5,2) + mrSges(6,3)) * t69 + t129 * t66 + t131 * t42 + t109;
t70 = mrSges(4,1) * t90 - mrSges(4,3) * t121;
t71 = -mrSges(4,2) * t90 + mrSges(4,3) * t120;
t10 = m(4) * (-t47 * t94 + t62) + t74 * mrSges(4,2) + t75 * mrSges(4,1) + t100 * t13 + t137 * t12 + (t101 * t70 - t104 * t71) * t125;
t72 = (-mrSges(4,1) * t104 + mrSges(4,2) * t101) * t125;
t11 = m(4) * t141 + t89 * mrSges(4,1) - t74 * mrSges(4,3) - t72 * t121 + t90 * t71 - t140;
t9 = m(4) * t123 - t89 * mrSges(4,2) - t75 * mrSges(4,3) - t100 * t12 + t72 * t120 + t137 * t13 - t90 * t70;
t115 = t101 * t9 + t104 * t11;
t4 = m(3) * t118 + qJDD(2) * mrSges(3,1) - t106 * mrSges(3,2) - t94 * t10 + t115 * t97;
t6 = m(3) * t68 + t10 * t97 + t115 * t94;
t8 = m(3) * t122 - t106 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t101 * t11 + t104 * t9;
t119 = m(2) * t92 + t4 * t126 + t8 * t127 + t98 * t6;
t2 = m(2) * t85 - t102 * t4 + t105 * t8;
t1 = m(2) * t84 - t6 * t95 + (t102 * t8 + t105 * t4) * t98;
t3 = [-m(1) * g(1) - t1 * t93 + t2 * t96, t2, t8, t9, t13, -t42 * mrSges(6,2) - t66 * t56 - t114, t16; -m(1) * g(2) + t1 * t96 + t2 * t93, t1, t4, t11, t12, t42 * mrSges(6,1) - t69 * mrSges(6,3) + t66 * t51 + t82 * t57 - t109, t15; -m(1) * g(3) + t119, t119, t6, t10, t140, t43 * mrSges(6,1) + t69 * mrSges(6,2) + t67 * t51 - t82 * t56 - t112, t111;];
f_new  = t3;
