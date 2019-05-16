% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRR4
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
% Datum: 2019-05-05 11:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:23:57
% EndTime: 2019-05-05 11:24:12
% DurationCPUTime: 7.62s
% Computational Cost: add. (145071->184), mult. (299662->258), div. (0->0), fcn. (243946->16), ass. (0->109)
t101 = cos(pkin(7));
t113 = qJD(2) ^ 2;
t107 = sin(qJ(2));
t112 = cos(qJ(2));
t99 = sin(pkin(6));
t131 = t112 * t99;
t102 = cos(pkin(6));
t100 = cos(pkin(13));
t97 = sin(pkin(13));
t87 = g(1) * t97 - g(2) * t100;
t133 = t102 * t87;
t88 = -g(1) * t100 - g(2) * t97;
t96 = -g(3) + qJDD(1);
t121 = -t107 * t88 + t112 * t133 + t96 * t131;
t98 = sin(pkin(7));
t60 = pkin(9) * t113 * t98 + qJDD(2) * pkin(2) + t121;
t72 = t102 * t96 - t87 * t99;
t139 = t101 * t60 + t72 * t98;
t106 = sin(qJ(3));
t111 = cos(qJ(3));
t132 = t107 * t99;
t126 = t107 * t133 + t112 * t88 + t96 * t132;
t128 = qJDD(2) * t98;
t61 = -pkin(2) * t113 + pkin(9) * t128 + t126;
t138 = -t106 * t61 + t139 * t111;
t104 = sin(qJ(5));
t109 = cos(qJ(5));
t105 = sin(qJ(4));
t110 = cos(qJ(4));
t129 = qJD(2) * t111;
t124 = t98 * t129;
t127 = t139 * t106 + t111 * t61;
t130 = qJD(2) * t98;
t78 = (-pkin(3) * t111 - pkin(10) * t106) * t130;
t94 = qJD(2) * t101 + qJD(3);
t92 = t94 ^ 2;
t93 = qJDD(2) * t101 + qJDD(3);
t33 = -pkin(3) * t92 + pkin(10) * t93 + t78 * t124 + t127;
t68 = t101 * t72;
t79 = (qJD(3) * t129 + qJDD(2) * t106) * t98;
t125 = t106 * t130;
t80 = -qJD(3) * t125 + t111 * t128;
t40 = -t80 * pkin(3) - t79 * pkin(10) + t68 + (-t60 + (pkin(3) * t106 - pkin(10) * t111) * t94 * qJD(2)) * t98;
t122 = -t105 * t33 + t110 * t40;
t70 = -t105 * t125 + t110 * t94;
t53 = qJD(4) * t70 + t105 * t93 + t110 * t79;
t71 = t105 * t94 + t110 * t125;
t74 = qJDD(4) - t80;
t86 = qJD(4) - t124;
t24 = (t70 * t86 - t53) * pkin(11) + (t70 * t71 + t74) * pkin(4) + t122;
t135 = t105 * t40 + t110 * t33;
t52 = -qJD(4) * t71 - t105 * t79 + t110 * t93;
t65 = pkin(4) * t86 - pkin(11) * t71;
t69 = t70 ^ 2;
t26 = -t69 * pkin(4) + t52 * pkin(11) - t65 * t86 + t135;
t136 = t104 * t24 + t109 * t26;
t103 = sin(qJ(6));
t108 = cos(qJ(6));
t58 = -t104 * t71 + t109 * t70;
t59 = t104 * t70 + t109 * t71;
t46 = -pkin(5) * t58 - pkin(12) * t59;
t73 = qJDD(5) + t74;
t85 = qJD(5) + t86;
t84 = t85 ^ 2;
t21 = -pkin(5) * t84 + pkin(12) * t73 + t58 * t46 + t136;
t32 = -t93 * pkin(3) - t92 * pkin(10) + t78 * t125 - t138;
t115 = -t52 * pkin(4) - t69 * pkin(11) + t71 * t65 + t32;
t36 = -t59 * qJD(5) - t104 * t53 + t109 * t52;
t37 = t58 * qJD(5) + t104 * t52 + t109 * t53;
t22 = (-t58 * t85 - t37) * pkin(12) + (t59 * t85 - t36) * pkin(5) + t115;
t47 = -t103 * t59 + t108 * t85;
t30 = t47 * qJD(6) + t103 * t73 + t108 * t37;
t35 = qJDD(6) - t36;
t48 = t103 * t85 + t108 * t59;
t41 = -mrSges(7,1) * t47 + mrSges(7,2) * t48;
t55 = qJD(6) - t58;
t42 = -mrSges(7,2) * t55 + mrSges(7,3) * t47;
t18 = m(7) * (-t103 * t21 + t108 * t22) - t30 * mrSges(7,3) + t35 * mrSges(7,1) - t48 * t41 + t55 * t42;
t29 = -t48 * qJD(6) - t103 * t37 + t108 * t73;
t43 = mrSges(7,1) * t55 - mrSges(7,3) * t48;
t19 = m(7) * (t103 * t22 + t108 * t21) + t29 * mrSges(7,3) - t35 * mrSges(7,2) + t47 * t41 - t55 * t43;
t45 = -mrSges(6,1) * t58 + mrSges(6,2) * t59;
t50 = mrSges(6,1) * t85 - t59 * mrSges(6,3);
t14 = m(6) * t136 - t73 * mrSges(6,2) + t36 * mrSges(6,3) - t103 * t18 + t108 * t19 + t58 * t45 - t85 * t50;
t118 = -t104 * t26 + t109 * t24;
t116 = m(7) * (-pkin(5) * t73 - pkin(12) * t84 + t59 * t46 - t118) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t42 + t48 * t43;
t49 = -mrSges(6,2) * t85 + t58 * mrSges(6,3);
t15 = m(6) * t118 + t73 * mrSges(6,1) - t37 * mrSges(6,3) - t59 * t45 + t85 * t49 - t116;
t62 = -mrSges(5,1) * t70 + mrSges(5,2) * t71;
t63 = -mrSges(5,2) * t86 + mrSges(5,3) * t70;
t11 = m(5) * t122 + t74 * mrSges(5,1) - t53 * mrSges(5,3) + t104 * t14 + t109 * t15 - t71 * t62 + t86 * t63;
t64 = mrSges(5,1) * t86 - mrSges(5,3) * t71;
t12 = m(5) * t135 - t74 * mrSges(5,2) + t52 * mrSges(5,3) - t104 * t15 + t109 * t14 + t70 * t62 - t86 * t64;
t75 = mrSges(4,1) * t94 - mrSges(4,3) * t125;
t76 = -mrSges(4,2) * t94 + mrSges(4,3) * t124;
t10 = m(4) * (-t98 * t60 + t68) + t79 * mrSges(4,2) - t80 * mrSges(4,1) + t105 * t12 + t110 * t11 + (t106 * t75 - t111 * t76) * t130;
t117 = -m(6) * t115 + t36 * mrSges(6,1) - t37 * mrSges(6,2) - t103 * t19 - t108 * t18 + t58 * t49 - t59 * t50;
t114 = m(5) * t32 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t70 * t63 + t71 * t64 - t117;
t77 = (-mrSges(4,1) * t111 + mrSges(4,2) * t106) * t130;
t13 = m(4) * t138 + t93 * mrSges(4,1) - t79 * mrSges(4,3) - t77 * t125 + t94 * t76 - t114;
t9 = m(4) * t127 - t93 * mrSges(4,2) + t80 * mrSges(4,3) - t105 * t11 + t110 * t12 + t77 * t124 - t94 * t75;
t119 = t106 * t9 + t111 * t13;
t4 = m(3) * t121 + qJDD(2) * mrSges(3,1) - t113 * mrSges(3,2) - t98 * t10 + t119 * t101;
t6 = m(3) * t72 + t101 * t10 + t119 * t98;
t8 = m(3) * t126 - t113 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t106 * t13 + t111 * t9;
t123 = m(2) * t96 + t102 * t6 + t4 * t131 + t8 * t132;
t2 = m(2) * t88 - t107 * t4 + t112 * t8;
t1 = m(2) * t87 - t99 * t6 + (t107 * t8 + t112 * t4) * t102;
t3 = [-m(1) * g(1) - t1 * t97 + t100 * t2, t2, t8, t9, t12, t14, t19; -m(1) * g(2) + t1 * t100 + t2 * t97, t1, t4, t13, t11, t15, t18; -m(1) * g(3) + t123, t123, t6, t10, t114, -t117, t116;];
f_new  = t3;
