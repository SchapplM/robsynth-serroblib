% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:21:50
% EndTime: 2019-05-08 04:22:04
% DurationCPUTime: 4.11s
% Computational Cost: add. (54975->205), mult. (118794->260), div. (0->0), fcn. (88305->10), ass. (0->98)
t102 = sin(qJ(2));
t107 = cos(qJ(2));
t109 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t108 = cos(qJ(1));
t123 = t103 * g(1) - t108 * g(2);
t116 = -qJDD(1) * pkin(1) - t123;
t129 = qJD(1) * t102;
t127 = qJD(1) * qJD(2);
t88 = qJDD(1) * t107 - t102 * t127;
t91 = qJD(2) * pkin(2) - pkin(8) * t129;
t98 = t107 ^ 2;
t113 = -t88 * pkin(2) + t91 * t129 + (-pkin(8) * t98 - pkin(7)) * t109 + t116;
t104 = cos(qJ(5));
t100 = sin(qJ(4));
t105 = cos(qJ(4));
t101 = sin(qJ(3));
t106 = cos(qJ(3));
t118 = -g(1) * t108 - g(2) * t103;
t84 = -pkin(1) * t109 + qJDD(1) * pkin(7) + t118;
t130 = t102 * t84;
t136 = pkin(2) * t109;
t87 = qJDD(1) * t102 + t107 * t127;
t58 = qJDD(2) * pkin(2) - t87 * pkin(8) - t130 + (pkin(8) * t127 + t102 * t136 - g(3)) * t107;
t122 = -g(3) * t102 + t107 * t84;
t59 = pkin(8) * t88 - qJD(2) * t91 - t98 * t136 + t122;
t120 = -t101 * t59 + t106 * t58;
t81 = (-t101 * t102 + t106 * t107) * qJD(1);
t66 = qJD(3) * t81 + t101 * t88 + t106 * t87;
t82 = (t101 * t107 + t102 * t106) * qJD(1);
t96 = qJDD(2) + qJDD(3);
t97 = qJD(2) + qJD(3);
t26 = (t81 * t97 - t66) * pkin(9) + (t81 * t82 + t96) * pkin(3) + t120;
t131 = t101 * t58 + t106 * t59;
t65 = -qJD(3) * t82 - t101 * t87 + t106 * t88;
t78 = pkin(3) * t97 - pkin(9) * t82;
t80 = t81 ^ 2;
t31 = -pkin(3) * t80 + t65 * pkin(9) - t78 * t97 + t131;
t133 = t100 * t26 + t105 * t31;
t73 = -t100 * t82 + t105 * t81;
t74 = t100 * t81 + t105 * t82;
t54 = -pkin(4) * t73 - pkin(10) * t74;
t94 = qJD(4) + t97;
t92 = t94 ^ 2;
t93 = qJDD(4) + t96;
t21 = -pkin(4) * t92 + pkin(10) * t93 + t54 * t73 + t133;
t111 = -t65 * pkin(3) - t80 * pkin(9) + t82 * t78 + t113;
t41 = -qJD(4) * t74 - t100 * t66 + t105 * t65;
t42 = qJD(4) * t73 + t100 * t65 + t105 * t66;
t24 = (-t73 * t94 - t42) * pkin(10) + (t74 * t94 - t41) * pkin(4) + t111;
t99 = sin(qJ(5));
t121 = t104 * t24 - t99 * t21;
t61 = t104 * t94 - t74 * t99;
t33 = t61 * qJD(5) + t104 * t42 + t93 * t99;
t39 = qJDD(5) - t41;
t70 = qJD(5) - t73;
t47 = -mrSges(7,2) * t70 + t61 * mrSges(7,3);
t62 = t104 * t74 + t94 * t99;
t126 = m(7) * (-0.2e1 * qJD(6) * t62 + (t61 * t70 - t33) * qJ(6) + (t61 * t62 + t39) * pkin(5) + t121) + t39 * mrSges(7,1) + t70 * t47;
t45 = -mrSges(7,1) * t61 + mrSges(7,2) * t62;
t46 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t48 = -mrSges(6,2) * t70 + t61 * mrSges(6,3);
t11 = m(6) * t121 + t39 * mrSges(6,1) + t70 * t48 + (-t46 - t45) * t62 + (-mrSges(6,3) - mrSges(7,3)) * t33 + t126;
t134 = t104 * t21 + t99 * t24;
t32 = -t62 * qJD(5) + t104 * t93 - t42 * t99;
t49 = pkin(5) * t70 - t62 * qJ(6);
t60 = t61 ^ 2;
t125 = m(7) * (-t60 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t61 - t49 * t70 + t134) + t32 * mrSges(7,3) + t61 * t45;
t50 = mrSges(7,1) * t70 - t62 * mrSges(7,3);
t51 = mrSges(6,1) * t70 - t62 * mrSges(6,3);
t13 = m(6) * t134 + t32 * mrSges(6,3) + t61 * t46 + (-t51 - t50) * t70 + (-mrSges(6,2) - mrSges(7,2)) * t39 + t125;
t68 = -mrSges(5,2) * t94 + mrSges(5,3) * t73;
t69 = mrSges(5,1) * t94 - mrSges(5,3) * t74;
t115 = -m(5) * t111 + t41 * mrSges(5,1) - t42 * mrSges(5,2) - t104 * t11 - t99 * t13 + t73 * t68 - t74 * t69;
t76 = -mrSges(4,2) * t97 + mrSges(4,3) * t81;
t77 = mrSges(4,1) * t97 - mrSges(4,3) * t82;
t112 = -m(4) * t113 + t65 * mrSges(4,1) - t66 * mrSges(4,2) + t81 * t76 - t82 * t77 + t115;
t89 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129;
t128 = qJD(1) * t107;
t90 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t139 = (t102 * t89 - t107 * t90) * qJD(1) + m(3) * (-t109 * pkin(7) + t116) - t88 * mrSges(3,1) + t87 * mrSges(3,2) - t112;
t119 = -t100 * t31 + t105 * t26;
t20 = -pkin(4) * t93 - pkin(10) * t92 + t74 * t54 - t119;
t124 = m(7) * (-t32 * pkin(5) - t60 * qJ(6) + t62 * t49 + qJDD(6) + t20) + t33 * mrSges(7,2) + t62 * t50;
t138 = m(6) * t20 + t33 * mrSges(6,2) - (t48 + t47) * t61 - (mrSges(6,1) + mrSges(7,1)) * t32 + t62 * t51 + t124;
t53 = -mrSges(5,1) * t73 + mrSges(5,2) * t74;
t14 = m(5) * t119 + t93 * mrSges(5,1) - t42 * mrSges(5,3) - t74 * t53 + t94 * t68 - t138;
t75 = -mrSges(4,1) * t81 + mrSges(4,2) * t82;
t9 = m(5) * t133 - t93 * mrSges(5,2) + t41 * mrSges(5,3) + t104 * t13 - t99 * t11 + t73 * t53 - t94 * t69;
t6 = m(4) * t120 + t96 * mrSges(4,1) - t66 * mrSges(4,3) + t100 * t9 + t105 * t14 - t82 * t75 + t97 * t76;
t7 = m(4) * t131 - t96 * mrSges(4,2) + t65 * mrSges(4,3) - t100 * t14 + t105 * t9 + t81 * t75 - t97 * t77;
t86 = (-mrSges(3,1) * t107 + mrSges(3,2) * t102) * qJD(1);
t4 = m(3) * (-t107 * g(3) - t130) - t87 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t86 * t129 + qJD(2) * t90 + t101 * t7 + t106 * t6;
t5 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t88 * mrSges(3,3) - qJD(2) * t89 - t101 * t6 + t106 * t7 + t86 * t128;
t137 = t102 * t5 + t107 * t4;
t8 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) - t139;
t1 = m(2) * t118 - t109 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t102 * t4 + t107 * t5;
t2 = [-m(1) * g(1) + t1 * t108 - t103 * t8, t1, t5, t7, t9, t13, -t39 * mrSges(7,2) - t70 * t50 + t125; -m(1) * g(2) + t1 * t103 + t108 * t8, t8, t4, t6, t14, t11, -t33 * mrSges(7,3) - t62 * t45 + t126; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t139, -t112, -t115, t138, -t32 * mrSges(7,1) - t61 * t47 + t124;];
f_new  = t2;
