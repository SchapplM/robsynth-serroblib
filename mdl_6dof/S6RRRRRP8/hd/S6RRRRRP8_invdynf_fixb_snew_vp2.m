% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:37:48
% EndTime: 2019-05-08 05:38:05
% DurationCPUTime: 5.26s
% Computational Cost: add. (93394->212), mult. (198717->279), div. (0->0), fcn. (159964->12), ass. (0->108)
t106 = sin(qJ(4));
t110 = cos(qJ(4));
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t103 = sin(pkin(6));
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t131 = qJD(1) * t112;
t104 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t124 = t109 * g(1) - g(2) * t113;
t87 = pkin(8) * t103 * t114 + qJDD(1) * pkin(1) + t124;
t135 = t104 * t87;
t121 = -g(1) * t113 - g(2) * t109;
t129 = qJDD(1) * t103;
t88 = -pkin(1) * t114 + pkin(8) * t129 + t121;
t136 = t108 * t135 + t112 * t88;
t132 = qJD(1) * t103;
t90 = (-pkin(2) * t112 - pkin(9) * t108) * t132;
t100 = qJD(1) * t104 + qJD(2);
t98 = t100 ^ 2;
t99 = qJDD(1) * t104 + qJDD(2);
t58 = -t98 * pkin(2) + t99 * pkin(9) + (-g(3) * t108 + t90 * t131) * t103 + t136;
t144 = t104 * g(3);
t91 = (qJD(2) * t131 + qJDD(1) * t108) * t103;
t126 = t108 * t132;
t92 = -qJD(2) * t126 + t112 * t129;
t59 = -t92 * pkin(2) - t91 * pkin(9) - t144 + (-t87 + (pkin(2) * t108 - pkin(9) * t112) * t100 * qJD(1)) * t103;
t123 = -t107 * t58 + t111 * t59;
t79 = t100 * t111 - t107 * t126;
t66 = qJD(3) * t79 + t107 * t99 + t111 * t91;
t80 = t100 * t107 + t111 * t126;
t84 = qJDD(3) - t92;
t125 = t103 * t131;
t96 = qJD(3) - t125;
t27 = (t79 * t96 - t66) * pkin(10) + (t79 * t80 + t84) * pkin(3) + t123;
t137 = t107 * t59 + t111 * t58;
t65 = -qJD(3) * t80 - t107 * t91 + t111 * t99;
t74 = pkin(3) * t96 - pkin(10) * t80;
t78 = t79 ^ 2;
t30 = -pkin(3) * t78 + pkin(10) * t65 - t74 * t96 + t137;
t122 = -t106 * t30 + t110 * t27;
t69 = -t106 * t80 + t110 * t79;
t70 = t106 * t79 + t110 * t80;
t53 = -pkin(4) * t69 - pkin(11) * t70;
t83 = qJDD(4) + t84;
t95 = qJD(4) + t96;
t94 = t95 ^ 2;
t22 = -t83 * pkin(4) - t94 * pkin(11) + t70 * t53 - t122;
t105 = sin(qJ(5));
t145 = cos(qJ(5));
t42 = qJD(4) * t69 + t106 * t65 + t110 * t66;
t61 = t105 * t95 + t145 * t70;
t32 = qJD(5) * t61 + t105 * t42 - t145 * t83;
t60 = t105 * t70 - t145 * t95;
t33 = -t60 * qJD(5) + t105 * t83 + t145 * t42;
t68 = qJD(5) - t69;
t47 = -mrSges(7,2) * t60 + mrSges(7,3) * t68;
t127 = m(7) * (-0.2e1 * qJD(6) * t61 + (t60 * t68 - t33) * qJ(6) + (t61 * t68 + t32) * pkin(5) + t22) + t32 * mrSges(7,1) + t60 * t47;
t48 = -mrSges(6,2) * t68 - mrSges(6,3) * t60;
t49 = mrSges(6,1) * t68 - mrSges(6,3) * t61;
t50 = -mrSges(7,1) * t68 + mrSges(7,2) * t61;
t147 = m(6) * t22 + t32 * mrSges(6,1) + (t49 - t50) * t61 + (mrSges(6,2) - mrSges(7,3)) * t33 + t60 * t48 + t127;
t140 = t106 * t27 + t110 * t30;
t23 = -pkin(4) * t94 + pkin(11) * t83 + t53 * t69 + t140;
t133 = t103 * t112;
t120 = -g(3) * t133 - t108 * t88 + t112 * t135;
t57 = -t99 * pkin(2) - t98 * pkin(9) + t90 * t126 - t120;
t116 = -t65 * pkin(3) - t78 * pkin(10) + t80 * t74 + t57;
t41 = -qJD(4) * t70 - t106 * t66 + t110 * t65;
t25 = (-t69 * t95 - t42) * pkin(11) + (t70 * t95 - t41) * pkin(4) + t116;
t119 = -t105 * t23 + t145 * t25;
t40 = qJDD(5) - t41;
t44 = pkin(5) * t60 - qJ(6) * t61;
t67 = t68 ^ 2;
t146 = m(7) * (-t40 * pkin(5) - t67 * qJ(6) + t61 * t44 + qJDD(6) - t119);
t142 = -mrSges(6,3) - mrSges(7,2);
t141 = t105 * t25 + t145 * t23;
t45 = mrSges(7,1) * t60 - mrSges(7,3) * t61;
t139 = -mrSges(6,1) * t60 - mrSges(6,2) * t61 - t45;
t134 = t103 * t108;
t128 = m(7) * (-pkin(5) * t67 + qJ(6) * t40 + 0.2e1 * qJD(6) * t68 - t44 * t60 + t141) + t68 * t50 + t40 * mrSges(7,3);
t14 = m(6) * t141 - t40 * mrSges(6,2) + t139 * t60 + t142 * t32 - t68 * t49 + t128;
t16 = m(6) * t119 - t146 + (t48 + t47) * t68 + t139 * t61 + (mrSges(6,1) + mrSges(7,1)) * t40 + t142 * t33;
t62 = -mrSges(5,2) * t95 + mrSges(5,3) * t69;
t63 = mrSges(5,1) * t95 - mrSges(5,3) * t70;
t118 = -m(5) * t116 + t41 * mrSges(5,1) - t42 * mrSges(5,2) - t105 * t14 - t145 * t16 + t69 * t62 - t70 * t63;
t72 = -mrSges(4,2) * t96 + mrSges(4,3) * t79;
t73 = mrSges(4,1) * t96 - mrSges(4,3) * t80;
t115 = m(4) * t57 - t65 * mrSges(4,1) + t66 * mrSges(4,2) - t79 * t72 + t80 * t73 - t118;
t86 = -mrSges(3,2) * t100 + mrSges(3,3) * t125;
t89 = (-mrSges(3,1) * t112 + mrSges(3,2) * t108) * t132;
t10 = m(3) * t120 + t99 * mrSges(3,1) - t91 * mrSges(3,3) + t100 * t86 - t89 * t126 - t115;
t52 = -mrSges(5,1) * t69 + mrSges(5,2) * t70;
t11 = m(5) * t140 - t83 * mrSges(5,2) + t41 * mrSges(5,3) - t105 * t16 + t145 * t14 + t69 * t52 - t95 * t63;
t12 = m(5) * t122 + t83 * mrSges(5,1) - t42 * mrSges(5,3) - t70 * t52 + t95 * t62 - t147;
t71 = -mrSges(4,1) * t79 + mrSges(4,2) * t80;
t7 = m(4) * t123 + t84 * mrSges(4,1) - t66 * mrSges(4,3) + t106 * t11 + t110 * t12 - t80 * t71 + t96 * t72;
t8 = m(4) * t137 - t84 * mrSges(4,2) + t65 * mrSges(4,3) - t106 * t12 + t110 * t11 + t79 * t71 - t96 * t73;
t85 = mrSges(3,1) * t100 - mrSges(3,3) * t126;
t4 = m(3) * (-g(3) * t134 + t136) + t92 * mrSges(3,3) - t99 * mrSges(3,2) + t89 * t125 - t100 * t85 + t111 * t8 - t107 * t7;
t6 = m(3) * (-t103 * t87 - t144) + t91 * mrSges(3,2) - t92 * mrSges(3,1) + t107 * t8 + t111 * t7 + (t108 * t85 - t112 * t86) * t132;
t130 = t10 * t133 + t104 * t6 + t4 * t134;
t2 = m(2) * t121 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t10 + t112 * t4;
t1 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t103 * t6 + (t10 * t112 + t108 * t4) * t104;
t3 = [-m(1) * g(1) - t1 * t109 + t113 * t2, t2, t4, t8, t11, t14, -t32 * mrSges(7,2) - t60 * t45 + t128; -m(1) * g(2) + t1 * t113 + t109 * t2, t1, t10, t7, t12, t16, -t33 * mrSges(7,3) - t61 * t50 + t127; (-m(1) - m(2)) * g(3) + t130, -m(2) * g(3) + t130, t6, t115, -t118, t147, -t40 * mrSges(7,1) + t33 * mrSges(7,2) + t61 * t45 - t68 * t47 + t146;];
f_new  = t3;
