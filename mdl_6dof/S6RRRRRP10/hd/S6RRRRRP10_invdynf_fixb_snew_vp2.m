% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP10
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
% Datum: 2019-05-08 06:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:18:16
% EndTime: 2019-05-08 06:18:33
% DurationCPUTime: 5.29s
% Computational Cost: add. (98059->210), mult. (208189->280), div. (0->0), fcn. (165608->12), ass. (0->109)
t103 = sin(pkin(6));
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t130 = qJD(1) * qJD(2);
t94 = (-qJDD(1) * t112 + t108 * t130) * t103;
t105 = sin(qJ(5));
t145 = cos(qJ(5));
t106 = sin(qJ(4));
t110 = cos(qJ(4));
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t132 = qJD(1) * t112;
t104 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t126 = t109 * g(1) - t113 * g(2);
t144 = pkin(8) * t103;
t89 = qJDD(1) * pkin(1) + t114 * t144 + t126;
t136 = t104 * t89;
t123 = -t113 * g(1) - t109 * g(2);
t90 = -t114 * pkin(1) + qJDD(1) * t144 + t123;
t137 = t108 * t136 + t112 * t90;
t133 = qJD(1) * t103;
t92 = (-pkin(2) * t112 - pkin(9) * t108) * t133;
t100 = t104 * qJD(1) + qJD(2);
t98 = t100 ^ 2;
t99 = t104 * qJDD(1) + qJDD(2);
t53 = -t98 * pkin(2) + t99 * pkin(9) + (-g(3) * t108 + t92 * t132) * t103 + t137;
t143 = t104 * g(3);
t93 = (qJDD(1) * t108 + t112 * t130) * t103;
t54 = t94 * pkin(2) - t93 * pkin(9) - t143 + (-t89 + (pkin(2) * t108 - pkin(9) * t112) * t100 * qJD(1)) * t103;
t138 = t107 * t54 + t111 * t53;
t128 = t108 * t133;
t81 = t111 * t100 - t107 * t128;
t82 = t107 * t100 + t111 * t128;
t68 = -t81 * pkin(3) - t82 * pkin(10);
t86 = qJDD(3) + t94;
t127 = t103 * t132;
t97 = qJD(3) - t127;
t95 = t97 ^ 2;
t31 = -t95 * pkin(3) + t86 * pkin(10) + t81 * t68 + t138;
t134 = t103 * t112;
t122 = -g(3) * t134 - t108 * t90 + t112 * t136;
t52 = -t99 * pkin(2) - t98 * pkin(9) + t92 * t128 - t122;
t65 = -t82 * qJD(3) - t107 * t93 + t111 * t99;
t66 = t81 * qJD(3) + t107 * t99 + t111 * t93;
t34 = (-t81 * t97 - t66) * pkin(10) + (t82 * t97 - t65) * pkin(3) + t52;
t125 = -t106 * t31 + t110 * t34;
t71 = -t106 * t82 + t110 * t97;
t43 = t71 * qJD(4) + t106 * t86 + t110 * t66;
t64 = qJDD(4) - t65;
t72 = t106 * t97 + t110 * t82;
t80 = qJD(4) - t81;
t20 = (t71 * t80 - t43) * pkin(11) + (t71 * t72 + t64) * pkin(4) + t125;
t140 = t106 * t34 + t110 * t31;
t42 = -t72 * qJD(4) - t106 * t66 + t110 * t86;
t61 = t80 * pkin(4) - t72 * pkin(11);
t70 = t71 ^ 2;
t22 = -t70 * pkin(4) + t42 * pkin(11) - t80 * t61 + t140;
t121 = -t105 * t22 + t145 * t20;
t55 = t105 * t72 - t145 * t71;
t56 = t105 * t71 + t145 * t72;
t37 = t55 * pkin(5) - t56 * qJ(6);
t63 = qJDD(5) + t64;
t78 = qJD(5) + t80;
t76 = t78 ^ 2;
t146 = m(7) * (-t63 * pkin(5) - t76 * qJ(6) + t56 * t37 + qJDD(6) - t121);
t142 = -mrSges(6,3) - mrSges(7,2);
t141 = t105 * t20 + t145 * t22;
t38 = t55 * mrSges(7,1) - t56 * mrSges(7,3);
t139 = -t55 * mrSges(6,1) - t56 * mrSges(6,2) - t38;
t135 = t103 * t108;
t124 = -t107 * t53 + t111 * t54;
t30 = -t86 * pkin(3) - t95 * pkin(10) + t82 * t68 - t124;
t118 = -t42 * pkin(4) - t70 * pkin(11) + t72 * t61 + t30;
t27 = t56 * qJD(5) + t105 * t43 - t145 * t42;
t28 = -t55 * qJD(5) + t105 * t42 + t145 * t43;
t44 = -t55 * mrSges(7,2) + t78 * mrSges(7,3);
t47 = -t78 * mrSges(7,1) + t56 * mrSges(7,2);
t119 = t28 * mrSges(7,3) + t56 * t47 - m(7) * (-0.2e1 * qJD(6) * t56 + (t55 * t78 - t28) * qJ(6) + (t56 * t78 + t27) * pkin(5) + t118) - t27 * mrSges(7,1) - t55 * t44;
t45 = -t78 * mrSges(6,2) - t55 * mrSges(6,3);
t46 = t78 * mrSges(6,1) - t56 * mrSges(6,3);
t117 = m(6) * t118 + t27 * mrSges(6,1) + t28 * mrSges(6,2) + t55 * t45 + t56 * t46 - t119;
t59 = -t80 * mrSges(5,2) + t71 * mrSges(5,3);
t60 = t80 * mrSges(5,1) - t72 * mrSges(5,3);
t115 = m(5) * t30 - t42 * mrSges(5,1) + t43 * mrSges(5,2) - t71 * t59 + t72 * t60 + t117;
t67 = -t81 * mrSges(4,1) + t82 * mrSges(4,2);
t73 = -t97 * mrSges(4,2) + t81 * mrSges(4,3);
t14 = m(4) * t124 + t86 * mrSges(4,1) - t66 * mrSges(4,3) - t82 * t67 + t97 * t73 - t115;
t87 = t100 * mrSges(3,1) - mrSges(3,3) * t128;
t129 = m(7) * (-t76 * pkin(5) + t63 * qJ(6) + 0.2e1 * qJD(6) * t78 - t55 * t37 + t141) + t78 * t47 + t63 * mrSges(7,3);
t12 = m(6) * t141 - t63 * mrSges(6,2) + t139 * t55 + t142 * t27 - t78 * t46 + t129;
t13 = m(6) * t121 - t146 + (t45 + t44) * t78 + (mrSges(6,1) + mrSges(7,1)) * t63 + t139 * t56 + t142 * t28;
t57 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t10 = m(5) * t125 + t64 * mrSges(5,1) - t43 * mrSges(5,3) + t105 * t12 + t145 * t13 - t72 * t57 + t80 * t59;
t11 = m(5) * t140 - t64 * mrSges(5,2) + t42 * mrSges(5,3) - t105 * t13 + t145 * t12 + t71 * t57 - t80 * t60;
t74 = t97 * mrSges(4,1) - t82 * mrSges(4,3);
t9 = m(4) * t138 - t86 * mrSges(4,2) + t65 * mrSges(4,3) - t106 * t10 + t110 * t11 + t81 * t67 - t97 * t74;
t91 = (-mrSges(3,1) * t112 + mrSges(3,2) * t108) * t133;
t4 = m(3) * (-g(3) * t135 + t137) - t94 * mrSges(3,3) - t99 * mrSges(3,2) + t91 * t127 - t100 * t87 + t111 * t9 - t107 * t14;
t88 = -t100 * mrSges(3,2) + mrSges(3,3) * t127;
t6 = m(3) * (-t103 * t89 - t143) + t93 * mrSges(3,2) + t94 * mrSges(3,1) + t107 * t9 + t111 * t14 + (t108 * t87 - t112 * t88) * t133;
t116 = m(4) * t52 - t65 * mrSges(4,1) + t66 * mrSges(4,2) + t110 * t10 + t106 * t11 - t81 * t73 + t82 * t74;
t8 = m(3) * t122 + t99 * mrSges(3,1) - t93 * mrSges(3,3) + t100 * t88 - t91 * t128 - t116;
t131 = t104 * t6 + t8 * t134 + t4 * t135;
t2 = m(2) * t123 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t8 + t112 * t4;
t1 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t103 * t6 + (t108 * t4 + t112 * t8) * t104;
t3 = [-m(1) * g(1) - t109 * t1 + t113 * t2, t2, t4, t9, t11, t12, -t27 * mrSges(7,2) - t55 * t38 + t129; -m(1) * g(2) + t113 * t1 + t109 * t2, t1, t8, t14, t10, t13, -t119; (-m(1) - m(2)) * g(3) + t131, -m(2) * g(3) + t131, t6, t116, t115, t117, -t63 * mrSges(7,1) + t28 * mrSges(7,2) + t56 * t38 - t78 * t44 + t146;];
f_new  = t3;
