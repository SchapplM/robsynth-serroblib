% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 08:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:22:10
% EndTime: 2019-05-05 08:22:15
% DurationCPUTime: 2.04s
% Computational Cost: add. (26986->176), mult. (51424->224), div. (0->0), fcn. (34881->12), ass. (0->95)
t102 = cos(qJ(2));
t92 = sin(pkin(11));
t94 = cos(pkin(11));
t78 = t92 * g(1) - t94 * g(2);
t95 = cos(pkin(6));
t132 = t78 * t95;
t79 = -t94 * g(1) - t92 * g(2);
t91 = -g(3) + qJDD(1);
t93 = sin(pkin(6));
t99 = sin(qJ(2));
t140 = (t91 * t93 + t132) * t102 - t99 * t79;
t103 = qJD(3) ^ 2;
t98 = sin(qJ(3));
t123 = qJD(2) * t98;
t101 = cos(qJ(3));
t104 = qJD(2) ^ 2;
t131 = t93 * t99;
t120 = t102 * t79 + t91 * t131 + t99 * t132;
t40 = -t104 * pkin(2) + qJDD(2) * pkin(8) + t120;
t62 = -t93 * t78 + t95 * t91;
t127 = t101 * t62 - t98 * t40;
t75 = (-pkin(3) * t101 - pkin(9) * t98) * qJD(2);
t107 = qJDD(3) * pkin(3) + t103 * pkin(9) - t75 * t123 + t127;
t135 = cos(qJ(4));
t97 = sin(qJ(4));
t72 = -t135 * qJD(3) + t97 * t123;
t122 = t101 * qJD(2);
t86 = qJD(4) - t122;
t133 = t72 * t86;
t121 = qJD(2) * qJD(3);
t116 = t101 * t121;
t76 = t98 * qJDD(2) + t116;
t48 = -t72 * qJD(4) + t97 * qJDD(3) + t135 * t76;
t139 = (-t48 + t133) * qJ(5) - t107;
t100 = cos(qJ(6));
t126 = t101 * t40 + t98 * t62;
t31 = -t103 * pkin(3) + qJDD(3) * pkin(9) + t75 * t122 + t126;
t118 = t98 * t121;
t39 = -qJDD(2) * pkin(2) - t104 * pkin(8) - t140;
t77 = t101 * qJDD(2) - t118;
t33 = (-t76 - t116) * pkin(9) + (-t77 + t118) * pkin(3) + t39;
t128 = t135 * t31 + t97 * t33;
t136 = 2 * qJD(5);
t73 = t97 * qJD(3) + t135 * t123;
t52 = t72 * pkin(4) - t73 * qJ(5);
t69 = qJDD(4) - t77;
t85 = t86 ^ 2;
t110 = -t85 * pkin(4) + t69 * qJ(5) + t86 * t136 - t72 * t52 + t128;
t114 = t135 * t33 - t97 * t31;
t21 = -t69 * pkin(4) - t85 * qJ(5) + t73 * t52 + qJDD(5) - t114;
t16 = (-t48 - t133) * pkin(10) + (t72 * t73 - t69) * pkin(5) + t21;
t47 = t73 * qJD(4) - t135 * qJDD(3) + t97 * t76;
t61 = -t86 * pkin(5) - t73 * pkin(10);
t68 = t72 ^ 2;
t17 = -t68 * pkin(5) + t47 * pkin(10) + t86 * t61 + t110;
t96 = sin(qJ(6));
t49 = t100 * t72 - t96 * t73;
t27 = t49 * qJD(6) + t100 * t48 + t96 * t47;
t50 = t100 * t73 + t96 * t72;
t36 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t84 = qJD(6) - t86;
t41 = -t84 * mrSges(7,2) + t49 * mrSges(7,3);
t67 = qJDD(6) - t69;
t14 = m(7) * (t100 * t16 - t96 * t17) - t27 * mrSges(7,3) + t67 * mrSges(7,1) - t50 * t36 + t84 * t41;
t26 = -t50 * qJD(6) + t100 * t47 - t96 * t48;
t42 = t84 * mrSges(7,1) - t50 * mrSges(7,3);
t15 = m(7) * (t100 * t17 + t96 * t16) + t26 * mrSges(7,3) - t67 * mrSges(7,2) + t49 * t36 - t84 * t42;
t59 = -t86 * mrSges(6,1) + t73 * mrSges(6,2);
t111 = m(6) * t110 + t69 * mrSges(6,3) + t100 * t15 - t96 * t14 + t86 * t59;
t53 = t72 * mrSges(6,1) - t73 * mrSges(6,3);
t125 = -t72 * mrSges(5,1) - t73 * mrSges(5,2) - t53;
t129 = -mrSges(5,3) - mrSges(6,2);
t58 = t86 * mrSges(5,1) - t73 * mrSges(5,3);
t10 = m(5) * t128 - t69 * mrSges(5,2) + t125 * t72 + t129 * t47 - t86 * t58 + t111;
t108 = -m(6) * t21 - t100 * t14 - t96 * t15;
t57 = -t86 * mrSges(5,2) - t72 * mrSges(5,3);
t60 = -t72 * mrSges(6,2) + t86 * mrSges(6,3);
t11 = m(5) * t114 + (t57 + t60) * t86 + t125 * t73 + (mrSges(5,1) + mrSges(6,1)) * t69 + t129 * t48 + t108;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t123;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t122;
t138 = m(4) * t39 - t77 * mrSges(4,1) + t76 * mrSges(4,2) - (t101 * t81 - t80 * t98) * qJD(2) + t97 * t10 + t135 * t11;
t115 = m(7) * (-t68 * pkin(10) + (-pkin(4) - pkin(5)) * t47 + (-pkin(4) * t86 + t136 + t61) * t73 - t139) + t27 * mrSges(7,2) - t26 * mrSges(7,1) + t50 * t42 - t49 * t41;
t109 = m(6) * (-0.2e1 * qJD(5) * t73 + (t73 * t86 + t47) * pkin(4) + t139) + t47 * mrSges(6,1) + t72 * t60 - t115;
t137 = -m(5) * t107 + t47 * mrSges(5,1) + (t58 - t59) * t73 + (mrSges(5,2) - mrSges(6,3)) * t48 + t72 * t57 + t109;
t8 = m(3) * t140 + qJDD(2) * mrSges(3,1) - t104 * mrSges(3,2) - t138;
t134 = t102 * t8;
t74 = (-mrSges(4,1) * t101 + mrSges(4,2) * t98) * qJD(2);
t12 = m(4) * t127 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t81 - t123 * t74 - t137;
t9 = m(4) * t126 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t80 + t10 * t135 - t97 * t11 + t122 * t74;
t4 = m(3) * t120 - t104 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t101 * t9 - t98 * t12;
t6 = m(3) * t62 + t101 * t12 + t98 * t9;
t119 = m(2) * t91 + t4 * t131 + t93 * t134 + t95 * t6;
t2 = m(2) * t79 + t102 * t4 - t99 * t8;
t1 = m(2) * t78 - t93 * t6 + (t4 * t99 + t134) * t95;
t3 = [-m(1) * g(1) - t92 * t1 + t94 * t2, t2, t4, t9, t10, -t47 * mrSges(6,2) - t72 * t53 + t111, t15; -m(1) * g(2) + t94 * t1 + t92 * t2, t1, t8, t12, t11, -t48 * mrSges(6,3) - t73 * t59 + t109, t14; -m(1) * g(3) + t119, t119, t6, t138, t137, -t69 * mrSges(6,1) + t48 * mrSges(6,2) + t73 * t53 - t86 * t60 - t108, t115;];
f_new  = t3;
