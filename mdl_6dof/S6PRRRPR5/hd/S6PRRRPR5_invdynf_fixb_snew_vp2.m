% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:00:05
% EndTime: 2019-05-05 08:00:18
% DurationCPUTime: 7.28s
% Computational Cost: add. (136732->184), mult. (289510->260), div. (0->0), fcn. (234130->16), ass. (0->108)
t101 = cos(pkin(7));
t111 = qJD(2) ^ 2;
t106 = sin(qJ(2));
t110 = cos(qJ(2));
t98 = sin(pkin(6));
t130 = t110 * t98;
t102 = cos(pkin(6));
t100 = cos(pkin(12));
t96 = sin(pkin(12));
t86 = t96 * g(1) - t100 * g(2);
t132 = t102 * t86;
t87 = -t100 * g(1) - t96 * g(2);
t94 = -g(3) + qJDD(1);
t120 = -t106 * t87 + t110 * t132 + t94 * t130;
t97 = sin(pkin(7));
t136 = pkin(9) * t97;
t61 = qJDD(2) * pkin(2) + t111 * t136 + t120;
t73 = t102 * t94 - t98 * t86;
t139 = t101 * t61 + t73 * t97;
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t128 = qJD(2) * qJD(3);
t80 = (-qJDD(2) * t109 + t105 * t128) * t97;
t131 = t106 * t98;
t125 = t106 * t132 + t110 * t87 + t94 * t131;
t62 = -t111 * pkin(2) + qJDD(2) * t136 + t125;
t138 = -t105 * t62 + t139 * t109;
t137 = 2 * qJD(5);
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t129 = qJD(2) * t97;
t123 = t109 * t129;
t126 = t139 * t105 + t109 * t62;
t78 = (-pkin(3) * t109 - pkin(10) * t105) * t129;
t92 = t101 * qJD(2) + qJD(3);
t90 = t92 ^ 2;
t91 = t101 * qJDD(2) + qJDD(3);
t33 = -t90 * pkin(3) + t91 * pkin(10) + t78 * t123 + t126;
t69 = t101 * t73;
t79 = (qJDD(2) * t105 + t109 * t128) * t97;
t36 = t80 * pkin(3) - t79 * pkin(10) + t69 + (-t61 + (pkin(3) * t105 - pkin(10) * t109) * t92 * qJD(2)) * t97;
t134 = t104 * t36 + t108 * t33;
t121 = -t104 * t33 + t108 * t36;
t124 = t105 * t129;
t71 = -t104 * t124 + t108 * t92;
t53 = t71 * qJD(4) + t104 * t91 + t108 * t79;
t72 = t104 * t92 + t108 * t124;
t74 = qJDD(4) + t80;
t85 = qJD(4) - t123;
t24 = (t71 * t85 - t53) * qJ(5) + (t71 * t72 + t74) * pkin(4) + t121;
t52 = -t72 * qJD(4) - t104 * t79 + t108 * t91;
t65 = t85 * pkin(4) - t72 * qJ(5);
t70 = t71 ^ 2;
t26 = -t70 * pkin(4) + t52 * qJ(5) - t85 * t65 + t134;
t95 = sin(pkin(13));
t99 = cos(pkin(13));
t59 = t99 * t71 - t95 * t72;
t127 = t59 * t137 + t95 * t24 + t99 * t26;
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t60 = t95 * t71 + t99 * t72;
t46 = -t59 * pkin(5) - t60 * pkin(11);
t84 = t85 ^ 2;
t21 = -t84 * pkin(5) + t74 * pkin(11) + t59 * t46 + t127;
t32 = -t91 * pkin(3) - t90 * pkin(10) + t78 * t124 - t138;
t113 = -t52 * pkin(4) - t70 * qJ(5) + t72 * t65 + qJDD(5) + t32;
t42 = t99 * t52 - t95 * t53;
t43 = t95 * t52 + t99 * t53;
t22 = (-t59 * t85 - t43) * pkin(11) + (t60 * t85 - t42) * pkin(5) + t113;
t47 = -t103 * t60 + t107 * t85;
t30 = t47 * qJD(6) + t103 * t74 + t107 * t43;
t48 = t103 * t85 + t107 * t60;
t37 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t56 = qJD(6) - t59;
t38 = -t56 * mrSges(7,2) + t47 * mrSges(7,3);
t41 = qJDD(6) - t42;
t18 = m(7) * (-t103 * t21 + t107 * t22) - t30 * mrSges(7,3) + t41 * mrSges(7,1) - t48 * t37 + t56 * t38;
t29 = -t48 * qJD(6) - t103 * t43 + t107 * t74;
t39 = t56 * mrSges(7,1) - t48 * mrSges(7,3);
t19 = m(7) * (t103 * t22 + t107 * t21) + t29 * mrSges(7,3) - t41 * mrSges(7,2) + t47 * t37 - t56 * t39;
t45 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t50 = t85 * mrSges(6,1) - t60 * mrSges(6,3);
t14 = m(6) * t127 - t74 * mrSges(6,2) + t42 * mrSges(6,3) - t103 * t18 + t107 * t19 + t59 * t45 - t85 * t50;
t119 = -t99 * t24 + t95 * t26;
t114 = m(7) * (-t74 * pkin(5) - t84 * pkin(11) + (t137 + t46) * t60 + t119) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t38 + t48 * t39;
t49 = -t85 * mrSges(6,2) + t59 * mrSges(6,3);
t15 = m(6) * (-0.2e1 * qJD(5) * t60 - t119) - t43 * mrSges(6,3) + t74 * mrSges(6,1) - t60 * t45 + t85 * t49 - t114;
t63 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t64 = -t85 * mrSges(5,2) + t71 * mrSges(5,3);
t11 = m(5) * t121 + t74 * mrSges(5,1) - t53 * mrSges(5,3) + t95 * t14 + t99 * t15 - t72 * t63 + t85 * t64;
t66 = t85 * mrSges(5,1) - t72 * mrSges(5,3);
t12 = m(5) * t134 - t74 * mrSges(5,2) + t52 * mrSges(5,3) + t99 * t14 - t95 * t15 + t71 * t63 - t85 * t66;
t75 = t92 * mrSges(4,1) - mrSges(4,3) * t124;
t76 = -t92 * mrSges(4,2) + mrSges(4,3) * t123;
t10 = m(4) * (-t97 * t61 + t69) + t79 * mrSges(4,2) + t80 * mrSges(4,1) + t104 * t12 + t108 * t11 + (t105 * t75 - t109 * t76) * t129;
t115 = -m(6) * t113 + t42 * mrSges(6,1) - t43 * mrSges(6,2) - t103 * t19 - t107 * t18 + t59 * t49 - t60 * t50;
t112 = m(5) * t32 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t71 * t64 + t72 * t66 - t115;
t77 = (-mrSges(4,1) * t109 + mrSges(4,2) * t105) * t129;
t13 = m(4) * t138 + t91 * mrSges(4,1) - t79 * mrSges(4,3) - t77 * t124 + t92 * t76 - t112;
t9 = m(4) * t126 - t91 * mrSges(4,2) - t80 * mrSges(4,3) - t104 * t11 + t108 * t12 + t77 * t123 - t92 * t75;
t117 = t105 * t9 + t109 * t13;
t4 = m(3) * t120 + qJDD(2) * mrSges(3,1) - t111 * mrSges(3,2) - t97 * t10 + t117 * t101;
t6 = m(3) * t73 + t101 * t10 + t117 * t97;
t8 = m(3) * t125 - t111 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t105 * t13 + t109 * t9;
t122 = m(2) * t94 + t102 * t6 + t4 * t130 + t8 * t131;
t2 = m(2) * t87 - t106 * t4 + t110 * t8;
t1 = m(2) * t86 - t98 * t6 + (t106 * t8 + t110 * t4) * t102;
t3 = [-m(1) * g(1) - t96 * t1 + t100 * t2, t2, t8, t9, t12, t14, t19; -m(1) * g(2) + t100 * t1 + t96 * t2, t1, t4, t13, t11, t15, t18; -m(1) * g(3) + t122, t122, t6, t10, t112, -t115, t114;];
f_new  = t3;
