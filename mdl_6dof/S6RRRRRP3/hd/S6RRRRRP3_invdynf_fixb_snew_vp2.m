% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP3
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
% Datum: 2019-05-08 04:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:37:42
% EndTime: 2019-05-08 04:37:58
% DurationCPUTime: 4.20s
% Computational Cost: add. (56308->204), mult. (112319->259), div. (0->0), fcn. (81510->10), ass. (0->98)
t104 = sin(qJ(3));
t109 = cos(qJ(3));
t105 = sin(qJ(2));
t110 = cos(qJ(2));
t132 = qJD(1) * qJD(2);
t112 = qJD(1) ^ 2;
t106 = sin(qJ(1));
t111 = cos(qJ(1));
t121 = -t111 * g(1) - t106 * g(2);
t87 = -t112 * pkin(1) + qJDD(1) * pkin(7) + t121;
t135 = t105 * t87;
t141 = pkin(2) * t112;
t91 = t105 * qJDD(1) + t110 * t132;
t58 = qJDD(2) * pkin(2) - t91 * pkin(8) - t135 + (pkin(8) * t132 + t105 * t141 - g(3)) * t110;
t101 = t110 ^ 2;
t125 = -t105 * g(3) + t110 * t87;
t92 = t110 * qJDD(1) - t105 * t132;
t134 = qJD(1) * t105;
t95 = qJD(2) * pkin(2) - pkin(8) * t134;
t59 = t92 * pkin(8) - qJD(2) * t95 - t101 * t141 + t125;
t122 = -t104 * t59 + t109 * t58;
t133 = qJD(1) * t110;
t84 = -t104 * t134 + t109 * t133;
t85 = (t104 * t110 + t105 * t109) * qJD(1);
t74 = -t84 * pkin(3) - t85 * pkin(9);
t100 = qJD(2) + qJD(3);
t98 = t100 ^ 2;
t99 = qJDD(2) + qJDD(3);
t36 = -t99 * pkin(3) - t98 * pkin(9) + t85 * t74 - t122;
t103 = sin(qJ(4));
t108 = cos(qJ(4));
t67 = t84 * qJD(3) + t104 * t92 + t109 * t91;
t77 = t103 * t100 + t108 * t85;
t43 = -t77 * qJD(4) - t103 * t67 + t108 * t99;
t83 = qJD(4) - t84;
t71 = t83 * pkin(4) - t77 * pkin(10);
t76 = t108 * t100 - t103 * t85;
t75 = t76 ^ 2;
t116 = -t43 * pkin(4) - t75 * pkin(10) + t77 * t71 + t36;
t102 = sin(qJ(5));
t107 = cos(qJ(5));
t44 = t76 * qJD(4) + t103 * t99 + t108 * t67;
t53 = t102 * t76 + t107 * t77;
t29 = -t53 * qJD(5) - t102 * t44 + t107 * t43;
t52 = -t102 * t77 + t107 * t76;
t30 = t52 * qJD(5) + t102 * t43 + t107 * t44;
t81 = qJD(5) + t83;
t47 = t81 * pkin(5) - t53 * qJ(6);
t48 = t81 * mrSges(7,1) - t53 * mrSges(7,3);
t51 = t52 ^ 2;
t129 = m(7) * (-t29 * pkin(5) - t51 * qJ(6) + t53 * t47 + qJDD(6) + t116) + t30 * mrSges(7,2) + t53 * t48;
t45 = -t81 * mrSges(7,2) + t52 * mrSges(7,3);
t46 = -t81 * mrSges(6,2) + t52 * mrSges(6,3);
t49 = t81 * mrSges(6,1) - t53 * mrSges(6,3);
t146 = m(6) * t116 + t30 * mrSges(6,2) + t53 * t49 + t129 - (t46 + t45) * t52 - (mrSges(6,1) + mrSges(7,1)) * t29;
t69 = -t83 * mrSges(5,2) + t76 * mrSges(5,3);
t70 = t83 * mrSges(5,1) - t77 * mrSges(5,3);
t145 = m(5) * t36 - t43 * mrSges(5,1) + t44 * mrSges(5,2) - t76 * t69 + t77 * t70 + t146;
t126 = t106 * g(1) - t111 * g(2);
t118 = -qJDD(1) * pkin(1) - t126;
t115 = -t92 * pkin(2) + t95 * t134 + (-pkin(8) * t101 - pkin(7)) * t112 + t118;
t66 = -t85 * qJD(3) - t104 * t91 + t109 * t92;
t33 = (-t100 * t84 - t67) * pkin(9) + (t100 * t85 - t66) * pkin(3) + t115;
t136 = t104 * t58 + t109 * t59;
t37 = -t98 * pkin(3) + t99 * pkin(9) + t84 * t74 + t136;
t123 = -t103 * t37 + t108 * t33;
t64 = qJDD(4) - t66;
t21 = (t76 * t83 - t44) * pkin(10) + (t76 * t77 + t64) * pkin(4) + t123;
t138 = t103 * t33 + t108 * t37;
t23 = -t75 * pkin(4) + t43 * pkin(10) - t83 * t71 + t138;
t124 = -t102 * t23 + t107 * t21;
t61 = qJDD(5) + t64;
t131 = m(7) * (-0.2e1 * qJD(6) * t53 + (t52 * t81 - t30) * qJ(6) + (t52 * t53 + t61) * pkin(5) + t124) + t81 * t45 + t61 * mrSges(7,1);
t40 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t41 = -t52 * mrSges(6,1) + t53 * mrSges(6,2);
t12 = m(6) * t124 + t61 * mrSges(6,1) + t81 * t46 + (-t41 - t40) * t53 + (-mrSges(6,3) - mrSges(7,3)) * t30 + t131;
t139 = t102 * t21 + t107 * t23;
t130 = m(7) * (-t51 * pkin(5) + t29 * qJ(6) + 0.2e1 * qJD(6) * t52 - t81 * t47 + t139) + t29 * mrSges(7,3) + t52 * t40;
t13 = m(6) * t139 + t29 * mrSges(6,3) + t52 * t41 + (-t49 - t48) * t81 + (-mrSges(6,2) - mrSges(7,2)) * t61 + t130;
t57 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t10 = m(5) * t123 + t64 * mrSges(5,1) - t44 * mrSges(5,3) + t102 * t13 + t107 * t12 - t77 * t57 + t83 * t69;
t11 = m(5) * t138 - t64 * mrSges(5,2) + t43 * mrSges(5,3) - t102 * t12 + t107 * t13 + t76 * t57 - t83 * t70;
t78 = -t100 * mrSges(4,2) + t84 * mrSges(4,3);
t79 = t100 * mrSges(4,1) - t85 * mrSges(4,3);
t117 = -m(4) * t115 + t66 * mrSges(4,1) - t67 * mrSges(4,2) - t108 * t10 - t103 * t11 + t84 * t78 - t85 * t79;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t143 = (t105 * t93 - t110 * t94) * qJD(1) + m(3) * (-t112 * pkin(7) + t118) - t92 * mrSges(3,1) + t91 * mrSges(3,2) - t117;
t73 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t14 = m(4) * t122 + t99 * mrSges(4,1) - t67 * mrSges(4,3) + t100 * t78 - t85 * t73 - t145;
t7 = m(4) * t136 - t99 * mrSges(4,2) + t66 * mrSges(4,3) - t103 * t10 - t100 * t79 + t108 * t11 + t84 * t73;
t90 = (-mrSges(3,1) * t110 + mrSges(3,2) * t105) * qJD(1);
t4 = m(3) * (-t110 * g(3) - t135) - t91 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t90 * t134 + qJD(2) * t94 + t104 * t7 + t109 * t14;
t5 = m(3) * t125 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t93 - t104 * t14 + t109 * t7 + t90 * t133;
t142 = t105 * t5 + t110 * t4;
t6 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t143;
t1 = m(2) * t121 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t4 + t110 * t5;
t2 = [-m(1) * g(1) + t111 * t1 - t106 * t6, t1, t5, t7, t11, t13, -t61 * mrSges(7,2) - t81 * t48 + t130; -m(1) * g(2) + t106 * t1 + t111 * t6, t6, t4, t14, t10, t12, -t30 * mrSges(7,3) - t53 * t40 + t131; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t117, t145, t146, -t29 * mrSges(7,1) - t52 * t45 + t129;];
f_new  = t2;
