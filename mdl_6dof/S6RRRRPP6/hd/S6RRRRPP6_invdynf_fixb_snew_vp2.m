% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:31:51
% EndTime: 2019-05-07 18:31:56
% DurationCPUTime: 2.04s
% Computational Cost: add. (22434->203), mult. (44298->238), div. (0->0), fcn. (29978->8), ass. (0->94)
t102 = qJD(2) ^ 2;
t97 = sin(qJ(2));
t127 = qJD(1) * t97;
t100 = cos(qJ(2));
t103 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t115 = -t101 * g(1) - t98 * g(2);
t76 = -t103 * pkin(1) + qJDD(1) * pkin(7) + t115;
t128 = -t100 * g(3) - t97 * t76;
t83 = (-pkin(2) * t100 - pkin(8) * t97) * qJD(1);
t48 = -qJDD(2) * pkin(2) - t102 * pkin(8) + t83 * t127 - t128;
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t81 = t96 * qJD(2) + t99 * t127;
t125 = qJD(1) * qJD(2);
t117 = t100 * t125;
t84 = t97 * qJDD(1) + t117;
t61 = -t81 * qJD(3) + t99 * qJDD(2) - t96 * t84;
t126 = t100 * qJD(1);
t91 = qJD(3) - t126;
t69 = t91 * pkin(3) - t81 * pkin(9);
t80 = t99 * qJD(2) - t96 * t127;
t78 = t80 ^ 2;
t106 = -t61 * pkin(3) - t78 * pkin(9) + t81 * t69 + t48;
t137 = cos(qJ(4));
t95 = sin(qJ(4));
t64 = -t137 * t80 + t95 * t81;
t89 = -qJD(4) - t91;
t136 = t64 * t89;
t140 = -2 * qJD(5);
t62 = t80 * qJD(3) + t96 * qJDD(2) + t99 * t84;
t33 = -t64 * qJD(4) + t137 * t62 + t95 * t61;
t65 = t137 * t81 + t95 * t80;
t105 = (-t33 - t136) * qJ(5) + t106 + (-t89 * pkin(4) + t140) * t65;
t139 = 2 * qJD(6);
t32 = t65 * qJD(4) - t137 * t61 + t95 * t62;
t50 = t65 * pkin(5) + t89 * qJ(6);
t51 = t65 * mrSges(7,1) + t89 * mrSges(7,3);
t53 = -t64 * mrSges(7,1) - t89 * mrSges(7,2);
t63 = t64 ^ 2;
t111 = t33 * mrSges(7,2) + t65 * t51 - m(7) * (t105 + (pkin(4) + qJ(6)) * t32 + t64 * t139 - t65 * t50 - t63 * pkin(5)) - t32 * mrSges(7,3) - t64 * t53;
t54 = t65 * mrSges(6,1) - t89 * mrSges(6,2);
t109 = t111 - m(6) * (t32 * pkin(4) + t105) + t33 * mrSges(6,3) + t65 * t54;
t52 = t64 * mrSges(6,1) + t89 * mrSges(6,3);
t129 = -t89 * mrSges(5,2) + t64 * mrSges(5,3) + t52;
t56 = -t89 * mrSges(5,1) - t65 * mrSges(5,3);
t144 = m(5) * t106 + t33 * mrSges(5,2) + t65 * t56 - t109 - t129 * t64 + (mrSges(5,1) - mrSges(6,2)) * t32;
t67 = -t91 * mrSges(4,2) + t80 * mrSges(4,3);
t68 = t91 * mrSges(4,1) - t81 * mrSges(4,3);
t143 = m(4) * t48 - t61 * mrSges(4,1) + t62 * mrSges(4,2) - t80 * t67 + t81 * t68 + t144;
t122 = t98 * g(1) - t101 * g(2);
t75 = -qJDD(1) * pkin(1) - t103 * pkin(7) - t122;
t92 = t97 * t125;
t85 = t100 * qJDD(1) - t92;
t45 = (-t84 - t117) * pkin(8) + (-t85 + t92) * pkin(2) + t75;
t121 = -t97 * g(3) + t100 * t76;
t49 = -t102 * pkin(2) + qJDD(2) * pkin(8) + t83 * t126 + t121;
t118 = t99 * t45 - t96 * t49;
t79 = qJDD(3) - t85;
t22 = (t80 * t91 - t62) * pkin(9) + (t80 * t81 + t79) * pkin(3) + t118;
t130 = t96 * t45 + t99 * t49;
t25 = -t78 * pkin(3) + t61 * pkin(9) - t91 * t69 + t130;
t132 = t137 * t25 + t95 * t22;
t40 = t64 * pkin(4) - t65 * qJ(5);
t77 = qJDD(4) + t79;
t88 = t89 ^ 2;
t110 = -t88 * pkin(4) + t77 * qJ(5) - t64 * t40 + t132;
t123 = t89 * t51 - m(7) * (-t32 * pkin(5) - t63 * qJ(6) + qJDD(6) + (t140 - t50) * t89 + t110) - t77 * mrSges(7,2);
t112 = m(6) * (0.2e1 * qJD(5) * t89 - t110) + t123;
t42 = -t64 * mrSges(6,2) - t65 * mrSges(6,3);
t131 = -t64 * mrSges(5,1) - t65 * mrSges(5,2) - t42;
t133 = -mrSges(5,3) - mrSges(6,1);
t39 = -t65 * mrSges(7,2) + t64 * mrSges(7,3);
t10 = m(5) * t132 + (t56 - t54) * t89 + (-mrSges(5,2) + mrSges(6,3)) * t77 + (-t39 + t131) * t64 + (-mrSges(7,1) + t133) * t32 - t112;
t66 = -t80 * mrSges(4,1) + t81 * mrSges(4,2);
t116 = t137 * t22 - t95 * t25;
t18 = -t77 * pkin(4) - t88 * qJ(5) + t65 * t40 + qJDD(5) - t116;
t124 = m(7) * (t89 * t139 + (t64 * t65 - t77) * qJ(6) + (t33 - t136) * pkin(5) + t18) + t33 * mrSges(7,1) + t65 * t39;
t113 = m(6) * t18 + t124;
t134 = mrSges(6,2) - mrSges(7,3);
t9 = m(5) * t116 + t131 * t65 + t133 * t33 + (-t53 + t129) * t89 + (mrSges(5,1) - t134) * t77 - t113;
t5 = m(4) * t118 + t79 * mrSges(4,1) - t62 * mrSges(4,3) + t95 * t10 + t137 * t9 - t81 * t66 + t91 * t67;
t6 = m(4) * t130 - t79 * mrSges(4,2) + t61 * mrSges(4,3) + t137 * t10 + t80 * t66 - t91 * t68 - t95 * t9;
t86 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t127;
t87 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t126;
t141 = m(3) * t75 - t85 * mrSges(3,1) + t84 * mrSges(3,2) - (t100 * t87 - t97 * t86) * qJD(1) + t99 * t5 + t96 * t6;
t82 = (-mrSges(3,1) * t100 + mrSges(3,2) * t97) * qJD(1);
t4 = m(3) * t121 - qJDD(2) * mrSges(3,2) + t85 * mrSges(3,3) - qJD(2) * t86 + t82 * t126 - t96 * t5 + t99 * t6;
t8 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t84 * mrSges(3,3) + qJD(2) * t87 - t82 * t127 - t143;
t138 = t100 * t8 + t97 * t4;
t2 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t103 * mrSges(2,2) - t141;
t1 = m(2) * t115 - t103 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t100 * t4 - t97 * t8;
t3 = [-m(1) * g(1) + t101 * t1 - t98 * t2, t1, t4, t6, t10, -t32 * mrSges(6,2) - t64 * t52 - t109, -t111; -m(1) * g(2) + t98 * t1 + t101 * t2, t2, t8, t5, t9, -t77 * mrSges(6,3) + t89 * t54 + (t39 + t42) * t64 + (mrSges(6,1) + mrSges(7,1)) * t32 + t112, -t77 * mrSges(7,3) + t89 * t53 + t124; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t141, t143, t144, t33 * mrSges(6,1) + t65 * t42 + (-t52 + t53) * t89 + t134 * t77 + t113, -t32 * mrSges(7,1) - t64 * t39 - t123;];
f_new  = t3;
