% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:20:39
% EndTime: 2019-05-05 21:20:42
% DurationCPUTime: 1.02s
% Computational Cost: add. (10067->175), mult. (18994->205), div. (0->0), fcn. (11027->8), ass. (0->81)
t125 = cos(qJ(4));
t113 = qJD(1) * qJD(3);
t91 = cos(qJ(3));
t108 = t91 * t113;
t90 = sin(qJ(1));
t92 = cos(qJ(1));
t109 = t90 * g(1) - t92 * g(2);
t66 = qJDD(1) * pkin(1) + t109;
t104 = -t92 * g(1) - t90 * g(2);
t94 = qJD(1) ^ 2;
t68 = -t94 * pkin(1) + t104;
t86 = sin(pkin(9));
t87 = cos(pkin(9));
t106 = t87 * t66 - t86 * t68;
t31 = -qJDD(1) * pkin(2) - t94 * pkin(7) - t106;
t89 = sin(qJ(3));
t70 = t89 * qJDD(1) + t108;
t80 = t89 * t113;
t71 = t91 * qJDD(1) - t80;
t21 = (-t70 - t108) * pkin(8) + (-t71 + t80) * pkin(3) + t31;
t114 = t91 * qJD(1);
t116 = t86 * t66 + t87 * t68;
t32 = -t94 * pkin(2) + qJDD(1) * pkin(7) + t116;
t85 = -g(3) + qJDD(2);
t119 = t91 * t32 + t89 * t85;
t69 = (-pkin(3) * t91 - pkin(8) * t89) * qJD(1);
t93 = qJD(3) ^ 2;
t25 = -t93 * pkin(3) + qJDD(3) * pkin(8) + t69 * t114 + t119;
t88 = sin(qJ(4));
t121 = t125 * t25 + t88 * t21;
t126 = -2 * qJD(5);
t115 = qJD(1) * t89;
t64 = -t125 * qJD(3) + t88 * t115;
t65 = t88 * qJD(3) + t125 * t115;
t44 = t64 * pkin(4) - t65 * qJ(5);
t63 = -qJDD(4) + t71;
t76 = -qJD(4) + t114;
t75 = t76 ^ 2;
t100 = -t75 * pkin(4) - t63 * qJ(5) + t76 * t126 - t64 * t44 + t121;
t53 = t76 * mrSges(6,1) + t65 * mrSges(6,2);
t130 = m(6) * t100 - t63 * mrSges(6,3) - t76 * t53;
t124 = t64 * t76;
t41 = -t64 * qJD(4) + t88 * qJDD(3) + t125 * t70;
t120 = -t89 * t32 + t91 * t85;
t97 = qJDD(3) * pkin(3) + t93 * pkin(8) - t69 * t115 + t120;
t129 = -(t41 + t124) * qJ(5) - t97;
t40 = t65 * qJD(4) - t125 * qJDD(3) + t88 * t70;
t46 = -t64 * mrSges(7,1) + t65 * mrSges(7,2);
t50 = t76 * pkin(5) - t65 * qJ(6);
t62 = t64 ^ 2;
t110 = m(7) * (-t62 * pkin(5) + t40 * qJ(6) + 0.2e1 * qJD(6) * t64 - t76 * t50 + t100) + t64 * t46 + t40 * mrSges(7,3);
t45 = t64 * mrSges(6,1) - t65 * mrSges(6,3);
t118 = -t64 * mrSges(5,1) - t65 * mrSges(5,2) - t45;
t122 = -mrSges(5,3) - mrSges(6,2);
t51 = t76 * mrSges(7,1) - t65 * mrSges(7,3);
t52 = -t76 * mrSges(5,1) - t65 * mrSges(5,3);
t10 = m(5) * t121 + (t52 - t51) * t76 + t118 * t64 + (mrSges(5,2) - mrSges(7,2)) * t63 + t122 * t40 + t110 + t130;
t72 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t115;
t73 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t114;
t103 = t125 * t21 - t88 * t25;
t17 = t63 * pkin(4) - t75 * qJ(5) + t65 * t44 + qJDD(5) - t103;
t48 = -t76 * mrSges(7,2) + t64 * mrSges(7,3);
t111 = m(7) * (-0.2e1 * qJD(6) * t65 + (-t41 + t124) * qJ(6) + (t64 * t65 + t63) * pkin(5) + t17) + t76 * t48 + t63 * mrSges(7,1);
t101 = m(6) * t17 + t111;
t49 = t76 * mrSges(5,2) - t64 * mrSges(5,3);
t54 = -t64 * mrSges(6,2) - t76 * mrSges(6,3);
t9 = m(5) * t103 + (-t49 - t54) * t76 + (-mrSges(5,1) - mrSges(6,1)) * t63 + (t46 + t118) * t65 + (mrSges(7,3) + t122) * t41 - t101;
t128 = m(4) * t31 - t71 * mrSges(4,1) + t70 * mrSges(4,2) + (t72 * t89 - t73 * t91) * qJD(1) + t88 * t10 + t125 * t9;
t105 = m(7) * (-t62 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t40 + (pkin(4) * t76 + (2 * qJD(5)) + t50) * t65 - t129) + t41 * mrSges(7,2) - t40 * mrSges(7,1) + t65 * t51 - t64 * t48;
t99 = m(6) * (t65 * t126 + (-t65 * t76 + t40) * pkin(4) + t129) + t40 * mrSges(6,1) + t64 * t54 - t105;
t127 = -m(5) * t97 + t40 * mrSges(5,1) + (t52 - t53) * t65 + (mrSges(5,2) - mrSges(6,3)) * t41 + t64 * t49 + t99;
t67 = (-mrSges(4,1) * t91 + mrSges(4,2) * t89) * qJD(1);
t6 = m(4) * t119 - qJDD(3) * mrSges(4,2) + t71 * mrSges(4,3) - qJD(3) * t72 + t125 * t10 + t67 * t114 - t88 * t9;
t8 = m(4) * t120 + qJDD(3) * mrSges(4,1) - t70 * mrSges(4,3) + qJD(3) * t73 - t67 * t115 - t127;
t112 = m(3) * t85 + t89 * t6 + t91 * t8;
t98 = -t63 * mrSges(7,2) - t76 * t51 + t110;
t4 = m(3) * t106 + qJDD(1) * mrSges(3,1) - t94 * mrSges(3,2) - t128;
t3 = m(3) * t116 - t94 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t91 * t6 - t89 * t8;
t2 = m(2) * t104 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t87 * t3 - t86 * t4;
t1 = m(2) * t109 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) + t86 * t3 + t87 * t4;
t5 = [-m(1) * g(1) - t90 * t1 + t92 * t2, t2, t3, t6, t10, -t40 * mrSges(6,2) - t64 * t45 + t130 + t98, t98; -m(1) * g(2) + t92 * t1 + t90 * t2, t1, t4, t8, t9, -t41 * mrSges(6,3) - t65 * t53 + t99, -t41 * mrSges(7,3) - t65 * t46 + t111; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t112, t128, t127, t63 * mrSges(6,1) + t76 * t54 + (t45 - t46) * t65 + (mrSges(6,2) - mrSges(7,3)) * t41 + t101, t105;];
f_new  = t5;
