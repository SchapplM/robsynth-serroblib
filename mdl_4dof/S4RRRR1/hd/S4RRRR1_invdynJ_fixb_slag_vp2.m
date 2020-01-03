% Calculate vector of inverse dynamics joint torques for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:10
% DurationCPUTime: 1.03s
% Computational Cost: add. (1414->189), mult. (2184->265), div. (0->0), fcn. (1052->12), ass. (0->104)
t84 = sin(qJ(4));
t144 = t84 / 0.2e1;
t150 = mrSges(4,2) - mrSges(5,3);
t88 = cos(qJ(4));
t53 = -t88 * mrSges(5,1) + t84 * mrSges(5,2);
t149 = -mrSges(4,1) + t53;
t110 = mrSges(5,1) * t84 + mrSges(5,2) * t88;
t123 = pkin(1) * qJD(1);
t86 = sin(qJ(2));
t116 = t86 * t123;
t85 = sin(qJ(3));
t111 = t85 * t116;
t90 = cos(qJ(2));
t115 = t90 * t123;
t80 = qJD(1) + qJD(2);
t50 = pkin(2) * t80 + t115;
t89 = cos(qJ(3));
t24 = t50 * t89 - t111;
t74 = qJD(3) + t80;
t18 = -pkin(3) * t74 - t24;
t148 = t18 * t110 + qJD(4) * (Ifges(5,5) * t88 - Ifges(5,6) * t84) / 0.2e1;
t141 = pkin(1) * t90;
t45 = -qJD(2) * t116 + qJDD(1) * t141;
t79 = qJDD(1) + qJDD(2);
t147 = t45 * mrSges(3,1) + Ifges(3,3) * t79;
t119 = qJD(4) * t84;
t25 = t89 * t116 + t50 * t85;
t19 = pkin(7) * t74 + t25;
t73 = qJDD(3) + t79;
t120 = qJD(3) * t89;
t34 = pkin(2) * t79 + t45;
t122 = qJD(2) * t90;
t46 = (qJD(1) * t122 + qJDD(1) * t86) * pkin(1);
t8 = -qJD(3) * t111 + t50 * t120 + t85 * t34 + t89 * t46;
t5 = pkin(7) * t73 + t8;
t2 = -t19 * t119 + t5 * t88;
t138 = t2 * t88;
t118 = qJD(4) * t88;
t3 = -t19 * t118 - t5 * t84;
t146 = -t3 * t84 + t138;
t114 = (t84 ^ 2 + t88 ^ 2) * t19;
t9 = -t25 * qJD(3) + t34 * t89 - t46 * t85;
t145 = m(5) * pkin(3);
t142 = m(4) + m(5);
t83 = qJ(1) + qJ(2);
t75 = sin(t83);
t140 = g(1) * t75;
t87 = sin(qJ(1));
t139 = g(1) * t87;
t136 = Ifges(5,1) * t84;
t135 = Ifges(5,4) * t84;
t134 = Ifges(5,4) * t88;
t133 = Ifges(5,2) * t88;
t130 = t74 * t88;
t128 = t84 * mrSges(5,3);
t127 = t85 * t86;
t126 = t86 * t89;
t77 = qJ(3) + t83;
t66 = sin(t77);
t67 = cos(t77);
t124 = t67 * pkin(3) + t66 * pkin(7);
t70 = pkin(2) + t141;
t41 = pkin(1) * t126 + t85 * t70;
t121 = qJD(3) * t85;
t76 = cos(t83);
t65 = pkin(2) * t76;
t117 = t65 + t124;
t31 = t53 * t74;
t112 = t74 * mrSges(4,1) - t31;
t109 = t133 + t135;
t47 = qJD(4) * mrSges(5,1) - t74 * t128;
t48 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t130;
t106 = t47 * t88 + t48 * t84;
t105 = -t84 * t47 + t88 * t48;
t104 = t85 * t90 + t126;
t103 = t89 * t90 - t127;
t40 = -pkin(1) * t127 + t70 * t89;
t101 = t84 * (Ifges(5,1) * t88 - t135);
t100 = t74 * mrSges(4,2) - t105;
t99 = t149 * t67 + t150 * t66;
t98 = -t76 * mrSges(3,1) + t75 * mrSges(3,2) + t99;
t32 = -t74 * t119 + t73 * t88;
t33 = t74 * t118 + t73 * t84;
t97 = -t106 * qJD(4) + t88 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t32) - t84 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t33);
t96 = (t145 - t149) * t66 + (-m(5) * pkin(7) + t150) * t67;
t26 = Ifges(5,6) * qJD(4) + t109 * t74;
t54 = Ifges(5,4) * t130;
t27 = Ifges(5,5) * qJD(4) + t74 * t136 + t54;
t6 = -pkin(3) * t73 - t9;
t95 = t9 * mrSges(4,1) + mrSges(5,3) * t138 - t3 * t128 + t6 * t53 + (Ifges(5,1) * t33 + Ifges(5,4) * t32) * t144 + t88 * (Ifges(5,4) * t33 + Ifges(5,2) * t32) / 0.2e1 + t32 * t109 / 0.2e1 + t33 * (t134 + t136) / 0.2e1 - t26 * t119 / 0.2e1 + Ifges(4,3) * t73 + (t27 + t74 * (-Ifges(5,2) * t84 + t134)) * t118 / 0.2e1 + (0.2e1 * Ifges(5,5) * t144 + Ifges(5,6) * t88) * qJDD(4) + (t101 * t74 / 0.2e1 + t148) * qJD(4);
t94 = t76 * mrSges(3,2) + t96;
t93 = m(5) * t146 + t97;
t92 = -t8 * mrSges(4,2) + t95;
t91 = cos(qJ(1));
t78 = t91 * pkin(1);
t69 = -pkin(2) * t89 - pkin(3);
t68 = pkin(2) * t85 + pkin(7);
t38 = t103 * t123;
t37 = t104 * t123;
t35 = -pkin(3) - t40;
t15 = t70 * t121 + (t104 * qJD(2) + t86 * t120) * pkin(1);
t14 = t70 * t120 + (t103 * qJD(2) - t86 * t121) * pkin(1);
t13 = -mrSges(5,1) * t32 + mrSges(5,2) * t33;
t1 = [(-t14 * t74 - t41 * t73 - t8) * mrSges(4,2) + (-t15 * t74 + t40 * t73) * mrSges(4,1) + (t142 * t139 + (-t80 * t122 - t86 * t79) * mrSges(3,2) + (-qJD(2) * t86 * t80 + t90 * t79) * mrSges(3,1) + (-g(2) * t91 + t45 * t90 + t46 * t86 + t139) * m(3)) * pkin(1) + (-t91 * mrSges(2,1) + t87 * mrSges(2,2) - m(5) * (t78 + t117) - m(4) * (t65 + t78) + t98) * g(2) + m(4) * (t14 * t25 - t15 * t24 + t40 * t9 + t41 * t8) + t95 - t46 * mrSges(3,2) + t15 * t31 + t35 * t13 + t93 * (pkin(7) + t41) + (t87 * mrSges(2,1) + t91 * mrSges(2,2) + (t142 * pkin(2) + mrSges(3,1)) * t75 + t94) * g(1) + m(5) * (t14 * t114 + t15 * t18 + t35 * t6) + t105 * t14 + Ifges(2,3) * qJDD(1) + t147; ((t89 * mrSges(4,1) - t85 * mrSges(4,2)) * t73 + (-g(2) * t76 + t8 * t85 + t89 * t9 + t140) * m(4) + (-t112 * t85 - t100 * t89 + m(4) * (-t24 * t85 + t25 * t89)) * qJD(3)) * pkin(2) + t97 * t68 + (t75 * mrSges(3,1) + t94) * g(1) + (t80 * t115 - t46) * mrSges(3,2) + t100 * t38 + t112 * t37 + t98 * g(2) + t92 - m(4) * (-t24 * t37 + t25 * t38) + t69 * t13 + t80 * mrSges(3,1) * t116 + ((t140 + (t114 * t89 + t18 * t85) * qJD(3)) * pkin(2) - t38 * t114 - t18 * t37 + t6 * t69 - t117 * g(2) + t146 * t68) * m(5) + t147; t93 * pkin(7) - m(5) * (t24 * t114 + t18 * t25) + t96 * g(1) + t112 * t25 + t100 * t24 - t6 * t145 + (-m(5) * t124 + t99) * g(2) + t92 - pkin(3) * t13; t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t33 + Ifges(5,6) * t32 + Ifges(5,3) * qJDD(4) + g(3) * t53 + t106 * t19 + (t26 * t144 + (-t101 / 0.2e1 + t133 * t144) * t74 - (t27 + t54) * t88 / 0.2e1 - t148) * t74 + (g(1) * t67 + g(2) * t66) * t110;];
tau = t1;
