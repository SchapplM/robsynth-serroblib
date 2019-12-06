% Calculate vector of inverse dynamics joint torques for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:40
% DurationCPUTime: 1.54s
% Computational Cost: add. (1574->216), mult. (3235->303), div. (0->0), fcn. (2442->14), ass. (0->108)
t88 = sin(qJ(5));
t155 = mrSges(6,2) * t88;
t82 = pkin(9) + qJ(3);
t80 = qJ(4) + t82;
t73 = sin(t80);
t74 = cos(t80);
t179 = -t73 * t155 + t74 * (-m(6) * pkin(7) - mrSges(6,3));
t167 = t88 / 0.2e1;
t128 = qJD(5) * t88;
t81 = qJDD(3) + qJDD(4);
t83 = qJD(3) + qJD(4);
t91 = cos(qJ(5));
t49 = -t83 * t128 + t81 * t91;
t36 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t49;
t127 = qJD(5) * t91;
t50 = t83 * t127 + t81 * t88;
t37 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t50;
t105 = t91 * t36 - t88 * t37;
t138 = t88 * mrSges(6,3);
t59 = qJD(5) * mrSges(6,1) - t83 * t138;
t144 = t83 * t91;
t60 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t144;
t178 = -t59 * t127 - t60 * t128 + t105;
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t90 = sin(qJ(3));
t93 = cos(qJ(3));
t56 = -t90 * t84 + t86 * t93;
t157 = mrSges(6,1) * t91;
t112 = t155 - t157;
t111 = mrSges(6,1) * t88 + mrSges(6,2) * t91;
t57 = t84 * t93 + t90 * t86;
t48 = t57 * qJD(1);
t89 = sin(qJ(4));
t150 = t48 * t89;
t47 = t56 * qJD(1);
t44 = qJD(3) * pkin(3) + t47;
t92 = cos(qJ(4));
t19 = t44 * t92 - t150;
t15 = -pkin(4) * t83 - t19;
t176 = t15 * t111 + qJD(5) * (Ifges(6,5) * t91 - Ifges(6,6) * t88) / 0.2e1;
t46 = t112 * t83;
t174 = t83 * mrSges(5,1) - t46;
t173 = (-mrSges(5,1) + t112) * t74 + (mrSges(5,2) - mrSges(6,3)) * t73;
t132 = t92 * t48;
t20 = t44 * t89 + t132;
t16 = pkin(7) * t83 + t20;
t13 = qJD(2) * t91 - t16 * t88;
t14 = qJD(2) * t88 + t16 * t91;
t171 = -t13 * t127 - t14 * t128;
t169 = -t83 * mrSges(5,2) - t88 * t59 + t91 * t60;
t52 = t57 * qJD(3);
t28 = -qJD(1) * t52 + t56 * qJDD(1);
t25 = qJDD(3) * pkin(3) + t28;
t51 = t56 * qJD(3);
t27 = qJD(1) * t51 + t57 * qJDD(1);
t9 = -qJD(4) * t20 + t25 * t92 - t27 * t89;
t168 = m(6) * pkin(4);
t78 = sin(t82);
t165 = pkin(3) * t78;
t79 = cos(t82);
t164 = pkin(3) * t79;
t163 = pkin(3) * t89;
t162 = pkin(3) * t92;
t129 = qJD(4) * t92;
t130 = qJD(4) * t89;
t8 = t44 * t129 - t48 * t130 + t89 * t25 + t92 * t27;
t5 = pkin(7) * t81 + t8;
t2 = qJD(5) * t13 + qJDD(2) * t88 + t5 * t91;
t159 = t2 * t91;
t3 = -qJD(5) * t14 + qJDD(2) * t91 - t5 * t88;
t158 = t3 * t88;
t156 = mrSges(5,2) * t74;
t154 = Ifges(6,1) * t88;
t153 = Ifges(6,4) * t88;
t152 = Ifges(6,4) * t91;
t151 = Ifges(6,2) * t91;
t147 = t81 * mrSges(5,2);
t85 = sin(pkin(8));
t143 = t85 * t88;
t142 = t85 * t91;
t87 = cos(pkin(8));
t140 = t87 * t88;
t139 = t87 * t91;
t131 = t74 * pkin(4) + t73 * pkin(7);
t126 = m(3) + m(4) + m(5);
t118 = m(6) + t126;
t115 = t179 * t85;
t114 = t179 * t87;
t110 = t151 + t153;
t108 = t13 * t91 + t14 * t88;
t107 = -t13 * t88 + t14 * t91;
t104 = t92 * t56 - t57 * t89;
t30 = t56 * t89 + t57 * t92;
t101 = t88 * (Ifges(6,1) * t91 - t153);
t98 = -t108 * qJD(5) - t158;
t97 = m(6) * (-pkin(4) * t73 - t165) - t73 * t157;
t96 = t156 + (mrSges(5,1) + t157 + t168) * t73;
t95 = t98 + t159;
t40 = Ifges(6,6) * qJD(5) + t110 * t83;
t67 = Ifges(6,4) * t144;
t41 = Ifges(6,5) * qJD(5) + t83 * t154 + t67;
t6 = -pkin(4) * t81 - t9;
t94 = t9 * mrSges(5,1) - t8 * mrSges(5,2) + mrSges(6,3) * t159 + t6 * t112 + (Ifges(6,1) * t50 + Ifges(6,4) * t49) * t167 + t91 * (Ifges(6,4) * t50 + Ifges(6,2) * t49) / 0.2e1 + t49 * t110 / 0.2e1 + t50 * (t152 + t154) / 0.2e1 - t40 * t128 / 0.2e1 + Ifges(5,3) * t81 + (t41 + t83 * (-Ifges(6,2) * t88 + t152)) * t127 / 0.2e1 + (0.2e1 * Ifges(6,5) * t167 + Ifges(6,6) * t91) * qJDD(5) + (t101 * t83 / 0.2e1 + t176) * qJD(5);
t26 = -mrSges(6,1) * t49 + mrSges(6,2) * t50;
t11 = t30 * qJD(4) + t51 * t89 + t92 * t52;
t10 = t104 * qJD(4) + t51 * t92 - t52 * t89;
t1 = [t11 * t46 - t104 * t26 + (-qJD(3) * t51 - qJDD(3) * t57) * mrSges(4,2) + (t104 * t81 - t11 * t83) * mrSges(5,1) + (-qJD(3) * t52 + qJDD(3) * t56) * mrSges(4,1) + t169 * t10 + (-t147 + (-t91 * t59 - t88 * t60) * qJD(5) + t105) * t30 + (-m(2) - t118) * g(3) + m(6) * (t107 * t10 - t104 * t6 + t11 * t15 + t95 * t30) + m(4) * (t27 * t57 + t28 * t56 - t47 * t52 + t48 * t51) + m(5) * (t10 * t20 + t104 * t9 - t11 * t19 + t30 * t8) + (m(2) + m(3) * (t84 ^ 2 + t86 ^ 2)) * qJDD(1); t60 * t127 - t59 * t128 + t88 * t36 + t91 * t37 + m(6) * (t107 * qJD(5) + t2 * t88 + t3 * t91) + t126 * qJDD(2) + (-t85 * g(1) + t87 * g(2)) * t118; -t3 * t138 - g(1) * (t97 * t87 - t114) - g(2) * (t97 * t85 - t115) + t94 - t147 * t163 + t81 * mrSges(5,1) * t162 + Ifges(4,3) * qJDD(3) + t171 * mrSges(6,3) + (t47 * qJD(3) - t27) * mrSges(4,2) + (t48 * qJD(3) + t28) * mrSges(4,1) + (m(5) * t19 - m(6) * t15 + t174) * (t47 * t89 + t132) + (-m(5) * t20 - m(6) * t107 - t169) * (t47 * t92 - t150) + (-m(5) * t164 - m(6) * (t131 + t164) - mrSges(4,1) * t79 + mrSges(4,2) * t78 + t173) * g(3) + (m(6) * t95 + t178) * (pkin(7) + t163) + (m(6) * t6 + t26) * (-pkin(4) - t162) + (m(5) * t165 + mrSges(4,1) * t78 + mrSges(5,1) * t73 + mrSges(4,2) * t79 + t156) * (g(1) * t87 + g(2) * t85) + (m(6) * (t107 * t92 + t15 * t89) * qJD(4) + m(5) * (t8 * t89 + t9 * t92 + (-t19 * t89 + t20 * t92) * qJD(4)) - t174 * t130 + t169 * t129) * pkin(3); (-m(6) * t131 + t173) * g(3) + t94 + (m(6) * (-t158 + t159 + t171) + t178) * pkin(7) + t98 * mrSges(6,3) + (t96 * t87 + t114) * g(1) + (t96 * t85 + t115) * g(2) - pkin(4) * t26 - m(6) * (t107 * t19 + t15 * t20) - t6 * t168 - t169 * t19 + t174 * t20; Ifges(6,5) * t50 + Ifges(6,6) * t49 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t13 * t60 + t14 * t59 - g(1) * ((-t74 * t140 + t142) * mrSges(6,1) + (-t74 * t139 - t143) * mrSges(6,2)) - g(2) * ((-t74 * t143 - t139) * mrSges(6,1) + (-t74 * t142 + t140) * mrSges(6,2)) + g(3) * t111 * t73 + (t40 * t167 + (-t101 / 0.2e1 + t151 * t167) * t83 + t108 * mrSges(6,3) - (t41 + t67) * t91 / 0.2e1 - t176) * t83;];
tau = t1;
