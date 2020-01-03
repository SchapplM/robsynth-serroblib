% Calculate vector of inverse dynamics joint torques for
% S4RPRR7
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:53:51
% DurationCPUTime: 5.09s
% Computational Cost: add. (2075->325), mult. (5101->455), div. (0->0), fcn. (3536->10), ass. (0->146)
t132 = qJD(1) * qJD(2);
t80 = qJDD(1) * qJ(2) + t132;
t90 = sin(pkin(7));
t91 = cos(pkin(7));
t138 = t90 ^ 2 + t91 ^ 2;
t94 = sin(qJ(3));
t97 = cos(qJ(3));
t72 = t90 * t97 + t91 * t94;
t67 = t72 * qJD(1);
t153 = t67 * mrSges(4,3);
t93 = sin(qJ(4));
t96 = cos(qJ(4));
t51 = qJD(3) * t96 - t67 * t93;
t52 = qJD(3) * t93 + t67 * t96;
t178 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t51 - mrSges(5,2) * t52 - t153;
t139 = pkin(5) + qJ(2);
t77 = t139 * t90;
t73 = qJD(1) * t77;
t78 = t139 * t91;
t74 = qJD(1) * t78;
t46 = -t73 * t97 - t74 * t94;
t39 = -qJD(3) * pkin(3) - t46;
t189 = -m(5) * t39 + t178;
t188 = -m(4) - m(5);
t186 = t51 * Ifges(5,6);
t71 = t90 * t94 - t97 * t91;
t66 = t71 * qJD(1);
t61 = qJD(4) + t66;
t185 = t61 * Ifges(5,3);
t184 = Ifges(4,5) * qJD(3);
t183 = Ifges(4,6) * qJD(3);
t69 = t72 * qJD(3);
t45 = -qJD(1) * t69 - qJDD(1) * t71;
t123 = m(3) * qJ(2) + mrSges(3,3);
t182 = mrSges(2,2) - mrSges(4,3) - t123;
t89 = pkin(7) + qJ(3);
t85 = sin(t89);
t86 = cos(t89);
t116 = mrSges(4,1) * t86 - mrSges(4,2) * t85;
t117 = -mrSges(3,1) * t91 + mrSges(3,2) * t90;
t180 = m(3) * pkin(1) + t85 * mrSges(5,3) + mrSges(2,1) + t116 - t117;
t83 = pkin(2) * t91 + pkin(1);
t76 = -qJD(1) * t83 + qJD(2);
t31 = pkin(3) * t66 - pkin(6) * t67 + t76;
t47 = -t73 * t94 + t74 * t97;
t40 = qJD(3) * pkin(6) + t47;
t10 = t31 * t96 - t40 * t93;
t68 = t71 * qJD(3);
t44 = -qJD(1) * t68 + qJDD(1) * t72;
t75 = -qJDD(1) * t83 + qJDD(2);
t12 = -pkin(3) * t45 - pkin(6) * t44 + t75;
t120 = pkin(5) * qJDD(1) + t80;
t57 = t120 * t90;
t58 = t120 * t91;
t15 = t46 * qJD(3) - t94 * t57 + t97 * t58;
t13 = qJDD(3) * pkin(6) + t15;
t1 = qJD(4) * t10 + t12 * t93 + t13 * t96;
t11 = t31 * t93 + t40 * t96;
t2 = -qJD(4) * t11 + t12 * t96 - t13 * t93;
t179 = t2 * mrSges(5,1) - t1 * mrSges(5,2);
t16 = -qJD(3) * t47 - t57 * t97 - t58 * t94;
t26 = qJD(4) * t51 + qJDD(3) * t93 + t44 * t96;
t174 = t26 / 0.2e1;
t27 = -qJD(4) * t52 + qJDD(3) * t96 - t44 * t93;
t173 = t27 / 0.2e1;
t38 = qJDD(4) - t45;
t172 = t38 / 0.2e1;
t171 = -t51 / 0.2e1;
t170 = -t52 / 0.2e1;
t169 = t52 / 0.2e1;
t168 = -t61 / 0.2e1;
t166 = -t66 / 0.2e1;
t164 = t67 / 0.2e1;
t163 = t96 / 0.2e1;
t159 = mrSges(5,3) * t93;
t158 = Ifges(4,4) * t67;
t157 = Ifges(5,4) * t52;
t156 = Ifges(5,4) * t93;
t155 = Ifges(5,4) * t96;
t154 = t46 * mrSges(4,3);
t152 = t72 * t93;
t20 = Ifges(5,2) * t51 + Ifges(5,6) * t61 + t157;
t147 = t93 * t20;
t95 = sin(qJ(1));
t146 = t93 * t95;
t98 = cos(qJ(1));
t145 = t93 * t98;
t144 = t95 * t96;
t143 = t96 * mrSges(5,3);
t50 = Ifges(5,4) * t51;
t21 = Ifges(5,1) * t52 + Ifges(5,5) * t61 + t50;
t142 = t96 * t21;
t141 = t96 * t98;
t137 = qJD(4) * t93;
t136 = qJD(4) * t96;
t135 = qJDD(1) * t90;
t134 = qJDD(1) * t91;
t130 = Ifges(5,5) * t26 + Ifges(5,6) * t27 + Ifges(5,3) * t38;
t129 = t72 * t136;
t128 = t142 / 0.2e1;
t127 = m(5) * pkin(6) + mrSges(5,3);
t124 = -t137 / 0.2e1;
t122 = -t45 * mrSges(4,1) + t44 * mrSges(4,2);
t119 = pkin(3) * t86 + pkin(6) * t85;
t118 = -mrSges(3,1) * t134 + mrSges(3,2) * t135;
t114 = -mrSges(5,1) * t96 + mrSges(5,2) * t93;
t113 = mrSges(5,1) * t93 + mrSges(5,2) * t96;
t112 = Ifges(5,1) * t96 - t156;
t111 = -Ifges(5,2) * t93 + t155;
t110 = Ifges(5,5) * t96 - Ifges(5,6) * t93;
t29 = -mrSges(5,2) * t61 + mrSges(5,3) * t51;
t30 = mrSges(5,1) * t61 - mrSges(5,3) * t52;
t108 = t29 * t96 - t30 * t93;
t43 = pkin(3) * t71 - pkin(6) * t72 - t83;
t49 = -t77 * t94 + t78 * t97;
t24 = t43 * t96 - t49 * t93;
t25 = t43 * t93 + t49 * t96;
t106 = -t97 * t77 - t78 * t94;
t105 = -t68 * t93 + t129;
t104 = t137 * t72 + t68 * t96;
t103 = t39 * t113;
t101 = m(5) * pkin(3) - t114;
t84 = -qJDD(1) * pkin(1) + qJDD(2);
t65 = t141 * t86 + t146;
t64 = -t145 * t86 + t144;
t63 = -t144 * t86 + t145;
t62 = t146 * t86 + t141;
t59 = Ifges(4,4) * t66;
t53 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t66;
t42 = pkin(3) * t69 + pkin(6) * t68;
t41 = pkin(3) * t67 + pkin(6) * t66;
t35 = t67 * Ifges(4,1) + t184 - t59;
t34 = -t66 * Ifges(4,2) + t158 + t183;
t32 = -qJD(2) * t71 + qJD(3) * t106;
t19 = t52 * Ifges(5,5) + t185 + t186;
t18 = t41 * t93 + t46 * t96;
t17 = t41 * t96 - t46 * t93;
t14 = -qJDD(3) * pkin(3) - t16;
t9 = -mrSges(5,2) * t38 + mrSges(5,3) * t27;
t8 = mrSges(5,1) * t38 - mrSges(5,3) * t26;
t7 = -mrSges(5,1) * t27 + mrSges(5,2) * t26;
t6 = -qJD(4) * t25 - t32 * t93 + t42 * t96;
t5 = qJD(4) * t24 + t32 * t96 + t42 * t93;
t4 = t26 * Ifges(5,1) + t27 * Ifges(5,4) + t38 * Ifges(5,5);
t3 = t26 * Ifges(5,4) + t27 * Ifges(5,2) + t38 * Ifges(5,6);
t22 = [(-Ifges(4,1) * t164 - Ifges(4,4) * t166 - t76 * mrSges(4,2) - t184 / 0.2e1 - t35 / 0.2e1 + t154 + t147 / 0.2e1 - t128) * t68 + (t130 / 0.2e1 + t75 * mrSges(4,1) - t44 * Ifges(4,4) - t45 * Ifges(4,2) - Ifges(4,6) * qJDD(3) - t15 * mrSges(4,3) + Ifges(5,3) * t172 + Ifges(5,6) * t173 + Ifges(5,5) * t174 + t179) * t71 + (m(4) * t16 - m(5) * t14 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t44 - t7) * t106 + m(4) * (t15 * t49 + t32 * t47 - t75 * t83) + m(5) * (t1 * t25 + t10 * t6 + t11 * t5 + t2 * t24) + 0.2e1 * t138 * t80 * mrSges(3,3) - t20 * t129 / 0.2e1 + (-Ifges(5,1) * t104 - Ifges(5,4) * t105) * t169 + (-t1 * t152 + t10 * t104 - t11 * t105) * mrSges(5,3) + (Ifges(3,4) * t90 + Ifges(3,2) * t91) * t134 + (Ifges(3,1) * t90 + Ifges(3,4) * t91) * t135 + t32 * t53 + t49 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t45) + t5 * t29 + t6 * t30 + t84 * t117 - pkin(1) * t118 + (-m(4) * t46 - t189) * (qJD(2) * t72 + qJD(3) * t49) + t51 * (-Ifges(5,4) * t104 - Ifges(5,2) * t105) / 0.2e1 + t24 * t8 + t25 * t9 + t39 * (mrSges(5,1) * t105 - mrSges(5,2) * t104) - t83 * t122 + (t75 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,1) * t44 + Ifges(4,4) * t45 + Ifges(4,5) * qJDD(3) + t110 * t172 + t111 * t173 + t112 * t174 + t113 * t14 + t124 * t21 - t143 * t2 + t163 * t4) * t72 + t61 * (-Ifges(5,5) * t104 - Ifges(5,6) * t105) / 0.2e1 + m(3) * (-pkin(1) * t84 + (t132 + t80) * qJ(2) * t138) - t3 * t152 / 0.2e1 + (-Ifges(4,4) * t164 - Ifges(4,2) * t166 + Ifges(5,5) * t169 + t76 * mrSges(4,1) + t19 / 0.2e1 - t183 / 0.2e1 - t34 / 0.2e1 - t47 * mrSges(4,3) + t10 * mrSges(5,1) + t185 / 0.2e1 + t186 / 0.2e1 - t11 * mrSges(5,2)) * t69 + (-t63 * mrSges(5,1) - t62 * mrSges(5,2) + (t188 * t139 + t182) * t98 + (m(4) * t83 - m(5) * (-t119 - t83) + t180) * t95) * g(1) + (-t65 * mrSges(5,1) - t64 * mrSges(5,2) + t188 * (t139 * t95 + t98 * t83) + t182 * t95 + (-m(5) * t119 - t180) * t98) * g(2) + Ifges(2,3) * qJDD(1); t96 * t8 + t93 * t9 + t178 * t67 + t108 * qJD(4) - (-t108 - t53) * t66 + m(3) * t84 + t118 + t122 - t123 * qJD(1) ^ 2 * t138 + (-g(1) * t95 + g(2) * t98) * (m(3) - t188) + (t1 * t93 + t2 * t96 - t39 * t67 + t61 * (-t10 * t93 + t11 * t96)) * m(5) + (t46 * t67 + t47 * t66 + t75) * m(4); (t103 + t128) * qJD(4) - t66 * t154 + t3 * t163 + t34 * t164 + t147 * t166 + ((mrSges(4,2) - t127) * t86 + (mrSges(4,1) + t101) * t85) * (g(1) * t98 + g(2) * t95) + (-Ifges(4,2) * t67 + t142 + t35 - t59) * t66 / 0.2e1 + (t110 * t61 + t111 * t51 + t112 * t52) * qJD(4) / 0.2e1 + (Ifges(5,3) * t67 - t110 * t66) * t168 + (Ifges(5,5) * t67 - t112 * t66) * t170 + (Ifges(5,6) * t67 - t111 * t66) * t171 - (-Ifges(4,1) * t66 - t158 + t19) * t67 / 0.2e1 - t76 * (mrSges(4,1) * t67 - mrSges(4,2) * t66) - qJD(3) * (-Ifges(4,5) * t66 - Ifges(4,6) * t67) / 0.2e1 - t10 * (mrSges(5,1) * t67 + t143 * t66) - t11 * (-mrSges(5,2) * t67 + t159 * t66) + (-t10 * t136 - t11 * t137) * mrSges(5,3) + t1 * t143 + (-pkin(3) * t14 - t10 * t17 - t11 * t18) * m(5) + t93 * t4 / 0.2e1 - t46 * t53 + Ifges(4,5) * t44 + Ifges(4,6) * t45 - t18 * t29 - t17 * t30 + (t153 + t189) * t47 + t14 * t114 + (-t101 * t86 - t127 * t85 - t116) * g(3) + (m(5) * (t1 * t96 - t2 * t93 + (-t10 * t96 - t11 * t93) * qJD(4)) + t96 * t9 - t93 * t8 - t30 * t136 - t29 * t137) * pkin(6) - t15 * mrSges(4,2) + t16 * mrSges(4,1) - pkin(3) * t7 + (Ifges(5,5) * t93 + Ifges(5,6) * t96) * t172 + (Ifges(5,2) * t96 + t156) * t173 + (Ifges(5,1) * t93 + t155) * t174 + t66 * t103 - t2 * t159 + t20 * t124 + Ifges(4,3) * qJDD(3); -t39 * (mrSges(5,1) * t52 + mrSges(5,2) * t51) + (Ifges(5,1) * t51 - t157) * t170 + t20 * t169 + (Ifges(5,5) * t51 - Ifges(5,6) * t52) * t168 - t10 * t29 + t11 * t30 - g(1) * (mrSges(5,1) * t64 - mrSges(5,2) * t65) - g(2) * (-mrSges(5,1) * t62 + mrSges(5,2) * t63) + g(3) * t113 * t85 + (t10 * t51 + t11 * t52) * mrSges(5,3) + t130 + (-Ifges(5,2) * t52 + t21 + t50) * t171 + t179;];
tau = t22;
