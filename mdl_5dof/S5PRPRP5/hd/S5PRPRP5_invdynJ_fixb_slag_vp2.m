% Calculate vector of inverse dynamics joint torques for
% S5PRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:50
% DurationCPUTime: 5.62s
% Computational Cost: add. (1563->291), mult. (3666->373), div. (0->0), fcn. (2438->10), ass. (0->131)
t190 = mrSges(5,1) + mrSges(6,1);
t189 = mrSges(5,2) - mrSges(6,3);
t101 = cos(qJ(2));
t120 = -qJD(1) * t101 + qJD(3);
t94 = sin(pkin(8));
t96 = cos(pkin(8));
t145 = t94 ^ 2 + t96 ^ 2;
t179 = mrSges(4,3) * t145;
t191 = -m(5) - m(6);
t188 = Ifges(5,1) + Ifges(6,1);
t186 = Ifges(6,4) + Ifges(5,5);
t161 = cos(qJ(4));
t130 = t161 * t96;
t99 = sin(qJ(4));
t151 = t99 * t94;
t109 = t130 - t151;
t106 = t101 * t109;
t148 = pkin(6) + qJ(3);
t74 = t148 * t94;
t75 = t148 * t96;
t110 = -t161 * t74 - t99 * t75;
t184 = -qJD(1) * t106 + qJD(3) * t109 + qJD(4) * t110;
t40 = t161 * t75 - t99 * t74;
t69 = t161 * t94 + t99 * t96;
t183 = qJD(4) * t40 + t120 * t69;
t195 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t95 = sin(pkin(7));
t97 = cos(pkin(7));
t194 = g(1) * t97 + g(2) * t95;
t193 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t117 = -mrSges(4,1) * t96 + mrSges(4,2) * t94;
t93 = pkin(8) + qJ(4);
t88 = sin(t93);
t89 = cos(t93);
t192 = m(4) * pkin(2) - t189 * t88 + t190 * t89 + mrSges(3,1) - t117;
t187 = -Ifges(5,4) + Ifges(6,5);
t164 = t69 / 0.2e1;
t185 = Ifges(6,6) - Ifges(5,6);
t63 = t109 * qJD(4);
t33 = qJD(2) * t63 + qJDD(2) * t69;
t19 = -qJDD(4) * mrSges(6,1) + t33 * mrSges(6,2);
t182 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t33 + t19;
t137 = qJDD(2) * t94;
t64 = t69 * qJD(4);
t34 = qJD(2) * t64 - qJDD(2) * t130 + t137 * t99;
t21 = -mrSges(6,2) * t34 + qJDD(4) * mrSges(6,3);
t181 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t34 + t21;
t61 = t109 * qJD(2);
t159 = Ifges(6,5) * t61;
t59 = Ifges(5,4) * t61;
t62 = t69 * qJD(2);
t180 = t186 * qJD(4) + t188 * t62 - t159 + t59;
t46 = mrSges(6,2) * t61 + qJD(4) * mrSges(6,3);
t147 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t61 + t46;
t153 = t62 * mrSges(5,3);
t146 = -mrSges(6,2) * t62 + t190 * qJD(4) - t153;
t100 = sin(qJ(2));
t141 = qJD(1) * t100;
t78 = qJD(2) * qJ(3) + t141;
t123 = pkin(6) * qJD(2) + t78;
t54 = t123 * t96;
t152 = t99 * t54;
t53 = t123 * t94;
t14 = -t161 * t53 - t152;
t178 = -t14 + qJD(5);
t6 = t34 * mrSges(6,1) - t33 * mrSges(6,3);
t136 = qJDD(2) * t96;
t67 = -mrSges(4,1) * t136 + mrSges(4,2) * t137;
t7 = t34 * mrSges(5,1) + t33 * mrSges(5,2);
t177 = t6 + t67 + t7;
t27 = -mrSges(6,1) * t61 - mrSges(6,3) * t62;
t175 = -mrSges(5,1) * t61 + mrSges(5,2) * t62 + t117 * qJD(2) + t27;
t15 = t161 * t54 - t99 * t53;
t134 = qJD(1) * qJD(2);
t87 = t101 * t134;
t72 = t100 * qJDD(1) + t87;
t55 = t72 + t195;
t121 = pkin(6) * qJDD(2) + t55;
t37 = t121 * t94;
t38 = t121 * t96;
t5 = -qJD(4) * t15 - t161 * t37 - t99 * t38;
t11 = qJD(4) * qJ(5) + t15;
t85 = pkin(3) * t96 + pkin(2);
t65 = -qJD(2) * t85 + t120;
t9 = -pkin(4) * t61 - qJ(5) * t62 + t65;
t172 = t65 * mrSges(5,1) + t9 * mrSges(6,1) - t11 * mrSges(6,2);
t10 = -qJD(4) * pkin(4) + t178;
t171 = mrSges(5,2) * t65 + t10 * mrSges(6,2) - t14 * mrSges(5,3) - mrSges(6,3) * t9;
t169 = t61 / 0.2e1;
t168 = -t61 / 0.2e1;
t166 = t62 / 0.2e1;
t160 = Ifges(5,4) * t62;
t158 = g(3) * t100;
t149 = qJD(4) / 0.2e1;
t144 = t101 * t95;
t143 = t101 * t97;
t139 = qJD(2) * t100;
t138 = qJD(4) * t100;
t135 = m(4) - t191;
t124 = qJD(4) * t161;
t131 = -t53 * t124 + t161 * t38 - t99 * t37;
t129 = t145 * t55;
t128 = t145 * t78;
t86 = t100 * t134;
t71 = qJDD(1) * t101 - t86;
t118 = -mrSges(3,2) + t179;
t114 = pkin(4) * t89 + qJ(5) * t88;
t112 = qJDD(3) - t71;
t107 = t100 * (-qJD(2) * pkin(2) + t120) + t101 * t128;
t42 = -qJDD(2) * t85 + t112;
t102 = qJD(2) ^ 2;
t60 = -qJDD(2) * pkin(2) + t112;
t58 = Ifges(6,5) * t62;
t57 = t109 * t100;
t56 = t69 * t100;
t52 = t143 * t89 + t88 * t95;
t51 = t143 * t88 - t95 * t89;
t50 = t144 * t89 - t97 * t88;
t49 = t144 * t88 + t89 * t97;
t32 = -pkin(4) * t109 - qJ(5) * t69 - t85;
t26 = pkin(4) * t62 - qJ(5) * t61;
t23 = t61 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t160;
t22 = Ifges(6,6) * qJD(4) - t61 * Ifges(6,3) + t58;
t17 = t100 * t124 * t96 + t101 * t62 - t138 * t151;
t16 = qJD(2) * t106 - t138 * t69;
t8 = pkin(4) * t64 - qJ(5) * t63 - qJD(5) * t69;
t4 = -qJD(4) * t152 + t131;
t3 = pkin(4) * t34 - qJ(5) * t33 - qJD(5) * t62 + t42;
t2 = -qJDD(4) * pkin(4) + qJDD(5) - t5;
t1 = qJDD(4) * qJ(5) + (qJD(5) - t152) * qJD(4) + t131;
t12 = [m(2) * qJDD(1) + t181 * t57 + t182 * t56 - t146 * t17 + t147 * t16 + (-m(2) - m(3) - t135) * g(3) + (qJDD(2) * mrSges(3,1) + t102 * t118 - t177) * t101 + (-t102 * mrSges(3,1) + qJD(2) * t175 + t118 * qJDD(2)) * t100 + m(4) * (qJD(2) * t107 + t100 * t129 - t101 * t60) + m(3) * (t100 * t72 + t101 * t71) + m(5) * (-t101 * t42 + t139 * t65 - t14 * t17 + t15 * t16 + t4 * t57 - t5 * t56) + m(6) * (t1 * t57 + t10 * t17 - t101 * t3 + t11 * t16 + t139 * t9 + t2 * t56); (-t72 + t87) * mrSges(3,2) - t183 * t146 + t184 * t147 + (t71 + t86) * mrSges(3,1) + (-Ifges(5,2) * t169 + Ifges(6,3) * t168 + t22 / 0.2e1 - t23 / 0.2e1 - t15 * mrSges(5,3) + t187 * t166 + t185 * t149 + t172) * t64 + (t186 * qJDD(4) + t188 * t33) * t164 + t181 * t40 + (-m(5) * t65 - m(6) * t9 - t175) * t141 + (Ifges(4,4) * t94 + Ifges(4,2) * t96) * t136 + (Ifges(4,1) * t94 + Ifges(4,4) * t96) * t137 + (Ifges(5,4) * t169 + Ifges(6,5) * t168 + t186 * t149 + t188 * t166 + t171 + t180 / 0.2e1) * t63 - t85 * t7 - pkin(2) * t67 + t32 * t6 + t8 * t27 + (t191 * (t100 * t148 + t101 * t85) + (-m(6) * t114 - t192) * t101 + t193 * t100) * g(3) + t60 * t117 + t194 * ((t191 * t148 + t193) * t101 + (m(5) * t85 - m(6) * (-t114 - t85) + t192) * t100) + (-(Ifges(6,3) + Ifges(5,2)) * t109 + 0.2e1 * t187 * t164) * t34 + (-pkin(2) * t60 + qJ(3) * t129 - qJD(1) * t107 + qJD(3) * t128) * m(4) - t182 * t110 + (t1 * t40 + t10 * t183 + t11 * t184 - t110 * t2 + t3 * t32 + t8 * t9) * m(6) + (t110 * t5 - t14 * t183 + t15 * t184 + t4 * t40 - t42 * t85) * m(5) + t109 * (Ifges(5,4) * t33 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t3 * (-mrSges(6,1) * t109 - mrSges(6,3) * t69) + t42 * (-mrSges(5,1) * t109 + mrSges(5,2) * t69) + (t109 * t4 - t5 * t69) * mrSges(5,3) + (t1 * t109 + t2 * t69) * mrSges(6,2) + (-t109 * t185 + t186 * t69) * qJDD(4) / 0.2e1 + (-t109 * t187 + t188 * t69) * t33 / 0.2e1 - t109 * (Ifges(6,5) * t33 + Ifges(6,6) * qJDD(4)) / 0.2e1 + Ifges(3,3) * qJDD(2) + (t55 - t87 + t195) * t179; -t102 * t179 + t146 * t62 - t147 * t61 + (-t10 * t62 - t11 * t61 + t3) * m(6) + (t14 * t62 - t15 * t61 + t42) * m(5) + (-qJD(2) * t128 + t60) * m(4) + (t101 * g(3) - t194 * t100) * t135 + t177; t5 * mrSges(5,1) - t2 * mrSges(6,1) - t4 * mrSges(5,2) + t1 * mrSges(6,3) - pkin(4) * t19 + qJ(5) * t21 + qJD(5) * t46 + t159 * t169 + t23 * t166 - t26 * t27 + (Ifges(6,3) * t169 - t172) * t62 - t171 * t61 + t185 * t34 + t186 * t33 + (t189 * t89 + t190 * t88) * t158 - (t185 * t62 + t186 * t61) * qJD(4) / 0.2e1 - t147 * t14 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + (t189 * t50 + t190 * t49) * g(2) + (t189 * t52 + t190 * t51) * g(1) + (-t26 * t9 - pkin(4) * t2 + qJ(5) * t1 - g(2) * (-pkin(4) * t49 + qJ(5) * t50) - g(1) * (-pkin(4) * t51 + qJ(5) * t52) - (-pkin(4) * t88 + qJ(5) * t89) * t158 + t178 * t11) * m(6) + (-Ifges(5,2) * t62 + t180 + t59) * t168 - (t188 * t61 - t160 + t22 + t58) * t62 / 0.2e1 + (-m(6) * t10 + t146 + t153) * t15; -qJD(4) * t46 + t62 * t27 + (-g(1) * t51 - g(2) * t49 - t11 * qJD(4) - t158 * t88 + t9 * t62 + t2) * m(6) + t19;];
tau = t12;
