% Calculate vector of inverse dynamics joint torques for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:08:50
% DurationCPUTime: 6.98s
% Computational Cost: add. (1788->349), mult. (3778->442), div. (0->0), fcn. (2265->12), ass. (0->150)
t202 = m(6) + m(5);
t201 = mrSges(5,1) + mrSges(6,1);
t200 = mrSges(5,2) - mrSges(6,3);
t194 = Ifges(5,1) + Ifges(6,1);
t192 = Ifges(6,4) + Ifges(5,5);
t108 = sin(qJ(3));
t159 = Ifges(4,4) * t108;
t193 = -Ifges(5,4) + Ifges(6,5);
t105 = sin(pkin(7));
t87 = pkin(1) * t105 + pkin(6);
t77 = t87 * qJDD(1);
t199 = qJD(2) * qJD(3) + t77;
t110 = cos(qJ(3));
t83 = -t110 * mrSges(4,1) + t108 * mrSges(4,2);
t102 = qJ(3) + pkin(8);
t93 = sin(t102);
t95 = cos(t102);
t198 = t200 * t93 - t201 * t95 + t83;
t103 = qJ(1) + pkin(7);
t94 = sin(t103);
t96 = cos(t103);
t197 = g(1) * t96 + g(2) * t94;
t104 = sin(pkin(8));
t152 = cos(pkin(8));
t133 = t152 * t108;
t71 = t104 * t110 + t133;
t175 = t71 / 0.2e1;
t196 = -m(4) - m(3);
t162 = qJD(3) / 0.2e1;
t191 = Ifges(6,6) - Ifges(5,6);
t132 = t152 * t110;
t149 = qJD(1) * t108;
t62 = -qJD(1) * t132 + t104 * t149;
t169 = Ifges(6,5) * t62;
t58 = Ifges(5,4) * t62;
t64 = t71 * qJD(1);
t190 = t192 * qJD(3) + t194 * t64 + t169 - t58;
t144 = qJD(1) * qJD(3);
t73 = qJDD(1) * t110 - t108 * t144;
t74 = qJDD(1) * t108 + t110 * t144;
t39 = t104 * t74 - t152 * t73;
t30 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t39;
t33 = -mrSges(6,2) * t39 + qJDD(3) * mrSges(6,3);
t189 = t30 + t33;
t40 = t104 * t73 + t152 * t74;
t31 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t40;
t32 = -qJDD(3) * mrSges(6,1) + t40 * mrSges(6,2);
t188 = -t31 + t32;
t164 = t62 * mrSges(5,3);
t172 = mrSges(6,2) * t62;
t53 = qJD(3) * mrSges(6,3) - t172;
t161 = -qJD(3) * mrSges(5,2) - t164 + t53;
t171 = mrSges(6,2) * t64;
t160 = -mrSges(5,3) * t64 + qJD(3) * t201 - t171;
t100 = t110 * pkin(3);
t106 = cos(pkin(7));
t89 = -pkin(1) * t106 - pkin(2);
t76 = -t100 + t89;
t146 = qJD(3) * t108;
t79 = t87 * qJD(1);
t26 = t108 * qJDD(2) + t110 * t199 - t146 * t79;
t147 = qJD(2) * t108;
t55 = t110 * t79 + t147;
t98 = t110 * qJDD(2);
t27 = -t55 * qJD(3) - t108 * t77 + t98;
t186 = -t108 * t27 + t110 * t26;
t183 = m(4) * pkin(2) + mrSges(3,1) - t198;
t182 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t180 = -t62 / 0.2e1;
t179 = t62 / 0.2e1;
t177 = t64 / 0.2e1;
t143 = qJD(1) * qJD(4);
t145 = qJD(3) * t110;
t10 = -t79 * t145 + qJDD(3) * pkin(3) - qJ(4) * t74 + t98 + (-t143 - t199) * t108;
t13 = qJ(4) * t73 + t110 * t143 + t26;
t4 = t104 * t10 + t152 * t13;
t170 = Ifges(5,4) * t64;
t109 = sin(qJ(1));
t168 = pkin(1) * t109;
t167 = pkin(3) * t104;
t111 = cos(qJ(1));
t101 = t111 * pkin(1);
t129 = qJ(4) * qJD(1) + t79;
t48 = t110 * t129 + t147;
t42 = t152 * t48;
t99 = t110 * qJD(2);
t47 = -t108 * t129 + t99;
t44 = qJD(3) * pkin(3) + t47;
t12 = t104 * t44 + t42;
t165 = t12 * mrSges(5,3);
t158 = Ifges(4,4) * t110;
t157 = t104 * t48;
t153 = qJ(4) + t87;
t150 = t104 * t108;
t148 = qJD(1) * t110;
t141 = pkin(3) * t149;
t140 = pkin(3) * t146;
t139 = mrSges(4,3) * t149;
t138 = mrSges(4,3) * t148;
t137 = t152 * pkin(3);
t136 = t39 * mrSges(5,1) + t40 * mrSges(5,2);
t135 = t39 * mrSges(6,1) - t40 * mrSges(6,3);
t130 = qJD(3) * t153;
t128 = -g(1) * t94 + g(2) * t96;
t124 = pkin(4) * t95 + qJ(5) * t93;
t78 = t89 * qJDD(1);
t122 = mrSges(4,1) * t108 + mrSges(4,2) * t110;
t121 = t110 * Ifges(4,2) + t159;
t120 = Ifges(4,5) * t110 - Ifges(4,6) * t108;
t119 = t89 * qJD(1) * t122;
t3 = t10 * t152 - t104 * t13;
t11 = t152 * t44 - t157;
t118 = t108 * (Ifges(4,1) * t110 - t159);
t117 = t132 - t150;
t114 = -qJD(4) * t108 - t110 * t130;
t61 = qJD(1) * t76 + qJD(4);
t46 = -pkin(3) * t73 + qJDD(4) + t78;
t107 = -qJ(4) - pkin(6);
t92 = Ifges(4,4) * t148;
t91 = t100 + pkin(2);
t88 = -t137 - pkin(4);
t85 = qJ(5) + t167;
t82 = -qJD(3) * mrSges(4,2) + t138;
t80 = qJD(3) * mrSges(4,1) - t139;
t69 = t153 * t110;
t67 = Ifges(4,1) * t149 + Ifges(4,5) * qJD(3) + t92;
t66 = Ifges(4,6) * qJD(3) + qJD(1) * t121;
t65 = t117 * qJD(3);
t63 = t71 * qJD(3);
t60 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t74;
t59 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t73;
t57 = Ifges(6,5) * t64;
t54 = -t108 * t79 + t99;
t49 = qJD(4) * t110 - t108 * t130;
t35 = mrSges(5,1) * t62 + mrSges(5,2) * t64;
t34 = mrSges(6,1) * t62 - mrSges(6,3) * t64;
t23 = -t62 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t170;
t22 = Ifges(6,6) * qJD(3) + t62 * Ifges(6,3) + t57;
t21 = -pkin(4) * t117 - qJ(5) * t71 + t76;
t20 = pkin(4) * t64 + qJ(5) * t62 + t141;
t19 = pkin(4) * t62 - qJ(5) * t64 + t61;
t16 = t152 * t47 - t157;
t15 = t104 * t47 + t42;
t14 = pkin(4) * t63 - qJ(5) * t65 - qJD(5) * t71 + t140;
t8 = qJD(3) * qJ(5) + t12;
t7 = -qJD(3) * pkin(4) + qJD(5) - t11;
t5 = pkin(4) * t39 - qJ(5) * t40 - qJD(5) * t64 + t46;
t2 = -qJDD(3) * pkin(4) + qJDD(5) - t3;
t1 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t4;
t6 = [(t19 * mrSges(6,1) + t61 * mrSges(5,1) + t22 / 0.2e1 - t23 / 0.2e1 + Ifges(6,3) * t179 - Ifges(5,2) * t180 - t165 + t193 * t177 + t191 * t162) * t63 + (t192 * qJDD(3) + t194 * t40) * t175 + (t120 * t162 + t119) * qJD(3) + (-m(5) * t11 + m(6) * t7 - t160) * (t104 * t49 - t114 * t152) + (m(5) * t12 + m(6) * t8 + t161) * (t104 * t114 + t152 * t49) + (-m(5) * t3 + m(6) * t2 + t188) * (t104 * t69 + t133 * t153) + (m(5) * t4 + m(6) * t1 + t189) * (-t150 * t153 + t152 * t69) + (-t108 * t60 + t110 * t59 - t80 * t145 - t82 * t146 + m(4) * ((-t55 * t108 - t54 * t110) * qJD(3) + t186)) * t87 + (-t145 * t54 - t146 * t55 + t186) * mrSges(4,3) + t74 * t108 * Ifges(4,1) + (m(4) * t89 + t83) * t78 + (t1 * t117 + t2 * t71 - t63 * t8) * mrSges(6,2) + (t159 + t121) * t73 / 0.2e1 + (-t117 * t191 + t192 * t71) * qJDD(3) / 0.2e1 + (m(3) * (t105 ^ 2 + t106 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3) + 0.2e1 * (mrSges(3,1) * t106 - mrSges(3,2) * t105) * pkin(1)) * qJDD(1) + m(6) * (t14 * t19 + t21 * t5) + m(5) * (t140 * t61 + t46 * t76) + (-mrSges(2,1) * t111 + mrSges(2,2) * t109 - t202 * (-t107 * t94 + t96 * t91 + t101) + t196 * t101 + (-m(6) * t124 - t183) * t96 + t182 * t94) * g(2) + (mrSges(2,1) * t109 + mrSges(2,2) * t111 - t196 * t168 - t202 * (-t107 * t96 - t168) + t182 * t96 + (-m(6) * (-t124 - t91) + m(5) * t91 + t183) * t94) * g(1) + (-t117 * t193 + t194 * t71) * t40 / 0.2e1 + t117 * (Ifges(5,4) * t40 + Ifges(5,6) * qJDD(3)) / 0.2e1 - t117 * (Ifges(6,5) * t40 + Ifges(6,6) * qJDD(3)) / 0.2e1 + t5 * (-mrSges(6,1) * t117 - mrSges(6,3) * t71) + t46 * (-mrSges(5,1) * t117 + mrSges(5,2) * t71) + t110 * (Ifges(4,4) * t74 + Ifges(4,2) * t73) / 0.2e1 + t108 * Ifges(4,5) * qJDD(3) + t74 * t158 / 0.2e1 + t89 * (-mrSges(4,1) * t73 + mrSges(4,2) * t74) + t14 * t34 + t110 * qJDD(3) * Ifges(4,6) + (-(Ifges(6,3) + Ifges(5,2)) * t117 + 0.2e1 * t193 * t175) * t39 + (t117 * t4 - t3 * t71) * mrSges(5,3) + (t118 + t110 * (-Ifges(4,2) * t108 + t158)) * t144 / 0.2e1 + t35 * t140 + t21 * t135 + t76 * t136 + t67 * t145 / 0.2e1 - t66 * t146 / 0.2e1 + (mrSges(5,2) * t61 - mrSges(6,3) * t19 + Ifges(5,4) * t180 + Ifges(6,5) * t179 + t192 * t162 + t194 * t177 + t190 / 0.2e1 + t7 * mrSges(6,2) - t11 * mrSges(5,3)) * t65; m(3) * qJDD(2) + t108 * t59 + t110 * t60 + t189 * t71 - t188 * t117 + t161 * t65 - t160 * t63 + (-t108 * t80 + t110 * t82) * qJD(3) + m(4) * (t108 * t26 + t110 * t27 + (-t108 * t54 + t110 * t55) * qJD(3)) + m(5) * (-t11 * t63 + t117 * t3 + t12 * t65 + t4 * t71) + m(6) * (t1 * t71 - t117 * t2 + t63 * t7 + t65 * t8) + (-t202 + t196) * g(3); t191 * t39 + t192 * t40 - (t191 * t64 - t192 * t62) * qJD(3) / 0.2e1 - (-t194 * t62 - t170 + t22 + t57) * t64 / 0.2e1 + t160 * t15 - t161 * t16 + (-Ifges(5,2) * t64 + t190 - t58) * t179 + (Ifges(6,2) + Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (-t119 - t118 * qJD(1) / 0.2e1) * qJD(1) + t197 * (t122 + t202 * pkin(3) * t108 + (-m(6) * qJ(5) + t200) * t95 + (m(6) * pkin(4) + t201) * t93) + ((t104 * t4 + t152 * t3) * pkin(3) + t11 * t15 - t12 * t16 - t61 * t141) * m(5) + (t138 - t82) * t54 + t85 * t33 + t88 * t32 + Ifges(4,6) * t73 + Ifges(4,5) * t74 - t19 * (mrSges(6,1) * t64 + mrSges(6,3) * t62) - t61 * (mrSges(5,1) * t64 - mrSges(5,2) * t62) + qJD(5) * t53 - t20 * t34 - t26 * mrSges(4,2) + t27 * mrSges(4,1) + t1 * mrSges(6,3) - t2 * mrSges(6,1) + t3 * mrSges(5,1) - t4 * mrSges(5,2) + (t139 + t80) * t55 + (-m(6) * (t100 + t124) - m(5) * t100 + t198) * g(3) - (-Ifges(4,2) * t149 + t67 + t92) * t148 / 0.2e1 + t23 * t177 + (Ifges(6,3) * t64 - t169) * t180 + t64 * t165 + t30 * t167 + t8 * t171 + t7 * t172 + (t1 * t85 - t15 * t7 - t19 * t20 + t2 * t88 + (-t16 + qJD(5)) * t8) * m(6) - t35 * t141 - t120 * t144 / 0.2e1 + t66 * t149 / 0.2e1 + t31 * t137 - t11 * t164; t160 * t64 + t161 * t62 + t135 + t136 + (t62 * t8 - t64 * t7 + t128 + t5) * m(6) + (t11 * t64 + t12 * t62 + t128 + t46) * m(5); -qJD(3) * t53 + t64 * t34 + (g(3) * t95 - t8 * qJD(3) + t19 * t64 - t197 * t93 + t2) * m(6) + t32;];
tau = t6;
