% Calculate vector of inverse dynamics joint torques for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:44
% DurationCPUTime: 4.68s
% Computational Cost: add. (1710->309), mult. (3947->424), div. (0->0), fcn. (2522->8), ass. (0->152)
t150 = qJD(1) * qJD(2);
t149 = qJDD(1) * qJ(2);
t90 = t149 + t150;
t203 = t90 + t150;
t115 = sin(pkin(7));
t112 = t115 ^ 2;
t202 = t112 * t90;
t201 = pkin(3) + qJ(2);
t117 = cos(pkin(7));
t113 = t117 ^ 2;
t200 = t203 * t113;
t114 = sin(pkin(8));
t129 = mrSges(5,3) * t114 * t117 + mrSges(5,1) * t115;
t118 = sin(qJ(5));
t162 = t117 * t118;
t120 = cos(qJ(5));
t167 = t115 * t120;
t61 = t114 * t162 + t167;
t51 = t61 * qJD(1);
t161 = t117 * t120;
t144 = t114 * t161;
t156 = qJD(1) * t115;
t55 = qJD(1) * t144 - t118 * t156;
t176 = -mrSges(6,1) * t51 - mrSges(6,2) * t55 - t129 * qJD(1);
t116 = cos(pkin(8));
t137 = -qJ(3) * t115 - pkin(1);
t177 = -pkin(2) - qJ(4);
t67 = t177 * t117 + t137;
t42 = qJD(1) * t67 + qJD(2);
t88 = qJ(2) * t156 + qJD(3);
t70 = pkin(3) * t156 + t88;
t20 = -t114 * t42 + t116 * t70;
t18 = -pkin(4) * t156 - t20;
t199 = -m(5) * t20 + m(6) * t18 + t176;
t119 = sin(qJ(1));
t121 = cos(qJ(1));
t192 = -g(1) * t121 - g(2) * t119;
t198 = m(5) + m(6);
t52 = t61 * qJD(5);
t168 = t115 * t118;
t62 = t144 - t168;
t26 = qJD(1) * t52 - qJDD(1) * t62;
t197 = Ifges(6,5) * t26;
t53 = t62 * qJD(5);
t27 = qJD(1) * t53 + qJDD(1) * t61;
t196 = Ifges(6,6) * t27;
t147 = qJDD(1) * t117;
t140 = t116 * t147;
t81 = qJDD(5) + t140;
t195 = Ifges(6,3) * t81;
t194 = Ifges(3,4) + Ifges(4,6);
t21 = t114 * t70 + t116 * t42;
t19 = pkin(6) * t156 + t21;
t127 = (pkin(4) * t116 + pkin(6) * t114) * t117;
t155 = qJD(1) * t117;
t71 = t201 * t155 + qJD(4);
t36 = qJD(1) * t127 + t71;
t5 = -t118 * t19 + t120 * t36;
t6 = t118 * t36 + t120 * t19;
t193 = -t118 * t5 + t120 * t6;
t171 = t115 * mrSges(4,3);
t132 = t117 * mrSges(4,2) - t171;
t174 = mrSges(3,1) * t117;
t133 = -mrSges(3,2) * t115 + t174;
t191 = -t133 + t132 - mrSges(2,1);
t190 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t48 = pkin(3) * t147 + t117 * t90 + qJDD(4);
t30 = qJDD(1) * t127 + t48;
t152 = qJD(3) * t115;
t80 = -qJD(4) * t117 - t152;
t31 = qJD(1) * t80 + qJDD(1) * t67 + qJDD(2);
t148 = qJDD(1) * t115;
t73 = t115 * t90 + qJDD(3);
t47 = pkin(3) * t148 + t73;
t13 = t114 * t47 + t116 * t31;
t9 = pkin(6) * t148 + t13;
t1 = qJD(5) * t5 + t118 * t30 + t120 * t9;
t2 = -qJD(5) * t6 - t118 * t9 + t120 * t30;
t188 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t187 = t113 * t149 + t192 + t200 + t202;
t185 = -t55 / 0.2e1;
t59 = t129 * qJDD(1);
t7 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t183 = t59 - t7;
t182 = Ifges(6,4) * t55;
t83 = t201 * t115;
t33 = t114 * t83 + t116 * t67;
t175 = t200 * qJ(2);
t173 = mrSges(5,2) * t114;
t122 = qJD(1) ^ 2;
t170 = qJ(2) * t122;
t169 = t114 * t119;
t166 = t115 * t121;
t165 = t116 * t117;
t164 = t116 * t119;
t163 = t116 * t121;
t160 = t117 * t121;
t84 = t201 * t117;
t158 = t121 * pkin(1) + t119 * qJ(2);
t154 = qJD(2) * t115;
t153 = qJD(2) * t117;
t151 = m(4) + t198;
t146 = t195 + t196 + t197;
t141 = t114 * t147;
t136 = pkin(2) * t160 + qJ(3) * t166 + t158;
t135 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t12 = -t114 * t31 + t116 * t47;
t32 = -t114 * t67 + t116 * t83;
t29 = pkin(6) * t115 + t33;
t40 = t127 + t84;
t11 = t118 * t40 + t120 * t29;
t10 = -t118 * t29 + t120 * t40;
t82 = t116 * t155 + qJD(5);
t34 = -mrSges(6,2) * t82 + mrSges(6,3) * t51;
t35 = mrSges(6,1) * t82 + mrSges(6,3) * t55;
t131 = -t118 * t35 + t120 * t34;
t130 = t119 * pkin(3) + qJ(4) * t160 + t136;
t79 = -pkin(2) * t117 + t137;
t128 = -mrSges(5,2) * t115 - mrSges(5,3) * t165;
t126 = (mrSges(5,1) * t116 - t173) * t117;
t125 = t1 * t120 - t118 * t2 + (-t118 * t6 - t120 * t5) * qJD(5);
t14 = mrSges(6,1) * t81 - mrSges(6,3) * t26;
t15 = -mrSges(6,2) * t81 + mrSges(6,3) * t27;
t60 = t128 * qJDD(1);
t124 = -t118 * t14 + t120 * t15 + t60 + (-t118 * t34 - t120 * t35) * qJD(5);
t69 = t128 * qJD(1);
t123 = t199 * t114 + (m(5) * t21 + m(6) * t193 + t131 + t69) * t116;
t108 = t121 * qJ(2);
t104 = -qJDD(1) * pkin(1) + qJDD(2);
t99 = mrSges(4,2) * t147;
t98 = t113 * t170;
t96 = mrSges(3,2) * t148;
t85 = mrSges(5,1) * t140;
t74 = t132 * qJD(1);
t66 = -t115 * t169 + t163;
t64 = t114 * t166 + t164;
t58 = qJD(1) * t79 + qJD(2);
t56 = qJD(1) * t126;
t54 = (t114 * t167 + t162) * qJD(1);
t50 = (-t114 * t168 + t161) * qJD(1);
t46 = Ifges(6,4) * t51;
t45 = t114 * t154 + t116 * t80;
t41 = -qJD(1) * t152 + qJDD(1) * t79 + qJDD(2);
t38 = t118 * t160 + t120 * t64;
t37 = -t118 * t64 + t120 * t160;
t28 = -pkin(4) * t115 - t32;
t17 = -Ifges(6,1) * t55 + Ifges(6,5) * t82 + t46;
t16 = Ifges(6,2) * t51 + Ifges(6,6) * t82 - t182;
t8 = -pkin(4) * t148 - t12;
t4 = -qJD(5) * t11 - t118 * t45 + t120 * t153;
t3 = qJD(5) * t10 + t118 * t153 + t120 * t45;
t22 = [(mrSges(2,2) * t121 + (-m(6) * pkin(4) - mrSges(6,1) * t120 + mrSges(6,2) * t118 - mrSges(5,1)) * t66 - t198 * (t121 * pkin(3) + t108) + t190 * (t114 * t121 + t115 * t164) + (-m(3) - m(4)) * t108 + (m(3) * pkin(1) - m(4) * t79 - t198 * t137 + (t118 * mrSges(6,1) + t120 * mrSges(6,2) - t198 * t177 + mrSges(5,3)) * t117 - t191) * t119) * g(1) + (t146 / 0.2e1 + t195 / 0.2e1 + t196 / 0.2e1 + t197 / 0.2e1 + t188) * t165 + (-t79 * mrSges(4,3) + (-Ifges(5,5) * t114 - Ifges(5,6) * t116 + t194) * t117 + (Ifges(5,3) + Ifges(4,2) + Ifges(3,1)) * t115) * t148 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3)) * t117 + t194 * t115) * t147 + (-m(5) * t130 - t64 * mrSges(5,1) - mrSges(5,3) * t160 - m(6) * (pkin(4) * t64 + t130) - t38 * mrSges(6,1) - t37 * mrSges(6,2) - m(3) * t158 - m(4) * t136 + mrSges(2,2) * t119 + t190 * (-t115 * t163 + t169) + t191 * t121) * g(2) + (t73 * t115 + t187) * mrSges(4,1) + m(3) * (t203 * t112 * qJ(2) - pkin(1) * t104 + t175) + m(4) * (t41 * t79 + (qJ(2) * t73 + qJD(2) * t88 - qJD(3) * t58) * t115 + t175) - t74 * t152 - (Ifges(5,6) * t115 + (-Ifges(5,4) * t114 - Ifges(5,2) * t116) * t117) * t140 - (Ifges(5,5) * t115 + (-Ifges(5,1) * t114 - Ifges(5,4) * t116) * t117) * t141 + t84 * (-mrSges(5,2) * t141 + t85) - t8 * t135 + t41 * t132 - t104 * t133 + t13 * t128 + t12 * t129 + (mrSges(6,3) * t2 - Ifges(6,1) * t26 - Ifges(6,4) * t27 - Ifges(6,5) * t81) * t62 + m(6) * (t1 * t11 + t10 * t2 + t28 * t8 + t3 * t6 + t4 * t5) + (mrSges(6,3) * t1 + Ifges(6,4) * t26 + Ifges(6,2) * t27 + Ifges(6,6) * t81) * t61 + (-t5 * t52 + t6 * t53) * mrSges(6,3) + t199 * (t114 * t80 - t116 * t154) + (Ifges(6,1) * t52 + Ifges(6,4) * t53) * t185 + t56 * t153 + t48 * t126 - pkin(1) * t96 + m(5) * (t12 * t32 + t13 * t33 + t71 * t153 + t21 * t45 + t48 * t84) + (t187 + t202) * mrSges(3,3) + t45 * t69 + t82 * (Ifges(6,5) * t52 + Ifges(6,6) * t53) / 0.2e1 + t32 * t59 + t33 * t60 + t52 * t17 / 0.2e1 + t51 * (Ifges(6,4) * t52 + Ifges(6,2) * t53) / 0.2e1 + t18 * (-mrSges(6,1) * t53 + mrSges(6,2) * t52) + t53 * t16 / 0.2e1 + t28 * t7 + t3 * t34 + t4 * t35 + t11 * t15 + t10 * t14 + t79 * t99 + Ifges(2,3) * qJDD(1); -t56 * t155 - t54 * t34 - t50 * t35 + t96 + t99 + (-t171 - t174) * qJDD(1) + (-t156 * t69 - t183) * t114 + (t156 * t176 + t124) * t116 + (mrSges(4,1) + mrSges(3,3)) * t122 * (-t112 - t113) + (-g(1) * t119 + g(2) * t121) * (m(3) + t151) + (t114 * t8 - t5 * t50 - t54 * t6 + (t156 * t18 + t125) * t116) * m(6) + (-(t117 * t71 + (t114 * t21 + t116 * t20) * t115) * qJD(1) - t114 * t12 + t116 * t13) * m(5) + (-t156 * t88 + t41 - t98) * m(4) + (-t112 * t170 + t104 - t98) * m(3); t183 * t116 + t151 * t117 * g(3) + t124 * t114 + m(6) * (t114 * t125 - t116 * t8) + m(5) * (t114 * t13 + t116 * t12) + m(4) * t73 + (qJDD(1) * mrSges(4,1) + (m(4) * t58 + t123 + t74) * qJD(1) + t192 * t151) * t115; t118 * t15 + t120 * t14 + t85 + t131 * qJD(5) - t198 * t115 * g(3) + m(5) * t48 + m(6) * (t193 * qJD(5) + t1 * t118 + t120 * t2) + (qJD(1) * t123 - qJDD(1) * t173 + t198 * t192) * t117; -t18 * (-mrSges(6,1) * t55 + mrSges(6,2) * t51) + t55 * (Ifges(6,1) * t51 + t182) / 0.2e1 + t16 * t185 - t82 * (Ifges(6,5) * t51 + Ifges(6,6) * t55) / 0.2e1 - t5 * t34 + t6 * t35 - g(1) * (mrSges(6,1) * t37 - mrSges(6,2) * t38) - g(2) * ((t118 * t66 + t119 * t161) * mrSges(6,1) + (-t119 * t162 + t120 * t66) * mrSges(6,2)) - g(3) * t135 + (t5 * t51 - t55 * t6) * mrSges(6,3) + t146 - (Ifges(6,2) * t55 + t17 + t46) * t51 / 0.2e1 + t188;];
tau = t22;
