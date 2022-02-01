% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:20
% EndTime: 2022-01-20 09:51:25
% DurationCPUTime: 1.97s
% Computational Cost: add. (2407->267), mult. (3901->355), div. (0->0), fcn. (2356->16), ass. (0->134)
t149 = sin(pkin(8));
t154 = sin(qJ(2));
t191 = pkin(1) * qJD(1);
t177 = t154 * t191;
t106 = t149 * t177;
t151 = cos(pkin(8));
t157 = cos(qJ(2));
t176 = t157 * t191;
t80 = t151 * t176 - t106;
t204 = qJD(4) - t80;
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t119 = pkin(2) * t149 + qJ(4);
t148 = sin(pkin(9));
t90 = (-pkin(7) - t119) * t148;
t150 = cos(pkin(9));
t139 = t150 * pkin(7);
t91 = t119 * t150 + t139;
t54 = t153 * t90 + t156 * t91;
t97 = t148 * t156 + t150 * t153;
t206 = -t54 * qJD(5) - t204 * t97;
t53 = -t153 * t91 + t156 * t90;
t96 = -t148 * t153 + t150 * t156;
t205 = t53 * qJD(5) + t204 * t96;
t131 = t150 * qJDD(3);
t142 = qJDD(1) + qJDD(2);
t146 = qJD(1) + qJD(2);
t190 = pkin(1) * qJD(2);
t175 = t154 * t190;
t195 = pkin(1) * t157;
t94 = -qJD(1) * t175 + qJDD(1) * t195;
t75 = pkin(2) * t142 + t94;
t179 = qJD(2) * t157;
t95 = (qJD(1) * t179 + qJDD(1) * t154) * pkin(1);
t48 = t149 * t75 + t151 * t95;
t32 = qJ(4) * t142 + qJD(4) * t146 + t48;
t23 = -t148 * t32 + t131;
t24 = t148 * qJDD(3) + t150 * t32;
t203 = -t23 * t148 + t150 * t24;
t167 = -t150 * mrSges(5,1) + t148 * mrSges(5,2);
t145 = pkin(9) + qJ(5);
t134 = sin(t145);
t135 = cos(t145);
t202 = t135 * mrSges(6,1) - t134 * mrSges(6,2);
t201 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t200 = -mrSges(4,1) - t202 + t167;
t173 = mrSges(5,3) * (t148 ^ 2 + t150 ^ 2);
t199 = m(3) * pkin(1);
t71 = t97 * t146;
t197 = t71 / 0.2e1;
t196 = Ifges(6,4) * t71;
t147 = qJ(1) + qJ(2);
t138 = cos(t147);
t125 = pkin(2) * t138;
t194 = pkin(2) * t151;
t193 = pkin(4) * t150;
t70 = t96 * t146;
t192 = -mrSges(6,1) * t70 + mrSges(6,2) * t71 + t167 * t146;
t101 = pkin(2) * t146 + t176;
t64 = t149 * t101 + t151 * t177;
t58 = qJ(4) * t146 + t64;
t51 = t148 * qJD(3) + t150 * t58;
t183 = t142 * t148;
t182 = t142 * t150;
t181 = t151 * t154;
t126 = pkin(2) + t195;
t84 = pkin(1) * t181 + t149 * t126;
t178 = m(4) + m(5) + m(6);
t136 = pkin(8) + t147;
t121 = sin(t136);
t122 = cos(t136);
t174 = t122 * pkin(3) + t121 * qJ(4) + t125;
t123 = pkin(3) + t193;
t88 = t96 * qJD(5);
t44 = t97 * t142 + t146 * t88;
t89 = t97 * qJD(5);
t45 = t96 * t142 - t146 * t89;
t11 = -t45 * mrSges(6,1) + t44 * mrSges(6,2);
t47 = -t149 * t95 + t151 * t75;
t63 = t101 * t151 - t106;
t83 = -t149 * t154 * pkin(1) + t126 * t151;
t77 = -pkin(3) - t83;
t82 = -mrSges(5,1) * t182 + mrSges(5,2) * t183;
t152 = -pkin(7) - qJ(4);
t172 = -t121 * t152 + t122 * t123 + t125;
t171 = -g(1) * t121 + g(2) * t122;
t170 = qJDD(4) - t47;
t169 = qJD(4) - t63;
t168 = mrSges(3,1) * t154 + mrSges(3,2) * t157;
t133 = t150 * qJD(3);
t50 = -t148 * t58 + t133;
t165 = -t148 * t50 + t150 * t51;
t40 = t133 + (-pkin(7) * t146 - t58) * t148;
t41 = t146 * t139 + t51;
t9 = -t153 * t41 + t156 * t40;
t10 = t153 * t40 + t156 * t41;
t76 = qJ(4) + t84;
t59 = (-pkin(7) - t76) * t148;
t60 = t150 * t76 + t139;
t21 = -t153 * t60 + t156 * t59;
t22 = t153 * t59 + t156 * t60;
t81 = t151 * pkin(1) * t179 - t149 * t175;
t164 = pkin(1) * (t149 * t157 + t181);
t137 = sin(t147);
t161 = -t138 * mrSges(3,1) + t137 * mrSges(3,2) + t201 * t121 + t200 * t122;
t160 = t138 * mrSges(3,2) + (t178 * pkin(2) + mrSges(3,1)) * t137 + (m(5) * pkin(3) + m(6) * t123 - t200) * t121 + (-m(5) * qJ(4) + m(6) * t152 + t201) * t122;
t16 = t131 + (-pkin(7) * t142 - t32) * t148;
t17 = pkin(7) * t182 + t24;
t2 = t9 * qJD(5) + t153 * t16 + t156 * t17;
t3 = -t10 * qJD(5) - t153 * t17 + t156 * t16;
t30 = -t123 * t142 + t170;
t37 = Ifges(6,2) * t70 + Ifges(6,6) * qJD(5) + t196;
t69 = Ifges(6,4) * t70;
t38 = Ifges(6,1) * t71 + Ifges(6,5) * qJD(5) + t69;
t39 = -pkin(3) * t142 + t170;
t52 = -t123 * t146 + t169;
t159 = (Ifges(6,1) * t88 - Ifges(6,4) * t89) * t197 + qJD(5) * (Ifges(6,5) * t88 - Ifges(6,6) * t89) / 0.2e1 + t70 * (Ifges(6,4) * t88 - Ifges(6,2) * t89) / 0.2e1 + t52 * (mrSges(6,1) * t89 + mrSges(6,2) * t88) + t39 * t167 + t94 * mrSges(3,1) - t95 * mrSges(3,2) + t88 * t38 / 0.2e1 - t89 * t37 / 0.2e1 + t47 * mrSges(4,1) - t48 * mrSges(4,2) + (Ifges(5,4) * t148 + Ifges(5,2) * t150) * t182 + (Ifges(5,1) * t148 + Ifges(5,4) * t150) * t183 + (Ifges(3,3) + Ifges(4,3)) * t142 + (-t10 * t89 - t9 * t88) * mrSges(6,3) + t203 * mrSges(5,3) + (t30 * mrSges(6,2) - t3 * mrSges(6,3) + Ifges(6,1) * t44 + Ifges(6,4) * t45 + Ifges(6,5) * qJDD(5)) * t97 + (-t30 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * qJDD(5)) * t96;
t158 = cos(qJ(1));
t155 = sin(qJ(1));
t140 = t158 * pkin(1);
t124 = -pkin(3) - t194;
t104 = -t123 - t194;
t79 = qJD(2) * t164;
t78 = qJD(1) * t164;
t72 = qJD(4) + t81;
t67 = t77 - t193;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t61 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t57 = -pkin(3) * t146 + t169;
t34 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t45;
t33 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t44;
t6 = -t22 * qJD(5) - t97 * t72;
t5 = t21 * qJD(5) + t96 * t72;
t1 = [(t158 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t178) * pkin(1)) * t155 + t160) * g(1) + (-m(6) * (t140 + t172) + t155 * mrSges(2,2) - m(5) * (t140 + t174) - m(4) * (t125 + t140) + (-mrSges(2,1) - t199) * t158 + t161) * g(2) + (t154 * t95 + t157 * t94) * t199 + t192 * t79 + m(5) * (t39 * t77 + t57 * t79 + (t24 * t76 + t51 * t72) * t150 + (-t23 * t76 - t50 * t72) * t148) + t159 + t77 * t82 + t5 * t61 + t6 * t62 + t67 * t11 + t21 * t33 + t22 * t34 + (mrSges(4,1) * t83 - mrSges(4,2) * t84 + (mrSges(3,1) * t157 - mrSges(3,2) * t154) * pkin(1) + t76 * t173) * t142 + (-t79 * mrSges(4,1) - t81 * mrSges(4,2) - t168 * t190 + t72 * t173) * t146 + Ifges(2,3) * qJDD(1) + m(6) * (t10 * t5 + t2 * t22 + t21 * t3 + t30 * t67 + t52 * t79 + t6 * t9) + m(4) * (t47 * t83 + t48 * t84 - t63 * t79 + t64 * t81); t124 * t82 + t159 - t192 * t78 + t104 * t11 + t53 * t33 + t54 * t34 + (t78 * mrSges(4,1) + t80 * mrSges(4,2) + t168 * t191 + t204 * t173) * t146 + t160 * g(1) + t161 * g(2) + ((mrSges(4,1) * t151 - mrSges(4,2) * t149) * pkin(2) + t119 * t173) * t142 + t205 * t61 + t206 * t62 + (-t172 * g(2) + t205 * t10 + t104 * t30 + t2 * t54 + t206 * t9 + t3 * t53 - t52 * t78) * m(6) + (-t174 * g(2) + t203 * t119 + t124 * t39 + t204 * t165 - t57 * t78) * m(5) + ((t149 * t48 + t151 * t47) * pkin(2) + t63 * t78 - t64 * t80 - t125 * g(2)) * m(4); m(4) * qJDD(3) + t96 * t33 + t97 * t34 + t88 * t61 - t89 * t62 + m(5) * (t148 * t24 + t150 * t23) + m(6) * (t10 * t88 + t2 * t97 + t3 * t96 - t89 * t9) - t178 * g(3); -t146 ^ 2 * t173 - t70 * t61 + t71 * t62 + t11 + t82 + (-t10 * t70 + t9 * t71 + t171 + t30) * m(6) + (-t165 * t146 + t171 + t39) * m(5); Ifges(6,5) * t44 + Ifges(6,6) * t45 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t52 * (mrSges(6,1) * t71 + mrSges(6,2) * t70) - t71 * (Ifges(6,1) * t70 - t196) / 0.2e1 + t37 * t197 - qJD(5) * (Ifges(6,5) * t70 - Ifges(6,6) * t71) / 0.2e1 - t9 * t61 + t10 * t62 - g(3) * t202 + (t10 * t71 + t70 * t9) * mrSges(6,3) - (-Ifges(6,2) * t71 + t38 + t69) * t70 / 0.2e1 + (g(1) * t122 + g(2) * t121) * (mrSges(6,1) * t134 + mrSges(6,2) * t135);];
tau = t1;
