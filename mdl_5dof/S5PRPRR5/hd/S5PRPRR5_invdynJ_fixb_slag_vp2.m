% Calculate vector of inverse dynamics joint torques for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:27
% EndTime: 2019-12-05 15:53:43
% DurationCPUTime: 5.88s
% Computational Cost: add. (3000->357), mult. (7133->489), div. (0->0), fcn. (5247->14), ass. (0->161)
t133 = sin(pkin(9));
t195 = pkin(6) + qJ(3);
t111 = t195 * t133;
t135 = cos(pkin(9));
t112 = t195 * t135;
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t186 = t135 * t142;
t157 = t133 * t139 - t186;
t180 = qJD(4) * t142;
t143 = cos(qJ(2));
t182 = qJD(1) * t143;
t222 = -t111 * t180 + qJD(3) * t186 + (-qJD(3) * t133 - qJD(4) * t112) * t139 + t157 * t182;
t106 = t133 * t142 + t135 * t139;
t152 = t106 * t143;
t74 = -t139 * t111 + t142 * t112;
t221 = qJD(1) * t152 - t106 * qJD(3) - qJD(4) * t74;
t98 = t106 * qJD(4);
t230 = -pkin(7) * t98 + t222;
t97 = t157 * qJD(4);
t229 = pkin(7) * t97 + t221;
t184 = t133 ^ 2 + t135 ^ 2;
t220 = mrSges(4,3) * t184;
t228 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t131 = pkin(9) + qJ(4);
t126 = qJ(5) + t131;
t117 = sin(t126);
t118 = cos(t126);
t119 = pkin(3) * t135 + pkin(2);
t123 = sin(t131);
t124 = cos(t131);
t159 = -mrSges(4,1) * t135 + mrSges(4,2) * t133;
t227 = mrSges(3,1) + m(6) * (pkin(4) * t124 + t119) + mrSges(6,1) * t118 - mrSges(6,2) * t117 + m(5) * t119 + mrSges(5,1) * t124 - mrSges(5,2) * t123 + m(4) * pkin(2) - t159;
t226 = mrSges(3,2) + m(6) * (-pkin(7) - t195) - mrSges(6,3) - m(5) * t195 - mrSges(5,3) - m(4) * qJ(3) - mrSges(4,3);
t134 = sin(pkin(8));
t136 = cos(pkin(8));
t225 = g(1) * t136 + g(2) * t134;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t73 = -t142 * t111 - t112 * t139;
t48 = -pkin(7) * t106 + t73;
t49 = -pkin(7) * t157 + t74;
t25 = t138 * t48 + t141 * t49;
t224 = -qJD(5) * t25 - t230 * t138 + t229 * t141;
t24 = -t138 * t49 + t141 * t48;
t223 = qJD(5) * t24 + t229 * t138 + t230 * t141;
t214 = m(6) * pkin(4);
t219 = mrSges(5,1) + t214;
t95 = t157 * qJD(2);
t96 = t106 * qJD(2);
t168 = -t138 * t96 - t141 * t95;
t174 = qJDD(2) * t135;
t175 = qJDD(2) * t133;
t104 = -mrSges(4,1) * t174 + mrSges(4,2) * t175;
t63 = -qJD(2) * t97 + qJDD(2) * t106;
t64 = -qJD(2) * t98 - qJDD(2) * t157;
t31 = -t64 * mrSges(5,1) + t63 * mrSges(5,2);
t20 = qJD(5) * t168 + t138 * t64 + t141 * t63;
t60 = -t138 * t95 + t141 * t96;
t21 = -qJD(5) * t60 - t138 * t63 + t141 * t64;
t4 = -t21 * mrSges(6,1) + t20 * mrSges(6,2);
t218 = t104 + t31 + t4;
t28 = -mrSges(6,1) * t168 + mrSges(6,2) * t60;
t217 = mrSges(5,1) * t95 + mrSges(5,2) * t96 + t159 * qJD(2) + t28;
t213 = -t168 / 0.2e1;
t212 = -t60 / 0.2e1;
t211 = t60 / 0.2e1;
t209 = t96 / 0.2e1;
t208 = pkin(4) * t96;
t207 = pkin(4) * t98;
t132 = qJD(4) + qJD(5);
t206 = -t132 / 0.2e1;
t205 = Ifges(5,4) * t96;
t204 = Ifges(6,4) * t60;
t140 = sin(qJ(2));
t201 = g(3) * t140;
t183 = qJD(1) * t140;
t114 = qJD(2) * qJ(3) + t183;
t166 = pkin(6) * qJD(2) + t114;
t88 = t166 * t133;
t89 = t166 * t135;
t45 = -t139 * t88 + t142 * t89;
t36 = -pkin(7) * t95 + t45;
t191 = t138 * t36;
t189 = t139 * t89;
t44 = -t142 * t88 - t189;
t35 = -pkin(7) * t96 + t44;
t32 = qJD(4) * pkin(4) + t35;
t11 = t141 * t32 - t191;
t200 = t11 * mrSges(6,3);
t188 = t141 * t36;
t12 = t138 * t32 + t188;
t199 = t12 * mrSges(6,3);
t198 = t44 * mrSges(5,3);
t197 = t45 * mrSges(5,3);
t187 = t134 * t143;
t194 = (-t117 * t187 - t118 * t136) * mrSges(6,1) + (t117 * t136 - t118 * t187) * mrSges(6,2);
t185 = t136 * t143;
t193 = (-t117 * t185 + t118 * t134) * mrSges(6,1) + (-t117 * t134 - t118 * t185) * mrSges(6,2);
t178 = qJD(1) * qJD(2);
t121 = t143 * t178;
t110 = t140 * qJDD(1) + t121;
t181 = qJD(2) * t140;
t179 = m(4) + m(5) + m(6);
t120 = t140 * t178;
t90 = t110 + t228;
t167 = t184 * t90;
t165 = pkin(6) * qJDD(2) + t90;
t164 = t184 * t114;
t163 = qJD(3) - t182;
t109 = qJDD(1) * t143 - t120;
t161 = -mrSges(3,2) + t220;
t158 = -mrSges(6,1) * t117 - mrSges(6,2) * t118;
t91 = t106 * t140;
t92 = t157 * t140;
t46 = t138 * t92 - t141 * t91;
t47 = -t138 * t91 - t141 * t92;
t65 = -t106 * t138 - t141 * t157;
t66 = t106 * t141 - t138 * t157;
t156 = qJDD(3) - t109;
t127 = qJDD(4) + qJDD(5);
t71 = t165 * t133;
t72 = t165 * t135;
t23 = -qJD(4) * t45 - t139 * t72 - t142 * t71;
t7 = qJDD(4) * pkin(4) - pkin(7) * t63 + t23;
t22 = -qJD(4) * t189 - t139 * t71 + t142 * t72 - t88 * t180;
t8 = pkin(7) * t64 + t22;
t2 = qJD(5) * t11 + t138 * t7 + t141 * t8;
t3 = -qJD(5) * t12 - t138 * t8 + t141 * t7;
t155 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t20 + Ifges(6,6) * t21 + Ifges(6,3) * t127;
t148 = (-qJD(2) * pkin(2) + t163) * t140 + t143 * t164;
t99 = -qJD(2) * t119 + t163;
t80 = -qJDD(2) * t119 + t156;
t144 = qJD(2) ^ 2;
t94 = -qJDD(2) * pkin(2) + t156;
t93 = Ifges(5,4) * t95;
t86 = pkin(4) * t157 - t119;
t83 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t96;
t82 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t95;
t67 = pkin(4) * t95 + t99;
t56 = t96 * Ifges(5,1) + Ifges(5,5) * qJD(4) - t93;
t55 = -t95 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t205;
t54 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t64;
t53 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t63;
t52 = Ifges(6,4) * t168;
t51 = -qJD(2) * t152 + t140 * t97;
t50 = -t140 * t98 - t143 * t95;
t41 = mrSges(6,1) * t132 - mrSges(6,3) * t60;
t40 = -mrSges(6,2) * t132 + mrSges(6,3) * t168;
t37 = -pkin(4) * t64 + t80;
t30 = -qJD(5) * t66 + t138 * t97 - t141 * t98;
t29 = qJD(5) * t65 - t138 * t98 - t141 * t97;
t27 = Ifges(6,1) * t60 + Ifges(6,5) * t132 + t52;
t26 = Ifges(6,2) * t168 + Ifges(6,6) * t132 + t204;
t16 = -mrSges(6,2) * t127 + mrSges(6,3) * t21;
t15 = mrSges(6,1) * t127 - mrSges(6,3) * t20;
t14 = t141 * t35 - t191;
t13 = -t138 * t35 - t188;
t10 = -qJD(5) * t47 - t138 * t50 + t141 * t51;
t9 = qJD(5) * t46 + t138 * t51 + t141 * t50;
t1 = [m(2) * qJDD(1) + t10 * t41 + t46 * t15 + t47 * t16 + t9 * t40 + t50 * t82 + t51 * t83 - t91 * t53 - t92 * t54 + (-m(2) - m(3) - t179) * g(3) + (qJDD(2) * mrSges(3,1) + t144 * t161 - t218) * t143 + (-t144 * mrSges(3,1) + qJD(2) * t217 + t161 * qJDD(2)) * t140 + m(4) * (qJD(2) * t148 + t140 * t167 - t143 * t94) + m(5) * (-t143 * t80 + t181 * t99 - t22 * t92 - t23 * t91 + t44 * t51 + t45 * t50) + m(3) * (t109 * t143 + t110 * t140) + m(6) * (t10 * t11 + t12 * t9 - t143 * t37 + t181 * t67 + t2 * t47 + t3 * t46); (t226 * g(3) + t225 * t227) * t140 + (-t227 * g(3) + t225 * t226) * t143 + (Ifges(6,1) * t29 + Ifges(6,4) * t30) * t211 - t98 * t197 + t30 * t199 + (Ifges(4,4) * t133 + Ifges(4,2) * t135) * t174 + (Ifges(4,1) * t133 + Ifges(4,4) * t135) * t175 + t132 * (Ifges(6,5) * t29 + Ifges(6,6) * t30) / 0.2e1 + t223 * t40 + (t11 * t224 + t223 * t12 + t2 * t25 + t207 * t67 + t24 * t3 + t37 * t86) * m(6) + t224 * t41 - t119 * t31 - pkin(2) * t104 - t97 * t56 / 0.2e1 - t98 * t55 / 0.2e1 + t86 * t4 + t67 * (-mrSges(6,1) * t30 + mrSges(6,2) * t29) + t73 * t53 + t74 * t54 + t30 * t26 / 0.2e1 + t24 * t15 + t25 * t16 + t29 * t27 / 0.2e1 + t221 * t83 + (-t119 * t80 + t22 * t74 + t221 * t44 + t222 * t45 + t23 * t73) * m(5) + t222 * t82 + t168 * (Ifges(6,4) * t29 + Ifges(6,2) * t30) / 0.2e1 + (-Ifges(5,1) * t97 - Ifges(5,4) * t98) * t209 - t95 * (-Ifges(5,4) * t97 - Ifges(5,2) * t98) / 0.2e1 + qJD(4) * (-Ifges(5,5) * t97 - Ifges(5,6) * t98) / 0.2e1 + t99 * (mrSges(5,1) * t98 - mrSges(5,2) * t97) + (-m(5) * t99 - m(6) * t67 - t217) * t183 + (-t121 + t90 + t228) * t220 + (mrSges(5,2) * t80 - mrSges(5,3) * t23 + Ifges(5,1) * t63 + Ifges(5,4) * t64 + Ifges(5,5) * qJDD(4)) * t106 + t94 * t159 + t28 * t207 - t29 * t200 + t97 * t198 - (-mrSges(5,1) * t80 + mrSges(5,3) * t22 + Ifges(5,4) * t63 + Ifges(5,2) * t64 + Ifges(5,6) * qJDD(4)) * t157 + (-pkin(2) * t94 + qJ(3) * t167 - t148 * qJD(1) + qJD(3) * t164) * m(4) + (mrSges(6,2) * t37 - mrSges(6,3) * t3 + Ifges(6,1) * t20 + Ifges(6,4) * t21 + Ifges(6,5) * t127) * t66 + (-mrSges(6,1) * t37 + mrSges(6,3) * t2 + Ifges(6,4) * t20 + Ifges(6,2) * t21 + Ifges(6,6) * t127) * t65 + (-t110 + t121) * mrSges(3,2) + (t109 + t120) * mrSges(3,1) + Ifges(3,3) * qJDD(2); -t144 * t220 - t168 * t40 + t60 * t41 + t95 * t82 + t96 * t83 + (t11 * t60 - t12 * t168 + t37) * m(6) + (t44 * t96 + t45 * t95 + t80) * m(5) + (-qJD(2) * t164 + t94) * m(4) + (t143 * g(3) - t140 * t225) * t179 + t218; t55 * t209 + (t138 * t2 + t141 * t3 + (-t11 * t138 + t12 * t141) * qJD(5)) * t214 + t96 * t197 - t95 * t198 - t44 * t82 + t45 * t83 + Ifges(5,5) * t63 + Ifges(5,6) * t64 - t14 * t40 - t13 * t41 - t22 * mrSges(5,2) + t23 * mrSges(5,1) + t155 + (Ifges(6,5) * t206 + Ifges(6,1) * t212 + Ifges(6,4) * t213 + t200 - t67 * mrSges(6,2) - t27 / 0.2e1) * t168 - qJD(4) * (-Ifges(5,5) * t95 - Ifges(5,6) * t96) / 0.2e1 - t99 * (mrSges(5,1) * t96 - mrSges(5,2) * t95) - t96 * (-Ifges(5,1) * t95 - t205) / 0.2e1 + (mrSges(5,2) * t124 + t123 * t219 - t158) * t201 + (-(-t123 * t134 - t124 * t185) * mrSges(5,2) - t193 - t219 * (-t123 * t185 + t124 * t134)) * g(1) + (-(t123 * t136 - t124 * t187) * mrSges(5,2) - t194 - t219 * (-t123 * t187 - t124 * t136)) * g(2) - t28 * t208 - m(6) * (t11 * t13 + t12 * t14 + t208 * t67) - (Ifges(6,6) * t206 + Ifges(6,4) * t212 + Ifges(6,2) * t213 + t67 * mrSges(6,1) - t26 / 0.2e1 - t199) * t60 + (-Ifges(5,2) * t96 + t56 - t93) * t95 / 0.2e1 + Ifges(5,3) * qJDD(4) + ((-t138 * t41 + t141 * t40) * qJD(5) + t138 * t16 + t141 * t15) * pkin(4); -t67 * (mrSges(6,1) * t60 + mrSges(6,2) * t168) + (Ifges(6,1) * t168 - t204) * t212 + t26 * t211 + (Ifges(6,5) * t168 - Ifges(6,6) * t60) * t206 - t11 * t40 + t12 * t41 - g(1) * t193 - g(2) * t194 - t158 * t201 + (t11 * t168 + t12 * t60) * mrSges(6,3) + t155 + (-Ifges(6,2) * t60 + t27 + t52) * t213;];
tau = t1;
