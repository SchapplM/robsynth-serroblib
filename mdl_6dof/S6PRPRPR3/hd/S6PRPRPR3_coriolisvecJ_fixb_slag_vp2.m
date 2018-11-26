% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:56:30
% EndTime: 2018-11-23 14:56:34
% DurationCPUTime: 4.11s
% Computational Cost: add. (2610->362), mult. (6579->499), div. (0->0), fcn. (4396->10), ass. (0->188)
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t135 = pkin(9) * t115 - qJ(5) * t118;
t128 = t135 * qJD(4);
t168 = qJD(5) * t115;
t171 = qJD(4) * t115;
t149 = pkin(4) * t171 - t168;
t111 = sin(pkin(11));
t112 = sin(pkin(6));
t113 = cos(pkin(11));
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t66 = (t111 * t119 + t113 * t116) * t112;
t62 = qJD(1) * t66;
t245 = t128 + t149 - t62;
t244 = mrSges(6,2) - mrSges(5,1);
t243 = -mrSges(5,3) - mrSges(6,1);
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t172 = qJD(2) * t118;
t83 = -qJD(4) * t114 - t117 * t172;
t170 = qJD(4) * t117;
t84 = -t114 * t172 + t170;
t44 = -mrSges(7,1) * t83 + mrSges(7,2) * t84;
t94 = -mrSges(6,1) * t172 - qJD(4) * mrSges(6,3);
t192 = -t94 + t44;
t180 = cos(pkin(6));
t98 = qJD(1) * t180 + qJD(3);
t182 = t115 * t98;
t174 = qJD(1) * t112;
t158 = t116 * t174;
t157 = t119 * t174;
t89 = qJD(2) * pkin(2) + t157;
t51 = t111 * t89 + t113 * t158;
t48 = qJD(2) * pkin(8) + t51;
t35 = t118 * t48 + t182;
t31 = -qJD(4) * qJ(5) - t35;
t224 = -m(5) * t35 + m(6) * t31;
t107 = pkin(5) * t172;
t25 = t107 - t31;
t93 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t172;
t242 = m(7) * t25 + t192 - t224 + t93;
t142 = mrSges(7,1) * t117 - mrSges(7,2) * t114;
t120 = -pkin(4) - pkin(9);
t151 = pkin(5) * qJD(2) + t48;
t143 = t151 * t115;
t91 = t118 * t98;
t28 = t91 - t143;
t237 = qJD(5) - t28;
t20 = qJD(4) * t120 + t237;
t152 = -qJ(5) * t115 - pkin(3);
t126 = t118 * t120 + t152;
t90 = t111 * t158;
t50 = t113 * t89 - t90;
t32 = qJD(2) * t126 - t50;
t5 = -t114 * t32 + t117 * t20;
t6 = t114 * t20 + t117 * t32;
t144 = t114 * t5 - t117 * t6;
t205 = -t117 / 0.2e1;
t206 = -t114 / 0.2e1;
t173 = qJD(2) * t115;
t101 = qJD(6) + t173;
t204 = Ifges(7,4) * t84;
t40 = Ifges(7,2) * t83 + Ifges(7,6) * t101 + t204;
t80 = Ifges(7,4) * t83;
t41 = Ifges(7,1) * t84 + Ifges(7,5) * t101 + t80;
t241 = -t144 * mrSges(7,3) - t25 * t142 - t205 * t40 - t206 * t41;
t176 = t115 * t117;
t203 = pkin(2) * t113;
t67 = t126 - t203;
t103 = pkin(2) * t111 + pkin(8);
t196 = pkin(5) + t103;
t81 = t196 * t115;
t38 = t114 * t81 + t117 * t67;
t64 = t113 * t157 - t90;
t82 = t196 * t118;
t70 = qJD(4) * t82;
t240 = -qJD(6) * t38 - t114 * t245 + t117 * t70 - t176 * t64;
t177 = t114 * t115;
t37 = -t114 * t67 + t117 * t81;
t239 = qJD(6) * t37 + t114 * t70 + t117 * t245 - t177 * t64;
t228 = (t111 * t116 - t113 * t119) * t112;
t125 = qJD(2) * t228;
t56 = qJD(1) * t125;
t183 = t115 * t56;
t10 = -t183 + (t118 * t151 + t182) * qJD(4);
t166 = qJD(4) * qJD(2);
t154 = t115 * t166;
t100 = pkin(4) * t154;
t124 = t62 - t168;
t27 = t100 + (t128 + t124) * qJD(2);
t1 = qJD(6) * t5 + t10 * t114 + t117 * t27;
t2 = -qJD(6) * t6 + t10 * t117 - t114 * t27;
t225 = t1 * t114 + t117 * t2;
t238 = m(7) * t225;
t186 = Ifges(7,6) * t117;
t189 = Ifges(7,5) * t114;
t136 = t186 + t189;
t191 = Ifges(7,4) * t114;
t137 = Ifges(7,2) * t117 + t191;
t190 = Ifges(7,4) * t117;
t139 = Ifges(7,1) * t114 + t190;
t207 = t101 / 0.2e1;
t211 = t84 / 0.2e1;
t230 = qJD(4) / 0.2e1;
t231 = -qJD(4) / 0.2e1;
t232 = -qJD(2) / 0.2e1;
t131 = -pkin(4) * t118 + t152;
t36 = qJD(2) * t131 - t50;
t47 = -qJD(2) * pkin(3) - t50;
t236 = t31 * mrSges(6,1) + t47 * mrSges(5,1) + Ifges(5,6) * t231 + (Ifges(5,4) * t115 + t118 * Ifges(5,2)) * t232 + Ifges(6,5) * t230 + (-Ifges(6,6) * t115 - t118 * Ifges(6,3)) * qJD(2) / 0.2e1 + t136 * t207 - t35 * mrSges(5,3) - t36 * mrSges(6,2) + t83 * t137 / 0.2e1 + t139 * t211 + t241;
t194 = -qJD(4) * t244 + t173 * t243;
t184 = t115 * t48;
t34 = -t91 + t184;
t227 = -qJD(5) - t34;
t30 = -qJD(4) * pkin(4) - t227;
t217 = m(5) * t34 + m(6) * t30 - t194;
t85 = (mrSges(6,2) * t118 - mrSges(6,3) * t115) * qJD(2);
t235 = -t85 + (mrSges(5,1) * t118 - mrSges(5,2) * t115 + mrSges(4,1)) * qJD(2) - m(6) * t36;
t167 = qJD(6) * t118;
t155 = t117 * t167;
t165 = qJD(4) * qJD(6);
t57 = -t114 * t165 + (t114 * t171 - t155) * qJD(2);
t156 = t114 * t167;
t58 = -t117 * t165 + (t115 * t170 + t156) * qJD(2);
t223 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t57 + Ifges(7,6) * t58;
t59 = -mrSges(7,2) * t101 + mrSges(7,3) * t83;
t60 = mrSges(7,1) * t101 - mrSges(7,3) * t84;
t133 = t114 * t60 - t117 * t59;
t153 = t118 * t166;
t42 = mrSges(7,1) * t153 - mrSges(7,3) * t57;
t43 = -mrSges(7,2) * t153 + mrSges(7,3) * t58;
t134 = t114 * t43 + t117 * t42;
t220 = t133 * qJD(6) - t134;
t169 = qJD(4) * t118;
t195 = -t118 * t56 + t169 * t98;
t11 = (-qJD(5) + t184) * qJD(4) - t195;
t12 = -t171 * t48 + t195;
t219 = m(5) * t12 - m(6) * t11;
t215 = 0.2e1 * m(5);
t214 = -t64 / 0.2e1;
t213 = -t83 / 0.2e1;
t212 = -t84 / 0.2e1;
t208 = -t101 / 0.2e1;
t13 = qJD(4) * t35 - t183;
t150 = t180 * t118;
t45 = t115 * t66 - t150;
t200 = t13 * t45;
t63 = qJD(2) * t66;
t55 = qJD(1) * t63;
t199 = t55 * t228;
t193 = t93 - t94;
t188 = Ifges(7,5) * t117;
t187 = Ifges(7,6) * t114;
t163 = -0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(6,6);
t162 = Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t161 = -Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1;
t160 = t103 * t13 / 0.2e1;
t159 = qJ(5) * t169;
t141 = mrSges(7,1) * t114 + mrSges(7,2) * t117;
t140 = Ifges(7,1) * t117 - t191;
t138 = -Ifges(7,2) * t114 + t190;
t21 = -t114 * t228 + t117 * t45;
t22 = t114 * t45 + t117 * t228;
t46 = t115 * t180 + t118 * t66;
t105 = Ifges(5,4) * t172;
t123 = t30 * mrSges(6,1) + t34 * mrSges(5,3) + t47 * mrSges(5,2) + t5 * mrSges(7,1) + t101 * Ifges(7,3) + t84 * Ifges(7,5) + t83 * Ifges(7,6) + Ifges(5,1) * t173 / 0.2e1 + Ifges(5,5) * t230 + t105 / 0.2e1 + Ifges(6,4) * t231 + (-t115 * Ifges(6,2) - Ifges(6,6) * t118) * t232 - t36 * mrSges(6,3) - t6 * mrSges(7,2);
t106 = pkin(4) * t173;
t99 = Ifges(7,3) * t153;
t87 = -qJ(5) * t172 + t106;
t79 = t131 - t203;
t77 = (mrSges(5,1) * t115 + mrSges(5,2) * t118) * t166;
t76 = (-mrSges(6,2) * t115 - mrSges(6,3) * t118) * t166;
t71 = t149 - t159;
t69 = t196 * t171;
t68 = qJD(2) * t135 + t106;
t33 = t100 + (t124 - t159) * qJD(2);
t29 = t107 + t35;
t26 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t19 = qJD(4) * t46 - t115 * t125;
t17 = Ifges(7,1) * t57 + Ifges(7,4) * t58 + Ifges(7,5) * t153;
t16 = Ifges(7,4) * t57 + Ifges(7,2) * t58 + Ifges(7,6) * t153;
t15 = t114 * t29 + t117 * t68;
t14 = -t114 * t68 + t117 * t29;
t9 = (qJD(5) - t143) * qJD(4) + t195;
t4 = -qJD(6) * t22 - t114 * t63 + t117 * t19;
t3 = qJD(6) * t21 + t114 * t19 + t117 * t63;
t7 = [t3 * t59 + t4 * t60 + t22 * t43 + t46 * t26 + t21 * t42 + m(4) * (-t125 * t51 - t56 * t66 + t199) + m(5) * (t12 * t46 + t199 + t200) + m(6) * (-t11 * t46 + t200) + m(7) * (t1 * t22 + t2 * t21 + t3 * t6 + t4 * t5 + t46 * t9) + (m(6) * t33 + t76 + t77) * t228 + t217 * t19 + (-m(4) * t50 + m(5) * t47 - t235) * t63 - t242 * (-qJD(4) * t150 + t118 * t125 + t171 * t66) - t243 * (t153 * t45 - t154 * t46) + (mrSges(4,2) * t228 + (-mrSges(3,1) * t116 - mrSges(3,2) * t119) * t112) * qJD(2) ^ 2; -t55 * mrSges(4,1) + t82 * t26 + t37 * t42 + t38 * t43 - t69 * t44 + t71 * t85 + t79 * t76 + t240 * t60 + t239 * t59 + (qJD(2) * t64 + t56) * mrSges(4,2) + m(6) * (t33 * t79 + t36 * t71) + (t99 / 0.2e1 + t55 * mrSges(5,2) - t33 * mrSges(6,3) + t194 * t64 - t243 * t13 + 0.2e1 * (t214 * t30 + t160) * m(6) + (t214 * t34 + t160) * t215 + (t163 * t173 + (-t193 + t224) * t103 + t161 * qJD(4) + t236) * qJD(4) + t223) * t115 + (t12 * mrSges(5,3) - t11 * mrSges(6,1) - t57 * t139 / 0.2e1 - t58 * t137 / 0.2e1 + t9 * t142 - t55 * mrSges(5,1) + t33 * mrSges(6,2) + t16 * t205 + t17 * t206 + (-t1 * t117 + t2 * t114) * mrSges(7,3) + ((t187 - t188) * t207 + t138 * t213 + t140 * t212 - t25 * t141 + t41 * t205 + t114 * t40 / 0.2e1 + (t114 * t6 + t117 * t5) * mrSges(7,3)) * qJD(6) + (((-t189 / 0.2e1 - t186 / 0.2e1 - t163) * t118 + (-0.3e1 / 0.2e1 * Ifges(6,3) - 0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,2)) * t115) * qJD(2) + t123 + t162 * qJD(4)) * qJD(4) + (qJD(4) * t217 + t219) * t103 - t242 * t64) * t118 + (t77 + t55 * t215 / 0.2e1) * (-pkin(3) - t203) + (-t47 * t215 / 0.2e1 + t235) * t62 + (t1 * t38 + t2 * t37 + t239 * t6 + t240 * t5 - t25 * t69 + t82 * t9) * m(7) + (t50 * t62 - t51 * t64 + (-t111 * t56 - t113 * t55) * pkin(2)) * m(4); m(7) * (-t6 * t155 + t5 * t156) + (m(7) * (t176 * t5 + t177 * t6) - t243 * qJD(2) * (-t115 ^ 2 - t118 ^ 2)) * qJD(4) + (-m(6) - m(5)) * t118 * t13 + (t26 + m(7) * t9 + (t114 * t59 + t117 * t60 + t217) * qJD(4) + t219) * t115 + (t242 * qJD(4) + t220 - t238) * t118; t57 * t140 / 0.2e1 + t58 * t138 / 0.2e1 + t9 * t141 + t117 * t17 / 0.2e1 + t16 * t206 - t87 * t85 - t15 * t59 - t14 * t60 - t28 * t44 + qJ(5) * t26 - t11 * mrSges(6,3) - t12 * mrSges(5,2) + t194 * t35 + t193 * t34 + t244 * t13 + t192 * qJD(5) - t225 * mrSges(7,3) + (t136 * t208 + t137 * t213 + t139 * t212 - t241) * qJD(6) + (((t188 / 0.2e1 - t187 / 0.2e1 - pkin(4) * mrSges(6,1) + t162) * qJD(4) - t105 / 0.2e1 - t123 - Ifges(6,6) * t172 / 0.2e1) * t118 + ((-mrSges(6,1) * qJ(5) + t161) * qJD(4) + (Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t172 + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t173 - t236) * t115) * qJD(2) + (qJ(5) * t9 - t14 * t5 - t15 * t6 + t237 * t25) * m(7) + (t134 + t238 + (-m(7) * t144 - t133) * qJD(6)) * t120 + (-pkin(4) * t13 - qJ(5) * t11 + t227 * t31 - t30 * t35 - t36 * t87) * m(6); -t192 * qJD(4) + (mrSges(6,1) * t169 + (-t133 + t85) * t115) * qJD(2) + (-qJD(4) * t25 - t101 * t144 + t225) * m(7) + (qJD(4) * t31 + t173 * t36 + t13) * m(6) - t220; t99 - t25 * (mrSges(7,1) * t84 + mrSges(7,2) * t83) + (Ifges(7,1) * t83 - t204) * t212 + t40 * t211 + (Ifges(7,5) * t83 - Ifges(7,6) * t84) * t208 - t5 * t59 + t6 * t60 + (t5 * t83 + t6 * t84) * mrSges(7,3) + (-Ifges(7,2) * t84 + t41 + t80) * t213 + t223;];
tauc  = t7(:);
