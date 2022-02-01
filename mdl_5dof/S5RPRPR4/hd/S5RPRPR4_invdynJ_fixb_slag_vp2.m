% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:37
% DurationCPUTime: 6.46s
% Computational Cost: add. (3408->389), mult. (7405->535), div. (0->0), fcn. (5001->16), ass. (0->173)
t244 = qJD(3) / 0.2e1;
t150 = qJDD(3) + qJDD(5);
t152 = qJD(3) + qJD(5);
t155 = sin(pkin(9));
t157 = cos(pkin(9));
t161 = sin(qJ(3));
t164 = cos(qJ(3));
t112 = -t155 * t161 + t157 * t164;
t104 = t112 * qJD(1);
t187 = t164 * qJD(1);
t188 = t161 * qJD(1);
t105 = -t155 * t187 - t157 * t188;
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t178 = t163 * t104 + t105 * t160;
t186 = qJD(1) * qJD(3);
t116 = qJDD(1) * t164 - t161 * t186;
t117 = qJDD(1) * t161 + t164 * t186;
t68 = t116 * t157 - t117 * t155;
t69 = t116 * t155 + t117 * t157;
t19 = qJD(5) * t178 + t160 * t68 + t163 * t69;
t156 = sin(pkin(8));
t133 = pkin(1) * t156 + pkin(6);
t122 = t133 * qJD(1);
t145 = t164 * qJDD(2);
t185 = qJD(1) * qJD(4);
t189 = qJD(3) * t164;
t120 = t133 * qJDD(1);
t242 = qJD(2) * qJD(3) + t120;
t37 = -t122 * t189 + qJDD(3) * pkin(3) - qJ(4) * t117 + t145 + (-t185 - t242) * t161;
t190 = qJD(3) * t161;
t57 = t161 * qJDD(2) - t122 * t190 + t164 * t242;
t40 = qJ(4) * t116 + t164 * t185 + t57;
t12 = -t155 * t40 + t157 * t37;
t6 = qJDD(3) * pkin(4) - pkin(7) * t69 + t12;
t211 = pkin(7) * t105;
t176 = qJ(4) * qJD(1) + t122;
t191 = qJD(2) * t161;
t83 = t164 * t176 + t191;
t72 = t155 * t83;
t147 = t164 * qJD(2);
t82 = -t161 * t176 + t147;
t76 = qJD(3) * pkin(3) + t82;
t38 = t157 * t76 - t72;
t26 = qJD(3) * pkin(4) + t211 + t38;
t212 = pkin(7) * t104;
t201 = t157 * t83;
t39 = t155 * t76 + t201;
t27 = t39 + t212;
t7 = -t160 * t27 + t163 * t26;
t13 = t155 * t37 + t157 * t40;
t9 = pkin(7) * t68 + t13;
t2 = qJD(5) * t7 + t160 * t6 + t163 * t9;
t56 = t104 * t160 - t105 * t163;
t20 = -qJD(5) * t56 - t160 * t69 + t163 * t68;
t215 = Ifges(6,4) * t56;
t50 = Ifges(6,4) * t178;
t24 = Ifges(6,1) * t56 + Ifges(6,5) * t152 + t50;
t8 = t160 * t26 + t163 * t27;
t3 = -qJD(5) * t8 - t160 * t9 + t163 * t6;
t158 = cos(pkin(8));
t135 = -pkin(1) * t158 - pkin(2);
t148 = t164 * pkin(3);
t119 = t135 - t148;
t103 = qJD(1) * t119 + qJD(4);
t64 = -pkin(4) * t104 + t103;
t243 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t150 - (Ifges(6,5) * t178 - Ifges(6,6) * t56) * t152 / 0.2e1 + (t178 * t7 + t56 * t8) * mrSges(6,3) - (-Ifges(6,2) * t56 + t24 + t50) * t178 / 0.2e1 - t64 * (mrSges(6,1) * t56 + mrSges(6,2) * t178) - (Ifges(6,1) * t178 - t215) * t56 / 0.2e1;
t126 = -t164 * mrSges(4,1) + t161 * mrSges(4,2);
t153 = qJ(3) + pkin(9);
t140 = sin(t153);
t142 = cos(t153);
t146 = qJ(5) + t153;
t131 = sin(t146);
t132 = cos(t146);
t179 = t132 * mrSges(6,1) - t131 * mrSges(6,2);
t241 = -mrSges(5,1) * t142 + mrSges(5,2) * t140 + t126 - t179;
t240 = Ifges(4,4) * t187 / 0.2e1 + Ifges(4,5) * t244;
t225 = m(3) + m(4) + m(6) + m(5);
t239 = pkin(1) * t225 + mrSges(2,1);
t23 = Ifges(6,2) * t178 + Ifges(6,6) * t152 + t215;
t237 = t23 / 0.2e1;
t134 = pkin(3) * t157 + pkin(4);
t213 = pkin(3) * t155;
t101 = t134 * t163 - t160 * t213;
t42 = -t155 * t82 - t201;
t30 = t42 - t212;
t43 = t157 * t82 - t72;
t31 = t43 + t211;
t235 = t101 * qJD(5) - t160 * t30 - t163 * t31;
t102 = t134 * t160 + t163 * t213;
t234 = -t102 * qJD(5) + t160 * t31 - t163 * t30;
t198 = qJDD(3) / 0.2e1;
t94 = t122 * t164 + t191;
t58 = -t94 * qJD(3) - t120 * t161 + t145;
t230 = -t161 * t58 + t164 * t57;
t154 = qJ(1) + pkin(8);
t141 = sin(t154);
t143 = cos(t154);
t229 = g(1) * t143 + g(2) * t141;
t205 = Ifges(4,4) * t161;
t172 = t164 * Ifges(4,2) + t205;
t196 = Ifges(4,6) * qJD(3);
t228 = -qJD(1) * t172 / 0.2e1 - t196 / 0.2e1 - t94 * mrSges(4,3);
t93 = -t122 * t161 + t147;
t227 = Ifges(4,1) * t188 / 0.2e1 - t93 * mrSges(4,3) + t240;
t226 = 0.2e1 * t198;
t192 = pkin(4) * t142 + t148;
t224 = -m(4) * pkin(2) - m(6) * (pkin(2) + t192) - m(5) * (t148 + pkin(2)) - mrSges(3,1) + t241;
t159 = -qJ(4) - pkin(6);
t223 = -m(4) * pkin(6) + m(5) * t159 - m(6) * (pkin(7) - t159) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t222 = m(5) * pkin(3);
t219 = t56 / 0.2e1;
t217 = -t105 / 0.2e1;
t193 = qJ(4) + t133;
t177 = qJD(3) * t193;
t86 = qJD(4) * t164 - t161 * t177;
t87 = -qJD(4) * t161 - t164 * t177;
t45 = t155 * t87 + t157 * t86;
t110 = t193 * t161;
t111 = t193 * t164;
t60 = -t155 * t110 + t157 * t111;
t206 = mrSges(4,2) * t164;
t204 = Ifges(4,4) * t164;
t203 = Ifges(5,4) * t105;
t183 = pkin(3) * t190;
t182 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t181 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t44 = -t155 * t86 + t157 * t87;
t59 = -t157 * t110 - t111 * t155;
t175 = -g(1) * t141 + g(2) * t143;
t121 = t135 * qJDD(1);
t173 = mrSges(6,1) * t131 + mrSges(6,2) * t132;
t113 = t155 * t164 + t157 * t161;
t46 = -pkin(7) * t113 + t59;
t47 = pkin(7) * t112 + t60;
t21 = -t160 * t47 + t163 * t46;
t22 = t160 * t46 + t163 * t47;
t65 = t112 * t163 - t113 * t160;
t66 = t112 * t160 + t113 * t163;
t81 = -pkin(3) * t116 + qJDD(4) + t121;
t165 = cos(qJ(1));
t162 = sin(qJ(1));
t125 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t187;
t124 = t135 * qJD(1);
t123 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t188;
t107 = t112 * qJD(3);
t106 = t113 * qJD(3);
t100 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t117;
t99 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t116;
t98 = Ifges(5,4) * t104;
t89 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t105;
t88 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t104;
t85 = pkin(4) * t106 + t183;
t84 = pkin(3) * t188 - pkin(4) * t105;
t80 = -pkin(4) * t112 + t119;
t63 = -mrSges(5,1) * t104 - mrSges(5,2) * t105;
t62 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t69;
t61 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t68;
t52 = -t105 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t98;
t51 = t104 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t203;
t49 = mrSges(6,1) * t152 - mrSges(6,3) * t56;
t48 = -mrSges(6,2) * t152 + mrSges(6,3) * t178;
t41 = -pkin(4) * t68 + t81;
t33 = -pkin(7) * t106 + t45;
t32 = -pkin(7) * t107 + t44;
t29 = -qJD(5) * t66 - t106 * t163 - t107 * t160;
t28 = qJD(5) * t65 - t106 * t160 + t107 * t163;
t25 = -mrSges(6,1) * t178 + mrSges(6,2) * t56;
t15 = -mrSges(6,2) * t150 + mrSges(6,3) * t20;
t14 = mrSges(6,1) * t150 - mrSges(6,3) * t19;
t5 = -qJD(5) * t22 - t160 * t33 + t163 * t32;
t4 = qJD(5) * t21 + t160 * t32 + t163 * t33;
t1 = [(Ifges(6,1) * t28 + Ifges(6,4) * t29) * t219 + t164 * (Ifges(4,4) * t117 + Ifges(4,2) * t116 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t152 * (Ifges(6,5) * t28 + Ifges(6,6) * t29) / 0.2e1 + t135 * (-mrSges(4,1) * t116 + mrSges(4,2) * t117) + (Ifges(4,5) * t161 + Ifges(4,6) * t164) * t198 + (m(4) * t135 + t126) * t121 + t107 * t52 / 0.2e1 - t106 * t51 / 0.2e1 + t85 * t25 + t45 * t88 + t44 * t89 + t64 * (-mrSges(6,1) * t29 + mrSges(6,2) * t28) + t60 * t61 + t59 * t62 + t4 * t48 + t5 * t49 + t21 * t14 + t22 * t15 + t28 * t24 / 0.2e1 + (mrSges(6,2) * t41 - mrSges(6,3) * t3 + Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t150) * t66 + t161 * (Ifges(4,4) * t116 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t117 * t204 / 0.2e1 + t116 * t172 / 0.2e1 + t63 * t183 + (t164 * t99 - t161 * t100 + m(4) * ((-t94 * t161 - t93 * t164) * qJD(3) + t230) - t123 * t189 - t125 * t190) * t133 + t230 * mrSges(4,3) + (t81 * mrSges(5,2) - t12 * mrSges(5,3) + t69 * Ifges(5,1) + Ifges(5,4) * t68 + Ifges(5,5) * t226) * t113 + (-t81 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t69 + t68 * Ifges(5,2) + Ifges(5,6) * t226) * t112 + t227 * t189 + t228 * t190 + m(5) * (t103 * t183 + t119 * t81 + t12 * t59 + t13 * t60 + t38 * t44 + t39 * t45) + t80 * t181 + t119 * t182 + (t164 * (-Ifges(4,2) * t161 + t204) + t161 * (Ifges(4,1) * t164 - t205)) * t186 / 0.2e1 + qJD(3) ^ 2 * (Ifges(4,5) * t164 - Ifges(4,6) * t161) / 0.2e1 + (Ifges(5,5) * t107 - Ifges(5,6) * t106) * t244 + t178 * (Ifges(6,4) * t28 + Ifges(6,2) * t29) / 0.2e1 + (Ifges(5,1) * t107 - Ifges(5,4) * t106) * t217 + t104 * (Ifges(5,4) * t107 - Ifges(5,2) * t106) / 0.2e1 + t103 * (mrSges(5,1) * t106 + mrSges(5,2) * t107) + (-t106 * t39 - t107 * t38) * mrSges(5,3) + (-t28 * t7 + t29 * t8) * mrSges(6,3) + (mrSges(2,2) * t165 - t141 * t224 + t143 * t223 + t162 * t239) * g(1) + (mrSges(2,2) * t162 + t141 * t223 + t143 * t224 - t165 * t239) * g(2) + t124 * (mrSges(4,1) * t161 + t206) * qJD(3) + m(6) * (t2 * t22 + t21 * t3 + t4 * t8 + t41 * t80 + t5 * t7 + t64 * t85) + t117 * t161 * Ifges(4,1) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t158 - 0.2e1 * mrSges(3,2) * t156 + m(3) * (t156 ^ 2 + t158 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t29 * t237 + (-mrSges(6,1) * t41 + mrSges(6,3) * t2 + Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t150) * t65; m(3) * qJDD(2) + t164 * t100 - t106 * t89 + t107 * t88 + t112 * t62 + t113 * t61 + t65 * t14 + t66 * t15 + t161 * t99 + t28 * t48 + t29 * t49 + (-t123 * t161 + t125 * t164) * qJD(3) + m(4) * (t161 * t57 + t164 * t58 + (-t161 * t93 + t164 * t94) * qJD(3)) + m(5) * (-t106 * t38 + t107 * t39 + t112 * t12 + t113 * t13) + m(6) * (t2 * t66 + t28 * t8 + t29 * t7 + t3 * t65) - t225 * g(3); t51 * t217 + (t12 * t157 + t13 * t155) * t222 + t234 * t49 + (-t192 * g(3) + t101 * t3 + t102 * t2 + t234 * t7 + t235 * t8 - t64 * t84) * m(6) + t235 * t48 + Ifges(4,6) * t116 + Ifges(4,5) * t117 + t94 * t123 - t93 * t125 - qJD(3) * (Ifges(5,5) * t104 + Ifges(5,6) * t105) / 0.2e1 - t103 * (-mrSges(5,1) * t105 + mrSges(5,2) * t104) + t101 * t14 + t102 * t15 - t43 * t88 - t42 * t89 + Ifges(5,5) * t69 - t84 * t25 + Ifges(5,6) * t68 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t12 * mrSges(5,1) - t13 * mrSges(5,2) + t229 * (-m(6) * (-pkin(3) * t161 - pkin(4) * t140) + mrSges(5,1) * t140 + t206 + mrSges(5,2) * t142 + (mrSges(4,1) + t222) * t161 + t173) - m(5) * (t38 * t42 + t39 * t43) - (Ifges(5,2) * t105 + t52 + t98) * t104 / 0.2e1 + (-m(5) * t148 + t241) * g(3) + t243 + (t155 * t61 + t157 * t62) * pkin(3) + (t104 * t38 - t105 * t39) * mrSges(5,3) + ((-t124 * mrSges(4,2) - t227 - t240) * t164 + (-t124 * mrSges(4,1) + t196 / 0.2e1 + (t205 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t164) * qJD(1) + (-m(5) * t103 - t63) * pkin(3) - t228) * t161) * qJD(1) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t56 * t237 + t105 * (Ifges(5,1) * t104 + t203) / 0.2e1; -t104 * t88 - t105 * t89 - t178 * t48 + t56 * t49 + t181 + t182 + (-t178 * t8 + t56 * t7 + t175 + t41) * m(6) + (-t104 * t39 - t105 * t38 + t175 + t81) * m(5); -g(3) * t179 + t229 * t173 + t23 * t219 - t7 * t48 + t8 * t49 + t243;];
tau = t1;
