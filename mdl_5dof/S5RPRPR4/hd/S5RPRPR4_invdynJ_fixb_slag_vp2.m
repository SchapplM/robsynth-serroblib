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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:38:16
% EndTime: 2020-01-03 11:38:40
% DurationCPUTime: 7.09s
% Computational Cost: add. (3408->390), mult. (7405->537), div. (0->0), fcn. (5001->16), ass. (0->175)
t251 = -qJD(3) / 0.2e1;
t153 = qJDD(3) + qJDD(5);
t155 = qJD(3) + qJD(5);
t158 = sin(pkin(9));
t160 = cos(pkin(9));
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t112 = -t158 * t164 + t160 * t167;
t104 = t112 * qJD(1);
t191 = t167 * qJD(1);
t192 = t164 * qJD(1);
t105 = -t158 * t191 - t160 * t192;
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t181 = t166 * t104 + t105 * t163;
t190 = qJD(1) * qJD(3);
t118 = qJDD(1) * t167 - t164 * t190;
t119 = qJDD(1) * t164 + t167 * t190;
t68 = t118 * t160 - t119 * t158;
t69 = t118 * t158 + t119 * t160;
t19 = qJD(5) * t181 + t163 * t68 + t166 * t69;
t159 = sin(pkin(8));
t135 = pkin(1) * t159 + pkin(6);
t124 = t135 * qJD(1);
t147 = t167 * qJDD(2);
t189 = qJD(1) * qJD(4);
t193 = qJD(3) * t167;
t122 = t135 * qJDD(1);
t247 = qJD(2) * qJD(3) + t122;
t37 = -t124 * t193 + qJDD(3) * pkin(3) - qJ(4) * t119 + t147 + (-t189 - t247) * t164;
t194 = qJD(3) * t164;
t57 = t164 * qJDD(2) - t124 * t194 + t167 * t247;
t40 = qJ(4) * t118 + t167 * t189 + t57;
t12 = -t158 * t40 + t160 * t37;
t6 = qJDD(3) * pkin(4) - pkin(7) * t69 + t12;
t216 = pkin(7) * t105;
t179 = qJ(4) * qJD(1) + t124;
t195 = qJD(2) * t164;
t83 = t167 * t179 + t195;
t72 = t158 * t83;
t149 = t167 * qJD(2);
t82 = -t164 * t179 + t149;
t76 = qJD(3) * pkin(3) + t82;
t38 = t160 * t76 - t72;
t26 = qJD(3) * pkin(4) + t216 + t38;
t217 = pkin(7) * t104;
t205 = t160 * t83;
t39 = t158 * t76 + t205;
t27 = t39 + t217;
t7 = -t163 * t27 + t166 * t26;
t13 = t158 * t37 + t160 * t40;
t9 = pkin(7) * t68 + t13;
t2 = qJD(5) * t7 + t163 * t6 + t166 * t9;
t56 = t104 * t163 - t105 * t166;
t20 = -qJD(5) * t56 - t163 * t69 + t166 * t68;
t219 = Ifges(6,4) * t56;
t50 = Ifges(6,4) * t181;
t24 = Ifges(6,1) * t56 + Ifges(6,5) * t155 + t50;
t156 = qJ(3) + pkin(9);
t148 = qJ(5) + t156;
t133 = sin(t148);
t134 = cos(t148);
t246 = mrSges(6,1) * t133 + mrSges(6,2) * t134;
t157 = qJ(1) + pkin(8);
t145 = cos(t157);
t249 = t145 * g(3);
t8 = t163 * t26 + t166 * t27;
t3 = -qJD(5) * t8 - t163 * t9 + t166 * t6;
t161 = cos(pkin(8));
t137 = -pkin(1) * t161 - pkin(2);
t151 = t167 * pkin(3);
t121 = t137 - t151;
t103 = qJD(1) * t121 + qJD(4);
t64 = -pkin(4) * t104 + t103;
t250 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t153 - t246 * t249 - (Ifges(6,5) * t181 - Ifges(6,6) * t56) * t155 / 0.2e1 + (t181 * t7 + t56 * t8) * mrSges(6,3) - (-Ifges(6,2) * t56 + t24 + t50) * t181 / 0.2e1 - t64 * (mrSges(6,1) * t56 + mrSges(6,2) * t181) - (Ifges(6,1) * t181 - t219) * t56 / 0.2e1;
t231 = -m(3) - m(4) - m(5) - m(6);
t248 = pkin(1) * t231 - mrSges(2,1);
t128 = -t167 * mrSges(4,1) + t164 * mrSges(4,2);
t142 = sin(t156);
t144 = cos(t156);
t182 = t134 * mrSges(6,1) - t133 * mrSges(6,2);
t245 = -mrSges(5,1) * t144 + mrSges(5,2) * t142 + t128 - t182;
t244 = -Ifges(4,4) * t191 / 0.2e1 + Ifges(4,5) * t251;
t23 = Ifges(6,2) * t181 + Ifges(6,6) * t155 + t219;
t242 = t23 / 0.2e1;
t136 = pkin(3) * t160 + pkin(4);
t218 = pkin(3) * t158;
t101 = t136 * t166 - t163 * t218;
t42 = -t158 * t82 - t205;
t30 = t42 - t217;
t43 = t160 * t82 - t72;
t31 = t43 + t216;
t240 = t101 * qJD(5) - t163 * t30 - t166 * t31;
t102 = t136 * t163 + t166 * t218;
t239 = -t102 * qJD(5) + t163 * t31 - t166 * t30;
t202 = qJDD(3) / 0.2e1;
t94 = t124 * t167 + t195;
t58 = -qJD(3) * t94 - t122 * t164 + t147;
t235 = -t164 * t58 + t167 * t57;
t209 = Ifges(4,4) * t164;
t175 = t167 * Ifges(4,2) + t209;
t200 = Ifges(4,6) * qJD(3);
t234 = -t94 * mrSges(4,3) - qJD(1) * t175 / 0.2e1 - t200 / 0.2e1;
t93 = -t124 * t164 + t149;
t233 = -t93 * mrSges(4,3) + Ifges(4,1) * t192 / 0.2e1 - t244;
t232 = 0.2e1 * t202;
t211 = mrSges(4,2) * t167;
t226 = m(5) * pkin(3);
t230 = -t164 * (mrSges(4,1) + t226) + m(6) * (-pkin(3) * t164 - pkin(4) * t142) - mrSges(5,1) * t142 - mrSges(5,2) * t144 - t211;
t196 = pkin(4) * t144 + t151;
t229 = -m(6) * (pkin(2) + t196) - m(5) * (t151 + pkin(2)) - m(4) * pkin(2) - mrSges(3,1) + t245;
t162 = -qJ(4) - pkin(6);
t227 = m(4) * pkin(6) - m(5) * t162 - m(6) * (-pkin(7) + t162) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t223 = t56 / 0.2e1;
t221 = -t105 / 0.2e1;
t143 = sin(t157);
t215 = g(2) * t143;
t197 = qJ(4) + t135;
t180 = qJD(3) * t197;
t86 = qJD(4) * t167 - t164 * t180;
t87 = -qJD(4) * t164 - t167 * t180;
t45 = t158 * t87 + t160 * t86;
t110 = t197 * t164;
t111 = t197 * t167;
t60 = -t158 * t110 + t160 * t111;
t208 = Ifges(4,4) * t167;
t207 = Ifges(5,4) * t105;
t187 = pkin(3) * t194;
t185 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t184 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t44 = -t158 * t86 + t160 * t87;
t59 = -t160 * t110 - t111 * t158;
t178 = g(2) * t145 + g(3) * t143;
t123 = t137 * qJDD(1);
t113 = t158 * t167 + t160 * t164;
t46 = -pkin(7) * t113 + t59;
t47 = pkin(7) * t112 + t60;
t21 = -t163 * t47 + t166 * t46;
t22 = t163 * t46 + t166 * t47;
t65 = t112 * t166 - t113 * t163;
t66 = t112 * t163 + t113 * t166;
t81 = -pkin(3) * t118 + qJDD(4) + t123;
t168 = cos(qJ(1));
t165 = sin(qJ(1));
t127 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t191;
t126 = t137 * qJD(1);
t125 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t192;
t107 = t112 * qJD(3);
t106 = t113 * qJD(3);
t100 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t99 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t118;
t98 = Ifges(5,4) * t104;
t89 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t105;
t88 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t104;
t85 = pkin(4) * t106 + t187;
t84 = pkin(3) * t192 - pkin(4) * t105;
t80 = -pkin(4) * t112 + t121;
t63 = -mrSges(5,1) * t104 - mrSges(5,2) * t105;
t62 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t69;
t61 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t68;
t52 = -t105 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t98;
t51 = t104 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t207;
t49 = mrSges(6,1) * t155 - mrSges(6,3) * t56;
t48 = -mrSges(6,2) * t155 + mrSges(6,3) * t181;
t41 = -pkin(4) * t68 + t81;
t33 = -pkin(7) * t106 + t45;
t32 = -pkin(7) * t107 + t44;
t29 = -qJD(5) * t66 - t106 * t166 - t107 * t163;
t28 = qJD(5) * t65 - t106 * t163 + t107 * t166;
t25 = -mrSges(6,1) * t181 + mrSges(6,2) * t56;
t15 = -mrSges(6,2) * t153 + mrSges(6,3) * t20;
t14 = mrSges(6,1) * t153 - mrSges(6,3) * t19;
t5 = -qJD(5) * t22 - t163 * t33 + t166 * t32;
t4 = qJD(5) * t21 + t163 * t32 + t166 * t33;
t1 = [t80 * t184 + t121 * t185 + t126 * (mrSges(4,1) * t164 + t211) * qJD(3) + t167 * (Ifges(4,4) * t119 + Ifges(4,2) * t118 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t155 * (Ifges(6,5) * t28 + Ifges(6,6) * t29) / 0.2e1 + t137 * (-mrSges(4,1) * t118 + mrSges(4,2) * t119) - t106 * t51 / 0.2e1 + t107 * t52 / 0.2e1 + t44 * t89 + t85 * t25 + t45 * t88 + t60 * t61 + t59 * t62 + t64 * (-mrSges(6,1) * t29 + mrSges(6,2) * t28) + t4 * t48 + t5 * t49 + (t167 * (-Ifges(4,2) * t164 + t208) + t164 * (Ifges(4,1) * t167 - t209)) * t190 / 0.2e1 + t181 * (Ifges(6,4) * t28 + Ifges(6,2) * t29) / 0.2e1 + t103 * (mrSges(5,1) * t106 + mrSges(5,2) * t107) + t104 * (Ifges(5,4) * t107 - Ifges(5,2) * t106) / 0.2e1 + qJD(3) * (Ifges(5,5) * t107 - Ifges(5,6) * t106) / 0.2e1 + (Ifges(5,1) * t107 - Ifges(5,4) * t106) * t221 + (-t106 * t39 - t107 * t38) * mrSges(5,3) + (-mrSges(6,1) * t41 + mrSges(6,3) * t2 + Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t153) * t65 + t119 * t208 / 0.2e1 + t164 * (Ifges(4,4) * t118 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t119 * t164 * Ifges(4,1) + t28 * t24 / 0.2e1 + t21 * t14 + t22 * t15 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t161 - 0.2e1 * t159 * mrSges(3,2) + m(3) * (t159 ^ 2 + t161 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + m(5) * (t103 * t187 + t12 * t59 + t121 * t81 + t13 * t60 + t38 * t44 + t39 * t45) + (m(4) * t137 + t128) * t123 + t118 * t175 / 0.2e1 + (m(4) * ((-t94 * t164 - t93 * t167) * qJD(3) + t235) + t167 * t99 - t127 * t194 - t164 * t100 - t125 * t193) * t135 + t235 * mrSges(4,3) + (Ifges(6,1) * t28 + Ifges(6,4) * t29) * t223 + t29 * t242 + (Ifges(4,5) * t164 + Ifges(4,6) * t167) * t202 + (mrSges(6,2) * t41 - mrSges(6,3) * t3 + Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t153) * t66 + t63 * t187 + m(6) * (t2 * t22 + t21 * t3 + t4 * t8 + t41 * t80 + t5 * t7 + t64 * t85) + (-t28 * t7 + t29 * t8) * mrSges(6,3) + (-mrSges(2,2) * t168 + t143 * t229 + t145 * t227 + t165 * t248) * g(3) + (mrSges(2,2) * t165 - t143 * t227 + t145 * t229 + t168 * t248) * g(2) + (t81 * mrSges(5,2) - t12 * mrSges(5,3) + t69 * Ifges(5,1) + Ifges(5,4) * t68 + Ifges(5,5) * t232) * t113 + (-t81 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t69 + t68 * Ifges(5,2) + Ifges(5,6) * t232) * t112 + t233 * t193 + t234 * t194 + qJD(3) ^ 2 * (Ifges(4,5) * t167 - Ifges(4,6) * t164) / 0.2e1; m(3) * qJDD(2) + t167 * t100 - t106 * t89 + t107 * t88 + t112 * t62 + t113 * t61 + t65 * t14 + t66 * t15 + t164 * t99 + t28 * t48 + t29 * t49 + (-t125 * t164 + t127 * t167) * qJD(3) + m(4) * (t164 * t57 + t167 * t58 + (-t164 * t93 + t167 * t94) * qJD(3)) + m(5) * (-t106 * t38 + t107 * t39 + t112 * t12 + t113 * t13) + m(6) * (t2 * t66 + t28 * t8 + t29 * t7 + t3 * t65) + t231 * g(1); (t246 - t230) * t215 + t105 * (Ifges(5,1) * t104 + t207) / 0.2e1 + (t104 * t38 - t105 * t39) * mrSges(5,3) + Ifges(4,5) * t119 + t94 * t125 - t93 * t127 + Ifges(4,6) * t118 + t101 * t14 + t102 * t15 - t103 * (-mrSges(5,1) * t105 + mrSges(5,2) * t104) - t42 * t89 - t84 * t25 - t43 * t88 + Ifges(5,6) * t68 + Ifges(5,5) * t69 - t57 * mrSges(4,2) + t58 * mrSges(4,1) - m(5) * (t38 * t42 + t39 * t43) + t250 - (Ifges(5,2) * t105 + t52 + t98) * t104 / 0.2e1 + t239 * t49 + (-t196 * g(1) + t101 * t3 + t102 * t2 + t239 * t7 + t240 * t8 - t64 * t84) * m(6) + t240 * t48 + t12 * mrSges(5,1) - t13 * mrSges(5,2) + t56 * t242 + t51 * t221 + (t12 * t160 + t13 * t158) * t226 + (Ifges(5,5) * t104 + Ifges(5,6) * t105) * t251 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t230 * t249 + (t158 * t61 + t160 * t62) * pkin(3) + ((-t126 * mrSges(4,2) - t233 + t244) * t167 + (-t126 * mrSges(4,1) + t200 / 0.2e1 + (t209 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t167) * qJD(1) + (-m(5) * t103 - t63) * pkin(3) - t234) * t164) * qJD(1) + (-m(5) * t151 + t245) * g(1); -t104 * t88 - t105 * t89 - t181 * t48 + t56 * t49 + t184 + t185 + (-t181 * t8 + t56 * t7 + t178 + t41) * m(6) + (-t104 * t39 - t105 * t38 + t178 + t81) * m(5); -g(1) * t182 + t215 * t246 + t23 * t223 - t7 * t48 + t8 * t49 + t250;];
tau = t1;
