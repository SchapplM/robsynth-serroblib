% Calculate vector of inverse dynamics joint torques for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:40
% EndTime: 2019-12-05 15:50:04
% DurationCPUTime: 6.62s
% Computational Cost: add. (2515->408), mult. (6034->597), div. (0->0), fcn. (4684->12), ass. (0->199)
t128 = sin(qJ(4));
t131 = cos(qJ(4));
t167 = pkin(4) * t128 - pkin(8) * t131;
t123 = sin(pkin(5));
t121 = sin(pkin(10));
t124 = cos(pkin(10));
t129 = sin(qJ(2));
t132 = cos(qJ(2));
t152 = t121 * t132 + t124 * t129;
t80 = t152 * t123;
t71 = qJD(1) * t80;
t266 = t167 * qJD(4) - t71;
t127 = sin(qJ(5));
t191 = t128 * qJD(2);
t130 = cos(qJ(5));
t196 = qJD(4) * t130;
t95 = -t127 * t191 + t196;
t265 = -t95 / 0.2e1;
t198 = qJD(4) * t127;
t96 = t130 * t191 + t198;
t264 = -t96 / 0.2e1;
t190 = t131 * qJD(2);
t115 = qJD(5) - t190;
t263 = -t115 / 0.2e1;
t184 = mrSges(5,3) * t191;
t216 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t95 + mrSges(6,2) * t96 + t184;
t164 = -mrSges(6,1) * t130 + mrSges(6,2) * t127;
t143 = m(6) * pkin(4) - t164;
t165 = mrSges(5,1) * t131 - mrSges(5,2) * t128;
t185 = m(6) * pkin(8) + mrSges(6,3);
t262 = t128 * t185 + t131 * t143 + t165;
t126 = cos(pkin(5));
t113 = qJD(1) * t126 + qJD(3);
t199 = qJD(1) * t132;
t179 = t123 * t199;
t103 = qJD(2) * pkin(2) + t179;
t209 = t123 * t129;
t180 = qJD(1) * t209;
t65 = t121 * t103 + t124 * t180;
t62 = qJD(2) * pkin(7) + t65;
t39 = t113 * t128 + t131 * t62;
t35 = qJD(4) * pkin(8) + t39;
t150 = -pkin(4) * t131 - pkin(8) * t128 - pkin(3);
t105 = t121 * t180;
t64 = t103 * t124 - t105;
t42 = qJD(2) * t150 - t64;
t10 = t127 * t42 + t130 * t35;
t188 = qJD(2) * qJD(4);
t101 = qJDD(2) * t131 - t128 * t188;
t102 = qJDD(2) * t128 + t131 * t188;
t213 = qJDD(2) * pkin(2);
t170 = qJD(2) * t180;
t207 = t123 * t132;
t82 = qJDD(1) * t207 - t170;
t77 = t82 + t213;
t175 = qJD(2) * t199;
t83 = (qJDD(1) * t129 + t175) * t123;
t36 = -t121 * t83 + t124 * t77;
t32 = -qJDD(2) * pkin(3) - t36;
t13 = -pkin(4) * t101 - pkin(8) * t102 + t32;
t112 = qJDD(1) * t126 + qJDD(3);
t195 = qJD(4) * t131;
t197 = qJD(4) * t128;
t37 = t121 * t77 + t124 * t83;
t33 = qJDD(2) * pkin(7) + t37;
t7 = t128 * t112 + t113 * t195 + t131 * t33 - t197 * t62;
t5 = qJDD(4) * pkin(8) + t7;
t9 = -t127 * t35 + t130 * t42;
t1 = qJD(5) * t9 + t127 * t13 + t130 * t5;
t2 = -qJD(5) * t10 - t127 * t5 + t13 * t130;
t166 = t1 * t130 - t127 * t2;
t192 = qJD(5) * t130;
t194 = qJD(5) * t127;
t261 = -t10 * t194 - t9 * t192 + t166;
t122 = sin(pkin(9));
t125 = cos(pkin(9));
t205 = t126 * t132;
t260 = -t122 * t129 + t125 * t205;
t259 = -m(6) - m(5);
t258 = t101 / 0.2e1;
t257 = t102 / 0.2e1;
t118 = pkin(2) * t121 + pkin(7);
t178 = t118 * t197;
t201 = t130 * t131;
t203 = t127 * t131;
t233 = pkin(2) * t124;
t90 = t150 - t233;
t56 = -t118 * t203 + t130 * t90;
t74 = t124 * t179 - t105;
t256 = qJD(5) * t56 + t266 * t127 - t130 * t178 - t201 * t74;
t57 = t118 * t201 + t127 * t90;
t255 = -qJD(5) * t57 + t127 * t178 + t266 * t130 + t203 * t74;
t54 = qJD(5) * t95 + qJDD(4) * t127 + t102 * t130;
t55 = -qJD(5) * t96 + qJDD(4) * t130 - t102 * t127;
t14 = -mrSges(6,1) * t55 + mrSges(6,2) * t54;
t225 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t102 + t14;
t254 = mrSges(5,1) + t143;
t253 = mrSges(5,2) - t185;
t120 = Ifges(5,4) * t190;
t91 = Ifges(6,4) * t95;
t51 = Ifges(6,1) * t96 + Ifges(6,5) * t115 + t91;
t252 = Ifges(5,1) * t191 + Ifges(5,5) * qJD(4) + t130 * t51 + t120;
t200 = t132 * t124;
t79 = t121 * t209 - t123 * t200;
t251 = (mrSges(3,1) * t132 - mrSges(3,2) * t129) * t123 - t79 * mrSges(4,1) - t80 * mrSges(4,2);
t214 = qJD(4) * t39;
t8 = t112 * t131 - t128 * t33 - t214;
t92 = qJDD(5) - t101;
t29 = mrSges(6,1) * t92 - mrSges(6,3) * t54;
t30 = -mrSges(6,2) * t92 + mrSges(6,3) * t55;
t248 = -t127 * t29 + t130 * t30;
t247 = -t128 * t8 + t131 * t7;
t224 = Ifges(5,4) * t128;
t160 = Ifges(5,2) * t131 + t224;
t246 = Ifges(5,6) * qJD(4) / 0.2e1 + qJD(2) * t160 / 0.2e1 + Ifges(6,5) * t264 + Ifges(6,6) * t265 + Ifges(6,3) * t263;
t245 = -mrSges(4,1) - t262;
t163 = t127 * mrSges(6,1) + t130 * mrSges(6,2);
t142 = m(6) * pkin(7) + t163;
t244 = m(5) * pkin(7) - mrSges(4,2) + t142;
t243 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t242 = qJD(2) ^ 2;
t234 = Ifges(6,4) * t96;
t50 = Ifges(6,2) * t95 + Ifges(6,6) * t115 + t234;
t241 = -t50 / 0.2e1;
t240 = t54 / 0.2e1;
t239 = t55 / 0.2e1;
t238 = t92 / 0.2e1;
t236 = t96 / 0.2e1;
t6 = -qJDD(4) * pkin(4) - t8;
t230 = t128 * t6;
t223 = Ifges(5,4) * t131;
t222 = Ifges(6,4) * t127;
t221 = Ifges(6,4) * t130;
t219 = t128 * t74;
t114 = pkin(2) * t207;
t215 = -t79 * pkin(3) + t114;
t210 = t123 * t128;
t208 = t123 * t131;
t206 = t126 * t129;
t204 = t127 * t128;
t202 = t128 * t130;
t193 = qJD(5) * t128;
t189 = m(4) - t259;
t187 = Ifges(6,5) * t54 + Ifges(6,6) * t55 + Ifges(6,3) * t92;
t183 = mrSges(5,3) * t190;
t176 = t127 * t195;
t171 = t188 / 0.2e1;
t169 = t260 * pkin(2);
t162 = Ifges(6,1) * t130 - t222;
t161 = Ifges(6,1) * t127 + t221;
t159 = -Ifges(6,2) * t127 + t221;
t158 = Ifges(6,2) * t130 + t222;
t157 = Ifges(5,5) * t131 - Ifges(5,6) * t128;
t156 = Ifges(6,5) * t130 - Ifges(6,6) * t127;
t155 = Ifges(6,5) * t127 + Ifges(6,6) * t130;
t81 = t152 * t126;
t93 = t121 * t129 - t200;
t43 = t122 * t93 - t125 * t81;
t46 = t122 * t81 + t125 * t93;
t60 = t126 * t128 + t131 * t80;
t22 = t127 * t79 + t130 * t60;
t21 = -t127 * t60 + t130 * t79;
t38 = t113 * t131 - t128 * t62;
t153 = t126 * t131 - t128 * t80;
t149 = -t122 * t205 - t125 * t129;
t61 = -qJD(2) * pkin(3) - t64;
t147 = t61 * (mrSges(5,1) * t128 + mrSges(5,2) * t131);
t146 = t128 * (Ifges(5,1) * t131 - t224);
t144 = t93 * t126;
t141 = t149 * pkin(2);
t140 = -t127 * t193 + t130 * t195;
t139 = t128 * t192 + t176;
t137 = Ifges(6,5) * t128 + t131 * t162;
t136 = Ifges(6,6) * t128 + t131 * t159;
t135 = Ifges(6,3) * t128 + t131 * t156;
t119 = -pkin(3) - t233;
t108 = -qJD(4) * mrSges(5,2) + t183;
t98 = t167 * qJD(2);
t97 = t165 * qJD(2);
t84 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t101;
t73 = t93 * t123 * qJD(2);
t72 = qJD(2) * t80;
t68 = mrSges(6,1) * t115 - mrSges(6,3) * t96;
t67 = -mrSges(6,2) * t115 + mrSges(6,3) * t95;
t63 = -mrSges(5,1) * t101 + mrSges(5,2) * t102;
t47 = t122 * t144 - t125 * t152;
t44 = -t122 * t152 - t125 * t144;
t34 = -qJD(4) * pkin(4) - t38;
t28 = t122 * t210 - t131 * t46;
t26 = -t125 * t210 - t131 * t43;
t20 = qJD(4) * t153 - t131 * t73;
t19 = qJD(4) * t60 - t128 * t73;
t16 = t127 * t98 + t130 * t38;
t15 = -t127 * t38 + t130 * t98;
t12 = t54 * Ifges(6,1) + t55 * Ifges(6,4) + t92 * Ifges(6,5);
t11 = t54 * Ifges(6,4) + t55 * Ifges(6,2) + t92 * Ifges(6,6);
t4 = qJD(5) * t21 + t127 * t72 + t130 * t20;
t3 = -qJD(5) * t22 - t127 * t20 + t130 * t72;
t17 = [m(2) * qJDD(1) + t20 * t108 + t21 * t29 + t22 * t30 + t3 * t68 + t4 * t67 + t60 * t84 + t79 * t63 - t72 * t97 - t225 * t153 + t216 * t19 + (-mrSges(3,1) * t129 - mrSges(3,2) * t132) * t242 * t123 + (-mrSges(4,1) * t72 + mrSges(4,2) * t73) * qJD(2) + t251 * qJDD(2) + (-m(2) - m(3) - t189) * g(3) + m(3) * (qJDD(1) * t126 ^ 2 + (t129 * t83 + t132 * t82) * t123) + m(4) * (t112 * t126 - t36 * t79 + t37 * t80 - t64 * t72 - t65 * t73) + m(5) * (t153 * t8 - t19 * t38 + t20 * t39 + t32 * t79 + t60 * t7 + t61 * t72) + m(6) * (t1 * t22 + t10 * t4 - t153 * t6 + t19 * t34 + t2 * t21 + t3 * t9); t95 * (qJD(4) * t136 - t158 * t193) / 0.2e1 + t115 * (qJD(4) * t135 - t155 * t193) / 0.2e1 + t225 * t118 * t128 - (t127 * t51 + t130 * t50) * t193 / 0.2e1 + qJDD(4) * (Ifges(5,5) * t128 + Ifges(5,6) * t131) + t119 * t63 + t71 * t97 - t108 * t178 + (t82 + t170) * mrSges(3,1) + t216 * (t118 * t195 - t219) + ((t121 * t37 + t124 * t36) * pkin(2) + t64 * t71 - t65 * t74) * m(4) + (qJD(2) * t71 + t124 * t213 + t36) * mrSges(4,1) + (-m(4) * t169 - t260 * mrSges(3,1) - (-t122 * t132 - t125 * t206) * mrSges(3,2) + t259 * (t44 * pkin(3) + t169) + t244 * t43 + t245 * t44) * g(2) + (-m(4) * t114 - m(6) * t215 - t142 * t80 - m(5) * (pkin(7) * t80 + t215) + t262 * t79 - t251) * g(3) + t56 * t29 + t57 * t30 - t32 * t165 + qJD(4) ^ 2 * t157 / 0.2e1 + (-t1 * t204 - t10 * t139 - t140 * t9 - t2 * t202) * mrSges(6,3) + (qJD(2) * t74 - t121 * t213 - t37) * mrSges(4,2) + t252 * t195 / 0.2e1 + t34 * (mrSges(6,1) * t139 + mrSges(6,2) * t140) + (t9 * mrSges(6,1) - t10 * mrSges(6,2) - t246) * t197 + (g(1) * t46 + g(2) * t43 - g(3) * t80 - t195 * t38 - t197 * t39 + t247) * mrSges(5,3) + (t119 * t32 + ((-t39 * t128 - t38 * t131) * qJD(4) + t247) * t118 - t61 * t71 - (-t128 * t38 + t131 * t39) * t74) * m(5) - t11 * t204 / 0.2e1 + t12 * t202 / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + qJD(4) * t147 + (t123 * t175 - t83) * mrSges(3,2) + ((-Ifges(5,2) * t128 + t223) * t171 + Ifges(5,4) * t257 + Ifges(5,2) * t258 - t187 / 0.2e1 + t118 * t84 - t74 * t108 - Ifges(6,3) * t238 - Ifges(6,6) * t239 - Ifges(6,5) * t240 + t243) * t131 + (Ifges(5,1) * t102 + Ifges(5,4) * t258 + t156 * t238 + t159 * t239 + t162 * t240) * t128 + (-m(4) * t141 - t149 * mrSges(3,1) - (t122 * t206 - t125 * t132) * mrSges(3,2) + t259 * (t47 * pkin(3) + t141) + t244 * t46 + t245 * t47) * g(1) + t255 * t68 + t256 * t67 + (t1 * t57 + t2 * t56 + (t195 * t34 + t230) * t118 - t219 * t34 + t255 * t9 + t256 * t10) * m(6) + t146 * t171 + t163 * t230 + (qJD(4) * t137 - t161 * t193) * t236 + t176 * t241 + t223 * t257 + t160 * t258; m(4) * t112 + ((-t127 * t68 + t130 * t67 + t108) * qJD(4) + m(5) * (t8 + t214) + m(6) * (t10 * t196 - t198 * t9 - t6) - t225) * t131 + (t84 + (-t127 * t67 - t130 * t68) * qJD(5) + t216 * qJD(4) + m(5) * (-qJD(4) * t38 + t7) + m(6) * (qJD(4) * t34 + t261) + t248) * t128 + (-t126 * g(3) + (-g(1) * t122 + g(2) * t125) * t123) * t189; -t242 * t146 / 0.2e1 + t51 * t192 / 0.2e1 - t157 * t188 / 0.2e1 + t115 * t34 * t163 + (t115 * t156 + t159 * t95 + t162 * t96) * qJD(5) / 0.2e1 - (t115 * t135 + t95 * t136 + t137 * t96) * qJD(2) / 0.2e1 + (t190 * t50 + t12) * t127 / 0.2e1 + t130 * t11 / 0.2e1 + Ifges(5,6) * t101 + Ifges(5,5) * t102 - t16 * t67 - t15 * t68 + t261 * mrSges(6,3) - pkin(4) * t14 - t7 * mrSges(5,2) + t8 * mrSges(5,1) + t6 * t164 + (t253 * t26 - t254 * (-t125 * t208 + t128 * t43)) * g(2) + (t253 * t28 - t254 * (t122 * t208 + t128 * t46)) * g(1) + (-t153 * t254 + t253 * t60) * g(3) - (-Ifges(5,2) * t191 + t120 + t252) * t190 / 0.2e1 + (-m(6) * t34 + t184 - t216) * t39 + (-t9 * (mrSges(6,1) * t128 - mrSges(6,3) * t201) - t10 * (-mrSges(6,2) * t128 - mrSges(6,3) * t203) - t147) * qJD(2) + t246 * t191 + (-t68 * t192 - t67 * t194 + m(6) * ((-t10 * t127 - t9 * t130) * qJD(5) + t166) + t248) * pkin(8) + (t183 - t108) * t38 + (-pkin(4) * t6 - t10 * t16 - t15 * t9) * m(6) + t155 * t238 + t158 * t239 + t161 * t240 + t194 * t241 + Ifges(5,3) * qJDD(4); -t34 * (mrSges(6,1) * t96 + mrSges(6,2) * t95) + (Ifges(6,1) * t95 - t234) * t264 + t50 * t236 + (Ifges(6,5) * t95 - Ifges(6,6) * t96) * t263 - t9 * t67 + t10 * t68 - g(1) * ((-t127 * t28 - t130 * t47) * mrSges(6,1) + (t127 * t47 - t130 * t28) * mrSges(6,2)) - g(2) * ((-t127 * t26 - t130 * t44) * mrSges(6,1) + (t127 * t44 - t130 * t26) * mrSges(6,2)) - g(3) * (mrSges(6,1) * t21 - mrSges(6,2) * t22) + (t10 * t96 + t9 * t95) * mrSges(6,3) + t187 + (-Ifges(6,2) * t96 + t51 + t91) * t265 - t243;];
tau = t17;
