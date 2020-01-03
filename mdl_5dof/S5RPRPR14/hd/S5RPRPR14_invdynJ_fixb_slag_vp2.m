% Calculate vector of inverse dynamics joint torques for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:47
% DurationCPUTime: 8.89s
% Computational Cost: add. (3321->444), mult. (6563->617), div. (0->0), fcn. (4106->10), ass. (0->201)
t274 = m(6) + m(5);
t273 = -mrSges(5,2) + mrSges(6,3);
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t166 = -mrSges(6,1) * t134 + mrSges(6,2) * t131;
t272 = -m(6) * pkin(4) + t166;
t240 = cos(qJ(3));
t270 = t240 / 0.2e1;
t133 = sin(qJ(1));
t135 = cos(qJ(1));
t262 = -g(1) * t133 + g(2) * t135;
t127 = qJ(3) + pkin(8);
t120 = sin(t127);
t121 = cos(t127);
t132 = sin(qJ(3));
t153 = mrSges(4,1) * t132 + mrSges(4,2) * t240;
t269 = -mrSges(5,1) * t120 + t121 * t273 - t153;
t243 = -m(3) - m(4);
t129 = sin(pkin(8));
t208 = cos(pkin(8));
t192 = qJD(1) * qJD(3);
t96 = qJDD(1) * t240 - t132 * t192;
t176 = qJD(3) * t240;
t97 = -qJD(1) * t176 - qJDD(1) * t132;
t56 = t129 * t97 + t208 * t96;
t49 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t56;
t169 = t208 * t240;
t199 = qJD(1) * t132;
t83 = -qJD(1) * t169 + t129 * t199;
t62 = qJD(3) * t134 + t131 * t83;
t26 = qJD(5) * t62 + qJDD(3) * t131 + t134 * t56;
t63 = qJD(3) * t131 - t134 * t83;
t27 = -qJD(5) * t63 + qJDD(3) * t134 - t131 * t56;
t7 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t267 = t7 - t49;
t239 = mrSges(5,3) * t83;
t266 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t62 - mrSges(6,2) * t63 + t239;
t265 = -qJDD(3) / 0.2e1;
t144 = -t129 * t240 - t132 * t208;
t84 = t144 * qJD(1);
t264 = qJD(5) - t84;
t136 = -pkin(1) - pkin(6);
t104 = qJDD(1) * t136 + qJDD(2);
t105 = qJD(1) * t136 + qJD(2);
t198 = qJD(3) * t132;
t60 = t104 * t240 - t105 * t198;
t61 = t104 * t132 + t105 * t176;
t152 = t132 * t61 + t240 * t60;
t193 = qJD(1) * qJD(2);
t106 = qJDD(1) * qJ(2) + t193;
t55 = -t129 * t96 + t208 * t97;
t54 = qJDD(5) - t55;
t10 = mrSges(6,1) * t54 - mrSges(6,3) * t26;
t11 = -mrSges(6,2) * t54 + mrSges(6,3) * t27;
t263 = -t10 * t131 + t11 * t134;
t156 = mrSges(4,1) * t240 - mrSges(4,2) * t132;
t185 = Ifges(4,4) * t240;
t261 = qJ(2) * t156 + (-Ifges(4,1) * t132 - t185) * t270;
t260 = m(4) * t152 + t132 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t97);
t259 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t175 = t240 * qJD(4);
t33 = qJDD(3) * pkin(3) - qJ(4) * t96 - qJD(1) * t175 + t60;
t194 = t132 * qJD(4);
t39 = qJ(4) * t97 - qJD(1) * t194 + t61;
t14 = -t129 * t39 + t208 * t33;
t15 = t129 * t33 + t208 * t39;
t73 = (-qJ(4) * qJD(1) + t105) * t132;
t216 = t129 * t73;
t177 = qJD(1) * t240;
t98 = t240 * t105;
t74 = -qJ(4) * t177 + t98;
t71 = qJD(3) * pkin(3) + t74;
t34 = t208 * t71 - t216;
t69 = t208 * t73;
t35 = t129 * t71 + t69;
t85 = -qJD(3) * t169 + t129 * t198;
t86 = t144 * qJD(3);
t90 = t129 * t132 - t169;
t258 = t14 * t90 + t144 * t15 - t34 * t86 + t35 * t85;
t31 = qJD(3) * pkin(7) + t35;
t99 = pkin(3) * t199 + qJD(1) * qJ(2) + qJD(4);
t38 = -pkin(4) * t84 + pkin(7) * t83 + t99;
t12 = -t131 * t31 + t134 * t38;
t65 = -pkin(3) * t97 + qJDD(4) + t106;
t18 = -pkin(4) * t55 - pkin(7) * t56 + t65;
t9 = qJDD(3) * pkin(7) + t15;
t1 = qJD(5) * t12 + t131 * t18 + t134 * t9;
t13 = t131 * t38 + t134 * t31;
t2 = -qJD(5) * t13 - t131 * t9 + t134 * t18;
t257 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t234 = pkin(7) * t121;
t256 = mrSges(2,2) - m(6) * (pkin(4) * t120 - t234) - mrSges(3,3) + t269;
t254 = qJD(1) ^ 2;
t253 = t26 / 0.2e1;
t252 = t27 / 0.2e1;
t251 = t54 / 0.2e1;
t250 = -t62 / 0.2e1;
t249 = -t63 / 0.2e1;
t248 = t63 / 0.2e1;
t247 = -t264 / 0.2e1;
t246 = -t83 / 0.2e1;
t8 = -qJDD(3) * pkin(4) - t14;
t242 = t8 * t90;
t241 = t131 / 0.2e1;
t238 = mrSges(5,3) * t84;
t237 = Ifges(5,4) * t83;
t236 = Ifges(6,4) * t63;
t235 = pkin(3) * t129;
t125 = t132 * pkin(3);
t223 = mrSges(6,3) * t131;
t222 = mrSges(6,3) * t134;
t221 = Ifges(4,4) * t132;
t220 = Ifges(6,4) * t131;
t219 = Ifges(6,4) * t134;
t214 = t131 * t90;
t211 = t134 * t90;
t207 = t131 * t133;
t206 = t131 * t135;
t102 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t177;
t205 = t132 * t102;
t204 = t133 * t134;
t203 = t134 * t135;
t115 = qJ(2) + t125;
t202 = qJ(4) - t136;
t200 = pkin(1) * t135 + qJ(2) * t133;
t197 = qJD(5) * t131;
t196 = qJD(5) * t134;
t195 = qJDD(1) * mrSges(3,2);
t107 = pkin(3) * t176 + qJD(2);
t190 = Ifges(6,5) * t26 + Ifges(6,6) * t27 + Ifges(6,3) * t54;
t22 = Ifges(6,2) * t62 + Ifges(6,6) * t264 + t236;
t183 = t22 * t241;
t59 = Ifges(6,4) * t62;
t23 = Ifges(6,1) * t63 + Ifges(6,5) * t264 + t59;
t181 = -t134 * t23 / 0.2e1;
t180 = t240 * t136;
t179 = t208 * pkin(3);
t178 = -t55 * mrSges(5,1) + mrSges(5,2) * t56;
t174 = t196 / 0.2e1;
t124 = t135 * qJ(2);
t173 = -pkin(1) * t133 + t124;
t172 = -t192 / 0.2e1;
t171 = (t106 + t193) * qJ(2);
t170 = pkin(3) * t177;
t168 = -g(1) * t135 - g(2) * t133;
t165 = mrSges(6,1) * t131 + mrSges(6,2) * t134;
t164 = Ifges(6,1) * t134 - t220;
t163 = -Ifges(6,2) * t131 + t219;
t162 = Ifges(6,5) * t134 - Ifges(6,6) * t131;
t161 = t12 * t134 + t13 * t131;
t160 = t12 * t131 - t13 * t134;
t36 = -mrSges(6,2) * t264 + mrSges(6,3) * t62;
t37 = mrSges(6,1) * t264 - mrSges(6,3) * t63;
t159 = -t131 * t36 - t134 * t37;
t51 = -pkin(4) * t144 + pkin(7) * t90 + t115;
t100 = t202 * t132;
t147 = -qJ(4) * t240 + t180;
t58 = -t100 * t208 + t129 * t147;
t19 = -t131 * t58 + t134 * t51;
t20 = t131 * t51 + t134 * t58;
t155 = t240 * Ifges(4,1) - t221;
t154 = -Ifges(4,2) * t132 + t185;
t151 = -Ifges(4,5) * t132 - Ifges(4,6) * t240;
t150 = -t131 * t86 + t196 * t90;
t149 = t134 * t86 + t197 * t90;
t145 = t132 * (-Ifges(4,2) * t240 - t221);
t139 = t198 * t202 - t175;
t138 = -qJD(5) * t161 + t1 * t134 - t131 * t2;
t130 = -qJ(4) - pkin(6);
t119 = -pkin(1) * qJDD(1) + qJDD(2);
t114 = -t179 - pkin(4);
t101 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t199;
t93 = t153 * qJD(1);
t88 = Ifges(4,5) * qJD(3) + qJD(1) * t155;
t87 = Ifges(4,6) * qJD(3) + qJD(1) * t154;
t81 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t96;
t80 = t120 * t203 - t207;
t79 = t120 * t206 + t204;
t78 = t120 * t204 + t206;
t77 = -t120 * t207 + t203;
t75 = Ifges(5,4) * t84;
t72 = qJD(3) * t147 - t194;
t67 = -qJD(3) * mrSges(5,2) + t238;
t50 = -mrSges(5,1) * t84 - mrSges(5,2) * t83;
t48 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t55;
t47 = -t83 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t75;
t46 = t84 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t237;
t45 = -pkin(4) * t83 - pkin(7) * t84 + t170;
t44 = -pkin(4) * t85 - pkin(7) * t86 + t107;
t43 = t208 * t74 - t216;
t42 = t129 * t74 + t69;
t41 = t129 * t139 + t208 * t72;
t30 = -qJD(3) * pkin(4) - t34;
t21 = Ifges(6,5) * t63 + t62 * Ifges(6,6) + Ifges(6,3) * t264;
t17 = t131 * t45 + t134 * t43;
t16 = -t131 * t43 + t134 * t45;
t6 = -qJD(5) * t20 - t131 * t41 + t134 * t44;
t5 = qJD(5) * t19 + t131 * t44 + t134 * t41;
t4 = Ifges(6,1) * t26 + Ifges(6,4) * t27 + Ifges(6,5) * t54;
t3 = Ifges(6,4) * t26 + Ifges(6,2) * t27 + Ifges(6,6) * t54;
t24 = [(-m(5) * t34 + m(6) * t30 - t266) * (t129 * t72 - t139 * t208) + (qJD(5) * t23 + t3) * t214 / 0.2e1 + t99 * (-mrSges(5,1) * t85 + mrSges(5,2) * t86) + t84 * (Ifges(5,4) * t86 + Ifges(5,2) * t85) / 0.2e1 + qJD(3) * (Ifges(5,5) * t86 + Ifges(5,6) * t85) / 0.2e1 + (Ifges(5,1) * t86 + Ifges(5,4) * t85) * t246 + t145 * t172 + m(5) * (t107 * t99 + t115 * t65 + t15 * t58 + t35 * t41) + m(6) * (t1 * t20 + t12 * t6 + t13 * t5 + t19 * t2) + t119 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t97 + mrSges(4,2) * t96) + t107 * t50 + qJD(2) * t93 + t86 * t47 / 0.2e1 + t85 * t46 / 0.2e1 - t85 * t21 / 0.2e1 + t41 * t67 + t58 * t48 + (Ifges(4,1) * t96 + Ifges(4,4) * t97) * t270 + (t1 * t214 - t12 * t149 + t13 * t150 + t2 * t211) * mrSges(6,3) + (t101 * t176 - t102 * t198 + t260) * t136 - t132 * (Ifges(4,4) * t96 + Ifges(4,2) * t97) / 0.2e1 - t152 * mrSges(4,3) + t261 * t192 + t5 * t36 + t6 * t37 + t19 * t10 + t20 * t11 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t65 * mrSges(5,2) + Ifges(5,5) * t265 - t162 * t251 - t163 * t252 - t164 * t253 + t22 * t174) * t90 + t258 * mrSges(5,3) + t81 * t180 + (-m(5) * t14 + m(6) * t8 + t267) * (-t100 * t129 - t147 * t208) + (-m(3) * t173 - m(4) * t124 - t80 * mrSges(6,1) + t79 * mrSges(6,2) - t274 * (t125 * t135 + t130 * t133 + t173) + (-m(4) * t136 - t259) * t133 + t256 * t135) * g(1) + (-t78 * mrSges(6,1) - t77 * mrSges(6,2) + t243 * t200 - t274 * (t125 * t133 - t130 * t135 + t200) + (-m(4) * pkin(6) + t259) * t135 + t256 * t133) * g(2) + m(4) * t171 + t62 * (Ifges(6,4) * t149 + Ifges(6,2) * t150 - Ifges(6,6) * t85) / 0.2e1 + t30 * (-mrSges(6,1) * t150 + mrSges(6,2) * t149) + qJD(3) ^ 2 * t151 / 0.2e1 + t97 * t154 / 0.2e1 + t96 * t155 / 0.2e1 - t86 * t181 - t86 * t183 + t264 * (Ifges(6,5) * t149 + Ifges(6,6) * t150 - Ifges(6,3) * t85) / 0.2e1 + t56 * (-Ifges(5,1) * t90 + Ifges(5,4) * t144) + t55 * (-Ifges(5,4) * t90 + Ifges(5,2) * t144) - (t65 * mrSges(5,1) + Ifges(5,6) * t265 + Ifges(6,3) * t251 + Ifges(6,6) * t252 + Ifges(6,5) * t253 + t190 / 0.2e1 + t257) * t144 + (Ifges(4,5) * t240 - Ifges(5,5) * t90 / 0.2e1 - Ifges(4,6) * t132 + Ifges(5,6) * t144 / 0.2e1) * qJDD(3) + (t153 + 0.2e1 * mrSges(3,3)) * t106 + (Ifges(6,1) * t149 + Ifges(6,4) * t150 - Ifges(6,5) * t85) * t248 + m(3) * (-pkin(1) * t119 + t171) - t87 * t176 / 0.2e1 + t115 * t178 - t12 * mrSges(6,1) * t85 + t13 * mrSges(6,2) * t85 - pkin(1) * t195 - t88 * t198 / 0.2e1 - t4 * t211 / 0.2e1 - t165 * t242; t195 + t240 * t81 + t267 * t90 + t266 * t86 + (t101 * t240 - t205) * qJD(3) + (t131 * t37 - t134 * t36 - t67) * t85 + (qJ(2) * t243 - mrSges(3,3)) * t254 - (qJD(5) * t159 + t263 + t48) * t144 + m(6) * (-t138 * t144 + t160 * t85 - t30 * t86 + t242) + m(3) * t119 - m(5) * t258 + (-m(5) * t99 - m(6) * t161 + t159 - t50 - t93) * qJD(1) + t262 * (-t243 + t274) + t260; t266 * t42 - (Ifges(5,2) * t83 + t47 + t75) * t84 / 0.2e1 + (Ifges(5,1) * t84 + t21 + t237) * t83 / 0.2e1 + t151 * t172 + (m(5) * t125 - m(6) * (-t125 + t234) - t272 * t120 - t269) * g(3) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t134 * t3 / 0.2e1 + t114 * t7 + Ifges(4,5) * t96 + Ifges(4,6) * t97 - t99 * (-mrSges(5,1) * t83 + mrSges(5,2) * t84) - qJD(3) * (Ifges(5,5) * t84 + Ifges(5,6) * t83) / 0.2e1 - t43 * t67 + Ifges(5,5) * t56 + t60 * mrSges(4,1) - t61 * mrSges(4,2) + Ifges(5,6) * t55 + (t145 / 0.2e1 - t261) * t254 + (m(6) * t138 - t196 * t37 - t197 * t36 + t263) * (pkin(7) + t235) + (t156 + t274 * t240 * pkin(3) + (mrSges(5,1) - t272) * t121 + (m(6) * pkin(7) + t273) * t120) * t262 - t17 * t36 - t16 * t37 + t14 * mrSges(5,1) - t15 * mrSges(5,2) + (-t12 * t196 - t13 * t197) * mrSges(6,3) + t264 * t30 * t165 + t105 * t205 + t49 * t179 + (t114 * t8 - t12 * t16 - t13 * t17 - t30 * t42) * m(6) + ((t129 * t15 + t14 * t208) * pkin(3) - t99 * t170 + t34 * t42 - t35 * t43) * m(5) + t1 * t222 + t23 * t174 + t84 * t181 + t84 * t183 + t48 * t235 + t34 * t238 + t4 * t241 + t46 * t246 + (-Ifges(6,3) * t83 + t162 * t84) * t247 + (t162 * t264 + t163 * t62 + t164 * t63) * qJD(5) / 0.2e1 + (-Ifges(6,5) * t83 + t164 * t84) * t249 + (-Ifges(6,6) * t83 + t163 * t84) * t250 + (Ifges(6,5) * t131 + Ifges(6,6) * t134) * t251 + (Ifges(6,2) * t134 + t220) * t252 + (Ifges(6,1) * t131 + t219) * t253 + t8 * t166 - t50 * t170 + t87 * t177 / 0.2e1 - t101 * t98 - t22 * t197 / 0.2e1 + t88 * t199 / 0.2e1 - t12 * (-mrSges(6,1) * t83 - t222 * t84) - t13 * (mrSges(6,2) * t83 - t223 * t84) - t2 * t223 - t35 * t239; -t84 * t67 - t266 * t83 + (t264 * t36 + t10) * t134 + (-t264 * t37 + t11) * t131 + t178 + (t1 * t131 + t134 * t2 - t160 * t264 + t30 * t83 + t168) * m(6) + (-t34 * t83 - t35 * t84 + t168 + t65) * m(5); -t30 * (mrSges(6,1) * t63 + mrSges(6,2) * t62) + (Ifges(6,1) * t62 - t236) * t249 + t22 * t248 + (Ifges(6,5) * t62 - Ifges(6,6) * t63) * t247 - t12 * t36 + t13 * t37 - g(1) * (mrSges(6,1) * t77 - mrSges(6,2) * t78) - g(2) * (mrSges(6,1) * t79 + mrSges(6,2) * t80) + g(3) * t165 * t121 + (t12 * t62 + t13 * t63) * mrSges(6,3) + t190 + (-Ifges(6,2) * t63 + t23 + t59) * t250 + t257;];
tau = t24;
