% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:17
% EndTime: 2019-03-09 02:15:24
% DurationCPUTime: 4.58s
% Computational Cost: add. (4777->409), mult. (10898->518), div. (0->0), fcn. (7276->6), ass. (0->188)
t263 = -Ifges(5,1) / 0.2e1;
t261 = Ifges(6,1) + Ifges(7,1);
t258 = Ifges(6,5) + Ifges(7,4);
t262 = Ifges(6,6) - Ifges(7,6);
t130 = sin(pkin(9));
t131 = cos(pkin(9));
t227 = sin(qJ(4));
t228 = cos(qJ(4));
t106 = t227 * t130 - t228 * t131;
t107 = t130 * t228 + t131 * t227;
t101 = t107 * qJD(1);
t129 = qJD(1) * qJ(2);
t124 = qJD(3) + t129;
t126 = t130 * pkin(3);
t111 = qJD(1) * t126 + t124;
t133 = sin(qJ(5));
t134 = cos(qJ(5));
t132 = -pkin(1) - qJ(3);
t115 = qJD(1) * t132 + qJD(2);
t178 = -pkin(7) * qJD(1) + t115;
t96 = t178 * t130;
t97 = t178 * t131;
t70 = t227 * t97 + t228 * t96;
t62 = qJD(4) * pkin(8) + t70;
t102 = t106 * qJD(1);
t63 = pkin(4) * t101 + pkin(8) * t102 + t111;
t21 = -t133 * t62 + t134 * t63;
t22 = t133 * t63 + t134 * t62;
t148 = t22 * t133 + t21 * t134;
t248 = qJD(6) - t21;
t254 = qJD(5) + t101;
t15 = -pkin(5) * t254 + t248;
t16 = qJ(6) * t254 + t22;
t150 = t16 * t133 - t15 * t134;
t204 = Ifges(7,5) * t134;
t154 = Ifges(7,3) * t133 + t204;
t156 = Ifges(6,5) * t134 - Ifges(6,6) * t133;
t158 = Ifges(7,4) * t134 + Ifges(7,6) * t133;
t206 = Ifges(6,4) * t134;
t160 = -Ifges(6,2) * t133 + t206;
t205 = Ifges(7,5) * t133;
t162 = Ifges(7,1) * t134 + t205;
t207 = Ifges(6,4) * t133;
t164 = Ifges(6,1) * t134 - t207;
t165 = mrSges(7,1) * t133 - mrSges(7,3) * t134;
t167 = mrSges(6,1) * t133 + mrSges(6,2) * t134;
t229 = t134 / 0.2e1;
t143 = t134 * qJD(4) + t102 * t133;
t69 = -t227 * t96 + t228 * t97;
t61 = -qJD(4) * pkin(4) - t69;
t87 = qJD(4) * t133 - t102 * t134;
t23 = -pkin(5) * t143 - t87 * qJ(6) + t61;
t231 = t133 / 0.2e1;
t232 = -t133 / 0.2e1;
t233 = t254 / 0.2e1;
t236 = t87 / 0.2e1;
t238 = -t143 / 0.2e1;
t239 = t143 / 0.2e1;
t222 = Ifges(7,5) * t143;
t85 = Ifges(6,4) * t143;
t251 = t258 * t254 + t261 * t87 - t222 + t85;
t84 = Ifges(7,5) * t87;
t31 = Ifges(7,6) * t254 - Ifges(7,3) * t143 + t84;
t223 = Ifges(6,4) * t87;
t34 = Ifges(6,2) * t143 + Ifges(6,6) * t254 + t223;
t245 = t150 * mrSges(7,2) + t148 * mrSges(6,3) - t154 * t238 - t160 * t239 - t165 * t23 - t167 * t61 - t231 * t31 - t232 * t34 - (t162 + t164) * t236 - (t156 + t158) * t233 - t251 * t229;
t259 = t102 * t263;
t260 = t111 * mrSges(5,2) - t69 * mrSges(5,3) - Ifges(5,4) * t101 + Ifges(5,5) * qJD(4) - t245 + t259;
t257 = t106 * qJD(3);
t256 = t258 * t133 + t262 * t134;
t255 = t261 * t133 - t204 + t206;
t253 = -t87 / 0.2e1;
t234 = -t254 / 0.2e1;
t94 = qJD(4) * t101;
t55 = qJD(5) * t143 - t134 * t94;
t56 = qJD(5) * t87 - t133 * t94;
t95 = t102 * qJD(4);
t252 = -t258 * t95 + (-Ifges(6,4) + Ifges(7,5)) * t56 + t261 * t55;
t193 = t130 ^ 2 + t131 ^ 2;
t250 = mrSges(4,3) * t193;
t151 = pkin(5) * t133 - qJ(6) * t134;
t249 = -qJD(6) * t133 + t151 * t254 - t70;
t215 = -pkin(7) + t132;
t109 = t215 * t130;
t110 = t215 * t131;
t82 = t227 * t109 - t228 * t110;
t191 = qJD(5) * t134;
t192 = qJD(5) * t133;
t139 = t107 * qJD(3);
t43 = -qJD(1) * t139 + t69 * qJD(4);
t190 = qJD(1) * qJD(2);
t68 = -pkin(4) * t95 + pkin(8) * t94 + t190;
t3 = t133 * t68 + t134 * t43 + t63 * t191 - t192 * t62;
t4 = -qJD(5) * t22 - t133 * t43 + t134 * t68;
t171 = -t133 * t4 + t134 * t3;
t1 = -qJ(6) * t95 + qJD(6) * t254 + t3;
t2 = pkin(5) * t95 - t4;
t173 = t1 * t134 + t133 * t2;
t121 = qJ(2) + t126;
t79 = pkin(4) * t107 + pkin(8) * t106 + t121;
t83 = t109 * t228 + t110 * t227;
t209 = t133 * t79 + t134 * t83;
t57 = -t82 * qJD(4) - t139;
t179 = qJD(4) * t227;
t180 = qJD(4) * t228;
t103 = -t130 * t180 - t131 * t179;
t104 = -t130 * t179 + t131 * t180;
t77 = pkin(4) * t104 - pkin(8) * t103 + qJD(2);
t9 = -qJD(5) * t209 - t133 * t57 + t134 * t77;
t224 = mrSges(6,3) * t87;
t66 = mrSges(6,1) * t254 - t224;
t67 = -mrSges(7,1) * t254 + mrSges(7,2) * t87;
t210 = t66 - t67;
t64 = mrSges(7,2) * t143 + mrSges(7,3) * t254;
t225 = mrSges(6,3) * t143;
t65 = -mrSges(6,2) * t254 + t225;
t211 = t64 + t65;
t141 = -t133 * t211 - t134 * t210;
t246 = -m(6) * t148 - m(7) * t150 + t141;
t186 = -Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1;
t187 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t188 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t203 = t101 * Ifges(5,2);
t243 = t186 * t254 + t187 * t143 - t188 * t87 + Ifges(5,6) * qJD(4) + t70 * mrSges(5,3) - t111 * mrSges(5,1) - t16 * mrSges(7,3) + t15 * mrSges(7,1) + t22 * mrSges(6,2) - t21 * mrSges(6,1) + Ifges(6,6) * t238 + Ifges(7,6) * t239 - t102 * Ifges(5,4) - t203 / 0.2e1 + t258 * t253 + (Ifges(6,3) + Ifges(7,2)) * t234;
t241 = -t56 / 0.2e1;
t240 = t56 / 0.2e1;
t230 = -t134 / 0.2e1;
t226 = m(3) * qJ(2);
t44 = -qJD(1) * t257 + qJD(4) * t70;
t217 = t44 * t82;
t37 = -mrSges(6,1) * t95 - mrSges(6,3) * t55;
t38 = t95 * mrSges(7,1) + t55 * mrSges(7,2);
t214 = -t37 + t38;
t39 = mrSges(6,2) * t95 - mrSges(6,3) * t56;
t40 = -mrSges(7,2) * t56 - mrSges(7,3) * t95;
t213 = t39 + t40;
t212 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t143 + mrSges(6,2) * t87 - t102 * mrSges(5,3);
t80 = -pkin(4) * t102 + pkin(8) * t101;
t27 = t133 * t80 + t134 * t69;
t208 = m(4) * qJD(3);
t200 = t106 * t44;
t169 = mrSges(4,1) * t130 + mrSges(4,2) * t131;
t198 = mrSges(5,1) * t101 - mrSges(5,2) * t102 + qJD(1) * t169;
t195 = t104 * t133;
t194 = t104 * t134;
t46 = -mrSges(7,1) * t143 - mrSges(7,3) * t87;
t189 = -t46 - t212;
t181 = -t95 * mrSges(5,1) - t94 * mrSges(5,2);
t177 = qJD(1) * t193;
t172 = t1 * t133 - t134 * t2;
t170 = t133 * t3 + t134 * t4;
t168 = -mrSges(6,1) * t134 + mrSges(6,2) * t133;
t166 = -mrSges(7,1) * t134 - mrSges(7,3) * t133;
t159 = Ifges(6,2) * t134 + t207;
t153 = -Ifges(7,3) * t134 + t205;
t152 = -pkin(5) * t134 - qJ(6) * t133;
t149 = t133 * t15 + t134 * t16;
t147 = -t133 * t21 + t134 * t22;
t26 = -t133 * t69 + t134 * t80;
t29 = -t133 * t83 + t134 * t79;
t8 = t133 * t77 + t134 * t57 + t79 * t191 - t192 * t83;
t142 = t133 * t214 + t134 * t213;
t140 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t58 = qJD(4) * t83 - t257;
t135 = qJD(1) ^ 2;
t112 = -pkin(4) + t152;
t93 = Ifges(7,2) * t95;
t92 = Ifges(6,3) * t95;
t89 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t101;
t53 = Ifges(7,4) * t55;
t52 = Ifges(6,5) * t55;
t51 = Ifges(6,6) * t56;
t50 = Ifges(7,6) * t56;
t45 = pkin(5) * t87 - qJ(6) * t143;
t42 = -t106 * t151 + t82;
t25 = -pkin(5) * t107 - t29;
t24 = qJ(6) * t107 + t209;
t20 = mrSges(6,1) * t56 + mrSges(6,2) * t55;
t19 = mrSges(7,1) * t56 - mrSges(7,3) * t55;
t18 = pkin(5) * t102 - t26;
t17 = -qJ(6) * t102 + t27;
t12 = t55 * Ifges(6,4) - t56 * Ifges(6,2) - t95 * Ifges(6,6);
t11 = t55 * Ifges(7,5) - t95 * Ifges(7,6) + t56 * Ifges(7,3);
t10 = t151 * t103 + (qJD(5) * t152 + qJD(6) * t134) * t106 + t58;
t7 = -pkin(5) * t104 - t9;
t6 = t56 * pkin(5) - t55 * qJ(6) - t87 * qJD(6) + t44;
t5 = qJ(6) * t104 + qJD(6) * t107 + t8;
t13 = [(-mrSges(5,2) * t190 - t6 * t165 + t154 * t241 + t160 * t240 + t11 * t232 + Ifges(5,1) * t94 + (-mrSges(5,3) - t167) * t44 + t170 * mrSges(6,3) + t172 * mrSges(7,2) + (mrSges(7,2) * t149 + mrSges(6,3) * t147 + t153 * t238 + t159 * t239 + t166 * t23 + t168 * t61 + t229 * t34 + t233 * t256 + t236 * t255) * qJD(5) + (-t162 / 0.2e1 - t164 / 0.2e1) * t55 + (t251 * qJD(5) + t12) * t231 + (qJD(5) * t31 + t252) * t230 - (-t158 / 0.2e1 - t156 / 0.2e1 + Ifges(5,4)) * t95) * t106 + (-t82 * t94 + t83 * t95) * mrSges(5,3) + (t259 + t260) * t103 + (-t43 * mrSges(5,3) + mrSges(5,1) * t190 + t52 / 0.2e1 - t51 / 0.2e1 - t92 / 0.2e1 + t53 / 0.2e1 - t93 / 0.2e1 + t50 / 0.2e1 + Ifges(5,4) * t94 + t187 * t56 + t188 * t55 - (Ifges(5,2) - t186) * t95 + t140) * t107 + (m(5) * (qJD(1) * t121 + t111) + m(4) * (t124 + t129) + t198 + ((2 * mrSges(3,3)) + t169 + 0.2e1 * t226) * qJD(1)) * qJD(2) + m(7) * (t1 * t24 + t10 * t23 + t15 * t7 + t16 * t5 + t2 * t25 + t42 * t6) + (t203 / 0.2e1 - t243) * t104 + (-t115 * t193 - t132 * t177) * t208 + t209 * t39 + m(6) * (t209 * t3 + t21 * t9 + t22 * t8 + t29 * t4 + t58 * t61 + t217) + 0.2e1 * qJD(3) * t250 * qJD(1) + t29 * t37 + t25 * t38 + t24 * t40 + t42 * t19 + t10 * t46 + t121 * t181 + t5 * t64 + t8 * t65 + t9 * t66 + t7 * t67 + t82 * t20 + t57 * t89 + t212 * t58 + m(5) * (t43 * t83 + t57 * t70 - t58 * t69 + t217); (-mrSges(3,3) - t226) * t135 + (-t94 * mrSges(5,3) + t19 + t20) * t106 + t189 * t103 + (-t133 * t210 + t134 * t211 + t89) * t104 + m(5) * (t103 * t69 + t104 * t70 + t200) + m(7) * (-t103 * t23 + t106 * t6 + t15 * t195 + t16 * t194) + m(6) * (-t103 * t61 + t194 * t22 - t195 * t21 + t200) + (t95 * mrSges(5,3) + t141 * qJD(5) + m(5) * t43 + m(7) * (t15 * t191 - t16 * t192 + t173) + m(6) * (-t191 * t21 - t192 * t22 + t171) + t142) * t107 + (-m(4) * t124 - m(5) * t111 - t193 * t208 - t198 + t246) * qJD(1); t101 * t89 + (m(4) + m(5)) * t190 - t135 * t250 - t189 * t102 + (t211 * t254 - t214) * t134 + (-t210 * t254 + t213) * t133 - m(5) * (-t101 * t70 + t102 * t69) + m(4) * t115 * t177 + t181 + (t102 * t23 + t149 * t254 + t172) * m(7) + (t102 * t61 + t147 * t254 + t170) * m(6); t255 * t55 / 0.2e1 - (-(Ifges(5,2) / 0.2e1 + t263) * t102 - t260) * t101 - (t256 / 0.2e1 - Ifges(5,6)) * t95 - t245 * qJD(5) - t243 * t102 + (-pkin(4) * t44 - t21 * t26 - t22 * t27 - t61 * t70) * m(6) + t153 * t240 + t159 * t241 + t12 * t229 + t11 * t230 + (m(6) * t171 + m(7) * t173 + t246 * qJD(5) + t142) * pkin(8) + (t112 * t6 - t15 * t18 - t16 * t17 + t249 * t23) * m(7) + t249 * t46 + t252 * t231 - pkin(4) * t20 + t6 * t166 + (-mrSges(5,1) + t168) * t44 + t171 * mrSges(6,3) + t173 * mrSges(7,2) - t43 * mrSges(5,2) - t17 * t64 - t27 * t65 - t26 * t66 - t18 * t67 - t69 * t89 - Ifges(5,5) * t94 - t212 * t70 + t112 * t19; (-t143 * t15 + t16 * t87) * mrSges(7,2) - t93 - t92 + t53 + t52 - t51 + t50 + t34 * t236 + (Ifges(7,3) * t87 + t222) * t239 + t140 - pkin(5) * t38 + qJ(6) * t40 - t45 * t46 + qJD(6) * t64 - t23 * (mrSges(7,1) * t87 - mrSges(7,3) * t143) - t61 * (mrSges(6,1) * t87 + mrSges(6,2) * t143) + (t210 + t224) * t22 + (-t211 + t225) * t21 + (t258 * t143 - t262 * t87) * t234 + (-t2 * pkin(5) + t1 * qJ(6) - t15 * t22 + t248 * t16 - t23 * t45) * m(7) + (-Ifges(6,2) * t87 + t251 + t85) * t238 + (t143 * t261 - t223 + t31 + t84) * t253; t87 * t46 - t254 * t64 + 0.2e1 * (t2 / 0.2e1 + t16 * t234 + t23 * t236) * m(7) + t38;];
tauc  = t13(:);
