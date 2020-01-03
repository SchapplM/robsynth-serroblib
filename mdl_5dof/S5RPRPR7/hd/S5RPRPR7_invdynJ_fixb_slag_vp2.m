% Calculate vector of inverse dynamics joint torques for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:55
% EndTime: 2019-12-31 18:19:09
% DurationCPUTime: 8.24s
% Computational Cost: add. (3410->453), mult. (7303->617), div. (0->0), fcn. (4760->14), ass. (0->208)
t265 = m(6) + m(5);
t264 = mrSges(5,2) - mrSges(6,3);
t130 = sin(pkin(8));
t112 = pkin(1) * t130 + pkin(6);
t101 = t112 * qJDD(1);
t263 = qJD(2) * qJD(3) + t101;
t134 = sin(qJ(3));
t137 = cos(qJ(3));
t107 = -mrSges(4,1) * t137 + mrSges(4,2) * t134;
t127 = qJ(3) + pkin(9);
t118 = sin(t127);
t120 = cos(t127);
t262 = -mrSges(5,1) * t120 + t264 * t118 + t107;
t133 = sin(qJ(5));
t136 = cos(qJ(5));
t129 = sin(pkin(9));
t199 = cos(pkin(9));
t184 = qJD(1) * qJD(3);
t97 = qJDD(1) * t137 - t134 * t184;
t98 = qJDD(1) * t134 + t137 * t184;
t60 = t129 * t97 + t199 * t98;
t170 = t199 * t134;
t185 = t137 * qJD(1);
t88 = -qJD(1) * t170 - t129 * t185;
t68 = qJD(3) * t136 + t133 * t88;
t26 = qJD(5) * t68 + qJDD(3) * t133 + t136 * t60;
t240 = t26 / 0.2e1;
t69 = qJD(3) * t133 - t136 * t88;
t27 = -qJD(5) * t69 + qJDD(3) * t136 - t133 * t60;
t239 = t27 / 0.2e1;
t59 = -t129 * t98 + t199 * t97;
t58 = qJDD(5) - t59;
t238 = t58 / 0.2e1;
t261 = -m(4) - m(3);
t54 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t60;
t7 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t259 = -t54 + t7;
t146 = -t129 * t134 + t137 * t199;
t87 = t146 * qJD(1);
t82 = Ifges(5,4) * t87;
t103 = t112 * qJD(1);
t166 = qJ(4) * qJD(1) + t103;
t190 = qJD(2) * t134;
t71 = t137 * t166 + t190;
t62 = t199 * t71;
t124 = t137 * qJD(2);
t70 = -t134 * t166 + t124;
t64 = qJD(3) * pkin(3) + t70;
t34 = t129 * t64 + t62;
t31 = qJD(3) * pkin(7) + t34;
t131 = cos(pkin(8));
t114 = -pkin(1) * t131 - pkin(2);
t125 = t137 * pkin(3);
t100 = t114 - t125;
t86 = qJD(1) * t100 + qJD(4);
t41 = -pkin(4) * t87 + pkin(7) * t88 + t86;
t12 = -t133 * t31 + t136 * t41;
t258 = t12 * mrSges(6,1);
t13 = t133 * t41 + t136 * t31;
t257 = t13 * mrSges(6,2);
t256 = t68 * Ifges(6,6);
t252 = qJD(5) - t87;
t255 = t252 * Ifges(6,3);
t254 = t87 * Ifges(5,2);
t218 = t88 * mrSges(5,3);
t253 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t68 - mrSges(6,2) * t69 + t218;
t200 = qJDD(3) / 0.2e1;
t186 = qJD(5) * t136;
t90 = t146 * qJD(3);
t95 = t129 * t137 + t170;
t151 = t133 * t90 + t95 * t186;
t16 = mrSges(6,1) * t58 - mrSges(6,3) * t26;
t17 = -mrSges(6,2) * t58 + mrSges(6,3) * t27;
t250 = -t133 * t16 + t136 * t17;
t189 = qJD(3) * t134;
t49 = t134 * qJDD(2) - t103 * t189 + t137 * t263;
t123 = t137 * qJDD(2);
t80 = t103 * t137 + t190;
t50 = -t80 * qJD(3) - t101 * t134 + t123;
t249 = -t134 * t50 + t137 * t49;
t102 = t114 * qJDD(1);
t67 = -pkin(3) * t97 + qJDD(4) + t102;
t18 = -pkin(4) * t59 - pkin(7) * t60 + t67;
t183 = qJD(1) * qJD(4);
t188 = qJD(3) * t137;
t32 = -t103 * t188 + qJDD(3) * pkin(3) - qJ(4) * t98 + t123 + (-t183 - t263) * t134;
t35 = qJ(4) * t97 + t137 * t183 + t49;
t11 = t129 * t32 + t199 * t35;
t9 = qJDD(3) * pkin(7) + t11;
t1 = qJD(5) * t12 + t133 * t18 + t136 * t9;
t2 = -qJD(5) * t13 - t133 * t9 + t136 * t18;
t248 = t1 * t136 - t133 * t2;
t247 = 0.2e1 * t200;
t245 = m(4) * pkin(2) + mrSges(3,1) - t262;
t244 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t243 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t242 = Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t238;
t237 = -t68 / 0.2e1;
t236 = -t69 / 0.2e1;
t235 = t69 / 0.2e1;
t234 = -t252 / 0.2e1;
t232 = -t88 / 0.2e1;
t230 = t136 / 0.2e1;
t229 = Ifges(5,4) * t88;
t135 = sin(qJ(1));
t228 = pkin(1) * t135;
t227 = pkin(3) * t129;
t225 = pkin(7) * t118;
t138 = cos(qJ(1));
t126 = t138 * pkin(1);
t210 = t129 * t71;
t33 = t199 * t64 - t210;
t220 = t33 * mrSges(5,3);
t219 = t69 * Ifges(6,4);
t215 = Ifges(4,4) * t134;
t214 = Ifges(4,4) * t137;
t213 = Ifges(6,4) * t133;
t212 = Ifges(6,4) * t136;
t208 = t133 * t87;
t206 = t133 * t95;
t203 = t136 * t87;
t202 = t136 * t95;
t128 = qJ(1) + pkin(8);
t119 = sin(t128);
t197 = t119 * t133;
t196 = t119 * t136;
t121 = cos(t128);
t195 = t121 * t133;
t194 = t121 * t136;
t193 = qJ(4) + t112;
t191 = qJD(1) * t134;
t187 = qJD(5) * t133;
t181 = Ifges(6,5) * t26 + Ifges(6,6) * t27 + Ifges(6,3) * t58;
t180 = pkin(3) * t191;
t179 = pkin(3) * t189;
t178 = mrSges(4,3) * t191;
t177 = mrSges(4,3) * t185;
t65 = Ifges(6,4) * t68;
t23 = Ifges(6,1) * t69 + Ifges(6,5) * t252 + t65;
t175 = t23 * t230;
t174 = t199 * pkin(3);
t173 = -t59 * mrSges(5,1) + t60 * mrSges(5,2);
t172 = -t187 / 0.2e1;
t168 = t193 * t134;
t167 = qJD(3) * t193;
t164 = pkin(4) * t120 + t225;
t163 = -g(1) * t119 + g(2) * t121;
t162 = mrSges(4,1) * t134 + mrSges(4,2) * t137;
t160 = -mrSges(6,1) * t136 + mrSges(6,2) * t133;
t159 = mrSges(6,1) * t133 + mrSges(6,2) * t136;
t158 = Ifges(6,1) * t136 - t213;
t157 = t137 * Ifges(4,2) + t215;
t156 = -Ifges(6,2) * t133 + t212;
t155 = Ifges(4,5) * t137 - Ifges(4,6) * t134;
t154 = Ifges(6,5) * t136 - Ifges(6,6) * t133;
t152 = -t12 * t133 + t13 * t136;
t46 = -pkin(4) * t146 - pkin(7) * t95 + t100;
t93 = t193 * t137;
t52 = -t129 * t168 + t199 * t93;
t19 = -t133 * t52 + t136 * t46;
t20 = t133 * t46 + t136 * t52;
t150 = -t136 * t90 + t187 * t95;
t30 = -qJD(3) * pkin(4) - t33;
t149 = t30 * t159;
t10 = -t129 * t35 + t199 * t32;
t148 = t114 * qJD(1) * t162;
t147 = t134 * (Ifges(4,1) * t137 - t215);
t144 = m(6) * pkin(4) - t160;
t142 = -qJD(4) * t134 - t137 * t167;
t141 = (-t12 * t136 - t13 * t133) * qJD(5) + t248;
t132 = -qJ(4) - pkin(6);
t117 = Ifges(4,4) * t185;
t116 = t125 + pkin(2);
t113 = -t174 - pkin(4);
t106 = -qJD(3) * mrSges(4,2) + t177;
t104 = qJD(3) * mrSges(4,1) - t178;
t92 = Ifges(4,1) * t191 + Ifges(4,5) * qJD(3) + t117;
t91 = Ifges(4,6) * qJD(3) + qJD(1) * t157;
t89 = t95 * qJD(3);
t85 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t98;
t84 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t97;
t79 = -t103 * t134 + t124;
t77 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t87;
t76 = t120 * t194 + t197;
t75 = -t120 * t195 + t196;
t74 = -t120 * t196 + t195;
t73 = t120 * t197 + t194;
t72 = qJD(4) * t137 - t134 * t167;
t55 = -mrSges(5,1) * t87 - mrSges(5,2) * t88;
t53 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t59;
t48 = -t88 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t82;
t47 = Ifges(5,6) * qJD(3) - t229 + t254;
t45 = pkin(4) * t89 - pkin(7) * t90 + t179;
t44 = -pkin(4) * t88 - pkin(7) * t87 + t180;
t43 = mrSges(6,1) * t252 - mrSges(6,3) * t69;
t42 = -mrSges(6,2) * t252 + mrSges(6,3) * t68;
t40 = t129 * t142 + t199 * t72;
t37 = t199 * t70 - t210;
t36 = t129 * t70 + t62;
t22 = t68 * Ifges(6,2) + Ifges(6,6) * t252 + t219;
t21 = t69 * Ifges(6,5) + t255 + t256;
t15 = t133 * t44 + t136 * t37;
t14 = -t133 * t37 + t136 * t44;
t8 = -qJDD(3) * pkin(4) - t10;
t5 = t26 * Ifges(6,4) + t27 * Ifges(6,2) + t58 * Ifges(6,6);
t4 = -qJD(5) * t20 - t133 * t40 + t136 * t45;
t3 = qJD(5) * t19 + t133 * t45 + t136 * t40;
t6 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t131 - 0.2e1 * mrSges(3,2) * t130 + m(3) * (t130 ^ 2 + t131 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t252 * (-Ifges(6,5) * t150 - Ifges(6,6) * t151) / 0.2e1 - (-Ifges(5,4) * t60 + t67 * mrSges(5,1) - Ifges(5,2) * t59 - t11 * mrSges(5,3) + t181 / 0.2e1 + Ifges(6,3) * t238 + Ifges(6,6) * t239 + Ifges(6,5) * t240 - t247 * Ifges(5,6) + t244) * t146 + m(6) * (t1 * t20 + t12 * t4 + t13 * t3 + t19 * t2) + m(5) * (t100 * t67 + t11 * t52 + t179 * t86 + t34 * t40) + (-mrSges(2,1) * t138 - t76 * mrSges(6,1) + mrSges(2,2) * t135 - t75 * mrSges(6,2) - t265 * (t121 * t116 - t119 * t132 + t126) + t261 * t126 + t243 * t119 + (-m(6) * t164 - t245) * t121) * g(2) + (mrSges(2,1) * t135 - t74 * mrSges(6,1) + mrSges(2,2) * t138 - t73 * mrSges(6,2) - t261 * t228 - t265 * (-t121 * t132 - t228) + t243 * t121 + (-m(6) * (-t116 - t164) + m(5) * t116 + t245) * t119) * g(1) + t137 * (Ifges(4,4) * t98 + Ifges(4,2) * t97 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t114 * (-mrSges(4,1) * t97 + mrSges(4,2) * t98) + t40 * t77 + (t137 * (-Ifges(4,2) * t134 + t214) + t147) * t184 / 0.2e1 + t68 * (-Ifges(6,4) * t150 - Ifges(6,2) * t151) / 0.2e1 + (m(4) * t114 + t107) * t102 + t52 * t53 + t3 * t42 + t4 * t43 + t19 * t16 + t20 * t17 + t97 * t157 / 0.2e1 + t134 * t98 * Ifges(4,1) + t30 * (mrSges(6,1) * t151 - mrSges(6,2) * t150) + t92 * t188 / 0.2e1 - t91 * t189 / 0.2e1 + (-Ifges(6,1) * t150 - Ifges(6,4) * t151) * t235 + (Ifges(5,5) * t90 / 0.2e1 - Ifges(5,6) * t89 / 0.2e1 + t148 + t155 * qJD(3) / 0.2e1) * qJD(3) + t202 * t242 + (Ifges(4,5) * t134 + Ifges(4,6) * t137) * t200 + t100 * t173 + t134 * (Ifges(4,4) * t97 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t98 * t214 / 0.2e1 - t5 * t206 / 0.2e1 + (-t1 * t206 + t12 * t150 - t13 * t151 - t2 * t202) * mrSges(6,3) + t55 * t179 + (t67 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t60 + Ifges(5,4) * t59 + Ifges(5,5) * t247 + t154 * t238 + t156 * t239 + t158 * t240 + t8 * t159 + t23 * t172) * t95 + (t137 * t84 - t134 * t85 - t104 * t188 - t106 * t189 + m(4) * ((-t80 * t134 - t79 * t137) * qJD(3) + t249)) * t112 + (-t188 * t79 - t189 * t80 + t249) * mrSges(4,3) - t151 * t22 / 0.2e1 + (-m(5) * t33 + m(6) * t30 - t253) * (t129 * t72 - t142 * t199) + (t86 * mrSges(5,1) - t254 / 0.2e1 + t21 / 0.2e1 - t47 / 0.2e1 + t258 - t257 - t34 * mrSges(5,3) + t256 / 0.2e1 + t255 / 0.2e1 - Ifges(5,4) * t232 + Ifges(6,5) * t235) * t89 + (t86 * mrSges(5,2) + t82 / 0.2e1 + t48 / 0.2e1 + Ifges(5,1) * t232 - t220 + t175) * t90 + (-m(5) * t10 + m(6) * t8 + t259) * (t129 * t93 + t168 * t199); m(3) * qJDD(2) + t134 * t84 + t137 * t85 - t259 * t146 - t253 * t89 + (-t104 * t134 + t106 * t137) * qJD(3) + (-t133 * t43 + t136 * t42 + t77) * t90 + (t53 + (-t133 * t42 - t136 * t43) * qJD(5) + t250) * t95 + m(6) * (t141 * t95 - t146 * t8 + t152 * t90 + t30 * t89) + m(4) * (t134 * t49 + t137 * t50 + (-t134 * t79 + t137 * t80) * qJD(3)) + m(5) * (t10 * t146 + t11 * t95 - t33 * t89 + t34 * t90) + (-t265 + t261) * g(3); (-m(5) * t125 - m(6) * (t125 + t225) - t144 * t120 + t262) * g(3) + (t154 * t252 + t156 * t68 + t158 * t69) * qJD(5) / 0.2e1 + (t104 + t178) * t80 + t113 * t7 + Ifges(4,6) * t97 + Ifges(4,5) * t98 - qJD(3) * (Ifges(5,5) * t87 + Ifges(5,6) * t88) / 0.2e1 - t86 * (-mrSges(5,1) * t88 + mrSges(5,2) * t87) - t37 * t77 + (-t106 + t177) * t79 - (-Ifges(4,2) * t191 + t117 + t92) * t185 / 0.2e1 + (Ifges(5,1) * t87 + t21 + t229) * t88 / 0.2e1 - (Ifges(5,2) * t88 + t48 + t82) * t87 / 0.2e1 + (t208 / 0.2e1 + t172) * t22 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (t149 + t175) * qJD(5) + Ifges(5,5) * t60 + Ifges(5,6) * t59 - t15 * t42 - t14 * t43 - t49 * mrSges(4,2) + t50 * mrSges(4,1) + t10 * mrSges(5,1) - t11 * mrSges(5,2) - t87 * t149 + (-t148 - t147 * qJD(1) / 0.2e1) * qJD(1) + t54 * t174 + t91 * t191 / 0.2e1 - t155 * t184 / 0.2e1 - t55 * t180 + t5 * t230 + t47 * t232 + (-Ifges(6,3) * t88 + t154 * t87) * t234 + (-Ifges(6,5) * t88 + t158 * t87) * t236 + (-Ifges(6,6) * t88 + t156 * t87) * t237 + (Ifges(6,5) * t133 + Ifges(6,6) * t136) * t238 + (Ifges(6,2) * t136 + t213) * t239 + (Ifges(6,1) * t133 + t212) * t240 + t133 * t242 + t87 * t220 + t53 * t227 - t34 * t218 - t23 * t203 / 0.2e1 + (t113 * t8 - t12 * t14 - t13 * t15 - t30 * t36) * m(6) + ((-t187 + t208) * t13 + (-t186 + t203) * t12 + t248) * mrSges(6,3) + (m(6) * t141 - t186 * t43 - t187 * t42 + t250) * (pkin(7) + t227) + t8 * t160 + (g(1) * t121 + g(2) * t119) * (t162 + t265 * pkin(3) * t134 + (-m(6) * pkin(7) + t264) * t120 + (mrSges(5,1) + t144) * t118) + t253 * t36 - t88 * t257 + t88 * t258 + ((t10 * t199 + t11 * t129) * pkin(3) - t180 * t86 + t33 * t36 - t34 * t37) * m(5); -t87 * t77 - t253 * t88 + (t252 * t42 + t16) * t136 + (-t252 * t43 + t17) * t133 + t173 + (t1 * t133 + t2 * t136 + t152 * t252 + t30 * t88 + t163) * m(6) + (-t33 * t88 - t34 * t87 + t163 + t67) * m(5); -t30 * (mrSges(6,1) * t69 + mrSges(6,2) * t68) + (Ifges(6,1) * t68 - t219) * t236 + t22 * t235 + (Ifges(6,5) * t68 - Ifges(6,6) * t69) * t234 - t12 * t42 + t13 * t43 - g(1) * (mrSges(6,1) * t75 - mrSges(6,2) * t76) - g(2) * (-mrSges(6,1) * t73 + mrSges(6,2) * t74) + g(3) * t159 * t118 + (t12 * t68 + t13 * t69) * mrSges(6,3) + t181 + (-Ifges(6,2) * t69 + t23 + t65) * t237 + t244;];
tau = t6;
