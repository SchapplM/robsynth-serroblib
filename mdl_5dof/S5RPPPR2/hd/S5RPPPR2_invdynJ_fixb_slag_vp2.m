% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:20
% EndTime: 2020-01-03 11:22:39
% DurationCPUTime: 8.32s
% Computational Cost: add. (2714->374), mult. (6839->536), div. (0->0), fcn. (5079->10), ass. (0->173)
t196 = qJD(1) * qJD(2);
t136 = qJDD(1) * qJ(2) + t196;
t155 = sin(pkin(7));
t158 = cos(pkin(7));
t203 = t155 ^ 2 + t158 ^ 2;
t154 = sin(pkin(8));
t201 = qJD(1) * t155;
t187 = t154 * t201;
t153 = sin(pkin(9));
t156 = cos(pkin(9));
t200 = qJD(1) * t158;
t123 = -pkin(2) * t158 - qJ(3) * t155 - pkin(1);
t103 = qJD(1) * t123 + qJD(2);
t157 = cos(pkin(8));
t188 = qJ(2) * t200;
t65 = t154 * t103 + t157 * t188;
t51 = -qJ(4) * t200 + t65;
t130 = qJ(2) * t201 + qJD(3);
t173 = pkin(3) * t154 - qJ(4) * t157;
t77 = t173 * t201 + t130;
t25 = -t153 * t51 + t156 * t77;
t20 = -pkin(4) * t187 - t25;
t243 = m(6) * t20;
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t193 = qJDD(1) * t155;
t181 = t154 * t193;
t202 = qJD(1) * t154;
t185 = t161 * t202;
t182 = t156 * t201;
t183 = t153 * t200;
t95 = t157 * t182 - t183;
t53 = t155 * t185 - t159 * t95;
t212 = t157 * t155;
t105 = -t153 * t158 + t156 * t212;
t92 = t105 * qJDD(1);
t27 = qJD(5) * t53 + t159 * t181 + t161 * t92;
t242 = Ifges(6,5) * t27;
t186 = t159 * t202;
t55 = t155 * t186 + t161 * t95;
t28 = -qJD(5) * t55 - t159 * t92 + t161 * t181;
t241 = Ifges(6,6) * t28;
t104 = t153 * t212 + t156 * t158;
t91 = t104 * qJDD(1);
t87 = qJDD(5) + t91;
t240 = Ifges(6,3) * t87;
t239 = mrSges(2,2) - mrSges(3,3);
t238 = -mrSges(5,3) + mrSges(4,2);
t223 = -mrSges(5,1) * t187 - mrSges(6,1) * t53 + mrSges(6,2) * t55 + mrSges(5,3) * t95;
t174 = mrSges(3,1) * t158 - mrSges(3,2) * t155;
t237 = -mrSges(2,1) - t174;
t170 = -mrSges(4,1) * t158 - mrSges(4,3) * t212;
t93 = t104 * qJD(1);
t221 = -mrSges(5,1) * t93 - mrSges(5,2) * t95 + t170 * qJD(1);
t160 = sin(qJ(1));
t207 = t160 * t157;
t162 = cos(qJ(1));
t208 = t158 * t162;
t114 = t154 * t208 - t207;
t115 = t154 * t160 + t157 * t208;
t213 = t155 * t162;
t72 = t115 * t156 + t153 * t213;
t235 = -t114 * t161 + t159 * t72;
t234 = t114 * t159 + t161 * t72;
t233 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t195 = qJD(1) * qJD(4);
t211 = t157 * t158;
t198 = qJD(3) * t155;
t81 = -qJD(1) * t198 + qJDD(1) * t123 + qJDD(2);
t47 = t136 * t211 + t154 * t81;
t37 = (-qJ(4) * qJDD(1) - t195) * t158 + t47;
t121 = t155 * t136 + qJDD(3);
t48 = (qJDD(1) * t173 - t157 * t195) * t155 + t121;
t13 = t153 * t48 + t156 * t37;
t11 = pkin(6) * t181 + t13;
t192 = qJDD(1) * t158;
t210 = t158 * t154;
t46 = -t136 * t210 + t157 * t81;
t40 = pkin(3) * t192 + qJDD(4) - t46;
t18 = pkin(4) * t91 - pkin(6) * t92 + t40;
t64 = t103 * t157 - t154 * t188;
t50 = pkin(3) * t200 + qJD(4) - t64;
t19 = pkin(4) * t93 - pkin(6) * t95 + t50;
t26 = t153 * t77 + t156 * t51;
t21 = pkin(6) * t187 + t26;
t5 = -t159 * t21 + t161 * t19;
t1 = qJD(5) * t5 + t11 * t161 + t159 * t18;
t6 = t159 * t19 + t161 * t21;
t2 = -qJD(5) * t6 - t11 * t159 + t161 * t18;
t232 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t216 = t154 * t159;
t168 = t156 * t216 + t157 * t161;
t96 = (t153 * t155 + t156 * t211) * qJD(1);
t231 = -t168 * qJD(5) - t158 * t186 - t161 * t96;
t215 = t154 * t161;
t113 = t156 * t215 - t157 * t159;
t230 = -t113 * qJD(5) - t158 * t185 + t159 * t96;
t229 = qJD(1) ^ 2 * t203;
t228 = qJD(5) + t93;
t226 = t55 / 0.2e1;
t60 = mrSges(5,1) * t181 - mrSges(5,3) * t92;
t7 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t225 = t60 - t7;
t224 = Ifges(6,4) * t55;
t84 = qJ(2) * t211 + t154 * t123;
t74 = -qJ(4) * t158 + t84;
t89 = (qJ(2) + t173) * t155;
t36 = t153 * t89 + t156 * t74;
t218 = t153 * t154;
t217 = t154 * t155;
t214 = t155 * t160;
t209 = t158 * t160;
t180 = t157 * t193;
t205 = mrSges(4,1) * t181 + mrSges(4,2) * t180;
t204 = t162 * pkin(1) + t160 * qJ(2);
t199 = qJD(2) * t158;
t197 = m(4) + m(5) + m(6);
t191 = t240 + t241 + t242;
t42 = t91 * mrSges(5,1) + t92 * mrSges(5,2);
t179 = t160 * pkin(1) - qJ(2) * t162;
t83 = -qJ(2) * t210 + t123 * t157;
t178 = pkin(2) * t208 + qJ(3) * t213 + t204;
t66 = -t105 * t159 + t155 * t215;
t67 = t105 * t161 + t155 * t216;
t177 = mrSges(6,1) * t66 - mrSges(6,2) * t67;
t76 = t158 * pkin(3) - t83;
t175 = -mrSges(3,1) * t192 + mrSges(3,2) * t193;
t107 = -t154 * t198 + t157 * t199;
t12 = -t153 * t37 + t156 * t48;
t35 = -t153 * t74 + t156 * t89;
t29 = pkin(4) * t104 - pkin(6) * t105 + t76;
t31 = pkin(6) * t217 + t36;
t8 = -t159 * t31 + t161 * t29;
t9 = t159 * t29 + t161 * t31;
t172 = (-qJD(4) * t157 + qJD(2)) * t155;
t171 = pkin(2) * t209 + qJ(3) * t214 + t179;
t169 = mrSges(4,2) * t158 - mrSges(4,3) * t217;
t167 = (mrSges(4,1) * t154 + mrSges(4,2) * t157) * t155;
t166 = t115 * pkin(3) + qJ(4) * t114 + t178;
t111 = t154 * t209 + t157 * t162;
t165 = -g(1) * t217 - g(2) * t111 + g(3) * t114;
t112 = -t154 * t162 + t158 * t207;
t164 = t112 * pkin(3) + qJ(4) * t111 + t171;
t146 = -qJDD(1) * pkin(1) + qJDD(2);
t116 = t169 * qJD(1);
t109 = t170 * qJDD(1);
t108 = t169 * qJDD(1);
t106 = t154 * t199 + t157 * t198;
t99 = qJD(1) * t167;
t94 = t157 * t183 - t182;
t85 = -qJD(4) * t158 + t107;
t79 = t113 * t201;
t78 = t168 * t201;
t70 = t112 * t156 + t153 * t214;
t62 = -mrSges(5,2) * t187 - mrSges(5,3) * t93;
t59 = -mrSges(5,2) * t181 - mrSges(5,3) * t91;
t58 = t67 * qJD(5);
t57 = t66 * qJD(5);
t52 = Ifges(6,4) * t53;
t44 = t153 * t172 + t156 * t85;
t39 = t111 * t159 + t161 * t70;
t38 = t111 * t161 - t159 * t70;
t33 = mrSges(6,1) * t228 - mrSges(6,3) * t55;
t32 = -mrSges(6,2) * t228 + mrSges(6,3) * t53;
t30 = -pkin(4) * t217 - t35;
t17 = Ifges(6,1) * t55 + Ifges(6,5) * t228 + t52;
t16 = Ifges(6,2) * t53 + Ifges(6,6) * t228 + t224;
t15 = -mrSges(6,2) * t87 + mrSges(6,3) * t28;
t14 = mrSges(6,1) * t87 - mrSges(6,3) * t27;
t10 = -pkin(4) * t181 - t12;
t4 = -qJD(5) * t9 + t106 * t161 - t159 * t44;
t3 = qJD(5) * t8 + t106 * t159 + t161 * t44;
t22 = [0.2e1 * t203 * t136 * mrSges(3,3) + (-mrSges(6,3) * t2 + Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t87) * t67 + (-Ifges(5,6) * t104 + Ifges(5,5) * t105 + Ifges(5,3) * t217 + Ifges(4,6) * t158 - (Ifges(4,4) * t157 - Ifges(4,2) * t154) * t155) * t181 + (mrSges(5,1) * t12 - mrSges(5,2) * t13 + Ifges(5,5) * t92 - t91 * Ifges(5,6)) * t217 + t107 * t116 + t84 * t108 + t83 * t109 + t76 * t42 + t57 * t17 / 0.2e1 - t58 * t16 / 0.2e1 + t36 * t59 + t35 * t60 + t44 * t62 + t3 * t32 + t4 * t33 + t30 * t7 + t8 * t14 + t9 * t15 + m(6) * (t1 * t9 + t10 * t30 + t2 * t8 + t3 * t6 + t4 * t5) + (Ifges(3,4) * t193 - Ifges(4,5) * t180 + (Ifges(3,2) + Ifges(4,3)) * t192) * t158 + (-m(4) * t64 + m(5) * t50 - t221) * t106 + t121 * t167 + (qJ(2) * t205 + m(4) * (qJ(2) * t121 + qJD(2) * t130) + qJD(2) * t99 + Ifges(3,1) * t193 + (Ifges(4,1) * t157 - Ifges(4,4) * t154) * t180 + (-Ifges(4,5) * t157 + Ifges(4,6) * t154 + Ifges(3,4)) * t192) * t155 + (-t5 * t57 - t58 * t6) * mrSges(6,3) + m(5) * (t12 * t35 + t13 * t36 + t26 * t44 + t40 * t76) + (Ifges(6,1) * t57 - Ifges(6,4) * t58) * t226 + t53 * (Ifges(6,4) * t57 - Ifges(6,2) * t58) / 0.2e1 + t20 * (mrSges(6,1) * t58 + mrSges(6,2) * t57) + (t240 / 0.2e1 + t242 / 0.2e1 + t241 / 0.2e1 + t40 * mrSges(5,1) - Ifges(5,4) * t92 + Ifges(5,2) * t91 + t191 / 0.2e1 - t13 * mrSges(5,3) + t232) * t104 + (-m(5) * t166 - t72 * mrSges(5,1) - m(3) * t204 - m(4) * t178 - t115 * mrSges(4,1) - mrSges(4,3) * t213 - m(6) * (pkin(4) * t72 + t166) - t234 * mrSges(6,1) + t235 * mrSges(6,2) + t233 * (t115 * t153 - t156 * t213) + t237 * t162 + t239 * t160 + t238 * t114) * g(2) + (-m(6) * (pkin(4) * t70 + t164) - t39 * mrSges(6,1) - t38 * mrSges(6,2) - m(5) * t164 - t70 * mrSges(5,1) - m(3) * t179 - m(4) * t171 - t112 * mrSges(4,1) - mrSges(4,3) * t214 + t233 * (t112 * t153 - t156 * t214) - t239 * t162 + t237 * t160 + t238 * t111) * g(3) + (mrSges(6,3) * t1 + Ifges(6,4) * t27 + Ifges(6,2) * t28 + Ifges(6,6) * t87) * t66 + t228 * (Ifges(6,5) * t57 - Ifges(6,6) * t58) / 0.2e1 + (-m(5) * t25 + t223 + t243) * (t153 * t85 - t156 * t172) + m(4) * (t107 * t65 + t46 * t83 + t47 * t84) + t47 * t169 + t46 * t170 - t146 * t174 - pkin(1) * t175 - t10 * t177 + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t146 + (t136 + t196) * qJ(2) * t203) + (mrSges(5,2) * t40 - mrSges(5,3) * t12 + Ifges(5,1) * t92 - Ifges(5,4) * t91) * t105; -t99 * t201 + t175 - t168 * t14 + t113 * t15 - t96 * t62 + (-t225 * t153 + t156 * t59 + t221 * t200 + t108) * t154 - t223 * t94 + (-t116 * t200 + t109 - t42) * t157 + t230 * t33 + t231 * t32 - mrSges(3,3) * t229 + (g(2) * t162 + g(3) * t160) * (m(3) + t197) + (t1 * t113 + t10 * t218 - t168 * t2 - t20 * t94 + t230 * t5 + t231 * t6) * m(6) + (-(t130 * t155 - t64 * t210 + t65 * t211) * qJD(1) + t47 * t154 + t46 * t157) * m(4) + (-qJ(2) * t229 + t146) * m(3) + (t25 * t94 - t26 * t96 - t50 * t210 * qJD(1) - t157 * t40 + (-t12 * t153 + t13 * t156) * t154) * m(5); t79 * t32 - t78 * t33 + t225 * t156 + t197 * t158 * g(1) + (-t159 * t14 + t161 * t15 + t59 + (-t159 * t32 - t161 * t33) * qJD(5)) * t153 + m(5) * (t12 * t156 + t13 * t153) + m(4) * t121 + ((t221 * t157 + (t153 * t223 + t156 * t62 + t116) * t154 - m(5) * (-t154 * t156 * t26 + t157 * t50 + t218 * t25) - m(4) * (-t154 * t65 - t157 * t64) + t218 * t243) * qJD(1) + (-g(2) * t160 + g(3) * t162) * t197) * t155 + t205 + (-t10 * t156 + (t1 * t161 - t159 * t2 + (-t159 * t6 - t161 * t5) * qJD(5)) * t153 - t5 * t78 + t6 * t79) * m(6); t93 * t62 - t223 * t95 + (t228 * t32 + t14) * t161 + (-t228 * t33 + t15) * t159 + t42 + (t1 * t159 + t2 * t161 - t20 * t95 + t165 + t228 * (-t159 * t5 + t161 * t6)) * m(6) + (t25 * t95 + t26 * t93 + t165 + t40) * m(5); -t20 * (mrSges(6,1) * t55 + mrSges(6,2) * t53) - t55 * (Ifges(6,1) * t53 - t224) / 0.2e1 + t16 * t226 - t228 * (Ifges(6,5) * t53 - Ifges(6,6) * t55) / 0.2e1 - t5 * t32 + t6 * t33 - g(1) * t177 - g(2) * (mrSges(6,1) * t38 - mrSges(6,2) * t39) - g(3) * (mrSges(6,1) * t235 + mrSges(6,2) * t234) + (t5 * t53 + t55 * t6) * mrSges(6,3) + t191 - (-Ifges(6,2) * t55 + t17 + t52) * t53 / 0.2e1 + t232;];
tau = t22;
