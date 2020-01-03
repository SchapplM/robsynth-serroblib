% Calculate vector of inverse dynamics joint torques for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:11
% DurationCPUTime: 7.06s
% Computational Cost: add. (2584->374), mult. (6200->505), div. (0->0), fcn. (4412->10), ass. (0->171)
t182 = qJDD(1) * qJ(2);
t183 = qJD(1) * qJD(2);
t116 = t182 + t183;
t148 = sin(pkin(8));
t144 = t148 ^ 2;
t236 = t116 * t144;
t235 = t116 + t183;
t184 = m(4) + m(5) + m(6);
t234 = m(3) + t184;
t233 = Ifges(3,4) - Ifges(4,5);
t149 = cos(pkin(8));
t135 = t148 * qJ(3);
t174 = pkin(1) + t135;
t232 = -t149 * pkin(2) - t174;
t145 = t149 ^ 2;
t231 = t235 * t145;
t204 = t148 * mrSges(4,3);
t168 = t149 * mrSges(4,1) + t204;
t169 = mrSges(3,1) * t149 - mrSges(3,2) * t148;
t230 = -mrSges(2,1) - t168 - t169;
t192 = qJD(1) * t149;
t193 = qJD(1) * t148;
t92 = -qJD(1) * pkin(1) - pkin(2) * t192 - qJ(3) * t193 + qJD(2);
t70 = pkin(3) * t192 - t92;
t155 = cos(qJ(4));
t190 = qJD(1) * t155;
t152 = sin(qJ(4));
t191 = qJD(1) * t152;
t95 = -t148 * t191 - t149 * t190;
t96 = t148 * t190 - t149 * t191;
t229 = m(5) * t70 - mrSges(5,1) * t95 + mrSges(5,2) * t96 + t168 * qJD(1);
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t228 = -g(1) * t156 - g(2) * t153;
t227 = m(5) * pkin(6) - m(6) * (-pkin(7) - pkin(6)) + mrSges(2,2) + mrSges(5,3) + mrSges(6,3);
t226 = t145 * t182 + t228 + t231 + t236;
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t173 = -t151 * t96 + t154 * t95;
t209 = -pkin(6) + qJ(2);
t113 = t209 * t148;
t114 = t209 * t149;
t63 = t152 * t113 + t155 * t114;
t146 = qJD(4) + qJD(5);
t225 = t149 * (pkin(2) + pkin(3)) + t174;
t224 = m(6) * pkin(4);
t223 = -t173 / 0.2e1;
t166 = t151 * t95 + t154 * t96;
t222 = -t166 / 0.2e1;
t221 = t166 / 0.2e1;
t219 = t96 / 0.2e1;
t217 = pkin(4) * t96;
t216 = -t146 / 0.2e1;
t215 = Ifges(5,4) * t96;
t214 = Ifges(6,4) * t166;
t107 = qJD(1) * t114;
t115 = qJ(2) * t193 + qJD(3);
t99 = -pkin(6) * t193 + t115;
t55 = t107 * t155 + t152 * t99;
t32 = pkin(7) * t95 + t55;
t203 = t151 * t32;
t54 = -t107 * t152 + t155 * t99;
t31 = -pkin(7) * t96 + t54;
t30 = qJD(4) * pkin(4) + t31;
t14 = t154 * t30 - t203;
t213 = t14 * mrSges(6,3);
t201 = t154 * t32;
t15 = t151 * t30 + t201;
t212 = t15 * mrSges(6,3);
t211 = t54 * mrSges(5,3);
t210 = t55 * mrSges(5,3);
t147 = qJ(4) + qJ(5);
t136 = sin(t147);
t137 = cos(t147);
t164 = t136 * t149 - t137 * t148;
t73 = t164 * t153;
t163 = t136 * t148 + t137 * t149;
t74 = t163 * t153;
t208 = -t73 * mrSges(6,1) - t74 * mrSges(6,2);
t75 = t164 * t156;
t76 = t163 * t156;
t207 = -t75 * mrSges(6,1) - t76 * mrSges(6,2);
t206 = -t163 * mrSges(6,1) + t164 * mrSges(6,2);
t158 = qJD(1) ^ 2;
t200 = qJ(2) * t158;
t198 = t148 * t152;
t197 = t149 * t156;
t196 = t231 * qJ(2);
t195 = t156 * pkin(1) + t153 * qJ(2);
t189 = qJD(2) * t152;
t188 = qJD(2) * t155;
t187 = qJD(3) * t148;
t186 = qJD(4) * t152;
t185 = qJD(4) * t155;
t181 = qJDD(1) * t148;
t100 = t148 * t116 + qJDD(3);
t62 = t155 * t113 - t114 * t152;
t97 = t149 * pkin(3) - t232;
t105 = t148 * t155 - t149 * t152;
t162 = t149 * t155 + t198;
t93 = t162 * qJD(4);
t52 = -qJD(1) * t93 + qJDD(1) * t105;
t94 = t105 * qJD(4);
t53 = -qJD(1) * t94 - qJDD(1) * t162;
t171 = -t53 * mrSges(5,1) + t52 * mrSges(5,2);
t12 = qJD(5) * t173 + t151 * t53 + t154 * t52;
t13 = -qJD(5) * t166 - t151 * t52 + t154 * t53;
t170 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t167 = -mrSges(5,1) * t162 - mrSges(5,2) * t105;
t37 = -pkin(7) * t105 + t62;
t38 = -pkin(7) * t162 + t63;
t18 = -t151 * t38 + t154 * t37;
t19 = t151 * t37 + t154 * t38;
t165 = -qJD(1) * t187 + qJDD(2);
t56 = -t105 * t151 - t154 * t162;
t57 = t105 * t154 - t151 * t162;
t109 = t151 * t155 + t152 * t154;
t108 = -t151 * t152 + t154 * t155;
t143 = qJDD(4) + qJDD(5);
t80 = -pkin(6) * t181 + t100;
t89 = (-pkin(6) * qJDD(1) + t116) * t149;
t23 = -t107 * t186 + t152 * t80 + t155 * t89 + t99 * t185;
t11 = pkin(7) * t53 + t23;
t24 = -qJD(4) * t55 - t152 * t89 + t155 * t80;
t8 = qJDD(4) * pkin(4) - pkin(7) * t52 + t24;
t2 = qJD(5) * t14 + t11 * t154 + t151 * t8;
t3 = -qJD(5) * t15 - t11 * t151 + t154 * t8;
t161 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t12 + Ifges(6,6) * t13 + Ifges(6,3) * t143;
t35 = t113 * t185 - t114 * t186 + t148 * t189 + t149 * t188;
t159 = t105 * t224;
t36 = -qJD(4) * t63 + t148 * t188 - t149 * t189;
t58 = qJDD(1) * t225 - t165;
t134 = -qJDD(1) * pkin(1) + qJDD(2);
t132 = pkin(4) * t155 + pkin(3);
t129 = t145 * t200;
t126 = mrSges(3,2) * t181;
t90 = Ifges(5,4) * t95;
t88 = t162 * t156;
t87 = t105 * t156;
t86 = t162 * t153;
t85 = t105 * t153;
t78 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t96;
t77 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t95;
t69 = pkin(4) * t94 + t187;
t64 = qJDD(1) * t232 + t165;
t61 = t146 * t109;
t60 = t146 * t108;
t59 = pkin(4) * t162 + t97;
t44 = t96 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t90;
t43 = t95 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t215;
t42 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t53;
t41 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t52;
t40 = -pkin(4) * t95 + t70;
t39 = Ifges(6,4) * t173;
t34 = mrSges(6,1) * t146 - mrSges(6,3) * t166;
t33 = -mrSges(6,2) * t146 + mrSges(6,3) * t173;
t29 = pkin(7) * t93 + t36;
t28 = -pkin(7) * t94 + t35;
t27 = -pkin(4) * t53 + t58;
t26 = -qJD(5) * t57 + t151 * t93 - t154 * t94;
t25 = qJD(5) * t56 - t151 * t94 - t154 * t93;
t22 = -mrSges(6,1) * t173 + mrSges(6,2) * t166;
t21 = Ifges(6,1) * t166 + Ifges(6,5) * t146 + t39;
t20 = Ifges(6,2) * t173 + Ifges(6,6) * t146 + t214;
t17 = t154 * t31 - t203;
t16 = -t151 * t31 - t201;
t7 = -mrSges(6,2) * t143 + mrSges(6,3) * t13;
t6 = mrSges(6,1) * t143 - mrSges(6,3) * t12;
t5 = -qJD(5) * t19 - t151 * t28 + t154 * t29;
t4 = qJD(5) * t18 + t151 * t29 + t154 * t28;
t1 = [m(3) * (qJ(2) * t144 * t235 - pkin(1) * t134 + t196) + (t226 + t236) * mrSges(3,3) + m(6) * (t14 * t5 + t15 * t4 + t18 * t3 + t19 * t2 + t27 * t59 + t40 * t69) + (-mrSges(5,3) * t23 - Ifges(5,4) * t52 - Ifges(5,2) * t53 - Ifges(5,6) * qJDD(4)) * t162 - t25 * t213 + t93 * t211 + t59 * t170 + t97 * t171 - t58 * t167 - t64 * t168 - t134 * t169 + (-mrSges(5,3) * t24 + Ifges(5,1) * t52 + Ifges(5,4) * t53 + Ifges(5,5) * qJDD(4)) * t105 + t146 * (Ifges(6,5) * t25 + Ifges(6,6) * t26) / 0.2e1 - t94 * t43 / 0.2e1 - t93 * t44 / 0.2e1 + t35 * t77 + t36 * t78 + t62 * t41 + t63 * t42 + t69 * t22 + t40 * (-mrSges(6,1) * t26 + mrSges(6,2) * t25) + t4 * t33 + t5 * t34 + t25 * t21 / 0.2e1 + t26 * t20 / 0.2e1 + t18 * t6 + t19 * t7 - pkin(1) * t126 + (-m(5) * pkin(3) * t197 - m(3) * t195 - t88 * mrSges(5,1) - t76 * mrSges(6,1) - t87 * mrSges(5,2) + t75 * mrSges(6,2) - t184 * (pkin(2) * t197 + t156 * t135 + t195) + (-m(6) * (pkin(4) * t198 + t132 * t149) + t230) * t156 + t227 * t153) * g(2) + t229 * t187 + (t100 * t148 + t226) * mrSges(4,2) + m(5) * (t23 * t63 + t24 * t62 + t35 * t55 + t36 * t54 + t58 * t97) + (-mrSges(6,1) * t27 + mrSges(6,3) * t2 + Ifges(6,4) * t12 + Ifges(6,2) * t13 + Ifges(6,6) * t143) * t56 + (Ifges(6,1) * t25 + Ifges(6,4) * t26) * t221 - t94 * t210 + t26 * t212 + (mrSges(6,2) * t27 - mrSges(6,3) * t3 + Ifges(6,1) * t12 + Ifges(6,4) * t13 + Ifges(6,5) * t143) * t57 + m(4) * (t232 * t64 + (qJ(2) * t100 + qJD(2) * t115 - qJD(3) * t92) * t148 + t196) + (t86 * mrSges(5,1) + t74 * mrSges(6,1) + t85 * mrSges(5,2) - t73 * mrSges(6,2) + (m(3) * pkin(1) - m(4) * t232 + t225 * m(5) - m(6) * (-pkin(1) + (-pkin(2) - t132) * t149 + (-pkin(4) * t152 - qJ(3)) * t148) - t230) * t153 + (-qJ(2) * t234 + t227) * t156) * g(1) + t173 * (Ifges(6,4) * t25 + Ifges(6,2) * t26) / 0.2e1 + qJD(4) * (-Ifges(5,5) * t93 - Ifges(5,6) * t94) / 0.2e1 + t95 * (-Ifges(5,4) * t93 - Ifges(5,2) * t94) / 0.2e1 + t70 * (mrSges(5,1) * t94 - mrSges(5,2) * t93) + (-Ifges(5,1) * t93 - Ifges(5,4) * t94) * t219 + (t233 * t149 + (Ifges(3,1) + Ifges(4,1)) * t148) * t181 + ((pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t149 + t233 * t148) * t149 - t168 * t232 + Ifges(2,3)) * qJDD(1); t173 * t33 - t166 * t34 + t95 * t77 - t96 * t78 + t126 + (-t204 + (-mrSges(3,1) - mrSges(4,1)) * t149) * qJDD(1) - t170 - t171 + (mrSges(4,2) + mrSges(3,3)) * t158 * (-t144 - t145) + (-g(1) * t153 + g(2) * t156) * t234 + (-t14 * t166 + t15 * t173 - t27) * m(6) + (-t54 * t96 + t55 * t95 - t58) * m(5) + (-t115 * t193 - t129 + t64) * m(4) + (-t144 * t200 - t129 + t134) * m(3); t108 * t6 + t109 * t7 + t152 * t42 + t155 * t41 + t60 * t33 - t61 * t34 + (-t152 * t78 + t155 * t77) * qJD(4) + t184 * t149 * g(3) + m(5) * (t152 * t23 + t155 * t24 + (-t152 * t54 + t155 * t55) * qJD(4)) + m(4) * t100 + m(6) * (t108 * t3 + t109 * t2 - t14 * t61 + t15 * t60) + (qJDD(1) * mrSges(4,2) + (m(4) * t92 - m(6) * t40 - t22 - t229) * qJD(1) + t228 * t184) * t148; -t22 * t217 - m(6) * (t14 * t16 + t15 * t17 + t217 * t40) - t96 * (Ifges(5,1) * t95 - t215) / 0.2e1 - t70 * (mrSges(5,1) * t96 + mrSges(5,2) * t95) - qJD(4) * (Ifges(5,5) * t95 - Ifges(5,6) * t96) / 0.2e1 - t54 * t77 + t55 * t78 + Ifges(5,6) * t53 + Ifges(5,5) * t52 - t17 * t33 - t16 * t34 - t23 * mrSges(5,2) + t24 * mrSges(5,1) + (t151 * t7 + t154 * t6 + (-t151 * t34 + t154 * t33) * qJD(5)) * pkin(4) + t161 + (t151 * t2 + t154 * t3 + (-t14 * t151 + t15 * t154) * qJD(5)) * t224 + (t162 * t224 - t167 - t206) * g(3) + t43 * t219 + t96 * t210 + t95 * t211 + (-mrSges(5,1) * t87 + mrSges(5,2) * t88 - t156 * t159 - t207) * g(1) + Ifges(5,3) * qJDD(4) - (-t212 - t20 / 0.2e1 + t40 * mrSges(6,1) + Ifges(6,2) * t223 + Ifges(6,6) * t216 + Ifges(6,4) * t222) * t166 + (-t21 / 0.2e1 - t40 * mrSges(6,2) + Ifges(6,4) * t223 + t213 + Ifges(6,5) * t216 + Ifges(6,1) * t222) * t173 - (-Ifges(5,2) * t96 + t44 + t90) * t95 / 0.2e1 + (-mrSges(5,1) * t85 + mrSges(5,2) * t86 - t153 * t159 - t208) * g(2); -t40 * (mrSges(6,1) * t166 + mrSges(6,2) * t173) + (Ifges(6,1) * t173 - t214) * t222 + t20 * t221 + (Ifges(6,5) * t173 - Ifges(6,6) * t166) * t216 - t14 * t33 + t15 * t34 - g(1) * t207 - g(2) * t208 - g(3) * t206 + (t14 * t173 + t15 * t166) * mrSges(6,3) + t161 + (-Ifges(6,2) * t166 + t21 + t39) * t223;];
tau = t1;
