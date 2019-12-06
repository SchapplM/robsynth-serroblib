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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:30:58
% EndTime: 2019-12-05 17:31:18
% DurationCPUTime: 8.76s
% Computational Cost: add. (2714->371), mult. (6839->531), div. (0->0), fcn. (5079->10), ass. (0->171)
t154 = cos(qJ(1));
t239 = g(3) * t154;
t189 = qJD(1) * qJD(2);
t134 = qJDD(1) * qJ(2) + t189;
t147 = sin(pkin(7));
t150 = cos(pkin(7));
t196 = t147 ^ 2 + t150 ^ 2;
t146 = sin(pkin(8));
t194 = qJD(1) * t147;
t180 = t146 * t194;
t145 = sin(pkin(9));
t148 = cos(pkin(9));
t193 = qJD(1) * t150;
t123 = -pkin(2) * t150 - qJ(3) * t147 - pkin(1);
t103 = qJD(1) * t123 + qJD(2);
t149 = cos(pkin(8));
t181 = qJ(2) * t193;
t65 = t146 * t103 + t149 * t181;
t51 = -qJ(4) * t193 + t65;
t128 = qJ(2) * t194 + qJD(3);
t167 = pkin(3) * t146 - qJ(4) * t149;
t77 = t167 * t194 + t128;
t25 = -t145 * t51 + t148 * t77;
t20 = -pkin(4) * t180 - t25;
t238 = m(6) * t20;
t151 = sin(qJ(5));
t153 = cos(qJ(5));
t186 = qJDD(1) * t147;
t174 = t146 * t186;
t195 = qJD(1) * t146;
t178 = t153 * t195;
t175 = t148 * t194;
t176 = t145 * t193;
t95 = t149 * t175 - t176;
t53 = t147 * t178 - t151 * t95;
t204 = t147 * t149;
t105 = -t145 * t150 + t148 * t204;
t92 = t105 * qJDD(1);
t27 = qJD(5) * t53 + t151 * t174 + t153 * t92;
t237 = Ifges(6,5) * t27;
t179 = t151 * t195;
t55 = t147 * t179 + t153 * t95;
t28 = -qJD(5) * t55 - t151 * t92 + t153 * t174;
t236 = Ifges(6,6) * t28;
t104 = t145 * t204 + t148 * t150;
t91 = t104 * qJDD(1);
t87 = qJDD(5) + t91;
t235 = Ifges(6,3) * t87;
t234 = mrSges(5,3) - mrSges(4,2);
t216 = -mrSges(5,1) * t180 - mrSges(6,1) * t53 + mrSges(6,2) * t55 + mrSges(5,3) * t95;
t165 = -mrSges(4,1) * t150 - mrSges(4,3) * t204;
t93 = t104 * qJD(1);
t214 = -mrSges(5,1) * t93 - mrSges(5,2) * t95 + t165 * qJD(1);
t168 = -mrSges(3,1) * t150 + mrSges(3,2) * t147;
t233 = m(3) * pkin(1) - m(4) * t123 + t147 * mrSges(4,3) + mrSges(2,1) - t168;
t152 = sin(qJ(1));
t199 = t152 * t149;
t200 = t150 * t154;
t114 = t146 * t200 - t199;
t206 = t146 * t152;
t115 = t149 * t200 + t206;
t202 = t147 * t154;
t71 = t115 * t148 + t145 * t202;
t231 = -t114 * t153 + t151 * t71;
t230 = -t114 * t151 - t153 * t71;
t229 = -pkin(6) * m(6) + mrSges(5,2) - mrSges(6,3);
t188 = qJD(1) * qJD(4);
t201 = t149 * t150;
t191 = qJD(3) * t147;
t81 = -qJD(1) * t191 + qJDD(1) * t123 + qJDD(2);
t47 = t134 * t201 + t146 * t81;
t37 = (-qJ(4) * qJDD(1) - t188) * t150 + t47;
t121 = t147 * t134 + qJDD(3);
t48 = (qJDD(1) * t167 - t149 * t188) * t147 + t121;
t13 = t145 * t48 + t148 * t37;
t11 = pkin(6) * t174 + t13;
t185 = qJDD(1) * t150;
t208 = t146 * t150;
t46 = -t134 * t208 + t149 * t81;
t40 = pkin(3) * t185 + qJDD(4) - t46;
t18 = pkin(4) * t91 - pkin(6) * t92 + t40;
t64 = t103 * t149 - t146 * t181;
t50 = pkin(3) * t193 + qJD(4) - t64;
t19 = pkin(4) * t93 - pkin(6) * t95 + t50;
t26 = t145 * t77 + t148 * t51;
t21 = pkin(6) * t180 + t26;
t5 = -t151 * t21 + t153 * t19;
t1 = qJD(5) * t5 + t11 * t153 + t151 * t18;
t6 = t151 * t19 + t153 * t21;
t2 = -qJD(5) * t6 - t11 * t151 + t153 * t18;
t228 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t207 = t146 * t151;
t163 = t148 * t207 + t149 * t153;
t96 = (t145 * t147 + t148 * t201) * qJD(1);
t227 = -t163 * qJD(5) - t150 * t179 - t153 * t96;
t205 = t146 * t153;
t113 = t148 * t205 - t149 * t151;
t226 = -t113 * qJD(5) - t150 * t178 + t151 * t96;
t225 = -m(3) * qJ(2) - mrSges(3,3);
t224 = qJD(5) + t93;
t213 = qJ(2) * t152;
t223 = -t115 * pkin(3) - qJ(4) * t114 + t123 * t154 - t213;
t111 = t149 * t154 + t150 * t206;
t112 = t146 * t154 - t150 * t199;
t142 = t154 * qJ(2);
t222 = t112 * pkin(3) - qJ(4) * t111 + t123 * t152 + t142;
t220 = t55 / 0.2e1;
t60 = mrSges(5,1) * t174 - mrSges(5,3) * t92;
t7 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t219 = t60 - t7;
t217 = Ifges(6,4) * t55;
t84 = qJ(2) * t201 + t146 * t123;
t74 = -qJ(4) * t150 + t84;
t89 = (qJ(2) + t167) * t147;
t36 = t145 * t89 + t148 * t74;
t210 = t145 * t146;
t209 = t146 * t147;
t203 = t147 * t152;
t173 = t149 * t186;
t197 = mrSges(4,1) * t174 + mrSges(4,2) * t173;
t192 = qJD(2) * t150;
t190 = m(4) + m(5) + m(6);
t184 = t235 + t236 + t237;
t42 = t91 * mrSges(5,1) + t92 * mrSges(5,2);
t83 = -qJ(2) * t208 + t123 * t149;
t66 = -t105 * t151 + t147 * t205;
t67 = t105 * t153 + t147 * t207;
t172 = mrSges(6,1) * t66 - mrSges(6,2) * t67;
t76 = t150 * pkin(3) - t83;
t169 = -mrSges(3,1) * t185 + mrSges(3,2) * t186;
t107 = -t146 * t191 + t149 * t192;
t12 = -t145 * t37 + t148 * t48;
t35 = -t145 * t74 + t148 * t89;
t29 = pkin(4) * t104 - pkin(6) * t105 + t76;
t31 = pkin(6) * t209 + t36;
t8 = -t151 * t31 + t153 * t29;
t9 = t151 * t29 + t153 * t31;
t166 = (-qJD(4) * t149 + qJD(2)) * t147;
t164 = mrSges(4,2) * t150 - mrSges(4,3) * t209;
t161 = (mrSges(4,1) * t146 + mrSges(4,2) * t149) * t147;
t157 = -g(1) * t209 + g(2) * t111 - g(3) * t114;
t140 = -qJDD(1) * pkin(1) + qJDD(2);
t116 = t164 * qJD(1);
t109 = t165 * qJDD(1);
t108 = t164 * qJDD(1);
t106 = t146 * t192 + t149 * t191;
t99 = qJD(1) * t161;
t94 = t149 * t176 - t175;
t85 = -qJD(4) * t150 + t107;
t79 = t113 * t194;
t78 = t163 * t194;
t70 = t112 * t148 - t145 * t203;
t62 = -mrSges(5,2) * t180 - mrSges(5,3) * t93;
t59 = -mrSges(5,2) * t174 - mrSges(5,3) * t91;
t58 = t67 * qJD(5);
t57 = t66 * qJD(5);
t52 = Ifges(6,4) * t53;
t44 = t145 * t166 + t148 * t85;
t39 = -t111 * t151 + t153 * t70;
t38 = -t111 * t153 - t151 * t70;
t33 = mrSges(6,1) * t224 - mrSges(6,3) * t55;
t32 = -mrSges(6,2) * t224 + mrSges(6,3) * t53;
t30 = -pkin(4) * t209 - t35;
t17 = Ifges(6,1) * t55 + Ifges(6,5) * t224 + t52;
t16 = Ifges(6,2) * t53 + Ifges(6,6) * t224 + t217;
t15 = -mrSges(6,2) * t87 + mrSges(6,3) * t28;
t14 = mrSges(6,1) * t87 - mrSges(6,3) * t27;
t10 = -pkin(4) * t174 - t12;
t4 = -qJD(5) * t9 + t106 * t153 - t151 * t44;
t3 = qJD(5) * t8 + t106 * t151 + t153 * t44;
t22 = [(-t112 * mrSges(4,1) - t39 * mrSges(6,1) - t38 * mrSges(6,2) - (pkin(4) * t70 + t222) * m(6) - m(5) * t222 - t70 * mrSges(5,1) + t229 * (t112 * t145 + t148 * t203) + mrSges(2,2) * t154 + (-m(3) - m(4)) * t142 + t234 * t111 + t233 * t152) * g(3) + (mrSges(6,3) * t1 + Ifges(6,4) * t27 + Ifges(6,2) * t28 + Ifges(6,6) * t87) * t66 + (mrSges(5,1) * t12 - mrSges(5,2) * t13 + Ifges(5,5) * t92 - Ifges(5,6) * t91) * t209 + (-mrSges(6,3) * t2 + Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t87) * t67 - t10 * t172 + t140 * t168 - pkin(1) * t169 + t47 * t164 + t46 * t165 + (-m(5) * t25 + t216 + t238) * (t145 * t85 - t148 * t166) + (t40 * mrSges(5,1) + t235 / 0.2e1 + t237 / 0.2e1 + t236 / 0.2e1 - Ifges(5,4) * t92 + Ifges(5,2) * t91 + t184 / 0.2e1 - t13 * mrSges(5,3) + t228) * t104 + (-t230 * mrSges(6,1) - t231 * mrSges(6,2) - (-pkin(4) * t71 + t223) * m(6) - m(5) * t223 + t71 * mrSges(5,1) + m(4) * t213 + t115 * mrSges(4,1) + t229 * (-t115 * t145 + t148 * t202) + (-mrSges(2,2) - t225) * t152 + t234 * t114 + t233 * t154) * g(2) + (-m(4) * t64 + m(5) * t50 - t214) * t106 + m(5) * (t12 * t35 + t13 * t36 + t26 * t44 + t40 * t76) + t107 * t116 + t84 * t108 + t83 * t109 + m(3) * (-pkin(1) * t140 + (t134 + t189) * qJ(2) * t196) + t76 * t42 - t58 * t16 / 0.2e1 + t36 * t59 + t35 * t60 + t44 * t62 + t57 * t17 / 0.2e1 + t30 * t7 + t3 * t32 + t4 * t33 + t8 * t14 + t9 * t15 + (-Ifges(5,6) * t104 + Ifges(5,5) * t105 + Ifges(5,3) * t209 + Ifges(4,6) * t150 - (Ifges(4,4) * t149 - Ifges(4,2) * t146) * t147) * t174 + t20 * (mrSges(6,1) * t58 + mrSges(6,2) * t57) + t224 * (Ifges(6,5) * t57 - Ifges(6,6) * t58) / 0.2e1 + (qJD(2) * t99 + m(4) * (qJ(2) * t121 + qJD(2) * t128) + (Ifges(4,1) * t149 - Ifges(4,4) * t146) * t173 + qJ(2) * t197 + Ifges(3,1) * t186 + (-Ifges(4,5) * t149 + Ifges(4,6) * t146 + Ifges(3,4)) * t185) * t147 + (Ifges(3,4) * t186 - Ifges(4,5) * t173 + (Ifges(3,2) + Ifges(4,3)) * t185) * t150 + t53 * (Ifges(6,4) * t57 - Ifges(6,2) * t58) / 0.2e1 + (Ifges(6,1) * t57 - Ifges(6,4) * t58) * t220 + m(6) * (t1 * t9 + t10 * t30 + t2 * t8 + t3 * t6 + t4 * t5) + (-t5 * t57 - t58 * t6) * mrSges(6,3) + m(4) * (t107 * t65 + t46 * t83 + t47 * t84) + (0.2e1 * t196 * t134 - t239) * mrSges(3,3) + t121 * t161 + (mrSges(5,2) * t40 - mrSges(5,3) * t12 + Ifges(5,1) * t92 - Ifges(5,4) * t91) * t105 + Ifges(2,3) * qJDD(1); (-t145 * t219 + t148 * t59 + t193 * t214 + t108) * t146 + t113 * t15 - t163 * t14 - t99 * t194 - t96 * t62 + t169 - t216 * t94 + (-t116 * t193 + t109 - t42) * t149 + m(3) * t140 + t226 * t33 + t227 * t32 + (g(2) * t154 + g(3) * t152) * (-m(3) - t190) + t225 * t196 * qJD(1) ^ 2 + (t1 * t113 + t10 * t210 - t163 * t2 - t20 * t94 + t226 * t5 + t227 * t6) * m(6) + (-(t128 * t147 + t201 * t65 - t208 * t64) * qJD(1) + t146 * t47 + t149 * t46) * m(4) + (t25 * t94 - t26 * t96 - t50 * t208 * qJD(1) - t149 * t40 + (-t12 * t145 + t13 * t148) * t146) * m(5); t79 * t32 - t78 * t33 + t219 * t148 + t190 * t150 * g(1) + (-t151 * t14 + t153 * t15 + t59 + (-t151 * t32 - t153 * t33) * qJD(5)) * t145 + m(5) * (t12 * t148 + t13 * t145) + m(4) * t121 + ((t214 * t149 + (t145 * t216 + t148 * t62 + t116) * t146 - m(5) * (-t146 * t148 * t26 + t149 * t50 + t210 * t25) - m(4) * (-t146 * t65 - t149 * t64) + t210 * t238) * qJD(1) + (g(2) * t152 - t239) * t190) * t147 + t197 + (-t10 * t148 + (t1 * t153 - t151 * t2 + (-t151 * t6 - t153 * t5) * qJD(5)) * t145 - t5 * t78 + t6 * t79) * m(6); t93 * t62 - t216 * t95 + (t224 * t32 + t14) * t153 + (-t224 * t33 + t15) * t151 + t42 + (t1 * t151 + t153 * t2 - t20 * t95 + t157 + t224 * (-t151 * t5 + t153 * t6)) * m(6) + (t25 * t95 + t26 * t93 + t157 + t40) * m(5); -t20 * (mrSges(6,1) * t55 + mrSges(6,2) * t53) - t55 * (Ifges(6,1) * t53 - t217) / 0.2e1 + t16 * t220 - t224 * (Ifges(6,5) * t53 - Ifges(6,6) * t55) / 0.2e1 - t5 * t32 + t6 * t33 - g(1) * t172 - g(2) * (mrSges(6,1) * t38 - mrSges(6,2) * t39) - g(3) * (-mrSges(6,1) * t231 + mrSges(6,2) * t230) + (t5 * t53 + t55 * t6) * mrSges(6,3) + t184 - (-Ifges(6,2) * t55 + t17 + t52) * t53 / 0.2e1 + t228;];
tau = t22;
