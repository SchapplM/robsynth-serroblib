% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:02
% EndTime: 2019-12-05 16:35:20
% DurationCPUTime: 4.89s
% Computational Cost: add. (2620->408), mult. (7024->615), div. (0->0), fcn. (4852->10), ass. (0->199)
t141 = sin(qJ(3));
t144 = cos(qJ(3));
t159 = pkin(3) * t141 - qJ(4) * t144;
t103 = qJD(3) * t159 - qJD(4) * t141;
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t142 = sin(qJ(2));
t137 = sin(pkin(5));
t191 = qJD(1) * t137;
t145 = cos(qJ(2));
t193 = t144 * t145;
t258 = t136 * t103 - (t136 * t142 + t138 * t193) * t191;
t186 = qJD(3) * t141;
t257 = (-pkin(7) * t138 + pkin(8)) * t186 + t258;
t162 = pkin(4) * t136 - pkin(8) * t138;
t157 = pkin(7) + t162;
t175 = t145 * t191;
t185 = qJD(3) * t144;
t256 = -t141 * t175 + t157 * t185;
t188 = qJD(2) * t141;
t118 = t159 * qJD(2);
t170 = t142 * t191;
t120 = qJD(2) * pkin(7) + t170;
t113 = t141 * t120;
t139 = cos(pkin(5));
t199 = t139 * t144;
t173 = qJD(1) * t199;
t87 = -t113 + t173;
t47 = t136 * t118 + t138 * t87;
t255 = -pkin(8) * t188 + qJD(4) * t138 - t47;
t226 = qJD(3) / 0.2e1;
t203 = t136 * t144;
t254 = t175 * t203 + (t103 - t170) * t138;
t140 = sin(qJ(5));
t143 = cos(qJ(5));
t194 = t143 * t144;
t198 = t140 * t141;
t155 = t138 * t194 + t198;
t152 = t155 * qJD(3);
t182 = qJD(5) * t144;
t116 = t136 * qJD(3) + t138 * t188;
t183 = qJD(5) * t116;
t39 = -t140 * t183 + (-t143 * t182 + t152) * qJD(2);
t241 = t39 / 0.2e1;
t196 = t141 * t143;
t197 = t140 * t144;
t154 = -t138 * t197 + t196;
t151 = t154 * qJD(3);
t40 = -t143 * t183 + (t140 * t182 + t151) * qJD(2);
t240 = t40 / 0.2e1;
t253 = t188 / 0.2e1;
t169 = -Ifges(4,6) * qJD(3) / 0.2e1;
t252 = Ifges(4,5) * t226;
t125 = -pkin(3) * t144 - qJ(4) * t141 - pkin(2);
t200 = t138 * t144;
t93 = pkin(7) * t200 + t136 * t125;
t81 = -pkin(8) * t144 + t93;
t98 = t157 * t141;
t33 = -t140 * t81 + t143 * t98;
t251 = qJD(5) * t33 + t256 * t140 + t257 * t143;
t34 = t140 * t98 + t143 * t81;
t250 = -qJD(5) * t34 - t257 * t140 + t256 * t143;
t190 = qJD(1) * t141;
t174 = t139 * t190;
t51 = t174 + (qJD(2) * t162 + t120) * t144;
t122 = -pkin(4) * t138 - pkin(8) * t136 - pkin(3);
t206 = qJ(4) * t138;
t90 = t122 * t143 - t140 * t206;
t249 = qJD(5) * t90 - t140 * t51 + t255 * t143;
t91 = t122 * t140 + t143 * t206;
t248 = -qJD(5) * t91 - t255 * t140 - t143 * t51;
t176 = pkin(7) * t136 + pkin(4);
t247 = -t176 * t186 - t254;
t179 = pkin(7) * t186;
t246 = t136 * t179 + t254;
t245 = -t138 * t179 + t258;
t189 = qJD(2) * t137;
t171 = t145 * t189;
t53 = t120 * t185 + (qJD(3) * t139 + t171) * t190;
t181 = qJD(2) * qJD(3);
t167 = t144 * t181;
t164 = t136 * t167;
t96 = t138 * mrSges(5,2) * t167 + mrSges(5,1) * t164;
t244 = -m(5) * t53 - t96;
t115 = t138 * qJD(3) - t136 * t188;
t112 = qJD(5) - t115;
t243 = Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t164 / 0.2e1;
t187 = qJD(2) * t144;
t83 = t116 * t143 - t140 * t187;
t229 = Ifges(6,4) * t83;
t82 = -t116 * t140 - t143 * t187;
t20 = Ifges(6,2) * t82 + Ifges(6,6) * t112 + t229;
t242 = -t20 / 0.2e1;
t239 = -t82 / 0.2e1;
t238 = t82 / 0.2e1;
t237 = -t83 / 0.2e1;
t236 = t83 / 0.2e1;
t75 = -qJD(3) * pkin(3) + qJD(4) - t87;
t234 = m(5) * t75;
t233 = -t112 / 0.2e1;
t232 = t112 / 0.2e1;
t231 = t136 / 0.2e1;
t230 = -t140 / 0.2e1;
t36 = Ifges(6,5) * t39;
t35 = Ifges(6,6) * t40;
t30 = -mrSges(6,1) * t82 + mrSges(6,2) * t83;
t95 = -mrSges(5,1) * t187 - mrSges(5,3) * t116;
t225 = t30 - t95;
t165 = qJD(1) * t171;
t192 = qJD(3) * t173 + t144 * t165;
t50 = (qJD(4) - t113) * qJD(3) + t192;
t65 = (t103 + t170) * qJD(2);
t18 = t136 * t65 + t138 * t50;
t88 = t120 * t144 + t174;
t79 = qJD(3) * qJ(4) + t88;
t89 = qJD(2) * t125 - t175;
t29 = t136 * t89 + t138 * t79;
t224 = Ifges(5,4) * t136;
t223 = Ifges(5,4) * t138;
t222 = Ifges(6,4) * t140;
t221 = Ifges(6,4) * t143;
t220 = Ifges(5,5) * t138;
t219 = Ifges(5,6) * t136;
t202 = t137 * t142;
t104 = t141 * t202 - t199;
t218 = t104 * t53;
t102 = (mrSges(5,1) * t141 - mrSges(5,3) * t200) * t181;
t12 = -mrSges(6,1) * t40 + mrSges(6,2) * t39;
t210 = t12 - t102;
t209 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t115 - mrSges(5,2) * t116 - mrSges(4,3) * t188;
t204 = t125 * t138;
t201 = t137 * t145;
t195 = t141 * t145;
t180 = pkin(7) * t141 * t53;
t9 = Ifges(6,3) * t164 + t35 + t36;
t178 = Ifges(5,5) * t187;
t177 = Ifges(5,6) * t187;
t172 = t142 * t189;
t168 = t141 * t181;
t16 = pkin(8) * t168 + t18;
t27 = qJD(3) * t51 + t141 * t165;
t23 = -pkin(8) * t187 + t29;
t26 = -pkin(4) * t115 - pkin(8) * t116 + t75;
t5 = -t140 * t23 + t143 * t26;
t1 = qJD(5) * t5 + t140 * t27 + t143 * t16;
t6 = t140 * t26 + t143 * t23;
t2 = -qJD(5) * t6 - t140 * t16 + t143 * t27;
t163 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t161 = -t1 * t140 - t143 * t2;
t160 = t140 * t5 - t143 * t6;
t17 = -t136 * t50 + t138 * t65;
t28 = -t136 * t79 + t138 * t89;
t44 = -mrSges(6,2) * t112 + mrSges(6,3) * t82;
t45 = mrSges(6,1) * t112 - mrSges(6,3) * t83;
t158 = -t140 * t45 + t143 * t44;
t105 = t139 * t141 + t144 * t202;
t70 = t105 * t138 - t136 * t201;
t31 = t104 * t143 - t140 * t70;
t32 = t104 * t140 + t143 * t70;
t46 = t118 * t138 - t136 * t87;
t156 = -t17 * mrSges(5,3) + (Ifges(5,5) * t141 + (Ifges(5,1) * t138 - t224) * t144) * t181 / 0.2e1 + t53 * mrSges(5,2);
t108 = t138 * t196 - t197;
t107 = -t138 * t198 - t194;
t121 = -qJD(2) * pkin(2) - t175;
t135 = Ifges(4,4) * t187;
t153 = t121 * mrSges(4,2) + Ifges(4,1) * t253 + t135 / 0.2e1 + t252 - t87 * mrSges(4,3);
t150 = -t178 / 0.2e1 + t116 * Ifges(5,1) + t115 * Ifges(5,4) - t28 * mrSges(5,3) + t75 * mrSges(5,2);
t149 = -t18 * mrSges(5,3) - (Ifges(5,6) * t141 + (-Ifges(5,2) * t136 + t223) * t144) * t181 / 0.2e1 + t9 / 0.2e1 + t53 * mrSges(5,1) + t35 / 0.2e1 + t36 / 0.2e1 + t163;
t148 = t121 * mrSges(4,1) + t28 * mrSges(5,1) + t169 - (Ifges(4,4) * t141 + t144 * Ifges(4,2)) * qJD(2) / 0.2e1 - Ifges(5,3) * t187 / 0.2e1 + t116 * Ifges(5,5) + t115 * Ifges(5,6) - t29 * mrSges(5,2) - t88 * mrSges(4,3);
t147 = t177 / 0.2e1 - t116 * Ifges(5,4) - t115 * Ifges(5,2) + t112 * Ifges(6,3) + t83 * Ifges(6,5) + t82 * Ifges(6,6) - t29 * mrSges(5,3) + t75 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2);
t146 = qJD(2) ^ 2;
t130 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t187;
t117 = (-mrSges(4,1) * t144 + mrSges(4,2) * t141) * qJD(2);
t111 = (mrSges(4,1) * t141 + mrSges(4,2) * t144) * t181;
t101 = (-mrSges(5,2) * t141 - mrSges(5,3) * t203) * t181;
t100 = t155 * qJD(2);
t99 = t154 * qJD(2);
t94 = mrSges(5,2) * t187 + mrSges(5,3) * t115;
t92 = -pkin(7) * t203 + t204;
t80 = t144 * t176 - t204;
t78 = Ifges(6,4) * t82;
t72 = -qJD(3) * t104 + t144 * t171;
t71 = qJD(3) * t105 + t141 * t171;
t69 = t105 * t136 + t138 * t201;
t58 = -qJD(5) * t108 + t151;
t57 = qJD(5) * t107 + t152;
t52 = -t120 * t186 + t192;
t43 = t136 * t172 + t138 * t72;
t42 = t136 * t72 - t138 * t172;
t37 = -pkin(4) * t188 - t46;
t25 = -mrSges(6,2) * t164 + mrSges(6,3) * t40;
t24 = mrSges(6,1) * t164 - mrSges(6,3) * t39;
t22 = pkin(4) * t187 - t28;
t21 = Ifges(6,1) * t83 + Ifges(6,5) * t112 + t78;
t15 = -pkin(4) * t168 - t17;
t10 = t39 * Ifges(6,4) + t40 * Ifges(6,2) + Ifges(6,6) * t164;
t4 = qJD(5) * t31 + t140 * t71 + t143 * t43;
t3 = -qJD(5) * t32 - t140 * t43 + t143 * t71;
t7 = [t70 * t101 + t104 * t96 + t72 * t130 + t31 * t24 + t32 * t25 + t3 * t45 + t4 * t44 + t43 * t94 - t209 * t71 + t210 * t69 + t225 * t42 + (t104 * t144 - t105 * t141) * mrSges(4,3) * t181 + ((-mrSges(3,2) * t146 - t111) * t145 + (-mrSges(3,1) * t146 + qJD(2) * t117) * t142) * t137 + m(4) * (t218 + t105 * t52 - t71 * t87 + t72 * t88 + (t121 - t175) * t172) + m(5) * (-t17 * t69 + t18 * t70 - t28 * t42 + t29 * t43 + t71 * t75 + t218) + m(6) * (t1 * t32 + t15 * t69 + t2 * t31 + t22 * t42 + t3 * t5 + t4 * t6); (t1 * t107 - t108 * t2 - t5 * t57 + t58 * t6) * mrSges(6,3) - pkin(2) * t111 + t92 * t102 + t107 * t10 / 0.2e1 + t15 * (-mrSges(6,1) * t107 + mrSges(6,2) * t108) + t93 * t101 + t80 * t12 + t57 * t21 / 0.2e1 + t22 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) + t58 * t20 / 0.2e1 + t33 * t24 + t34 * t25 + t250 * t45 + t251 * t44 + (t1 * t34 + t15 * t80 + t2 * t33 + t247 * t22 + t250 * t5 + t251 * t6) * m(6) + ((t169 + (-m(4) * t88 - t130) * pkin(7) + t148) * t141 + (t252 + (-m(4) * t87 - t209 + t234) * pkin(7) + t150 * t138 + t147 * t136 + t153) * t144) * qJD(3) + 0.2e1 * (-m(4) * (t121 * t142 + t193 * t88 - t195 * t87) / 0.2e1 - t195 * t234 / 0.2e1) * t191 + t245 * t94 + (t17 * t92 + t18 * t93 + t245 * t29 + t246 * t28 + t180) * m(5) + t246 * t95 + t247 * t30 + (-m(4) * pkin(2) * t170 + (mrSges(4,2) * t170 + (-0.3e1 / 0.2e1 * Ifges(4,4) * qJD(3) + (-t219 + t220) * t226) * t141) * t141 + (-mrSges(4,1) * t170 + ((Ifges(6,5) * t108 + Ifges(6,6) * t107) * t231 + (-0.3e1 / 0.2e1 * t220 + 0.3e1 / 0.2e1 * t219 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t144 + (-0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(5,1) * t138 ^ 2 / 0.2e1 + (-t223 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t136) * t136) * t141) * qJD(3)) * t144) * qJD(2) + (t53 * mrSges(4,3) + pkin(7) * t96 + t136 * t149 + t138 * t156 + t175 * t209) * t141 + m(4) * (pkin(7) * t144 * t52 + t180) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) + t52 * mrSges(4,3) - t130 * t175) * t144 - t117 * t170 + (Ifges(6,5) * t57 + Ifges(6,6) * t58) * t232 + (Ifges(6,1) * t57 + Ifges(6,4) * t58) * t236 + (Ifges(6,4) * t57 + Ifges(6,2) * t58) * t238 + (Ifges(6,4) * t108 + Ifges(6,2) * t107) * t240 + (Ifges(6,1) * t108 + Ifges(6,4) * t107) * t241 + t108 * t243; t248 * t45 + t249 * t44 + (t100 * t5 - t6 * t99) * mrSges(6,3) - t87 * t130 - t46 * t95 - t22 * (-mrSges(6,1) * t99 + mrSges(6,2) * t100) - t100 * t21 / 0.2e1 + t90 * t24 + t91 * t25 - t47 * t94 - t52 * mrSges(4,2) - t53 * mrSges(4,1) - t37 * t30 - m(5) * (t28 * t46 + t29 * t47 + t75 * t88) + ((-Ifges(6,2) * t140 + t221) * t240 + (Ifges(6,1) * t143 - t222) * t241 + t15 * (mrSges(6,1) * t140 + mrSges(6,2) * t143) + t143 * t243 + t10 * t230 + t225 * qJD(4) + t210 * qJ(4) + t161 * mrSges(6,3) + m(5) * (-qJ(4) * t17 - qJD(4) * t28) + ((-Ifges(6,2) * t143 - t222) * t238 + (-Ifges(6,1) * t140 - t221) * t236 + t22 * (mrSges(6,1) * t143 - mrSges(6,2) * t140) + (-Ifges(6,5) * t140 - Ifges(6,6) * t143) * t232 + t143 * t242 + t21 * t230 + t160 * mrSges(6,3)) * qJD(5) + t156) * t136 + (((Ifges(5,5) * t136 + Ifges(5,6) * t138) * t226 + t169 + Ifges(4,4) * t253 - t148) * t141 + (-t135 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t188 + (t178 / 0.2e1 - t150) * t138 + (-t177 / 0.2e1 - t147) * t136 + (Ifges(4,5) / 0.2e1 - t136 * (Ifges(5,2) * t138 + t224) / 0.2e1 + (-Ifges(6,3) * t138 + (Ifges(6,5) * t143 - Ifges(6,6) * t140) * t136) * t231 + t138 * (Ifges(5,1) * t136 + t223) / 0.2e1) * qJD(3) - t153) * t144) * qJD(2) + t209 * t88 + (Ifges(6,5) * t100 + Ifges(6,6) * t99) * t233 + (Ifges(6,1) * t100 + Ifges(6,4) * t99) * t237 + (Ifges(6,4) * t100 + Ifges(6,2) * t99) * t239 + t99 * t242 + (m(5) * (qJ(4) * t18 + qJD(4) * t29) + qJD(4) * t94 + qJ(4) * t101 - t149) * t138 + t244 * pkin(3) + (t1 * t91 + t2 * t90 - t22 * t37 + (qJ(4) * t15 + qJD(4) * t22) * t136 + t249 * t6 + t248 * t5) * m(6); t140 * t25 + t143 * t24 - t225 * t116 + t158 * qJD(5) + (-t158 - t94) * t115 - m(5) * (t115 * t29 - t116 * t28) + (-t112 * t160 - t116 * t22 - t161) * m(6) - t244; -t22 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + (Ifges(6,1) * t82 - t229) * t237 + t20 * t236 + (Ifges(6,5) * t82 - Ifges(6,6) * t83) * t233 - t5 * t44 + t6 * t45 + (t5 * t82 + t6 * t83) * mrSges(6,3) + t163 + t9 + (-Ifges(6,2) * t83 + t21 + t78) * t239;];
tauc = t7(:);
