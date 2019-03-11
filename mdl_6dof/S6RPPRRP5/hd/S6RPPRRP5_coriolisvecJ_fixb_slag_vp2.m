% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:41
% EndTime: 2019-03-09 02:07:50
% DurationCPUTime: 4.54s
% Computational Cost: add. (2472->391), mult. (5219->509), div. (0->0), fcn. (2612->4), ass. (0->180)
t248 = Ifges(6,4) + Ifges(7,4);
t249 = Ifges(6,1) + Ifges(7,1);
t239 = Ifges(6,5) + Ifges(7,5);
t247 = Ifges(6,2) + Ifges(7,2);
t238 = Ifges(7,6) + Ifges(6,6);
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t107 = pkin(1) + qJ(3);
t109 = sin(qJ(4));
t111 = cos(qJ(4));
t87 = pkin(4) * t109 - pkin(8) * t111 + t107;
t62 = qJD(1) * t87 - qJD(2);
t102 = qJD(1) * qJ(2) + qJD(3);
t94 = -qJD(1) * pkin(7) + t102;
t86 = t109 * t94;
t71 = qJD(4) * pkin(8) + t86;
t20 = -t108 * t71 + t110 * t62;
t21 = t108 * t62 + t110 * t71;
t125 = t108 * t21 + t110 * t20;
t138 = mrSges(7,1) * t108 + mrSges(7,2) * t110;
t140 = mrSges(6,1) * t108 + mrSges(6,2) * t110;
t171 = qJD(4) * t110;
t175 = qJD(1) * t111;
t82 = -t108 * t175 + t171;
t15 = qJ(6) * t82 + t21;
t173 = qJD(4) * t108;
t83 = t110 * t175 + t173;
t14 = -qJ(6) * t83 + t20;
t176 = qJD(1) * t109;
t99 = qJD(5) + t176;
t7 = pkin(5) * t99 + t14;
t142 = t108 * t15 + t110 * t7;
t188 = Ifges(7,6) * t108;
t189 = Ifges(6,6) * t108;
t190 = Ifges(7,5) * t110;
t191 = Ifges(6,5) * t110;
t209 = t110 / 0.2e1;
t212 = -t108 / 0.2e1;
t213 = t99 / 0.2e1;
t215 = t83 / 0.2e1;
t217 = t82 / 0.2e1;
t245 = t248 * t82;
t226 = t239 * t99 + t249 * t83 + t245;
t242 = t248 * t108;
t229 = t249 * t110 - t242;
t243 = t248 * t110;
t231 = -t247 * t108 + t243;
t241 = t248 * t83;
t235 = t238 * t99 + t247 * t82 + t241;
t72 = -qJD(4) * pkin(4) - t111 * t94;
t34 = -pkin(5) * t82 + qJD(6) + t72;
t246 = t235 * t212 + t226 * t209 + t34 * t138 + t72 * t140 + (-t189 + t191 - t188 + t190) * t213 + t231 * t217 + t229 * t215 - t125 * mrSges(6,3) - t142 * mrSges(7,3);
t186 = Ifges(5,5) * qJD(4);
t196 = Ifges(5,4) * t109;
t240 = qJD(1) / 0.2e1;
t244 = t186 / 0.2e1 + (t111 * Ifges(5,1) - t196) * t240 + t246;
t170 = qJD(4) * t111;
t152 = qJD(1) * t170;
t163 = qJD(4) * qJD(5);
t166 = qJD(5) * t111;
t46 = t110 * t163 + (-t108 * t166 - t109 * t171) * qJD(1);
t155 = t110 * t166;
t172 = qJD(4) * t109;
t120 = t108 * t172 - t155;
t47 = qJD(1) * t120 - t108 * t163;
t237 = t238 * t152 + t247 * t47 + t248 * t46;
t236 = t239 * t152 + t248 * t47 + t249 * t46;
t53 = mrSges(7,1) * t99 - mrSges(7,3) * t83;
t54 = mrSges(6,1) * t99 - mrSges(6,3) * t83;
t198 = t53 + t54;
t51 = -mrSges(7,2) * t99 + mrSges(7,3) * t82;
t52 = -mrSges(6,2) * t99 + mrSges(6,3) * t82;
t199 = t51 + t52;
t119 = t108 * t198 - t110 * t199;
t181 = qJD(4) * mrSges(5,2);
t89 = -mrSges(5,3) * t176 - t181;
t234 = -t119 + t89;
t233 = t239 * t108 + t238 * t110;
t232 = t247 * t110 + t242;
t230 = t249 * t108 + t243;
t106 = qJ(2) - pkin(7);
t228 = qJD(2) * t109 + t106 * t170;
t153 = -Ifges(5,6) * qJD(4) / 0.2e1;
t227 = m(6) * t125;
t167 = qJD(5) * t110;
t168 = qJD(5) * t108;
t165 = qJD(1) * qJD(2);
t61 = t109 * t165 + t170 * t94;
t148 = pkin(4) * t111 + pkin(8) * t109;
t81 = qJD(4) * t148 + qJD(3);
t63 = t81 * qJD(1);
t4 = t108 * t63 + t110 * t61 + t62 * t167 - t168 * t71;
t5 = -qJD(5) * t21 - t108 * t61 + t110 * t63;
t145 = -t108 * t5 + t110 * t4;
t224 = qJD(1) * t107;
t223 = (t109 ^ 2 + t111 ^ 2) * t94;
t159 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t160 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t161 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t95 = -qJD(2) + t224;
t222 = -t159 * t99 + t160 * t82 + t161 * t83 - t20 * mrSges(6,1) - t7 * mrSges(7,1) - t95 * mrSges(5,1) - t153 + (Ifges(5,4) * t111 - Ifges(5,2) * t109) * t240 + t15 * mrSges(7,2) + t21 * mrSges(6,2) - t238 * t217 - t239 * t215 - (Ifges(7,3) + Ifges(6,3)) * t213;
t220 = t46 / 0.2e1;
t219 = t47 / 0.2e1;
t218 = -t82 / 0.2e1;
t216 = -t83 / 0.2e1;
t214 = -t99 / 0.2e1;
t208 = pkin(5) * t108;
t202 = -qJ(6) - pkin(8);
t30 = mrSges(7,1) * t152 - mrSges(7,3) * t46;
t31 = mrSges(6,1) * t152 - mrSges(6,3) * t46;
t201 = -t30 - t31;
t32 = -mrSges(7,2) * t152 + mrSges(7,3) * t47;
t33 = -mrSges(6,2) * t152 + mrSges(6,3) * t47;
t200 = t32 + t33;
t177 = t110 * t111;
t85 = t148 * qJD(1);
t37 = t108 * t85 + t94 * t177;
t182 = qJD(4) * mrSges(5,1);
t197 = mrSges(6,1) * t82 - mrSges(6,2) * t83 - mrSges(5,3) * t175 + t182;
t178 = t109 * t110;
t50 = t106 * t178 + t108 * t87;
t60 = -t111 * t165 + t172 * t94;
t187 = t106 * t60;
t183 = qJD(3) * t95;
t180 = t108 * t109;
t179 = t108 * t111;
t169 = qJD(5) * t106;
t164 = qJD(1) * qJD(3);
t38 = -mrSges(7,1) * t82 + mrSges(7,2) * t83;
t162 = t38 - t197;
t158 = t108 * t176;
t156 = t109 * t169;
t154 = -t186 / 0.2e1;
t16 = -t47 * mrSges(7,1) + t46 * mrSges(7,2);
t151 = qJD(5) * t202;
t150 = t108 * t81 + t110 * t228 + t87 * t167;
t149 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t36 = t110 * t85 - t179 * t94;
t1 = pkin(5) * t152 - qJ(6) * t46 - qJD(6) * t83 + t5;
t2 = qJ(6) * t47 + qJD(6) * t82 + t4;
t147 = -t1 * t110 - t108 * t2;
t146 = -t1 * t108 + t110 * t2;
t144 = -t108 * t4 - t110 * t5;
t143 = t108 * t7 - t110 * t15;
t141 = mrSges(6,1) * t110 - mrSges(6,2) * t108;
t139 = mrSges(7,1) * t110 - mrSges(7,2) * t108;
t124 = t108 * t20 - t110 * t21;
t118 = -t108 * t199 - t110 * t198;
t117 = t5 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2);
t116 = -qJD(6) * t111 + (qJ(6) * qJD(4) - t169) * t109;
t112 = qJD(1) ^ 2;
t101 = -pkin(5) * t110 - pkin(4);
t98 = Ifges(6,3) * t152;
t97 = Ifges(7,3) * t152;
t92 = t202 * t110;
t91 = t202 * t108;
t84 = qJD(1) * (mrSges(5,1) * t109 + mrSges(5,2) * t111);
t80 = (-t106 + t208) * t111;
t76 = t110 * t87;
t67 = t110 * t81;
t65 = -qJD(6) * t108 + t110 * t151;
t64 = qJD(6) * t110 + t108 * t151;
t59 = -pkin(5) * t158 + t86;
t49 = -t106 * t180 + t76;
t45 = Ifges(6,5) * t46;
t44 = Ifges(7,5) * t46;
t43 = Ifges(6,6) * t47;
t42 = Ifges(7,6) * t47;
t40 = -pkin(5) * t120 - qJD(2) * t111 + t106 * t172;
t35 = -qJ(6) * t179 + t50;
t29 = -qJ(6) * t177 + t76 + (-t106 * t108 + pkin(5)) * t109;
t22 = qJ(6) * t158 + t37;
t19 = (pkin(5) * t111 + qJ(6) * t178) * qJD(1) + t36;
t18 = -pkin(5) * t47 + t60;
t17 = -mrSges(6,1) * t47 + mrSges(6,2) * t46;
t9 = -t110 * t156 + t67 + (-qJD(5) * t87 - t228) * t108;
t8 = -t108 * t156 + t150;
t6 = -qJ(6) * t155 + t108 * t116 + t150;
t3 = pkin(5) * t170 + t67 + t116 * t110 + ((qJ(6) * t111 - t87) * qJD(5) - t228) * t108;
t10 = [qJD(3) * t84 + t80 * t16 + t29 * t30 + t3 * t53 + t49 * t31 + t35 * t32 + t50 * t33 + t40 * t38 + t6 * t51 + t8 * t52 + t9 * t54 + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t149) * qJD(1) + m(4) * (qJD(2) * t102 + t183 + (qJ(2) * qJD(2) + qJD(3) * t107) * qJD(1)) + m(6) * (t20 * t9 + t21 * t8 + t4 * t50 + t49 * t5) + m(5) * (qJD(2) * t223 + t107 * t164 + t183) + m(7) * (t1 * t29 + t15 * t6 + t18 * t80 + t2 * t35 + t3 * t7 + t34 * t40) + (t97 / 0.2e1 + t98 / 0.2e1 + t45 / 0.2e1 + t42 / 0.2e1 + t43 / 0.2e1 + t44 / 0.2e1 + mrSges(5,1) * t164 + qJD(2) * t89 + (m(5) * t106 - mrSges(5,3)) * t61 - t160 * t47 - t161 * t46 + (0.3e1 / 0.2e1 * Ifges(5,4) * t176 + (m(6) * t72 - t197) * t106 + (-t95 - t224) * mrSges(5,2) + t154 - t244) * qJD(4) + t117) * t109 + (t18 * t138 - t106 * t17 + (mrSges(5,3) + t140) * t60 + t197 * qJD(2) + t147 * mrSges(7,3) + t144 * mrSges(6,3) - m(5) * t187 + m(6) * (-qJD(2) * t72 - t187) + (t106 * t89 + t153 - t222) * qJD(4) + (mrSges(6,3) * t124 + mrSges(7,3) * t143 + t139 * t34 + t141 * t72 + t232 * t218 + t230 * t216 + t233 * t214 - t235 * t110 / 0.2e1) * qJD(5) + (qJD(3) * mrSges(5,2) + (t107 * mrSges(5,1) + (t190 / 0.2e1 - t188 / 0.2e1 + t191 / 0.2e1 - t189 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t111 + (0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + t159) * t109) * qJD(4)) * qJD(1) + t229 * t220 + t231 * t219 + t236 * t209) * t111 + (qJD(5) * t226 + t237) * t111 * t212; t201 * t110 - t200 * t108 + t119 * qJD(5) + m(6) * (t124 * qJD(5) + t144) + m(7) * (t143 * qJD(5) + t147) - t149 * t112 + ((-t102 - qJD(3)) * m(4) + (t162 - t182) * t111 + (t181 - t234) * t109 - m(6) * (-t111 * t72 + t178 * t21 - t180 * t20) - m(7) * (-t111 * t34 + t15 * t178 - t180 * t7) + (-qJD(3) - t223) * m(5)) * qJD(1); -t112 * mrSges(4,3) + (-t16 - t17 + t234 * qJD(4) + m(6) * (t171 * t21 - t173 * t20 - t60) + m(7) * (t15 * t171 - t173 * t7 - t18) - m(5) * t60) * t111 + (t200 * t110 + t201 * t108 + t162 * qJD(4) + t118 * qJD(5) + m(6) * (qJD(4) * t72 - t167 * t20 - t168 * t21 + t145) + m(7) * (qJD(4) * t34 - t15 * t168 - t167 * t7 + t146) + m(5) * t61) * t109 + (-m(5) * t95 - t84 + (-t95 + qJD(2)) * m(4) - t227 - m(7) * t142 + t118) * qJD(1); (t65 - t19) * t53 + ((t153 + Ifges(5,4) * t175 / 0.2e1 + t233 * qJD(4) / 0.2e1 + t222) * t111 + ((-t196 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t111) * qJD(1) + t95 * mrSges(5,2) + t154 + t244) * t109) * qJD(1) + ((m(7) * t34 + t38) * t208 + (-t108 * t52 - t110 * t54 - t227) * pkin(8) + t246) * qJD(5) + t101 * t16 + (t64 - t22) * t51 - m(7) * (t15 * t22 + t19 * t7 + t34 * t59) + t91 * t30 - t92 * t32 - t37 * t52 - t36 * t54 - t59 * t38 - t61 * mrSges(5,2) - pkin(4) * t17 + m(7) * (t1 * t91 + t101 * t18 + t15 * t64 - t2 * t92 + t65 * t7) + t236 * t108 / 0.2e1 + t237 * t209 + t230 * t220 + t232 * t219 + (-t108 * t31 + t110 * t33) * pkin(8) + t145 * mrSges(6,3) + t146 * mrSges(7,3) - t18 * t139 + (-mrSges(5,1) - t141) * t60 + (t109 * t197 - t111 * t89) * t94 + (-pkin(4) * t60 + pkin(8) * t145 - t20 * t36 - t21 * t37 - t72 * t86) * m(6); (-t38 * t83 + t30) * pkin(5) + (t15 * t83 + t7 * t82) * mrSges(7,3) + (t20 * t82 + t21 * t83) * mrSges(6,3) + t97 + t98 + t45 + t42 + t43 + t44 - t72 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) - t34 * (mrSges(7,1) * t83 + mrSges(7,2) * t82) - t20 * t52 + t15 * t53 + t21 * t54 - t14 * t51 + t117 + (-(t14 - t7) * t15 + (-t34 * t83 + t1) * pkin(5)) * m(7) + (t249 * t82 - t241) * t216 + t235 * t215 + (-t238 * t83 + t239 * t82) * t214 + (-t247 * t83 + t226 + t245) * t218; -t82 * t51 + t83 * t53 + 0.2e1 * (t18 / 0.2e1 + t15 * t218 + t7 * t215) * m(7) + t16;];
tauc  = t10(:);
