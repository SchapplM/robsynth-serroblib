% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:16
% EndTime: 2019-12-31 18:56:25
% DurationCPUTime: 3.89s
% Computational Cost: add. (1938->337), mult. (4415->448), div. (0->0), fcn. (2291->4), ass. (0->166)
t231 = Ifges(5,4) + Ifges(6,4);
t228 = Ifges(5,1) + Ifges(6,1);
t219 = Ifges(5,5) + Ifges(6,5);
t230 = Ifges(5,2) + Ifges(6,2);
t218 = Ifges(6,6) + Ifges(5,6);
t100 = cos(qJ(3));
t98 = sin(qJ(3));
t85 = pkin(3) * t98 - pkin(7) * t100 + qJ(2);
t67 = t85 * qJD(1);
t101 = -pkin(1) - pkin(6);
t92 = qJD(1) * t101 + qJD(2);
t84 = t98 * t92;
t70 = qJD(3) * pkin(7) + t84;
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t27 = t99 * t67 - t70 * t97;
t28 = t67 * t97 + t70 * t99;
t114 = t27 * t99 + t28 * t97;
t127 = mrSges(6,1) * t97 + mrSges(6,2) * t99;
t129 = mrSges(5,1) * t97 + mrSges(5,2) * t99;
t156 = qJD(1) * t100;
t160 = qJD(3) * t99;
t80 = -t156 * t97 + t160;
t15 = qJ(5) * t80 + t28;
t81 = qJD(3) * t97 + t156 * t99;
t14 = -qJ(5) * t81 + t27;
t163 = qJD(1) * t98;
t95 = qJD(4) + t163;
t7 = pkin(4) * t95 + t14;
t131 = t15 * t97 + t7 * t99;
t179 = Ifges(6,6) * t97;
t180 = Ifges(5,6) * t97;
t181 = Ifges(6,5) * t99;
t182 = Ifges(5,5) * t99;
t192 = t99 / 0.2e1;
t195 = -t97 / 0.2e1;
t196 = t95 / 0.2e1;
t198 = t81 / 0.2e1;
t200 = t80 / 0.2e1;
t222 = t231 * t97;
t207 = t228 * t99 - t222;
t223 = t231 * t99;
t208 = -t230 * t97 + t223;
t226 = t231 * t80;
t210 = t219 * t95 + t228 * t81 + t226;
t224 = t231 * t81;
t216 = t218 * t95 + t230 * t80 + t224;
t168 = t100 * t92;
t71 = -qJD(3) * pkin(3) - t168;
t34 = -pkin(4) * t80 + qJD(5) + t71;
t229 = t216 * t195 + t210 * t192 + t34 * t127 + t71 * t129 + (-t180 + t182 - t179 + t181) * t196 + t208 * t200 + t207 * t198 - t114 * mrSges(5,3) - t131 * mrSges(6,3);
t166 = Ifges(4,5) * qJD(3);
t187 = Ifges(4,4) * t98;
t220 = qJD(1) / 0.2e1;
t225 = t166 / 0.2e1 + (t100 * Ifges(4,1) - t187) * t220 + t229;
t155 = qJD(3) * t100;
t137 = qJD(1) * t155;
t147 = t98 * t160;
t152 = qJD(3) * qJD(4);
t153 = qJD(4) * t100;
t45 = t99 * t152 + (-t153 * t97 - t147) * qJD(1);
t142 = t99 * t153;
t161 = qJD(3) * t98;
t110 = t161 * t97 - t142;
t46 = qJD(1) * t110 - t152 * t97;
t221 = t218 * t137 + t230 * t46 + t231 * t45;
t217 = t219 * t137 + t228 * t45 + t231 * t46;
t215 = t218 * t99 + t219 * t97;
t214 = t230 * t99 + t222;
t213 = t228 * t97 + t223;
t212 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t139 = -Ifges(4,6) * qJD(3) / 0.2e1;
t211 = m(5) * t114;
t167 = t101 * t98;
t53 = t99 * t167 + t97 * t85;
t143 = t92 * t155;
t158 = qJD(4) * t99;
t159 = qJD(4) * t97;
t133 = pkin(3) * t100 + pkin(7) * t98;
t78 = qJD(3) * t133 + qJD(2);
t59 = t78 * qJD(1);
t4 = t99 * t143 + t67 * t158 - t159 * t70 + t97 * t59;
t108 = -qJD(4) * t28 + t99 * t59;
t5 = -t143 * t97 + t108;
t134 = t4 * t99 - t5 * t97;
t144 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t145 = -Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t146 = -Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t171 = Ifges(4,4) * t100;
t205 = -t144 * t95 + t145 * t80 + t146 * t81 - t27 * mrSges(5,1) - t7 * mrSges(6,1) - t139 + (-t98 * Ifges(4,2) + t171) * t220 + t15 * mrSges(6,2) + t28 * mrSges(5,2) - t218 * t200 - t219 * t198 - (Ifges(6,3) + Ifges(5,3)) * t196;
t203 = t45 / 0.2e1;
t202 = t46 / 0.2e1;
t201 = -t80 / 0.2e1;
t199 = -t81 / 0.2e1;
t197 = -t95 / 0.2e1;
t190 = pkin(4) * t97;
t175 = -qJ(5) - pkin(7);
t48 = -mrSges(6,2) * t95 + mrSges(6,3) * t80;
t49 = -mrSges(5,2) * t95 + mrSges(5,3) * t80;
t174 = t48 + t49;
t50 = mrSges(6,1) * t95 - mrSges(6,3) * t81;
t51 = mrSges(5,1) * t95 - mrSges(5,3) * t81;
t173 = -t50 - t51;
t83 = t133 * qJD(1);
t37 = t99 * t168 + t97 * t83;
t172 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t80 - mrSges(5,2) * t81 - mrSges(4,3) * t156;
t170 = qJ(2) * mrSges(4,1);
t169 = qJ(2) * mrSges(4,2);
t164 = qJ(5) * t100;
t162 = qJD(3) * mrSges(4,2);
t157 = qJD(5) * t99;
t154 = qJD(3) * t101;
t141 = t100 * t154;
t151 = t99 * t141 + t85 * t158 + t97 * t78;
t150 = t97 * t167;
t149 = t97 * t163;
t148 = t92 * t161;
t140 = -t166 / 0.2e1;
t12 = -t46 * mrSges(6,1) + t45 * mrSges(6,2);
t138 = -t101 * t97 + pkin(4);
t136 = qJD(4) * t175;
t36 = -t168 * t97 + t99 * t83;
t1 = -qJ(5) * t45 - qJD(5) * t81 + (pkin(4) * qJD(1) - t92 * t97) * t155 + t108;
t2 = qJ(5) * t46 + qJD(5) * t80 + t4;
t135 = -t1 * t97 + t2 * t99;
t132 = t15 * t99 - t7 * t97;
t130 = mrSges(5,1) * t99 - mrSges(5,2) * t97;
t128 = mrSges(6,1) * t99 - mrSges(6,2) * t97;
t113 = t27 * t97 - t28 * t99;
t109 = t173 * t99 - t174 * t97;
t107 = t5 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2);
t106 = -qJD(4) * t53 + t99 * t78;
t96 = -pkin(4) * t99 - pkin(3);
t94 = Ifges(5,3) * t137;
t93 = Ifges(6,3) * t137;
t90 = t175 * t99;
t89 = t175 * t97;
t87 = -mrSges(4,3) * t163 - t162;
t82 = qJD(1) * (t98 * mrSges(4,1) + t100 * mrSges(4,2));
t79 = (-t101 + t190) * t100;
t77 = t99 * t85;
t61 = -qJD(5) * t97 + t136 * t99;
t60 = t136 * t97 + t157;
t56 = -pkin(4) * t149 + t84;
t52 = t77 - t150;
t47 = -pkin(4) * t110 + t154 * t98;
t44 = Ifges(5,5) * t45;
t43 = Ifges(6,5) * t45;
t42 = Ifges(5,6) * t46;
t41 = Ifges(6,6) * t46;
t38 = -mrSges(6,1) * t80 + mrSges(6,2) * t81;
t35 = -t164 * t97 + t53;
t33 = -mrSges(5,2) * t137 + mrSges(5,3) * t46;
t32 = -mrSges(6,2) * t137 + mrSges(6,3) * t46;
t31 = mrSges(5,1) * t137 - mrSges(5,3) * t45;
t30 = mrSges(6,1) * t137 - mrSges(6,3) * t45;
t29 = t138 * t98 - t164 * t99 + t77;
t20 = qJ(5) * t149 + t37;
t19 = -pkin(4) * t46 + t148;
t18 = (qJ(5) * t98 * t99 + pkin(4) * t100) * qJD(1) + t36;
t17 = -t141 * t97 + t106;
t16 = -qJD(4) * t150 + t151;
t13 = -mrSges(5,1) * t46 + mrSges(5,2) * t45;
t6 = -qJ(5) * t142 + (-qJD(5) * t100 + (qJ(5) * qJD(3) - qJD(4) * t101) * t98) * t97 + t151;
t3 = qJ(5) * t147 + (qJ(5) * t159 + qJD(3) * t138 - t157) * t100 + t106;
t8 = [t79 * t12 + t16 * t49 + t17 * t51 + t29 * t30 + t3 * t50 + t52 * t31 + t35 * t32 + t53 * t33 + t47 * t38 + t6 * t48 + m(5) * (t16 * t28 + t17 * t27 + t4 * t53 + t5 * t52) + m(6) * (t1 * t29 + t15 * t6 + t19 * t79 + t2 * t35 + t3 * t7 + t34 * t47) + (t82 + 0.2e1 * t212) * qJD(2) + (t94 / 0.2e1 + t93 / 0.2e1 + t44 / 0.2e1 + t42 / 0.2e1 + t43 / 0.2e1 + t41 / 0.2e1 + qJD(1) * qJD(2) * mrSges(4,1) - t145 * t46 - t146 * t45 + ((0.3e1 / 0.2e1 * t187 - 0.2e1 * t169) * qJD(1) + t140 + (m(5) * t71 - t172) * t101 - t225) * qJD(3) + t107) * t98 + (t19 * t127 - t101 * t13 + (-t1 * t99 - t2 * t97) * mrSges(6,3) + (-t4 * t97 - t5 * t99) * mrSges(5,3) + (t129 * t84 + t139 + (-m(5) * t84 + t87) * t101 - t205) * qJD(3) + (mrSges(5,3) * t113 - mrSges(6,3) * t132 + t128 * t34 + t130 * t71 + t214 * t201 + t213 * t199 + t215 * t197 - t216 * t99 / 0.2e1) * qJD(4) + (qJD(2) * mrSges(4,2) + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t181 / 0.2e1 - t179 / 0.2e1 + t182 / 0.2e1 - t180 / 0.2e1) * t100 + 0.2e1 * t170 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + t144) * t98) * qJD(3)) * qJD(1) + t207 * t203 + t208 * t202 + t217 * t192 + (qJD(4) * t210 + t221) * t195) * t100; (-t12 - t13 - m(6) * t19 + (-m(5) * t113 + m(6) * t132 + t173 * t97 + t174 * t99 + t87) * qJD(3)) * t100 + ((t32 + t33) * t99 + (-t30 - t31) * t97 + (t38 - t172) * qJD(3) + t109 * qJD(4) + m(6) * (qJD(3) * t34 - t15 * t159 - t158 * t7 + t135) + m(5) * (qJD(3) * t71 - t158 * t27 - t159 * t28 + t134 - t143)) * t98 + (-m(6) * t131 + t109 - t211 - t212 - t82) * qJD(1); ((m(6) * t34 + t38) * t190 + (-t97 * t49 - t99 * t51 - t211) * pkin(7) + t229) * qJD(4) + t134 * mrSges(5,3) + t135 * mrSges(6,3) + (-t97 * t31 + t99 * t33) * pkin(7) + (-t18 + t61) * t50 + m(6) * (t1 * t89 + t15 * t60 + t19 * t96 - t2 * t90 + t61 * t7) + ((-t87 - t162) * t100 + ((-mrSges(4,1) - t130) * qJD(3) + t172) * t98) * t92 - t19 * t128 + t96 * t12 + t89 * t30 - t90 * t32 - t37 * t49 - t36 * t51 - t56 * t38 - pkin(3) * t13 + (t60 - t20) * t48 + ((t139 + (-t170 + t171 / 0.2e1) * qJD(1) + t215 * qJD(3) / 0.2e1 + t205) * t100 + ((t169 - t187 / 0.2e1 + (-Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t100) * qJD(1) + t140 + t225) * t98) * qJD(1) - m(6) * (t15 * t20 + t18 * t7 + t34 * t56) + t213 * t203 + t214 * t202 + t217 * t97 / 0.2e1 + t221 * t192 + (-pkin(3) * t148 + pkin(7) * t134 - t27 * t36 - t28 * t37 - t71 * t84) * m(5); (-t38 * t81 + t30) * pkin(4) + (t15 * t81 + t7 * t80) * mrSges(6,3) + (t27 * t80 + t28 * t81) * mrSges(5,3) + (-(t14 - t7) * t15 + (-t34 * t81 + t1) * pkin(4)) * m(6) + t94 + t93 + t44 + t42 + t43 + t41 - t34 * (t81 * mrSges(6,1) + t80 * mrSges(6,2)) - t71 * (mrSges(5,1) * t81 + mrSges(5,2) * t80) - t27 * t49 + t15 * t50 + t28 * t51 - t14 * t48 + t107 + (t228 * t80 - t224) * t199 + t216 * t198 + (-t218 * t81 + t219 * t80) * t197 + (-t230 * t81 + t210 + t226) * t201; -t80 * t48 + t81 * t50 + 0.2e1 * (t19 / 0.2e1 + t15 * t201 + t7 * t198) * m(6) + t12;];
tauc = t8(:);
