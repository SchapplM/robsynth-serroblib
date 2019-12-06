% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:08
% EndTime: 2019-12-05 16:54:23
% DurationCPUTime: 5.16s
% Computational Cost: add. (2150->373), mult. (5766->507), div. (0->0), fcn. (3680->8), ass. (0->180)
t249 = Ifges(5,4) + Ifges(6,4);
t250 = Ifges(5,1) + Ifges(6,1);
t240 = Ifges(5,5) + Ifges(6,5);
t248 = Ifges(5,2) + Ifges(6,2);
t239 = Ifges(6,6) + Ifges(5,6);
t247 = qJD(3) / 0.2e1;
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t170 = qJD(3) * t120;
t118 = sin(qJ(3));
t174 = qJD(2) * t118;
t94 = -t117 * t174 + t170;
t246 = t249 * t94;
t245 = t249 * t120;
t244 = t249 * t117;
t95 = qJD(3) * t117 + t120 * t174;
t243 = t249 * t95;
t122 = cos(qJ(2));
t115 = sin(pkin(5));
t176 = qJD(1) * t115;
t157 = t122 * t176;
t188 = qJD(2) * pkin(2);
t101 = -t157 - t188;
t121 = cos(qJ(3));
t172 = qJD(2) * t121;
t114 = Ifges(4,4) * t172;
t153 = Ifges(4,5) * t247;
t119 = sin(qJ(2));
t158 = t119 * t176;
t100 = qJD(2) * pkin(7) + t158;
t116 = cos(pkin(5));
t175 = qJD(1) * t116;
t156 = t118 * t175;
t72 = t100 * t121 + t156;
t55 = qJD(3) * pkin(8) + t72;
t102 = -pkin(3) * t121 - pkin(8) * t118 - pkin(2);
t74 = qJD(2) * t102 - t157;
t19 = t117 * t74 + t120 * t55;
t11 = qJ(5) * t94 + t19;
t18 = -t117 * t55 + t120 * t74;
t132 = t19 * t117 + t18 * t120;
t143 = mrSges(6,1) * t117 + mrSges(6,2) * t120;
t145 = mrSges(5,1) * t117 + mrSges(5,2) * t120;
t206 = t120 / 0.2e1;
t209 = -t117 / 0.2e1;
t215 = t95 / 0.2e1;
t110 = qJD(4) - t172;
t226 = t240 * t110 + t250 * t95 + t246;
t227 = t239 * t110 + t248 * t94 + t243;
t228 = t110 / 0.2e1;
t229 = t94 / 0.2e1;
t231 = t250 * t120 - t244;
t233 = -t248 * t117 + t245;
t71 = -t118 * t100 + t121 * t175;
t54 = -qJD(3) * pkin(3) - t71;
t27 = -pkin(4) * t94 + qJD(5) + t54;
t10 = -qJ(5) * t95 + t18;
t5 = pkin(4) * t110 + t10;
t222 = t132 * mrSges(5,3) + (t11 * t117 + t5 * t120) * mrSges(6,3) - t27 * t143 - t54 * t145 - t233 * t229 - t231 * t215 - (-t239 * t117 + t240 * t120) * t228 - t227 * t209 - t226 * t206;
t241 = t174 / 0.2e1;
t242 = t101 * mrSges(4,2) - t71 * mrSges(4,3) + Ifges(4,1) * t241 + t114 / 0.2e1 + t153 - t222;
t165 = qJD(2) * qJD(3);
t150 = t118 * t165;
t164 = qJD(3) * qJD(4);
t168 = qJD(4) * t117;
t169 = qJD(3) * t121;
t62 = t120 * t164 + (-t118 * t168 + t120 * t169) * qJD(2);
t167 = qJD(4) * t120;
t129 = t117 * t169 + t118 * t167;
t63 = -qJD(2) * t129 - t117 * t164;
t238 = t239 * t150 + t248 * t63 + t249 * t62;
t237 = t240 * t150 + t249 * t63 + t250 * t62;
t196 = -qJ(5) - pkin(8);
t149 = qJD(4) * t196;
t166 = qJD(5) * t120;
t148 = pkin(3) * t118 - pkin(8) * t121;
t97 = t148 * qJD(2);
t26 = t117 * t97 + t120 * t71;
t236 = t166 - t26 + (qJ(5) * t172 + t149) * t117;
t178 = t120 * t121;
t131 = pkin(4) * t118 - qJ(5) * t178;
t25 = -t117 * t71 + t120 * t97;
t235 = -qJD(2) * t131 - qJD(5) * t117 + t120 * t149 - t25;
t234 = t248 * t120 + t244;
t232 = t250 * t117 + t245;
t230 = -m(4) * t71 + m(5) * t54;
t152 = -Ifges(4,6) * qJD(3) / 0.2e1;
t225 = t240 * t117 + t239 * t120;
t180 = t115 * t122;
t154 = qJD(2) * t180;
t128 = qJD(1) * (qJD(3) * t116 + t154);
t171 = qJD(3) * t118;
t35 = -t100 * t171 + t121 * t128;
t98 = t148 * qJD(3);
t66 = (t98 + t158) * qJD(2);
t3 = t117 * t66 + t120 * t35 + t74 * t167 - t168 * t55;
t4 = -qJD(4) * t19 - t117 * t35 + t120 * t66;
t224 = -t117 * t4 + t120 * t3;
t160 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t161 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t162 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t221 = -t160 * t110 - t161 * t94 - t162 * t95 - t101 * mrSges(4,1) - t18 * mrSges(5,1) - t5 * mrSges(6,1) - t152 + (Ifges(4,4) * t118 + t121 * Ifges(4,2)) * qJD(2) / 0.2e1 + t11 * mrSges(6,2) + t19 * mrSges(5,2) + t72 * mrSges(4,3) - t239 * t229 - (Ifges(6,3) + Ifges(5,3)) * t228 - t240 * t215;
t220 = t62 / 0.2e1;
t219 = t63 / 0.2e1;
t218 = -t94 / 0.2e1;
t216 = -t95 / 0.2e1;
t212 = m(6) * t27;
t211 = -t110 / 0.2e1;
t205 = pkin(4) * t117;
t204 = pkin(7) * t117;
t36 = t100 * t169 + t118 * t128;
t181 = t115 * t119;
t79 = -t116 * t121 + t118 * t181;
t201 = t36 * t79;
t67 = -mrSges(6,2) * t110 + mrSges(6,3) * t94;
t68 = -mrSges(5,2) * t110 + mrSges(5,3) * t94;
t195 = t68 + t67;
t69 = mrSges(6,1) * t110 - t95 * mrSges(6,3);
t70 = mrSges(5,1) * t110 - mrSges(5,3) * t95;
t194 = t70 + t69;
t193 = t102 * t167 + t117 * t98;
t186 = t120 * t98 + t171 * t204;
t111 = pkin(7) * t178;
t76 = t117 * t102 + t111;
t185 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t94 + mrSges(5,2) * t95 + mrSges(4,3) * t174;
t182 = qJ(5) * t118;
t179 = t118 * t120;
t177 = t121 * t122;
t173 = qJD(2) * t119;
t43 = -mrSges(6,1) * t94 + mrSges(6,2) * t95;
t163 = t43 + t185;
t159 = mrSges(4,3) * t172;
t155 = t115 * t173;
t104 = -qJD(3) * mrSges(4,2) + t159;
t151 = -m(4) * t72 - t104;
t20 = -t63 * mrSges(6,1) + t62 * mrSges(6,2);
t146 = mrSges(5,1) * t120 - mrSges(5,2) * t117;
t144 = mrSges(6,1) * t120 - mrSges(6,2) * t117;
t80 = t116 * t118 + t121 * t181;
t47 = -t117 * t80 - t120 * t180;
t130 = t117 * t180 - t120 * t80;
t1 = pkin(4) * t150 - qJ(5) * t62 - qJD(5) * t95 + t4;
t2 = qJ(5) * t63 + qJD(5) * t94 + t3;
t127 = -t4 * mrSges(5,1) - t1 * mrSges(6,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t123 = qJD(2) ^ 2;
t113 = -pkin(4) * t120 - pkin(3);
t109 = Ifges(5,3) * t150;
t108 = Ifges(6,3) * t150;
t106 = t196 * t120;
t105 = t196 * t117;
t99 = (pkin(7) + t205) * t118;
t96 = (-mrSges(4,1) * t121 + mrSges(4,2) * t118) * qJD(2);
t93 = t120 * t102;
t87 = (mrSges(4,1) * t118 + mrSges(4,2) * t121) * t165;
t75 = -t121 * t204 + t93;
t73 = pkin(4) * t129 + pkin(7) * t169;
t65 = (t117 * t119 + t120 * t177) * t176;
t64 = (-t117 * t177 + t119 * t120) * t176;
t60 = Ifges(5,5) * t62;
t59 = Ifges(6,5) * t62;
t58 = Ifges(5,6) * t63;
t57 = Ifges(6,6) * t63;
t49 = -t117 * t182 + t76;
t46 = -qJD(3) * t79 + t121 * t154;
t45 = qJD(3) * t80 + t118 * t154;
t42 = t156 + (qJD(2) * t205 + t100) * t121;
t41 = -qJ(5) * t179 + t93 + (-pkin(4) - t204) * t121;
t40 = -mrSges(5,2) * t150 + mrSges(5,3) * t63;
t39 = -mrSges(6,2) * t150 + mrSges(6,3) * t63;
t38 = mrSges(5,1) * t150 - mrSges(5,3) * t62;
t37 = mrSges(6,1) * t150 - mrSges(6,3) * t62;
t24 = -qJD(4) * t76 + t186;
t23 = (-t118 * t170 - t121 * t168) * pkin(7) + t193;
t21 = -mrSges(5,1) * t63 + mrSges(5,2) * t62;
t12 = -pkin(4) * t63 + t36;
t9 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t179 + (-qJD(5) * t118 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t121) * t117 + t193;
t8 = qJD(4) * t47 + t117 * t155 + t120 * t46;
t7 = qJD(4) * t130 - t117 * t46 + t120 * t155;
t6 = -t118 * t166 + t131 * qJD(3) + (-t111 + (-t102 + t182) * t117) * qJD(4) + t186;
t13 = [-t80 * mrSges(4,3) * t150 + t46 * t104 + t195 * t8 + t194 * t7 - (t40 + t39) * t130 + (t37 + t38) * t47 + (qJD(3) * t159 + t20 + t21) * t79 + t163 * t45 + ((-mrSges(3,2) * t123 - t87) * t122 + (-mrSges(3,1) * t123 + qJD(2) * t96) * t119) * t115 + m(4) * (t35 * t80 + t201 - t45 * t71 + t46 * t72 + (t101 - t157) * t155) + m(5) * (-t130 * t3 + t18 * t7 + t19 * t8 + t4 * t47 + t45 * t54 + t201) + m(6) * (t1 * t47 + t11 * t8 + t12 * t79 - t130 * t2 + t27 * t45 + t5 * t7); -t96 * t158 - pkin(2) * t87 + t99 * t20 + t23 * t68 + t24 * t70 + t41 * t37 + t75 * t38 + t49 * t39 + t76 * t40 + t73 * t43 + t6 * t69 + t9 * t67 - t195 * t65 - t194 * t64 + m(5) * (t18 * t24 + t19 * t23 + t3 * t76 + t4 * t75) - m(5) * (t18 * t64 + t19 * t65) + 0.2e1 * (-t101 / 0.2e1 - t188 / 0.2e1) * m(4) * t158 + (-t108 / 0.2e1 - t109 / 0.2e1 - t59 / 0.2e1 - t60 / 0.2e1 - t57 / 0.2e1 - t58 / 0.2e1 - t161 * t63 - t162 * t62 + (m(4) * pkin(7) + mrSges(4,3)) * t35 + (-mrSges(4,1) * t173 + t122 * t151) * t176 + (0.3e1 / 0.2e1 * t114 + (t185 + t230) * pkin(7) + t153 + t242) * qJD(3) + t127) * t121 + (t12 * t143 + (-t1 * t120 - t117 * t2) * mrSges(6,3) + (-t117 * t3 - t120 * t4) * mrSges(5,3) + (t152 - t221) * qJD(3) + (t27 * t144 + t54 * t146 + (-t11 * t120 + t117 * t5) * mrSges(6,3) + (t117 * t18 - t120 * t19) * mrSges(5,3) + t234 * t218 + t232 * t216 + t225 * t211 - t227 * t120 / 0.2e1) * qJD(4) + (mrSges(4,2) * t158 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t162 * t120 - t161 * t117) * t118 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t160) * t121) * qJD(3)) * qJD(2) + t231 * t220 + t233 * t219 + (t226 * qJD(4) + t238) * t209 + t237 * t206 + (t151 * qJD(3) + t21) * pkin(7) + (-t163 - t212 - t230) * t157 + (mrSges(4,3) + t145 + (m(4) + m(5)) * pkin(7)) * t36) * t118 + (t1 * t41 + t12 * t99 + t2 * t49 + t27 * t73 + (t6 - t64) * t5 + (-t65 + t9) * t11) * m(6); -m(5) * (t18 * t25 + t19 * t26 + t54 * t72) + t232 * t220 + t234 * t219 + t113 * t20 - t71 * t104 + t105 * t37 - t106 * t39 - t26 * t68 - t25 * t70 - t35 * mrSges(4,2) - t42 * t43 - pkin(3) * t21 + ((t43 + t212) * t205 + (-m(5) * t132 - t117 * t68 - t120 * t70) * pkin(8) - t222) * qJD(4) - t12 * t144 + (-mrSges(4,1) - t146) * t36 + t237 * t117 / 0.2e1 + t238 * t206 + (-t1 * t117 + t120 * t2) * mrSges(6,3) + ((Ifges(4,4) * t241 + t225 * t247 + t152 + t221) * t118 + ((Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t174 - t114 / 0.2e1 + t153 - t242) * t121) * qJD(2) + (-t117 * t38 + t120 * t40) * pkin(8) + m(5) * (-pkin(3) * t36 + t224 * pkin(8)) + t224 * mrSges(5,3) - t185 * t72 + t235 * t69 + t236 * t67 + (t1 * t105 - t106 * t2 + t236 * t11 + t113 * t12 + t235 * t5 - t27 * t42) * m(6); (-t43 * t95 + t37) * pkin(4) + (t11 * t95 + t5 * t94) * mrSges(6,3) + (t18 * t94 + t19 * t95) * mrSges(5,3) + t108 + t109 + t59 + t60 + t57 + t58 - t54 * (mrSges(5,1) * t95 + mrSges(5,2) * t94) - t27 * (mrSges(6,1) * t95 + mrSges(6,2) * t94) - t10 * t67 - t18 * t68 + t11 * t69 + t19 * t70 - t127 + (-(t10 - t5) * t11 + (-t27 * t95 + t1) * pkin(4)) * m(6) + (t250 * t94 - t243) * t216 + t227 * t215 + (-t239 * t95 + t240 * t94) * t211 + (-t248 * t95 + t226 + t246) * t218; -t94 * t67 + t95 * t69 + 0.2e1 * (t12 / 0.2e1 + t11 * t218 + t5 * t215) * m(6) + t20;];
tauc = t13(:);
