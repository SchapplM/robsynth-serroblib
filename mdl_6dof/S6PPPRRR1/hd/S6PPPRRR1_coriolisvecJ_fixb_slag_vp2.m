% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPPRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:11
% EndTime: 2019-03-08 18:39:17
% DurationCPUTime: 2.64s
% Computational Cost: add. (5587->336), mult. (14984->517), div. (0->0), fcn. (14104->16), ass. (0->181)
t210 = qJD(5) / 0.2e1;
t114 = sin(qJ(5));
t159 = qJD(4) * t114;
t209 = t159 / 0.2e1;
t142 = Ifges(6,5) * t210;
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t105 = sin(pkin(8));
t110 = cos(pkin(8));
t103 = sin(pkin(14));
t104 = sin(pkin(13));
t107 = sin(pkin(6));
t108 = cos(pkin(14));
t109 = cos(pkin(13));
t111 = cos(pkin(7));
t162 = t109 * t111;
t124 = (-t103 * t104 + t108 * t162) * t107;
t106 = sin(pkin(7));
t165 = t106 * t108;
t112 = cos(pkin(6));
t98 = qJD(1) * t112 + qJD(2);
t54 = qJD(1) * t124 + t165 * t98;
t145 = t106 * t107 * t109;
t76 = -qJD(1) * t145 + t111 * t98 + qJD(3);
t129 = t105 * t76 + t110 * t54;
t168 = t104 * t108;
t55 = t103 * t106 * t98 + (t103 * t162 + t168) * t107 * qJD(1);
t31 = -t115 * t55 + t129 * t118;
t27 = t31 * qJD(4);
t117 = cos(qJ(5));
t157 = qJD(4) * t117;
t102 = Ifges(6,4) * t157;
t100 = qJD(6) - t157;
t113 = sin(qJ(6));
t116 = cos(qJ(6));
t177 = Ifges(7,4) * t116;
t132 = -Ifges(7,2) * t113 + t177;
t178 = Ifges(7,4) * t113;
t134 = Ifges(7,1) * t116 - t178;
t135 = mrSges(7,1) * t113 + mrSges(7,2) * t116;
t32 = t115 * t129 + t118 * t55;
t30 = qJD(4) * pkin(10) + t32;
t43 = -t105 * t54 + t110 * t76;
t22 = t114 * t43 + t117 * t30;
t20 = qJD(5) * pkin(11) + t22;
t95 = -pkin(5) * t117 - pkin(11) * t114 - pkin(4);
t26 = qJD(4) * t95 - t31;
t5 = -t113 * t20 + t116 * t26;
t6 = t113 * t26 + t116 * t20;
t137 = t6 * t113 + t5 * t116;
t175 = Ifges(7,6) * t113;
t176 = Ifges(7,5) * t116;
t189 = t116 / 0.2e1;
t21 = -t114 * t30 + t117 * t43;
t19 = -qJD(5) * pkin(5) - t21;
t190 = -t113 / 0.2e1;
t155 = qJD(5) * t113;
t91 = t116 * t159 + t155;
t195 = t91 / 0.2e1;
t188 = Ifges(7,4) * t91;
t154 = qJD(5) * t116;
t90 = -t113 * t159 + t154;
t57 = Ifges(7,2) * t90 + Ifges(7,6) * t100 + t188;
t89 = Ifges(7,4) * t90;
t58 = Ifges(7,1) * t91 + Ifges(7,5) * t100 + t89;
t120 = -t137 * mrSges(7,3) + t57 * t190 + t58 * t189 + t19 * t135 + t90 * t132 / 0.2e1 + t134 * t195 + t100 * (-t175 + t176) / 0.2e1;
t29 = -qJD(4) * pkin(4) - t31;
t208 = t120 + t29 * mrSges(6,2) - t21 * mrSges(6,3) + Ifges(6,1) * t209 + t102 / 0.2e1 + t142;
t138 = pkin(5) * t114 - pkin(11) * t117;
t94 = t138 * qJD(5);
t25 = (t94 + t32) * qJD(4);
t7 = qJD(5) * t21 + t117 * t27;
t1 = qJD(6) * t5 + t113 * t25 + t116 * t7;
t2 = -qJD(6) * t6 - t113 * t7 + t116 * t25;
t207 = t1 * t116 - t113 * t2;
t148 = qJD(5) * qJD(6);
t152 = qJD(6) * t113;
t153 = qJD(5) * t117;
t74 = t116 * t148 + (-t114 * t152 + t116 * t153) * qJD(4);
t151 = qJD(6) * t116;
t75 = -t113 * t148 + (-t113 * t153 - t114 * t151) * qJD(4);
t206 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t74 - Ifges(7,6) * t75;
t163 = t108 * t110;
t166 = t105 * t118;
t204 = t106 * (-t103 * t115 + t118 * t163) + t111 * t166;
t164 = t106 * t112;
t63 = t108 * t164 + t124;
t83 = t111 * t112 - t145;
t128 = t105 * t83 + t110 * t63;
t64 = t107 * t168 + (t107 * t162 + t164) * t103;
t203 = -t115 * t64 + t118 * t128;
t202 = 2 * m(6);
t201 = pkin(10) / 0.2e1;
t200 = -t31 / 0.2e1;
t199 = t74 / 0.2e1;
t198 = t75 / 0.2e1;
t197 = -t90 / 0.2e1;
t196 = -t91 / 0.2e1;
t36 = t115 * t128 + t118 * t64;
t47 = -t105 * t63 + t110 * t83;
t23 = t114 * t36 - t47 * t117;
t8 = qJD(5) * t22 + t114 * t27;
t194 = t23 * t8;
t167 = t105 * t115;
t66 = t111 * t167 + (t103 * t118 + t115 * t163) * t106;
t82 = -t105 * t165 + t110 * t111;
t48 = t114 * t66 - t117 * t82;
t193 = t48 * t8;
t84 = -t117 * t110 + t114 * t167;
t192 = t8 * t84;
t191 = -t100 / 0.2e1;
t28 = t32 * qJD(4);
t184 = t28 * t203;
t183 = t28 * t204;
t179 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t90 + mrSges(7,2) * t91 + mrSges(6,3) * t159;
t172 = t118 * t28;
t170 = Ifges(6,6) * qJD(5);
t169 = qJD(5) * mrSges(6,3);
t161 = t113 * t117;
t160 = t116 * t117;
t158 = qJD(4) * t115;
t156 = qJD(4) * t118;
t150 = qJD(6) * t117;
t149 = qJD(4) * qJD(5);
t147 = t8 * t201;
t144 = t105 * t158;
t143 = t105 * t156;
t141 = -t170 / 0.2e1;
t140 = t114 * t149;
t136 = mrSges(7,1) * t116 - mrSges(7,2) * t113;
t133 = Ifges(7,1) * t113 + t177;
t131 = Ifges(7,2) * t116 + t178;
t130 = Ifges(7,5) * t113 + Ifges(7,6) * t116;
t24 = t114 * t47 + t117 * t36;
t12 = -t113 * t203 + t116 * t24;
t11 = -t113 * t24 - t116 * t203;
t49 = t114 * t82 + t117 * t66;
t40 = -t113 * t204 + t116 * t49;
t39 = -t113 * t49 - t116 * t204;
t85 = t110 * t114 + t117 * t167;
t70 = -t113 * t85 - t116 * t166;
t127 = t113 * t166 - t116 * t85;
t122 = t22 * mrSges(6,3) + t6 * mrSges(7,2) - t100 * Ifges(7,3) - t91 * Ifges(7,5) - t90 * Ifges(7,6) + t170 / 0.2e1 + (Ifges(6,4) * t114 + t117 * Ifges(6,2)) * qJD(4) / 0.2e1 - t29 * mrSges(6,1) - t5 * mrSges(7,1);
t119 = qJD(4) ^ 2;
t99 = Ifges(7,3) * t140;
t97 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t157;
t93 = t138 * qJD(4);
t92 = (-mrSges(6,1) * t117 + mrSges(6,2) * t114) * qJD(4);
t88 = (mrSges(6,1) * t114 + mrSges(6,2) * t117) * t149;
t80 = pkin(10) * t160 + t113 * t95;
t79 = -pkin(10) * t161 + t116 * t95;
t78 = mrSges(7,1) * t100 - mrSges(7,3) * t91;
t77 = -mrSges(7,2) * t100 + mrSges(7,3) * t90;
t69 = qJD(5) * t85 + t114 * t143;
t68 = -qJD(5) * t84 + t117 * t143;
t62 = -mrSges(7,2) * t140 + mrSges(7,3) * t75;
t61 = mrSges(7,1) * t140 - mrSges(7,3) * t74;
t60 = t66 * qJD(4);
t59 = t204 * qJD(4);
t52 = -t95 * t152 + t116 * t94 + (t114 * t155 - t116 * t150) * pkin(10);
t51 = t95 * t151 + t113 * t94 + (-t113 * t150 - t114 * t154) * pkin(10);
t50 = -mrSges(7,1) * t75 + mrSges(7,2) * t74;
t46 = t74 * Ifges(7,1) + t75 * Ifges(7,4) + Ifges(7,5) * t140;
t45 = t74 * Ifges(7,4) + t75 * Ifges(7,2) + Ifges(7,6) * t140;
t42 = qJD(6) * t127 - t113 * t68 + t116 * t144;
t41 = qJD(6) * t70 + t113 * t144 + t116 * t68;
t38 = qJD(5) * t49 + t114 * t59;
t37 = -qJD(5) * t48 + t117 * t59;
t34 = t36 * qJD(4);
t33 = t203 * qJD(4);
t18 = t113 * t93 + t116 * t21;
t17 = -t113 * t21 + t116 * t93;
t16 = -qJD(6) * t40 - t113 * t37 + t116 * t60;
t15 = qJD(6) * t39 + t113 * t60 + t116 * t37;
t14 = t113 * t32 + t160 * t31;
t13 = t116 * t32 - t161 * t31;
t10 = -qJD(5) * t23 + t117 * t33;
t9 = qJD(5) * t24 + t114 * t33;
t4 = qJD(6) * t11 + t10 * t116 + t113 * t34;
t3 = -qJD(6) * t12 - t10 * t113 + t116 * t34;
t35 = [t10 * t97 + t11 * t61 + t12 * t62 + t23 * t50 + t3 * t78 + t34 * t92 - t203 * t88 + t4 * t77 + t179 * t9 + m(7) * (t1 * t12 + t11 * t2 + t19 * t9 + t3 * t5 + t4 * t6 + t194) + m(5) * (t27 * t36 - t31 * t34 + t32 * t33 - t184) + m(6) * (t10 * t22 - t21 * t9 + t24 * t7 + t29 * t34 - t184 + t194) + (-t34 * mrSges(5,1) - t33 * mrSges(5,2) + (-t114 * t24 + t117 * t23) * t169) * qJD(4); t15 * t77 + t16 * t78 + t37 * t97 + t39 * t61 + t40 * t62 + t48 * t50 + t60 * t92 - t204 * t88 + t179 * t38 + m(7) * (t1 * t40 + t15 * t6 + t16 * t5 + t19 * t38 + t2 * t39 + t193) + m(5) * (t27 * t66 - t31 * t60 + t32 * t59 - t183) + m(6) * (-t21 * t38 + t22 * t37 + t29 * t60 + t49 * t7 - t183 + t193) + (-t60 * mrSges(5,1) - t59 * mrSges(5,2) + (-t114 * t49 + t117 * t48) * t169) * qJD(4); t41 * t77 + t42 * t78 + t84 * t50 + t70 * t61 - t127 * t62 + t68 * t97 + t179 * t69 + (-t114 * t85 + t117 * t84) * mrSges(6,3) * t149 + m(7) * (-t1 * t127 + t19 * t69 + t2 * t70 + t41 * t6 + t42 * t5 + t192) + m(6) * (-t21 * t69 + t22 * t68 + t7 * t85 + t192) + ((-t119 * mrSges(5,2) - t88) * t118 + (-t119 * mrSges(5,1) + qJD(4) * t92) * t115 + m(5) * (t115 * t27 + t156 * t32 - t158 * t31 - t172) + m(6) * (t158 * t29 - t172)) * t105; -pkin(4) * t88 - t32 * t92 + t79 * t61 + t80 * t62 + (-t13 + t52) * t78 + (-t14 + t51) * t77 - m(7) * (t13 * t5 + t14 * t6) + m(7) * (t1 * t80 + t2 * t79 + t5 * t52 + t51 * t6) + (-pkin(4) * t28 / 0.2e1 - t29 * t32 / 0.2e1) * t202 + (-t99 / 0.2e1 + t7 * mrSges(6,3) - t31 * t97 - t28 * mrSges(6,1) + (t200 * t22 + t201 * t7) * t202 + (0.3e1 / 0.2e1 * t102 + t142 + (-m(6) * t21 + m(7) * t19 + t179) * pkin(10) + t208) * qJD(5) + t206) * t117 + (pkin(10) * t50 + t45 * t190 + t46 * t189 + t28 * mrSges(6,2) + t134 * t199 + t132 * t198 + (mrSges(6,3) + t135) * t8 - t179 * t31 + (-t1 * t113 - t2 * t116) * mrSges(7,3) + 0.2e1 * (t19 * t200 + t147) * m(7) + (t147 + t21 * t31 / 0.2e1) * t202 + (-t116 * t57 / 0.2e1 + t58 * t190 + t19 * t136 + t131 * t197 + t133 * t196 + t130 * t191 + (t5 * t113 - t6 * t116) * mrSges(7,3)) * qJD(6) + (t141 + (-m(6) * t22 - t97) * pkin(10) + ((t176 / 0.2e1 - t175 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,4)) * t114 + (-Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(6,1)) * t117) * qJD(4) - t122) * qJD(5)) * t114; -t7 * mrSges(6,2) - pkin(5) * t50 - t18 * t77 - t17 * t78 - t21 * t97 + t113 * t46 / 0.2e1 + t45 * t189 + t133 * t199 + t131 * t198 + (-mrSges(6,1) - t136) * t8 - t179 * t22 + t207 * mrSges(7,3) + t120 * qJD(6) + ((Ifges(6,4) * t209 + t130 * t210 + t122 + t141) * t114 + (t142 - t102 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t159 - t208) * t117) * qJD(4) + (-pkin(5) * t8 - t17 * t5 - t18 * t6 - t19 * t22) * m(7) + (-t113 * t61 + t116 * t62 + m(7) * t207 + (-m(7) * t137 - t113 * t77 - t116 * t78) * qJD(6)) * pkin(11); t99 - t19 * (mrSges(7,1) * t91 + mrSges(7,2) * t90) + (Ifges(7,1) * t90 - t188) * t196 + t57 * t195 + (Ifges(7,5) * t90 - Ifges(7,6) * t91) * t191 - t5 * t77 + t6 * t78 + (t5 * t90 + t6 * t91) * mrSges(7,3) + (-Ifges(7,2) * t91 + t58 + t89) * t197 - t206;];
tauc  = t35(:);
