% Calculate joint inertia matrix for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:05
% EndTime: 2019-03-08 23:39:09
% DurationCPUTime: 1.48s
% Computational Cost: add. (2591->389), mult. (6311->571), div. (0->0), fcn. (7037->14), ass. (0->147)
t189 = 2 * pkin(10);
t140 = sin(pkin(13));
t143 = cos(pkin(13));
t146 = sin(qJ(6));
t150 = cos(qJ(6));
t101 = t140 * t150 + t143 * t146;
t147 = sin(qJ(4));
t83 = t101 * t147;
t100 = -t140 * t146 + t143 * t150;
t84 = t100 * t147;
t48 = t83 * mrSges(7,1) + t84 * mrSges(7,2);
t162 = t143 * t147;
t165 = t140 * t147;
t92 = mrSges(6,1) * t165 + mrSges(6,2) * t162;
t188 = t48 + t92;
t151 = cos(qJ(4));
t142 = sin(pkin(6));
t145 = cos(pkin(6));
t148 = sin(qJ(3));
t149 = sin(qJ(2));
t152 = cos(qJ(3));
t144 = cos(pkin(7));
t153 = cos(qJ(2));
t161 = t144 * t153;
t141 = sin(pkin(7));
t164 = t141 * t148;
t54 = t145 * t164 + (t148 * t161 + t149 * t152) * t142;
t88 = -t141 * t142 * t153 + t144 * t145;
t34 = t147 * t54 - t88 * t151;
t32 = t34 ^ 2;
t163 = t141 * t152;
t52 = -t145 * t163 + (t148 * t149 - t152 * t161) * t142;
t187 = t52 ^ 2;
t90 = t144 * t147 + t151 * t164;
t59 = -t140 * t90 - t143 * t163;
t60 = -t140 * t163 + t143 * t90;
t89 = -t151 * t144 + t147 * t164;
t20 = Ifges(6,1) * t60 + Ifges(6,4) * t59 + Ifges(6,5) * t89;
t186 = t20 / 0.2e1;
t26 = -t146 * t60 + t150 * t59;
t185 = t26 / 0.2e1;
t27 = t146 * t59 + t150 * t60;
t184 = t27 / 0.2e1;
t169 = Ifges(6,4) * t143;
t78 = -Ifges(6,6) * t151 + (-Ifges(6,2) * t140 + t169) * t147;
t183 = t78 / 0.2e1;
t170 = Ifges(6,4) * t140;
t79 = -Ifges(6,5) * t151 + (Ifges(6,1) * t143 - t170) * t147;
t182 = t79 / 0.2e1;
t181 = -t83 / 0.2e1;
t180 = t84 / 0.2e1;
t179 = t100 / 0.2e1;
t178 = t101 / 0.2e1;
t114 = Ifges(6,1) * t140 + t169;
t177 = t114 / 0.2e1;
t176 = pkin(2) * t152;
t175 = pkin(10) * t147;
t174 = pkin(10) * t151;
t173 = pkin(11) + qJ(5);
t123 = pkin(9) * t164;
t80 = t123 + (-pkin(3) - t176) * t144;
t33 = pkin(4) * t89 - qJ(5) * t90 + t80;
t94 = t144 * t148 * pkin(2) + pkin(9) * t163;
t81 = pkin(10) * t144 + t94;
t82 = (-pkin(3) * t152 - pkin(10) * t148 - pkin(2)) * t141;
t42 = t147 * t82 + t151 * t81;
t37 = -qJ(5) * t163 + t42;
t12 = t140 * t33 + t143 * t37;
t28 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t65 = -mrSges(5,1) * t163 - mrSges(5,3) * t90;
t172 = -t65 + t28;
t171 = -Ifges(5,5) * t90 + Ifges(5,6) * t89;
t56 = Ifges(7,5) * t101 + Ifges(7,6) * t100;
t168 = t147 * t34;
t36 = t147 * t88 + t151 * t54;
t167 = t151 * t36;
t110 = -t143 * mrSges(6,1) + t140 * mrSges(6,2);
t166 = -mrSges(5,1) + t110;
t108 = -pkin(4) * t151 - qJ(5) * t147 - pkin(3);
t72 = t140 * t108 + t143 * t174;
t160 = Ifges(5,5) * t147 + Ifges(5,6) * t151;
t159 = t140 ^ 2 + t143 ^ 2;
t6 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t89;
t158 = Ifges(6,5) * t140 / 0.2e1 + Ifges(6,6) * t143 / 0.2e1 + t56 / 0.2e1;
t157 = Ifges(4,5) * t164 + Ifges(4,6) * t163 + Ifges(4,3) * t144;
t10 = -t26 * mrSges(7,1) + t27 * mrSges(7,2);
t55 = -t100 * mrSges(7,1) + t101 * mrSges(7,2);
t11 = -t140 * t37 + t143 * t33;
t41 = -t147 * t81 + t151 * t82;
t43 = Ifges(7,5) * t84 - Ifges(7,6) * t83 - Ifges(7,3) * t151;
t38 = pkin(4) * t163 - t41;
t13 = -t140 * t36 + t143 * t52;
t14 = t140 * t52 + t143 * t36;
t156 = -t13 * t140 + t14 * t143;
t155 = pkin(10) ^ 2;
t139 = t151 ^ 2;
t138 = t147 ^ 2;
t134 = t138 * t155;
t129 = -pkin(5) * t143 - pkin(4);
t117 = Ifges(5,1) * t147 + Ifges(5,4) * t151;
t116 = Ifges(5,4) * t147 + Ifges(5,2) * t151;
t115 = -mrSges(5,1) * t151 + mrSges(5,2) * t147;
t113 = Ifges(6,2) * t143 + t170;
t111 = t173 * t143;
t109 = t173 * t140;
t106 = (pkin(5) * t140 + pkin(10)) * t147;
t105 = -mrSges(6,1) * t151 - mrSges(6,3) * t162;
t104 = mrSges(6,2) * t151 - mrSges(6,3) * t165;
t103 = -mrSges(4,2) * t144 + mrSges(4,3) * t163;
t102 = mrSges(4,1) * t144 - mrSges(4,3) * t164;
t99 = t143 * t108;
t93 = t144 * t176 - t123;
t91 = (-mrSges(4,1) * t152 + mrSges(4,2) * t148) * t141;
t77 = -Ifges(6,3) * t151 + (Ifges(6,5) * t143 - Ifges(6,6) * t140) * t147;
t71 = -t140 * t174 + t99;
t67 = -mrSges(7,1) * t151 - mrSges(7,3) * t84;
t66 = mrSges(7,2) * t151 - mrSges(7,3) * t83;
t64 = mrSges(5,2) * t163 - mrSges(5,3) * t89;
t63 = -t109 * t146 + t111 * t150;
t62 = -t109 * t150 - t111 * t146;
t61 = -pkin(11) * t165 + t72;
t58 = Ifges(7,1) * t101 + Ifges(7,4) * t100;
t57 = Ifges(7,4) * t101 + Ifges(7,2) * t100;
t50 = -pkin(11) * t162 + t99 + (-pkin(10) * t140 - pkin(5)) * t151;
t49 = mrSges(5,1) * t89 + mrSges(5,2) * t90;
t47 = Ifges(5,1) * t90 - Ifges(5,4) * t89 - Ifges(5,5) * t163;
t46 = Ifges(5,4) * t90 - Ifges(5,2) * t89 - Ifges(5,6) * t163;
t45 = Ifges(7,1) * t84 - Ifges(7,4) * t83 - Ifges(7,5) * t151;
t44 = Ifges(7,4) * t84 - Ifges(7,2) * t83 - Ifges(7,6) * t151;
t40 = mrSges(6,1) * t89 - mrSges(6,3) * t60;
t39 = -mrSges(6,2) * t89 + mrSges(6,3) * t59;
t22 = t146 * t50 + t150 * t61;
t21 = -t146 * t61 + t150 * t50;
t19 = Ifges(6,4) * t60 + Ifges(6,2) * t59 + Ifges(6,6) * t89;
t18 = Ifges(6,5) * t60 + Ifges(6,6) * t59 + Ifges(6,3) * t89;
t17 = -pkin(5) * t59 + t38;
t16 = mrSges(7,1) * t89 - mrSges(7,3) * t27;
t15 = -mrSges(7,2) * t89 + mrSges(7,3) * t26;
t9 = pkin(11) * t59 + t12;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t89;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t89;
t5 = pkin(5) * t89 - pkin(11) * t60 + t11;
t4 = t13 * t146 + t14 * t150;
t3 = t13 * t150 - t14 * t146;
t2 = t146 * t5 + t150 * t9;
t1 = -t146 * t9 + t150 * t5;
t23 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t32) + m(6) * (t13 ^ 2 + t14 ^ 2 + t32) + m(5) * (t36 ^ 2 + t187 + t32) + m(4) * (t54 ^ 2 + t88 ^ 2 + t187) + m(3) * (t145 ^ 2 + (t149 ^ 2 + t153 ^ 2) * t142 ^ 2); t54 * t103 + t13 * t40 + t14 * t39 + t4 * t15 + t3 * t16 + t36 * t64 + t88 * t91 + (-t102 + t49) * t52 + (mrSges(3,1) * t153 - mrSges(3,2) * t149) * t142 + (t10 + t172) * t34 + m(7) * (t1 * t3 + t17 * t34 + t2 * t4) + m(6) * (t11 * t13 + t12 * t14 + t34 * t38) + m(5) * (-t34 * t41 + t36 * t42 + t52 * t80) + m(4) * (-pkin(2) * t141 * t88 - t52 * t93 + t54 * t94); t144 * t157 + t90 * t47 + 0.2e1 * t93 * t102 + 0.2e1 * t94 * t103 + 0.2e1 * t80 * t49 + t59 * t19 + t60 * t20 + 0.2e1 * t42 * t64 + 0.2e1 * t41 * t65 + 0.2e1 * t38 * t28 + 0.2e1 * t12 * t39 + 0.2e1 * t11 * t40 + t26 * t7 + t27 * t8 + 0.2e1 * t2 * t15 + 0.2e1 * t1 * t16 + 0.2e1 * t17 * t10 + Ifges(3,3) + (t18 + t6 - t46) * t89 + (-0.2e1 * pkin(2) * t91 + (Ifges(4,1) * t164 + Ifges(4,5) * t144) * t148 + (Ifges(4,6) * t144 + (0.2e1 * Ifges(4,4) * t148 + (Ifges(4,2) + Ifges(5,3)) * t152) * t141 + t171) * t152) * t141 + m(4) * (pkin(2) ^ 2 * t141 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(5) * (t41 ^ 2 + t42 ^ 2 + t80 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2 + t38 ^ 2) + m(7) * (t1 ^ 2 + t17 ^ 2 + t2 ^ 2); mrSges(5,3) * t167 - t54 * mrSges(4,2) + t14 * t104 + t13 * t105 + t3 * t67 + t4 * t66 + (-mrSges(4,1) + t115) * t52 + (t147 * mrSges(5,3) + t188) * t34 + m(7) * (t106 * t34 + t21 * t3 + t22 * t4) + m(6) * (pkin(10) * t168 + t13 * t71 + t14 * t72) + m(5) * (-pkin(3) * t52 + (t167 + t168) * pkin(10)); -t160 * t163 / 0.2e1 + (pkin(10) * t64 + t42 * mrSges(5,3) + t46 / 0.2e1 - t18 / 0.2e1 - t6 / 0.2e1) * t151 + m(7) * (t1 * t21 + t106 * t17 + t2 * t22) + (-t116 / 0.2e1 + t77 / 0.2e1 + t43 / 0.2e1) * t89 + t12 * t104 + t11 * t105 + t106 * t10 + t80 * t115 + t90 * t117 / 0.2e1 + t38 * t92 + t93 * mrSges(4,1) - t94 * mrSges(4,2) + t71 * t40 + t72 * t39 + t2 * t66 + t1 * t67 + t17 * t48 - pkin(3) * t49 + (t143 * t186 - t140 * t19 / 0.2e1 - t41 * mrSges(5,3) + t47 / 0.2e1 + t172 * pkin(10)) * t147 + t21 * t16 + t22 * t15 + m(6) * (t11 * t71 + t12 * t72 + t38 * t175) + t157 + m(5) * (-pkin(3) * t80 + (-t41 * t147 + t42 * t151) * pkin(10)) + t8 * t180 + t7 * t181 + t60 * t182 + t59 * t183 + t45 * t184 + t44 * t185; -0.2e1 * pkin(3) * t115 + 0.2e1 * t72 * t104 + 0.2e1 * t71 * t105 + 0.2e1 * t106 * t48 + 0.2e1 * t21 * t67 + 0.2e1 * t22 * t66 - t83 * t44 + t84 * t45 + Ifges(4,3) + (t138 + t139) * mrSges(5,3) * t189 + (-t43 - t77 + t116) * t151 + m(7) * (t106 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t71 ^ 2 + t72 ^ 2 + t134) + m(5) * (pkin(3) ^ 2 + t139 * t155 + t134) + (-t140 * t78 + t143 * t79 + t92 * t189 + t117) * t147; -t36 * mrSges(5,2) + (t100 * t4 - t101 * t3) * mrSges(7,3) + t156 * mrSges(6,3) + (t55 + t166) * t34 + m(7) * (t129 * t34 + t3 * t62 + t4 * t63) + m(6) * (-pkin(4) * t34 + t156 * qJ(5)); (-t1 * t101 + t2 * t100) * mrSges(7,3) - Ifges(5,3) * t163 + t158 * t89 + t38 * t110 + t59 * t113 / 0.2e1 + t60 * t177 + t129 * t10 + t7 * t179 + t8 * t178 + t57 * t185 + t58 * t184 + t62 * t16 + t63 * t15 + t17 * t55 + t41 * mrSges(5,1) - t42 * mrSges(5,2) - pkin(4) * t28 - t171 + (qJ(5) * t39 + t12 * mrSges(6,3) + t19 / 0.2e1) * t143 + (-t11 * mrSges(6,3) - qJ(5) * t40 + t186) * t140 + m(7) * (t1 * t62 + t129 * t17 + t2 * t63) + m(6) * (-pkin(4) * t38 + (-t11 * t140 + t12 * t143) * qJ(5)); t106 * t55 + t129 * t48 - pkin(4) * t92 + t44 * t179 + t45 * t178 + t58 * t180 + t57 * t181 + t63 * t66 + t62 * t67 + t166 * t175 + (t22 * t100 - t21 * t101) * mrSges(7,3) + (t72 * mrSges(6,3) + qJ(5) * t104 + t147 * t177 + t183) * t143 + (-t147 * t113 / 0.2e1 - qJ(5) * t105 - t71 * mrSges(6,3) + t182) * t140 + m(6) * (-pkin(4) * t175 + (-t71 * t140 + t72 * t143) * qJ(5)) + m(7) * (t106 * t129 + t21 * t62 + t22 * t63) + (-pkin(10) * mrSges(5,2) - t158) * t151 + t160; -0.2e1 * pkin(4) * t110 + t100 * t57 + t101 * t58 + t143 * t113 + t140 * t114 + 0.2e1 * t129 * t55 + Ifges(5,3) + m(7) * (t129 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t159 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (t100 * t63 - t101 * t62) * mrSges(7,3) + 0.2e1 * t159 * qJ(5) * mrSges(6,3); 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t34; m(6) * t38 + m(7) * t17 + t10 + t28; m(6) * t175 + m(7) * t106 + t188; -m(6) * pkin(4) + m(7) * t129 + t110 + t55; m(6) + m(7); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t6; mrSges(7,1) * t21 - mrSges(7,2) * t22 + t43; mrSges(7,1) * t62 - t63 * mrSges(7,2) + t56; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
