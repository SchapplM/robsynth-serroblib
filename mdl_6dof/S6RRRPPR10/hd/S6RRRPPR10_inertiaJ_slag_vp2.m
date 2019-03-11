% Calculate joint inertia matrix for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:31
% EndTime: 2019-03-09 16:21:34
% DurationCPUTime: 1.80s
% Computational Cost: add. (2645->396), mult. (5761->547), div. (0->0), fcn. (6088->10), ass. (0->149)
t202 = pkin(4) + pkin(9);
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t201 = t152 ^ 2 + t155 ^ 2;
t146 = sin(pkin(11));
t191 = -t146 / 0.2e1;
t200 = -m(5) * pkin(3) + mrSges(5,2);
t151 = sin(qJ(6));
t154 = cos(qJ(6));
t148 = cos(pkin(11));
t147 = sin(pkin(6));
t156 = cos(qJ(2));
t179 = t147 * t156;
t149 = cos(pkin(6));
t153 = sin(qJ(2));
t180 = t147 * t153;
t88 = -t149 * t155 + t152 * t180;
t59 = t146 * t179 + t148 * t88;
t60 = t146 * t88 - t148 * t179;
t28 = -t151 * t60 + t154 * t59;
t199 = t28 / 0.2e1;
t29 = t151 * t59 + t154 * t60;
t198 = t29 / 0.2e1;
t197 = t59 / 0.2e1;
t182 = Ifges(6,4) * t148;
t75 = Ifges(6,5) * t152 + (-Ifges(6,1) * t146 - t182) * t155;
t196 = t75 / 0.2e1;
t160 = t146 * t151 - t148 * t154;
t79 = t160 * t155;
t195 = t79 / 0.2e1;
t98 = -t154 * t146 - t151 * t148;
t80 = t98 * t155;
t194 = t80 / 0.2e1;
t193 = t98 / 0.2e1;
t192 = -t160 / 0.2e1;
t190 = -t148 / 0.2e1;
t189 = pkin(1) * t156;
t123 = pkin(8) * t180;
t91 = t149 * t189 - t123;
t188 = t91 * mrSges(3,1);
t92 = t149 * t153 * pkin(1) + pkin(8) * t179;
t187 = t92 * mrSges(3,2);
t186 = Ifges(5,1) + Ifges(4,3);
t185 = pkin(3) + qJ(5);
t184 = -pkin(10) - t185;
t77 = pkin(9) * t149 + t92;
t78 = (-pkin(2) * t156 - pkin(9) * t153 - pkin(1)) * t147;
t36 = -t152 * t77 + t155 * t78;
t33 = pkin(3) * t179 - t36;
t89 = t149 * t152 + t155 * t180;
t21 = pkin(4) * t89 + qJ(5) * t179 + t33;
t76 = t123 + (-pkin(2) - t189) * t149;
t159 = -qJ(4) * t89 + t76;
t23 = t185 * t88 + t159;
t6 = t146 * t21 + t148 * t23;
t37 = t152 * t78 + t155 * t77;
t56 = -Ifges(7,5) * t160 + Ifges(7,6) * t98;
t183 = Ifges(6,4) * t146;
t118 = t202 * t152;
t169 = -qJ(4) * t152 - pkin(2);
t96 = -t185 * t155 + t169;
t51 = t146 * t118 + t148 * t96;
t181 = t146 * t155;
t178 = t148 * t155;
t108 = t146 * mrSges(6,1) + t148 * mrSges(6,2);
t177 = Ifges(4,5) * t152 + Ifges(4,6) * t155;
t176 = t201 * pkin(9) ^ 2;
t119 = t202 * t155;
t175 = t146 ^ 2 + t148 ^ 2;
t174 = t160 ^ 2 + t98 ^ 2;
t7 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t89;
t38 = Ifges(7,5) * t80 + Ifges(7,6) * t79 + Ifges(7,3) * t152;
t173 = Ifges(6,5) * t148 / 0.2e1 + Ifges(6,6) * t191 + t56 / 0.2e1;
t172 = Ifges(3,5) * t180 + Ifges(3,6) * t179 + Ifges(3,3) * t149;
t171 = t179 / 0.2e1;
t30 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t10 = -t28 * mrSges(7,1) + t29 * mrSges(7,2);
t46 = -t79 * mrSges(7,1) + t80 * mrSges(7,2);
t55 = -t98 * mrSges(7,1) - mrSges(7,2) * t160;
t170 = m(6) * t175;
t5 = -t146 * t23 + t148 * t21;
t168 = t175 * mrSges(6,3);
t167 = (Ifges(5,4) - Ifges(4,5)) * t89 + (-Ifges(5,5) + Ifges(4,6)) * t88;
t90 = mrSges(6,1) * t178 - mrSges(6,2) * t181;
t3 = pkin(5) * t89 - pkin(10) * t60 + t5;
t4 = pkin(10) * t59 + t6;
t1 = -t151 * t4 + t154 * t3;
t2 = t151 * t3 + t154 * t4;
t165 = -t1 * t160 - t2 * t98;
t164 = t146 * t6 + t148 * t5;
t102 = t148 * t118;
t41 = pkin(5) * t152 + t102 + (pkin(10) * t155 - t96) * t146;
t47 = -pkin(10) * t178 + t51;
t12 = -t151 * t47 + t154 * t41;
t13 = t151 * t41 + t154 * t47;
t163 = -t12 * t160 - t13 * t98;
t105 = t184 * t146;
t106 = t184 * t148;
t53 = -t105 * t151 + t106 * t154;
t54 = t105 * t154 + t106 * t151;
t162 = -t160 * t53 - t54 * t98;
t50 = -t146 * t96 + t102;
t161 = t146 * t51 + t148 * t50;
t32 = qJ(4) * t179 - t37;
t62 = t89 * mrSges(5,1) - mrSges(5,2) * t179;
t24 = -pkin(4) * t88 - t32;
t157 = qJ(4) ^ 2;
t127 = pkin(5) * t146 + qJ(4);
t117 = Ifges(4,1) * t152 + Ifges(4,4) * t155;
t116 = Ifges(4,4) * t152 + Ifges(4,2) * t155;
t115 = -Ifges(5,2) * t152 - Ifges(5,6) * t155;
t114 = -Ifges(5,6) * t152 - Ifges(5,3) * t155;
t113 = -mrSges(4,1) * t155 + mrSges(4,2) * t152;
t112 = mrSges(5,2) * t155 - mrSges(5,3) * t152;
t111 = Ifges(6,1) * t148 - t183;
t110 = -Ifges(6,2) * t146 + t182;
t107 = -pkin(3) * t155 + t169;
t104 = -mrSges(6,2) * t152 - mrSges(6,3) * t178;
t103 = mrSges(6,1) * t152 + mrSges(6,3) * t181;
t87 = pkin(5) * t178 + t119;
t74 = Ifges(6,6) * t152 + (-Ifges(6,2) * t148 - t183) * t155;
t73 = Ifges(6,3) * t152 + (-Ifges(6,5) * t146 - Ifges(6,6) * t148) * t155;
t66 = mrSges(7,1) * t152 - mrSges(7,3) * t80;
t65 = -mrSges(7,2) * t152 + mrSges(7,3) * t79;
t64 = -mrSges(4,1) * t179 - mrSges(4,3) * t89;
t63 = mrSges(4,2) * t179 - mrSges(4,3) * t88;
t61 = mrSges(5,1) * t88 + mrSges(5,3) * t179;
t58 = -Ifges(7,1) * t160 + Ifges(7,4) * t98;
t57 = -Ifges(7,4) * t160 + Ifges(7,2) * t98;
t49 = -mrSges(5,2) * t88 - mrSges(5,3) * t89;
t48 = mrSges(4,1) * t88 + mrSges(4,2) * t89;
t45 = Ifges(4,1) * t89 - Ifges(4,4) * t88 - Ifges(4,5) * t179;
t44 = Ifges(4,4) * t89 - Ifges(4,2) * t88 - Ifges(4,6) * t179;
t43 = -Ifges(5,4) * t179 - Ifges(5,2) * t89 + Ifges(5,6) * t88;
t42 = -Ifges(5,5) * t179 - Ifges(5,6) * t89 + Ifges(5,3) * t88;
t40 = Ifges(7,1) * t80 + Ifges(7,4) * t79 + Ifges(7,5) * t152;
t39 = Ifges(7,4) * t80 + Ifges(7,2) * t79 + Ifges(7,6) * t152;
t35 = mrSges(6,1) * t89 - mrSges(6,3) * t60;
t34 = -mrSges(6,2) * t89 + mrSges(6,3) * t59;
t31 = pkin(3) * t88 + t159;
t20 = Ifges(6,1) * t60 + Ifges(6,4) * t59 + Ifges(6,5) * t89;
t19 = Ifges(6,4) * t60 + Ifges(6,2) * t59 + Ifges(6,6) * t89;
t18 = Ifges(6,5) * t60 + Ifges(6,6) * t59 + Ifges(6,3) * t89;
t15 = mrSges(7,1) * t89 - mrSges(7,3) * t29;
t14 = -mrSges(7,2) * t89 + mrSges(7,3) * t28;
t11 = -pkin(5) * t59 + t24;
t9 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t89;
t8 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t89;
t16 = [0.2e1 * t76 * t48 + t59 * t19 + t60 * t20 + 0.2e1 * t32 * t61 + 0.2e1 * t33 * t62 + 0.2e1 * t37 * t63 + 0.2e1 * t36 * t64 + 0.2e1 * t31 * t49 + t28 * t8 + t29 * t9 + 0.2e1 * t24 * t30 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 + 0.2e1 * t1 * t15 + 0.2e1 * t11 * t10 + 0.2e1 * t2 * t14 + m(3) * (pkin(1) ^ 2 * t147 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2 + t76 ^ 2) + m(6) * (t24 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + (t45 + t18 + t7 - t43) * t89 + (t42 - t44) * t88 + Ifges(2,3) + ((-0.2e1 * t91 * mrSges(3,3) + Ifges(3,5) * t149 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t153) * t147) * t153 + (0.2e1 * t92 * mrSges(3,3) + Ifges(3,6) * t149 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t153 + (Ifges(3,2) + t186) * t156) * t147 + t167) * t156) * t147 + (t172 - 0.2e1 * t187 + 0.2e1 * t188) * t149; -t187 + t188 + t31 * t112 + t76 * t113 + t119 * t30 + t5 * t103 + t6 * t104 + t107 * t49 + t24 * t90 + t87 * t10 + t2 * t65 + t1 * t66 - pkin(2) * t48 + t50 * t35 + t51 * t34 + t11 * t46 + t12 * t15 + t13 * t14 + (-t43 / 0.2e1 + t45 / 0.2e1 + t18 / 0.2e1 + t7 / 0.2e1 + Ifges(5,4) * t171 - t36 * mrSges(4,3) + t33 * mrSges(5,1) + (t62 - t64) * pkin(9)) * t152 + m(6) * (t119 * t24 + t5 * t50 + t51 * t6) + m(7) * (t1 * t12 + t11 * t87 + t13 * t2) + (-t115 / 0.2e1 + t117 / 0.2e1 + t73 / 0.2e1 + t38 / 0.2e1) * t89 + (t114 / 0.2e1 - t116 / 0.2e1) * t88 + m(5) * (t107 * t31 + (t152 * t33 - t155 * t32) * pkin(9)) + m(4) * (-pkin(2) * t76 + (-t152 * t36 + t155 * t37) * pkin(9)) + t9 * t194 + t8 * t195 + t60 * t196 + t74 * t197 + t40 * t198 + t39 * t199 + t172 - t177 * t179 / 0.2e1 + (-t42 / 0.2e1 + t44 / 0.2e1 + t19 * t190 + t20 * t191 + t37 * mrSges(4,3) - t32 * mrSges(5,1) + Ifges(5,5) * t171 + (-t61 + t63) * pkin(9)) * t155; -0.2e1 * pkin(2) * t113 + 0.2e1 * t50 * t103 + 0.2e1 * t51 * t104 + 0.2e1 * t107 * t112 + 0.2e1 * t119 * t90 + 0.2e1 * t12 * t66 + 0.2e1 * t13 * t65 + t79 * t39 + t80 * t40 + 0.2e1 * t87 * t46 + Ifges(3,3) + (-t146 * t75 - t148 * t74 - t114 + t116) * t155 + (-t115 + t117 + t73 + t38) * t152 + m(5) * (t107 ^ 2 + t176) + m(4) * (pkin(2) ^ 2 + t176) + m(6) * (t119 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(7) * (t12 ^ 2 + t13 ^ 2 + t87 ^ 2) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * pkin(9) * t201; (-t61 + t30) * qJ(4) + m(7) * (t1 * t53 + t11 * t127 + t2 * t54) + m(5) * (-pkin(3) * t33 - qJ(4) * t32) + t60 * t111 / 0.2e1 + t127 * t10 + t24 * t108 - pkin(3) * t62 + t53 * t15 + t54 * t14 + t11 * t55 + t36 * mrSges(4,1) - t37 * mrSges(4,2) - t32 * mrSges(5,3) + t33 * mrSges(5,2) + t173 * t89 - t165 * mrSges(7,3) + (t20 / 0.2e1 - t185 * t35 - t5 * mrSges(6,3)) * t148 + (-t19 / 0.2e1 - t185 * t34 - t6 * mrSges(6,3)) * t146 + m(6) * (qJ(4) * t24 - t164 * t185) - t167 + t9 * t192 + t8 * t193 + t110 * t197 + t58 * t198 + t57 * t199 - t186 * t179; t119 * t108 + t127 * t46 + t40 * t192 + qJ(4) * t90 + t39 * t193 + t87 * t55 + t57 * t195 + t58 * t194 + t54 * t65 + t53 * t66 - t163 * mrSges(7,3) + (-t50 * mrSges(6,3) - t103 * t185 + t196) * t148 + (-t74 / 0.2e1 - t185 * t104 - t51 * mrSges(6,3)) * t146 + m(6) * (qJ(4) * t119 - t161 * t185) + m(7) * (t12 * t53 + t127 * t87 + t13 * t54) + (qJ(4) * mrSges(5,1) + t110 * t190 + t111 * t191 - Ifges(5,5)) * t155 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t173) * t152 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t155 + (-mrSges(4,1) + t200) * t152) * pkin(9) + t177; -0.2e1 * pkin(3) * mrSges(5,2) - t160 * t58 - t146 * t110 + t148 * t111 + 0.2e1 * t127 * t55 + t98 * t57 + m(7) * (t127 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t175 * t185 ^ 2 + t157) + m(5) * (pkin(3) ^ 2 + t157) + t186 + 0.2e1 * (mrSges(5,3) + t108) * qJ(4) - 0.2e1 * t162 * mrSges(7,3) + 0.2e1 * t185 * t168; m(5) * t33 + m(6) * t164 + m(7) * t165 - t98 * t14 + t146 * t34 + t148 * t35 - t15 * t160 + t62; -t160 * t66 + t148 * t103 + t146 * t104 - t98 * t65 + (m(5) * pkin(9) + mrSges(5,1)) * t152 + m(7) * t163 + m(6) * t161; m(7) * t162 - mrSges(7,3) * t174 - t170 * t185 - t168 + t200; m(7) * t174 + m(5) + t170; m(6) * t24 + m(7) * t11 + t10 + t30; m(6) * t119 + m(7) * t87 + t46 + t90; m(6) * qJ(4) + m(7) * t127 + t108 + t55; 0; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t12 - mrSges(7,2) * t13 + t38; mrSges(7,1) * t53 - mrSges(7,2) * t54 + t56; -mrSges(7,1) * t160 + mrSges(7,2) * t98; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
