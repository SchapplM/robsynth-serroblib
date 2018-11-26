% Calculate joint inertia matrix for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:11
% EndTime: 2018-11-23 17:29:12
% DurationCPUTime: 1.82s
% Computational Cost: add. (2878->353), mult. (6074->502), div. (0->0), fcn. (6435->10), ass. (0->139)
t142 = sin(qJ(6));
t146 = cos(qJ(6));
t175 = t142 ^ 2 + t146 ^ 2;
t214 = mrSges(7,3) * t175;
t213 = Ifges(4,1) + Ifges(3,3);
t106 = -mrSges(7,1) * t146 + mrSges(7,2) * t142;
t212 = -mrSges(6,1) + t106;
t143 = sin(qJ(5));
t144 = sin(qJ(4));
t147 = cos(qJ(5));
t148 = cos(qJ(4));
t102 = t143 * t144 - t147 * t148;
t103 = t143 * t148 + t144 * t147;
t123 = t144 * pkin(4) + qJ(3);
t58 = pkin(5) * t103 + pkin(11) * t102 + t123;
t150 = -pkin(2) - pkin(9);
t195 = -pkin(10) + t150;
t105 = t195 * t144;
t169 = t195 * t148;
t69 = t147 * t105 + t143 * t169;
t27 = -t142 * t69 + t146 * t58;
t28 = t142 * t58 + t146 * t69;
t161 = -t142 * t27 + t146 * t28;
t140 = sin(pkin(6));
t149 = cos(qJ(2));
t176 = t140 * t149;
t141 = cos(pkin(6));
t145 = sin(qJ(2));
t90 = t141 * t145 * pkin(1) + pkin(8) * t176;
t78 = -t141 * qJ(3) - t90;
t70 = pkin(3) * t176 - t78;
t87 = -t141 * t144 - t148 * t176;
t45 = -pkin(4) * t87 + t70;
t88 = t141 * t148 - t144 * t176;
t51 = t143 * t88 - t147 * t87;
t52 = t143 * t87 + t147 * t88;
t15 = pkin(5) * t51 - pkin(11) * t52 + t45;
t177 = t140 * t145;
t118 = pkin(8) * t177;
t200 = pkin(1) * t149;
t171 = -pkin(2) - t200;
t61 = pkin(3) * t177 + t118 + (-pkin(9) + t171) * t141;
t167 = -qJ(3) * t145 - pkin(1);
t71 = (t149 * t150 + t167) * t140;
t29 = -t144 * t71 + t148 * t61;
t20 = pkin(4) * t177 - pkin(10) * t88 + t29;
t30 = t144 * t61 + t148 * t71;
t24 = pkin(10) * t87 + t30;
t9 = t143 * t20 + t147 * t24;
t6 = pkin(11) * t177 + t9;
t2 = -t142 * t6 + t146 * t15;
t3 = t142 * t15 + t146 * t6;
t163 = -t142 * t2 + t146 * t3;
t181 = t102 * t142;
t59 = -mrSges(7,2) * t103 + mrSges(7,3) * t181;
t180 = t102 * t146;
t60 = mrSges(7,1) * t103 + mrSges(7,3) * t180;
t211 = -t142 * t60 + t146 * t59;
t37 = -t142 * t52 + t146 * t177;
t17 = -mrSges(7,2) * t51 + mrSges(7,3) * t37;
t38 = t142 * t177 + t146 * t52;
t18 = mrSges(7,1) * t51 - mrSges(7,3) * t38;
t210 = -t142 * t18 + t146 * t17;
t209 = Ifges(5,5) * t88 + Ifges(5,6) * t87 + Ifges(5,3) * t177;
t67 = t105 * t143 - t147 * t169;
t208 = t67 ^ 2;
t207 = -2 * mrSges(6,3);
t97 = t102 ^ 2;
t206 = t38 / 0.2e1;
t108 = Ifges(7,5) * t142 + Ifges(7,6) * t146;
t205 = t108 / 0.2e1;
t191 = Ifges(7,4) * t142;
t109 = Ifges(7,2) * t146 + t191;
t204 = t109 / 0.2e1;
t203 = t142 / 0.2e1;
t201 = t146 / 0.2e1;
t89 = t141 * t200 - t118;
t197 = t89 * mrSges(3,1);
t196 = t90 * mrSges(3,2);
t16 = -mrSges(7,1) * t37 + mrSges(7,2) * t38;
t44 = mrSges(6,1) * t177 - mrSges(6,3) * t52;
t194 = t16 - t44;
t193 = -Ifges(6,5) * t102 - Ifges(6,6) * t103;
t190 = Ifges(7,4) * t146;
t189 = t102 * t67;
t182 = t148 * mrSges(5,1);
t125 = pkin(4) * t143 + pkin(11);
t179 = t125 * t142;
t178 = t125 * t146;
t99 = mrSges(4,1) * t177 + t141 * mrSges(4,2);
t174 = t144 ^ 2 + t148 ^ 2;
t12 = Ifges(7,5) * t38 + Ifges(7,6) * t37 + Ifges(7,3) * t51;
t173 = Ifges(6,5) * t52 - Ifges(6,6) * t51 + Ifges(6,3) * t177;
t111 = Ifges(7,1) * t142 + t190;
t172 = t146 * t109 + t142 * t111 + Ifges(6,3);
t170 = m(5) * t174;
t168 = t175 * pkin(11);
t166 = t174 * mrSges(5,3);
t165 = t175 * t125;
t164 = Ifges(3,5) * t177 + Ifges(3,6) * t176 + t141 * t213;
t162 = -mrSges(7,1) * t142 - mrSges(7,2) * t146;
t8 = -t143 * t24 + t147 * t20;
t160 = t144 * t30 + t148 * t29;
t159 = 0.2e1 * t214;
t158 = t102 * t147 - t103 * t143;
t40 = -Ifges(7,5) * t180 + Ifges(7,6) * t181 + Ifges(7,3) * t103;
t157 = t212 * t102 + (-mrSges(6,2) + t214) * t103;
t156 = (mrSges(6,1) * t147 - mrSges(6,2) * t143) * pkin(4);
t41 = Ifges(7,6) * t103 + (Ifges(7,2) * t142 - t190) * t102;
t42 = Ifges(7,5) * t103 + (-Ifges(7,1) * t146 + t191) * t102;
t155 = -t69 * mrSges(6,2) + t181 * t204 + t41 * t201 + t42 * t203 + t193 - t111 * t180 / 0.2e1 + t103 * t205 + t212 * t67 + t161 * mrSges(7,3);
t13 = Ifges(7,4) * t38 + Ifges(7,2) * t37 + Ifges(7,6) * t51;
t14 = Ifges(7,1) * t38 + Ifges(7,4) * t37 + Ifges(7,5) * t51;
t5 = -pkin(5) * t177 - t8;
t154 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + mrSges(7,3) * t163 + t5 * t106 + t111 * t206 + t13 * t201 + t14 * t203 + t37 * t204 + t51 * t205 + t173;
t151 = qJ(3) ^ 2;
t133 = Ifges(5,5) * t148;
t126 = -pkin(4) * t147 - pkin(5);
t112 = Ifges(5,1) * t148 - Ifges(5,4) * t144;
t110 = Ifges(5,4) * t148 - Ifges(5,2) * t144;
t107 = mrSges(5,1) * t144 + mrSges(5,2) * t148;
t98 = -mrSges(4,1) * t176 - mrSges(4,3) * t141;
t96 = t103 ^ 2;
t80 = t141 * t171 + t118;
t79 = (-pkin(2) * t149 + t167) * t140;
t75 = mrSges(5,1) * t177 - mrSges(5,3) * t88;
t74 = -mrSges(5,2) * t177 + mrSges(5,3) * t87;
t66 = -Ifges(6,1) * t102 - Ifges(6,4) * t103;
t65 = -Ifges(6,4) * t102 - Ifges(6,2) * t103;
t64 = mrSges(6,1) * t103 - mrSges(6,2) * t102;
t55 = t162 * t102;
t54 = -mrSges(5,1) * t87 + mrSges(5,2) * t88;
t47 = Ifges(5,1) * t88 + Ifges(5,4) * t87 + Ifges(5,5) * t177;
t46 = Ifges(5,4) * t88 + Ifges(5,2) * t87 + Ifges(5,6) * t177;
t43 = -mrSges(6,2) * t177 - mrSges(6,3) * t51;
t25 = mrSges(6,1) * t51 + mrSges(6,2) * t52;
t22 = Ifges(6,1) * t52 - Ifges(6,4) * t51 + Ifges(6,5) * t177;
t21 = Ifges(6,4) * t52 - Ifges(6,2) * t51 + Ifges(6,6) * t177;
t1 = [(m(3) * pkin(1) ^ 2 * t140 + (0.2e1 * t79 * mrSges(4,2) + 0.2e1 * t90 * mrSges(3,3) + (-(2 * Ifges(4,5)) + Ifges(3,6)) * t141 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t149) * t140) * t149 + (-0.2e1 * t89 * mrSges(3,3) - 0.2e1 * t79 * mrSges(4,3) + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t141 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t145 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t149) * t140 + t173 + t209) * t145) * t140 + 0.2e1 * t78 * t98 + 0.2e1 * t80 * t99 + t87 * t46 + t88 * t47 + 0.2e1 * t70 * t54 + 0.2e1 * t30 * t74 + 0.2e1 * t29 * t75 + t52 * t22 + t37 * t13 + t38 * t14 + 0.2e1 * t9 * t43 + 0.2e1 * t8 * t44 + 0.2e1 * t45 * t25 + 0.2e1 * t5 * t16 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + m(3) * (t89 ^ 2 + t90 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t70 ^ 2) + m(4) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t45 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + (t12 - t21) * t51 + Ifges(2,3) + (t164 - 0.2e1 * t196 + 0.2e1 * t197) * t141; (-t98 + t54) * qJ(3) + m(6) * (t123 * t45 - t67 * t8 + t69 * t9) + m(4) * (-pkin(2) * t80 - qJ(3) * t78) + m(7) * (t2 * t27 + t28 * t3 + t5 * t67) + t70 * t107 + t87 * t110 / 0.2e1 + t88 * t112 / 0.2e1 + t123 * t25 - pkin(2) * t99 + t69 * t43 - t78 * mrSges(4,3) + t80 * mrSges(4,2) + t5 * t55 + t3 * t59 + t2 * t60 + t45 * t64 + t52 * t66 / 0.2e1 + t164 + (-t65 / 0.2e1 + t40 / 0.2e1) * t51 + t37 * t41 / 0.2e1 + t27 * t18 + t28 * t17 + (-t21 / 0.2e1 + t12 / 0.2e1 - t9 * mrSges(6,3)) * t103 + (t150 * t74 - t46 / 0.2e1 - t30 * mrSges(5,3)) * t144 + t197 + m(5) * (qJ(3) * t70 + t150 * t160) - t196 + (-Ifges(4,4) * t145 - Ifges(4,5) * t149 + (-Ifges(5,6) * t144 + t133 + t193) * t145 / 0.2e1) * t140 + (-t29 * mrSges(5,3) + t150 * t75 + t47 / 0.2e1) * t148 + t194 * t67 + (-t22 / 0.2e1 - t146 * t14 / 0.2e1 + t13 * t203 + t8 * mrSges(6,3)) * t102 + t42 * t206; -0.2e1 * pkin(2) * mrSges(4,2) - t144 * t110 + t148 * t112 + 0.2e1 * t123 * t64 + 0.2e1 * t27 * t60 + 0.2e1 * t28 * t59 + 0.2e1 * t67 * t55 + (t207 * t69 + t40 - t65) * t103 + (t142 * t41 - t146 * t42 + t207 * t67 - t66) * t102 + m(7) * (t27 ^ 2 + t28 ^ 2 + t208) + m(6) * (t123 ^ 2 + t69 ^ 2 + t208) + m(5) * (t150 ^ 2 * t174 + t151) + m(4) * (pkin(2) ^ 2 + t151) + 0.2e1 * (t107 + mrSges(4,3)) * qJ(3) - 0.2e1 * t150 * t166 + t213; t144 * t74 + t148 * t75 + t194 * t102 + (t43 + t210) * t103 + m(7) * (t102 * t5 + t103 * t163) + m(6) * (-t102 * t8 + t103 * t9) + m(5) * t160 + m(4) * t80 + t99; -m(4) * pkin(2) - t97 * mrSges(6,3) + t102 * t55 + mrSges(4,2) - t166 + (-mrSges(6,3) * t103 + t211) * t103 + m(7) * (t103 * t161 + t189) + m(6) * (t103 * t69 + t189) + t150 * t170; m(4) + t170 + m(6) * (t96 + t97) + m(7) * (t175 * t96 + t97); t126 * t16 + m(7) * (t125 * t163 + t126 * t5) + (m(6) * (t143 * t9 + t147 * t8) + t147 * t44 + t143 * t43) * pkin(4) + t17 * t178 - t18 * t179 + t29 * mrSges(5,1) - t30 * mrSges(5,2) + t154 + t209; t150 * t182 + t133 + t126 * t55 + m(7) * (t125 * t161 + t126 * t67) + (m(6) * (t143 * t69 - t147 * t67) + t158 * mrSges(6,3)) * pkin(4) + t59 * t178 - t60 * t179 + t155 + (-mrSges(5,2) * t150 - Ifges(5,6)) * t144; t182 - t144 * mrSges(5,2) + m(7) * (t126 * t102 + t103 * t165) - m(6) * t158 * pkin(4) + t157; 0.2e1 * t126 * t106 + Ifges(5,3) + 0.2e1 * t156 + t125 * t159 + m(7) * (t125 ^ 2 * t175 + t126 ^ 2) + m(6) * (t143 ^ 2 + t147 ^ 2) * pkin(4) ^ 2 + t172; t154 + (-m(7) * t5 - t16) * pkin(5) + (m(7) * t163 + t210) * pkin(11); t155 + (-m(7) * t67 - t55) * pkin(5) + (m(7) * t161 + t211) * pkin(11); m(7) * (-pkin(5) * t102 + t103 * t168) + t157; m(7) * (-pkin(5) * t126 + pkin(11) * t165) + (t126 - pkin(5)) * t106 + t156 + (t165 + t168) * mrSges(7,3) + t172; -0.2e1 * pkin(5) * t106 + m(7) * (pkin(11) ^ 2 * t175 + pkin(5) ^ 2) + pkin(11) * t159 + t172; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t27 - mrSges(7,2) * t28 + t40; t162 * t103; t125 * t162 + t108; pkin(11) * t162 + t108; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
