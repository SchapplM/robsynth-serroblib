% Calculate joint inertia matrix for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:02:47
% EndTime: 2018-11-23 17:02:49
% DurationCPUTime: 1.77s
% Computational Cost: add. (3739->386), mult. (8945->555), div. (0->0), fcn. (9933->12), ass. (0->150)
t193 = Ifges(3,3) + Ifges(4,3);
t141 = sin(pkin(11));
t128 = pkin(2) * t141 + pkin(9);
t192 = 0.2e1 * t128;
t144 = cos(pkin(11));
t142 = sin(pkin(6));
t151 = cos(qJ(2));
t164 = t142 * t151;
t148 = sin(qJ(2));
t165 = t142 * t148;
t88 = t141 * t165 - t144 * t164;
t89 = (t141 * t151 + t144 * t148) * t142;
t191 = Ifges(4,5) * t89 - Ifges(4,6) * t88;
t140 = sin(pkin(12));
t143 = cos(pkin(12));
t146 = sin(qJ(6));
t149 = cos(qJ(6));
t103 = t140 * t149 + t143 * t146;
t147 = sin(qJ(4));
t90 = t103 * t147;
t102 = -t140 * t146 + t143 * t149;
t91 = t102 * t147;
t55 = t90 * mrSges(7,1) + t91 * mrSges(7,2);
t94 = (pkin(5) * t140 + t128) * t147;
t163 = t143 * t147;
t166 = t140 * t147;
t95 = mrSges(6,1) * t166 + mrSges(6,2) * t163;
t190 = -m(7) * t94 - t55 - t95;
t106 = (-pkin(2) * t151 - pkin(1)) * t142;
t189 = 0.2e1 * t106;
t145 = cos(pkin(6));
t150 = cos(qJ(4));
t68 = t145 * t147 + t150 * t89;
t46 = -t140 * t68 + t143 * t88;
t47 = t140 * t88 + t143 * t68;
t67 = -t145 * t150 + t147 * t89;
t16 = Ifges(6,1) * t47 + Ifges(6,4) * t46 + Ifges(6,5) * t67;
t188 = t16 / 0.2e1;
t28 = -t146 * t47 + t149 * t46;
t187 = t28 / 0.2e1;
t29 = t146 * t46 + t149 * t47;
t186 = t29 / 0.2e1;
t170 = Ifges(6,4) * t143;
t86 = -Ifges(6,6) * t150 + (-Ifges(6,2) * t140 + t170) * t147;
t185 = t86 / 0.2e1;
t171 = Ifges(6,4) * t140;
t87 = -Ifges(6,5) * t150 + (Ifges(6,1) * t143 - t171) * t147;
t184 = t87 / 0.2e1;
t183 = -t90 / 0.2e1;
t182 = t91 / 0.2e1;
t181 = m(6) + m(7);
t180 = t102 / 0.2e1;
t179 = t103 / 0.2e1;
t113 = Ifges(6,1) * t140 + t170;
t178 = t113 / 0.2e1;
t177 = pkin(1) * t145;
t176 = pkin(4) * t150;
t125 = t151 * t177;
t96 = -pkin(8) * t165 + t125;
t175 = t96 * mrSges(3,1);
t97 = pkin(8) * t164 + t148 * t177;
t174 = t97 * mrSges(3,2);
t173 = pkin(10) + qJ(5);
t71 = pkin(2) * t145 + t125 + (-pkin(8) - qJ(3)) * t165;
t78 = qJ(3) * t164 + t97;
t44 = t141 * t71 + t144 * t78;
t41 = pkin(9) * t145 + t44;
t48 = pkin(3) * t88 - pkin(9) * t89 + t106;
t24 = t147 * t48 + t150 * t41;
t18 = qJ(5) * t88 + t24;
t43 = -t141 * t78 + t144 * t71;
t40 = -pkin(3) * t145 - t43;
t22 = pkin(4) * t67 - qJ(5) * t68 + t40;
t6 = t140 * t22 + t143 * t18;
t30 = -t46 * mrSges(6,1) + t47 * mrSges(6,2);
t50 = mrSges(5,1) * t88 - mrSges(5,3) * t68;
t172 = t30 - t50;
t109 = -t143 * mrSges(6,1) + t140 * mrSges(6,2);
t169 = -mrSges(5,1) + t109;
t63 = Ifges(7,5) * t103 + Ifges(7,6) * t102;
t129 = -pkin(2) * t144 - pkin(3);
t101 = -qJ(5) * t147 + t129 - t176;
t167 = t128 * t150;
t58 = t140 * t101 + t143 * t167;
t168 = t128 * t147;
t162 = Ifges(5,5) * t147 + Ifges(5,6) * t150;
t161 = t140 ^ 2 + t143 ^ 2;
t138 = t147 ^ 2;
t139 = t150 ^ 2;
t160 = t138 + t139;
t7 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t67;
t159 = Ifges(5,5) * t68 - Ifges(5,6) * t67 + Ifges(5,3) * t88;
t158 = t63 / 0.2e1 + Ifges(6,5) * t140 / 0.2e1 + Ifges(6,6) * t143 / 0.2e1;
t10 = -t28 * mrSges(7,1) + t29 * mrSges(7,2);
t62 = -t102 * mrSges(7,1) + t103 * mrSges(7,2);
t5 = -t140 * t18 + t143 * t22;
t23 = -t147 * t41 + t150 * t48;
t157 = qJ(5) * t161;
t51 = Ifges(7,5) * t91 - Ifges(7,6) * t90 - Ifges(7,3) * t150;
t156 = -t140 * t5 + t143 * t6;
t93 = t143 * t101;
t57 = -t140 * t167 + t93;
t155 = -t140 * t57 + t143 * t58;
t19 = -pkin(4) * t88 - t23;
t154 = Ifges(3,5) * t165 + Ifges(3,6) * t164 + t193 * t145 + t191;
t130 = -pkin(5) * t143 - pkin(4);
t126 = t128 ^ 2;
t117 = t138 * t126;
t116 = Ifges(5,1) * t147 + Ifges(5,4) * t150;
t115 = Ifges(5,4) * t147 + Ifges(5,2) * t150;
t114 = -t150 * mrSges(5,1) + t147 * mrSges(5,2);
t112 = Ifges(6,2) * t143 + t171;
t110 = t173 * t143;
t108 = t173 * t140;
t105 = -mrSges(6,1) * t150 - mrSges(6,3) * t163;
t104 = mrSges(6,2) * t150 - mrSges(6,3) * t166;
t85 = -Ifges(6,3) * t150 + (Ifges(6,5) * t143 - Ifges(6,6) * t140) * t147;
t79 = t89 * mrSges(4,2);
t77 = -mrSges(7,1) * t150 - mrSges(7,3) * t91;
t76 = mrSges(7,2) * t150 - mrSges(7,3) * t90;
t73 = mrSges(4,1) * t145 - mrSges(4,3) * t89;
t72 = -mrSges(4,2) * t145 - mrSges(4,3) * t88;
t70 = -t108 * t146 + t110 * t149;
t69 = -t108 * t149 - t110 * t146;
t65 = Ifges(7,1) * t103 + Ifges(7,4) * t102;
t64 = Ifges(7,4) * t103 + Ifges(7,2) * t102;
t56 = -pkin(10) * t166 + t58;
t54 = -pkin(10) * t163 + t93 + (-t128 * t140 - pkin(5)) * t150;
t53 = Ifges(7,1) * t91 - Ifges(7,4) * t90 - Ifges(7,5) * t150;
t52 = Ifges(7,4) * t91 - Ifges(7,2) * t90 - Ifges(7,6) * t150;
t49 = -mrSges(5,2) * t88 - mrSges(5,3) * t67;
t39 = mrSges(5,1) * t67 + mrSges(5,2) * t68;
t36 = Ifges(5,1) * t68 - Ifges(5,4) * t67 + Ifges(5,5) * t88;
t35 = Ifges(5,4) * t68 - Ifges(5,2) * t67 + Ifges(5,6) * t88;
t34 = t146 * t54 + t149 * t56;
t33 = -t146 * t56 + t149 * t54;
t32 = mrSges(6,1) * t67 - mrSges(6,3) * t47;
t31 = -mrSges(6,2) * t67 + mrSges(6,3) * t46;
t15 = Ifges(6,4) * t47 + Ifges(6,2) * t46 + Ifges(6,6) * t67;
t14 = Ifges(6,5) * t47 + Ifges(6,6) * t46 + Ifges(6,3) * t67;
t13 = mrSges(7,1) * t67 - mrSges(7,3) * t29;
t12 = -mrSges(7,2) * t67 + mrSges(7,3) * t28;
t11 = -pkin(5) * t46 + t19;
t9 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t67;
t8 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t67;
t4 = pkin(10) * t46 + t6;
t3 = pkin(5) * t67 - pkin(10) * t47 + t5;
t2 = t146 * t3 + t149 * t4;
t1 = -t146 * t4 + t149 * t3;
t17 = [(t154 - 0.2e1 * t174 + 0.2e1 * t175 + t191) * t145 + t47 * t16 + 0.2e1 * t24 * t49 + 0.2e1 * t23 * t50 + 0.2e1 * t40 * t39 + t46 * t15 + 0.2e1 * t19 * t30 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + t28 * t8 + t29 * t9 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + 0.2e1 * t11 * t10 + Ifges(4,1) * t89 ^ 2 + m(5) * (t23 ^ 2 + t24 ^ 2 + t40 ^ 2) + m(4) * (t106 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(6) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) + (t14 + t7 - t35) * t67 + ((Ifges(3,5) * t148 + Ifges(3,6) * t151) * t145 + 0.2e1 * (-t148 * t96 + t151 * t97) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t151 + mrSges(3,2) * t148) + t151 * (Ifges(3,4) * t148 + Ifges(3,2) * t151) + t148 * (Ifges(3,1) * t148 + Ifges(3,4) * t151)) * t142) * t142 + m(3) * (t96 ^ 2 + t97 ^ 2) + t68 * t36 + 0.2e1 * t44 * t72 + 0.2e1 * t43 * t73 + (mrSges(4,1) * t189 - 0.2e1 * Ifges(4,4) * t89 + Ifges(4,2) * t88 + t159) * t88 + Ifges(2,3) + t79 * t189; t154 + (t141 * t72 + t144 * t73 + m(4) * (t141 * t44 + t144 * t43)) * pkin(2) + t11 * t55 + t57 * t32 + t58 * t31 + t43 * mrSges(4,1) - t44 * mrSges(4,2) + t33 * t13 + t34 * t12 + (t24 * mrSges(5,3) + t128 * t49 + t35 / 0.2e1 - t14 / 0.2e1 - t7 / 0.2e1) * t150 + m(7) * (t1 * t33 + t11 * t94 + t2 * t34) + (t51 / 0.2e1 + t85 / 0.2e1 - t115 / 0.2e1) * t67 + m(6) * (t168 * t19 + t5 * t57 + t58 * t6) + t88 * t162 / 0.2e1 - t174 + t175 + t9 * t182 + t8 * t183 + t47 * t184 + t46 * t185 + t2 * t76 + t1 * t77 + (-t23 * mrSges(5,3) - t140 * t15 / 0.2e1 + t143 * t188 + t36 / 0.2e1 + t172 * t128) * t147 + t94 * t10 + t19 * t95 + t6 * t104 + t5 * t105 + t40 * t114 + t68 * t116 / 0.2e1 + t129 * t39 + t53 * t186 + t52 * t187 + m(5) * (t129 * t40 + (-t147 * t23 + t150 * t24) * t128); 0.2e1 * t58 * t104 + 0.2e1 * t57 * t105 + 0.2e1 * t129 * t114 + 0.2e1 * t33 * t77 + 0.2e1 * t34 * t76 - t90 * t52 + t91 * t53 + 0.2e1 * t94 * t55 + (-t51 - t85 + t115) * t150 + (-t140 * t86 + t143 * t87 + t192 * t95 + t116) * t147 + m(7) * (t33 ^ 2 + t34 ^ 2 + t94 ^ 2) + m(6) * (t57 ^ 2 + t58 ^ 2 + t117) + m(5) * (t126 * t139 + t129 ^ 2 + t117) + m(4) * (t141 ^ 2 + t144 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t144 - mrSges(4,2) * t141) * pkin(2) + t160 * mrSges(5,3) * t192 + t193; t88 * mrSges(4,1) + t91 * t12 - t90 * t13 + t79 + (-t10 - t172) * t150 + (-t140 * t32 + t143 * t31 + t49) * t147 + m(7) * (-t1 * t90 - t11 * t150 + t2 * t91) + m(6) * (t147 * t156 - t150 * t19) + m(5) * (t147 * t24 + t150 * t23) + m(4) * t106; m(7) * (-t33 * t90 + t34 * t91) + t91 * t76 - t90 * t77 + (-t140 * t105 + t143 * t104 + m(6) * (t155 - t167)) * t147 + t190 * t150; m(4) + m(5) * t160 + m(6) * (t138 * t161 + t139) + m(7) * (t90 ^ 2 + t91 ^ 2 + t139); (t6 * mrSges(6,3) + qJ(5) * t31 + t15 / 0.2e1) * t143 + (-t5 * mrSges(6,3) - qJ(5) * t32 + t188) * t140 + m(7) * (t1 * t69 + t11 * t130 + t2 * t70) + t11 * t62 + t64 * t187 + t65 * t186 - pkin(4) * t30 + t23 * mrSges(5,1) - t24 * mrSges(5,2) + (-t1 * t103 + t102 * t2) * mrSges(7,3) + t158 * t67 + m(6) * (-pkin(4) * t19 + qJ(5) * t156) + t159 + t69 * t13 + t70 * t12 + t8 * t180 + t9 * t179 + t19 * t109 + t46 * t112 / 0.2e1 + t47 * t178 + t130 * t10; t70 * t76 + t69 * t77 + t64 * t183 + t65 * t182 + t94 * t62 - pkin(4) * t95 + t52 * t180 + t53 * t179 + t130 * t55 + t169 * t168 + (t102 * t34 - t103 * t33) * mrSges(7,3) + (t58 * mrSges(6,3) + qJ(5) * t104 + t147 * t178 + t185) * t143 + (-t57 * mrSges(6,3) - qJ(5) * t105 - t147 * t112 / 0.2e1 + t184) * t140 + m(7) * (t130 * t94 + t33 * t69 + t34 * t70) + m(6) * (-pkin(4) * t168 + qJ(5) * t155) + (-t128 * mrSges(5,2) - t158) * t150 + t162; (t102 * t91 + t103 * t90) * mrSges(7,3) + (-t62 - t169) * t150 + (mrSges(6,3) * t161 - mrSges(5,2)) * t147 + m(6) * (t147 * t157 + t176) + m(7) * (-t130 * t150 - t69 * t90 + t70 * t91); -0.2e1 * pkin(4) * t109 + t102 * t64 + t103 * t65 + t143 * t112 + t140 * t113 + 0.2e1 * t130 * t62 + Ifges(5,3) + m(7) * (t130 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (qJ(5) ^ 2 * t161 + pkin(4) ^ 2) + 0.2e1 * (t102 * t70 - t103 * t69) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t157; m(6) * t19 + m(7) * t11 + t10 + t30; m(6) * t168 - t190; -t181 * t150; -m(6) * pkin(4) + m(7) * t130 + t109 + t62; t181; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t33 - mrSges(7,2) * t34 + t51; -t55; mrSges(7,1) * t69 - t70 * mrSges(7,2) + t63; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
