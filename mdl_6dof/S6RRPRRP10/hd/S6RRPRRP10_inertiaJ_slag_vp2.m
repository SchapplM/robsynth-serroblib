% Calculate joint inertia matrix for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:36
% EndTime: 2019-03-09 12:36:41
% DurationCPUTime: 1.86s
% Computational Cost: add. (3275->383), mult. (7304->530), div. (0->0), fcn. (8140->10), ass. (0->139)
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t194 = t144 ^ 2 + t147 ^ 2;
t193 = 0.2e1 * t194;
t141 = sin(pkin(6));
t148 = cos(qJ(2));
t170 = t141 * t148;
t145 = sin(qJ(4));
t187 = cos(qJ(4));
t140 = sin(pkin(11));
t142 = cos(pkin(11));
t143 = cos(pkin(6));
t146 = sin(qJ(2));
t171 = t141 * t146;
t95 = -t140 * t171 + t142 * t143;
t96 = t140 * t143 + t142 * t171;
t63 = t145 * t95 + t187 * t96;
t43 = t144 * t63 + t147 * t170;
t44 = -t144 * t170 + t147 * t63;
t62 = t145 * t96 - t187 * t95;
t10 = Ifges(6,5) * t44 - Ifges(6,6) * t43 + Ifges(6,3) * t62;
t11 = Ifges(7,4) * t44 + Ifges(7,2) * t62 + Ifges(7,6) * t43;
t192 = t10 + t11;
t105 = t140 * t145 - t187 * t142;
t106 = t187 * t140 + t145 * t142;
t172 = t106 * t147;
t173 = t106 * t144;
t46 = Ifges(6,5) * t172 - Ifges(6,6) * t173 + Ifges(6,3) * t105;
t47 = Ifges(7,4) * t172 + Ifges(7,2) * t105 + Ifges(7,6) * t173;
t191 = t46 + t47;
t183 = pkin(9) + qJ(3);
t109 = t183 * t142;
t157 = t183 * t140;
t78 = t109 * t145 + t187 * t157;
t190 = t78 ^ 2;
t189 = 0.2e1 * t78;
t123 = pkin(8) * t171;
t186 = pkin(1) * t148;
t89 = t123 + (-pkin(2) - t186) * t143;
t65 = -pkin(3) * t95 + t89;
t18 = pkin(4) * t62 - pkin(10) * t63 + t65;
t98 = t143 * t146 * pkin(1) + pkin(8) * t170;
t86 = qJ(3) * t143 + t98;
t87 = (-pkin(2) * t148 - qJ(3) * t146 - pkin(1)) * t141;
t53 = -t140 * t86 + t142 * t87;
t31 = -pkin(3) * t170 - pkin(9) * t96 + t53;
t54 = t140 * t87 + t142 * t86;
t36 = pkin(9) * t95 + t54;
t16 = t145 * t31 + t187 * t36;
t8 = -pkin(10) * t170 + t16;
t4 = t144 * t18 + t147 * t8;
t97 = t143 * t186 - t123;
t185 = t97 * mrSges(3,1);
t184 = t98 * mrSges(3,2);
t21 = -mrSges(6,2) * t62 - mrSges(6,3) * t43;
t22 = -mrSges(7,2) * t43 + mrSges(7,3) * t62;
t182 = t21 + t22;
t23 = mrSges(6,1) * t62 - mrSges(6,3) * t44;
t24 = -t62 * mrSges(7,1) + t44 * mrSges(7,2);
t181 = t23 - t24;
t180 = -Ifges(5,5) * t63 + Ifges(5,6) * t62;
t70 = -mrSges(6,2) * t105 - mrSges(6,3) * t173;
t73 = -mrSges(7,2) * t173 + mrSges(7,3) * t105;
t179 = t70 + t73;
t71 = mrSges(6,1) * t105 - mrSges(6,3) * t172;
t72 = -t105 * mrSges(7,1) + mrSges(7,2) * t172;
t178 = t71 - t72;
t128 = -pkin(3) * t142 - pkin(2);
t69 = pkin(4) * t105 - pkin(10) * t106 + t128;
t80 = t187 * t109 - t145 * t157;
t33 = t144 * t69 + t147 * t80;
t177 = Ifges(6,4) * t144;
t176 = Ifges(6,4) * t147;
t175 = Ifges(7,5) * t144;
t174 = Ifges(7,5) * t147;
t169 = Ifges(5,5) * t106 - Ifges(5,6) * t105;
t115 = Ifges(6,5) * t144 + Ifges(6,6) * t147;
t168 = t194 * pkin(10) ^ 2;
t167 = t140 ^ 2 + t142 ^ 2;
t12 = Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t62;
t9 = Ifges(7,5) * t44 + Ifges(7,6) * t62 + Ifges(7,3) * t43;
t165 = t9 / 0.2e1 - t12 / 0.2e1;
t13 = Ifges(7,1) * t44 + Ifges(7,4) * t62 + Ifges(7,5) * t43;
t14 = Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t62;
t164 = t13 / 0.2e1 + t14 / 0.2e1;
t45 = Ifges(7,6) * t105 + (Ifges(7,3) * t144 + t174) * t106;
t48 = Ifges(6,6) * t105 + (-Ifges(6,2) * t144 + t176) * t106;
t163 = -t48 / 0.2e1 + t45 / 0.2e1;
t49 = Ifges(7,4) * t105 + (Ifges(7,1) * t147 + t175) * t106;
t50 = Ifges(6,5) * t105 + (Ifges(6,1) * t147 - t177) * t106;
t162 = t49 / 0.2e1 + t50 / 0.2e1;
t161 = Ifges(3,5) * t171 + Ifges(3,6) * t170 + Ifges(3,3) * t143;
t114 = -Ifges(7,3) * t147 + t175;
t117 = Ifges(6,2) * t147 + t177;
t160 = t114 / 0.2e1 - t117 / 0.2e1;
t116 = Ifges(7,4) * t144 - Ifges(7,6) * t147;
t159 = t115 / 0.2e1 + t116 / 0.2e1;
t118 = Ifges(7,1) * t144 - t174;
t119 = Ifges(6,1) * t144 + t176;
t158 = t118 / 0.2e1 + t119 / 0.2e1;
t66 = -t95 * mrSges(4,1) + t96 * mrSges(4,2);
t29 = t62 * mrSges(5,1) + t63 * mrSges(5,2);
t75 = t105 * mrSges(5,1) + t106 * mrSges(5,2);
t108 = -t142 * mrSges(4,1) + t140 * mrSges(4,2);
t15 = -t145 * t36 + t187 * t31;
t3 = -t144 * t8 + t147 * t18;
t113 = -t147 * mrSges(6,1) + t144 * mrSges(6,2);
t155 = t144 * mrSges(6,1) + t147 * mrSges(6,2);
t112 = -t147 * mrSges(7,1) - t144 * mrSges(7,3);
t154 = t144 * mrSges(7,1) - t147 * mrSges(7,3);
t153 = pkin(5) * t147 + qJ(6) * t144;
t152 = pkin(5) * t144 - qJ(6) * t147;
t151 = -t53 * t140 + t54 * t142;
t32 = -t144 * t80 + t147 * t69;
t7 = pkin(4) * t170 - t15;
t111 = Ifges(4,1) * t140 + Ifges(4,4) * t142;
t110 = Ifges(4,4) * t140 + Ifges(4,2) * t142;
t107 = -pkin(4) - t153;
t82 = -mrSges(4,1) * t170 - mrSges(4,3) * t96;
t81 = mrSges(4,2) * t170 + mrSges(4,3) * t95;
t77 = Ifges(5,1) * t106 - Ifges(5,4) * t105;
t76 = Ifges(5,4) * t106 - Ifges(5,2) * t105;
t68 = t155 * t106;
t67 = t154 * t106;
t56 = Ifges(4,1) * t96 + Ifges(4,4) * t95 - Ifges(4,5) * t170;
t55 = Ifges(4,4) * t96 + Ifges(4,2) * t95 - Ifges(4,6) * t170;
t52 = -mrSges(5,1) * t170 - mrSges(5,3) * t63;
t51 = mrSges(5,2) * t170 - mrSges(5,3) * t62;
t37 = t106 * t152 + t78;
t28 = -pkin(5) * t105 - t32;
t27 = qJ(6) * t105 + t33;
t26 = Ifges(5,1) * t63 - Ifges(5,4) * t62 - Ifges(5,5) * t170;
t25 = Ifges(5,4) * t63 - Ifges(5,2) * t62 - Ifges(5,6) * t170;
t20 = mrSges(6,1) * t43 + mrSges(6,2) * t44;
t19 = mrSges(7,1) * t43 - mrSges(7,3) * t44;
t5 = t43 * pkin(5) - t44 * qJ(6) + t7;
t2 = -pkin(5) * t62 - t3;
t1 = qJ(6) * t62 + t4;
t6 = [((-0.2e1 * t97 * mrSges(3,3) + Ifges(3,5) * t143 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t146) * t141) * t146 + (0.2e1 * t98 * mrSges(3,3) - Ifges(4,5) * t96 + Ifges(3,6) * t143 - Ifges(4,6) * t95 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t146 + (Ifges(4,3) + Ifges(5,3) + Ifges(3,2)) * t148) * t141 + t180) * t148) * t141 + (t161 - 0.2e1 * t184 + 0.2e1 * t185) * t143 + 0.2e1 * t89 * t66 + t95 * t55 + t96 * t56 + 0.2e1 * t54 * t81 + 0.2e1 * t53 * t82 + t63 * t26 + 0.2e1 * t65 * t29 + 0.2e1 * t16 * t51 + 0.2e1 * t15 * t52 + 0.2e1 * t4 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t3 * t23 + 0.2e1 * t2 * t24 + 0.2e1 * t5 * t19 + 0.2e1 * t7 * t20 + (t13 + t14) * t44 + (t9 - t12) * t43 + (-t25 + t192) * t62 + m(4) * (t53 ^ 2 + t54 ^ 2 + t89 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t65 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (pkin(1) ^ 2 * t141 ^ 2 + t97 ^ 2 + t98 ^ 2) + Ifges(2,3); -t184 + t185 + (t20 - t52) * t78 + m(5) * (t128 * t65 - t15 * t78 + t16 * t80) + m(6) * (t3 * t32 + t33 * t4 + t7 * t78) + m(7) * (t1 * t27 + t2 * t28 + t37 * t5) + (t26 / 0.2e1 - t15 * mrSges(5,3) + t164 * t147 + t165 * t144) * t106 + t162 * t44 + t163 * t43 + t151 * mrSges(4,3) + m(4) * (-pkin(2) * t89 + qJ(3) * t151) + (-t76 / 0.2e1 + t46 / 0.2e1 + t47 / 0.2e1) * t62 + t128 * t29 + t89 * t108 + t95 * t110 / 0.2e1 + t96 * t111 / 0.2e1 + t4 * t70 + t3 * t71 + t2 * t72 + t1 * t73 + t65 * t75 + t63 * t77 / 0.2e1 + t80 * t51 - pkin(2) * t66 + t5 * t67 + t7 * t68 + t32 * t23 + t33 * t21 + t37 * t19 + t27 * t22 + t28 * t24 + (-t140 * t82 + t142 * t81) * qJ(3) + (t11 / 0.2e1 + t10 / 0.2e1 - t25 / 0.2e1 - t16 * mrSges(5,3)) * t105 - (Ifges(4,5) * t140 + Ifges(4,6) * t142 + t169) * t170 / 0.2e1 + t161 + t140 * t56 / 0.2e1 + t142 * t55 / 0.2e1; -0.2e1 * pkin(2) * t108 + t142 * t110 + t140 * t111 + 0.2e1 * t128 * t75 + 0.2e1 * t27 * t73 + 0.2e1 * t28 * t72 + 0.2e1 * t32 * t71 + 0.2e1 * t33 * t70 + 0.2e1 * t37 * t67 + t68 * t189 + Ifges(3,3) + 0.2e1 * t167 * qJ(3) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t80 + t191 - t76) * t105 + m(5) * (t128 ^ 2 + t80 ^ 2 + t190) + m(6) * (t32 ^ 2 + t33 ^ 2 + t190) + m(7) * (t27 ^ 2 + t28 ^ 2 + t37 ^ 2) + m(4) * (t167 * qJ(3) ^ 2 + pkin(2) ^ 2) + (mrSges(5,3) * t189 + t77 + (t49 + t50) * t147 + (t45 - t48) * t144) * t106; t181 * t147 + t182 * t144 + m(6) * (t144 * t4 + t147 * t3) + m(7) * (t1 * t144 - t147 * t2) + m(5) * t65 + m(4) * t89 + t29 + t66; -m(4) * pkin(2) + t178 * t147 + t179 * t144 + m(7) * (t144 * t27 - t147 * t28) + m(6) * (t144 * t33 + t147 * t32) + m(5) * t128 + t108 + t75; m(4) + m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t193; -Ifges(5,3) * t170 + t15 * mrSges(5,1) - t16 * mrSges(5,2) - pkin(4) * t20 + t107 * t19 + t5 * t112 + t7 * t113 + t159 * t62 + t158 * t44 + t160 * t43 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + t182 * pkin(10) - t165) * t147 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) - t181 * pkin(10) + t164) * t144 + m(7) * (t107 * t5 + (t1 * t147 + t144 * t2) * pkin(10)) + m(6) * (-pkin(4) * t7 + (-t144 * t3 + t147 * t4) * pkin(10)) - t180; -t80 * mrSges(5,2) - pkin(4) * t68 + t107 * t67 + t37 * t112 + (t113 - mrSges(5,1)) * t78 + t159 * t105 + (t27 * mrSges(7,2) + t33 * mrSges(6,3) + t179 * pkin(10) - t163) * t147 + (t28 * mrSges(7,2) - t32 * mrSges(6,3) - t178 * pkin(10) + t162) * t144 + m(7) * (t107 * t37 + (t144 * t28 + t147 * t27) * pkin(10)) + m(6) * (-pkin(4) * t78 + (-t144 * t32 + t147 * t33) * pkin(10)) + (t144 * t160 + t147 * t158) * t106 + t169; 0; -0.2e1 * pkin(4) * t113 + 0.2e1 * t107 * t112 + Ifges(5,3) + (t117 - t114) * t147 + (t119 + t118) * t144 + m(7) * (t107 ^ 2 + t168) + m(6) * (pkin(4) ^ 2 + t168) + (mrSges(7,2) + mrSges(6,3)) * pkin(10) * t193; t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) - pkin(5) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t22 + t1 * mrSges(7,3) + t192; -pkin(5) * t72 + m(7) * (-pkin(5) * t28 + qJ(6) * t27) + qJ(6) * t73 + t27 * mrSges(7,3) + t32 * mrSges(6,1) - t28 * mrSges(7,1) - t33 * mrSges(6,2) + t191; m(7) * t153 - t112 - t113; -t152 * mrSges(7,2) + (-m(7) * t152 - t154 - t155) * pkin(10) + t116 + t115; Ifges(7,2) + Ifges(6,3) + 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2); m(7) * t2 + t24; m(7) * t28 + t72; -m(7) * t147; (m(7) * pkin(10) + mrSges(7,2)) * t144; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
