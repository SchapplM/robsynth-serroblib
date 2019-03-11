% Calculate joint inertia matrix for
% S6RRPRRP9
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:34
% EndTime: 2019-03-09 12:27:38
% DurationCPUTime: 1.69s
% Computational Cost: add. (3234->387), mult. (7225->525), div. (0->0), fcn. (8074->10), ass. (0->140)
t195 = -2 * mrSges(7,3);
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t166 = t148 ^ 2 + t151 ^ 2;
t194 = 0.2e1 * t166;
t145 = sin(pkin(6));
t152 = cos(qJ(2));
t169 = t145 * t152;
t149 = sin(qJ(4));
t186 = cos(qJ(4));
t144 = sin(pkin(11));
t146 = cos(pkin(11));
t147 = cos(pkin(6));
t150 = sin(qJ(2));
t170 = t145 * t150;
t97 = -t144 * t170 + t146 * t147;
t98 = t144 * t147 + t146 * t170;
t64 = t149 * t97 + t186 * t98;
t43 = -t148 * t64 - t151 * t169;
t44 = -t148 * t169 + t151 * t64;
t63 = t149 * t98 - t186 * t97;
t10 = Ifges(6,5) * t44 + Ifges(6,6) * t43 + Ifges(6,3) * t63;
t9 = Ifges(7,5) * t44 + Ifges(7,6) * t43 + Ifges(7,3) * t63;
t193 = t10 + t9;
t180 = mrSges(6,2) * t148;
t115 = -mrSges(6,1) * t151 + t180;
t192 = -m(6) * pkin(4) + t115;
t183 = pkin(9) + qJ(3);
t110 = t183 * t146;
t157 = t183 * t144;
t80 = t110 * t149 + t186 * t157;
t191 = t80 ^ 2;
t190 = 0.2e1 * t80;
t188 = m(7) * pkin(5);
t126 = pkin(8) * t170;
t185 = pkin(1) * t152;
t91 = t126 + (-pkin(2) - t185) * t147;
t67 = -pkin(3) * t97 + t91;
t19 = pkin(4) * t63 - pkin(10) * t64 + t67;
t100 = t147 * t150 * pkin(1) + pkin(8) * t169;
t88 = qJ(3) * t147 + t100;
t89 = (-pkin(2) * t152 - qJ(3) * t150 - pkin(1)) * t145;
t53 = -t144 * t88 + t146 * t89;
t32 = -pkin(3) * t169 - pkin(9) * t98 + t53;
t54 = t144 * t89 + t146 * t88;
t37 = pkin(9) * t97 + t54;
t16 = t149 * t32 + t186 * t37;
t8 = -pkin(10) * t169 + t16;
t4 = t148 * t19 + t151 * t8;
t99 = t147 * t185 - t126;
t184 = t99 * mrSges(3,1);
t182 = -qJ(6) - pkin(10);
t181 = -Ifges(5,5) * t64 + Ifges(5,6) * t63;
t107 = t144 * t149 - t186 * t146;
t108 = t186 * t144 + t149 * t146;
t131 = -pkin(3) * t146 - pkin(2);
t71 = pkin(4) * t107 - pkin(10) * t108 + t131;
t82 = t186 * t110 - t149 * t157;
t34 = t148 * t71 + t151 * t82;
t171 = t108 * t151;
t172 = t108 * t148;
t69 = mrSges(7,1) * t172 + mrSges(7,2) * t171;
t179 = Ifges(6,4) * t148;
t178 = Ifges(6,4) * t151;
t177 = Ifges(7,4) * t148;
t176 = Ifges(7,4) * t151;
t175 = t100 * mrSges(3,2);
t174 = Ifges(7,5) * t171 + Ifges(7,3) * t107;
t173 = Ifges(6,5) * t171 + Ifges(6,3) * t107;
t168 = Ifges(5,5) * t108 - Ifges(5,6) * t107;
t117 = Ifges(7,5) * t148 + Ifges(7,6) * t151;
t118 = Ifges(6,5) * t148 + Ifges(6,6) * t151;
t167 = t144 ^ 2 + t146 ^ 2;
t11 = Ifges(7,4) * t44 + Ifges(7,2) * t43 + Ifges(7,6) * t63;
t12 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t63;
t165 = -t11 / 0.2e1 - t12 / 0.2e1;
t13 = Ifges(7,1) * t44 + Ifges(7,4) * t43 + Ifges(7,5) * t63;
t14 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t63;
t164 = t13 / 0.2e1 + t14 / 0.2e1;
t47 = Ifges(7,6) * t107 + (-Ifges(7,2) * t148 + t176) * t108;
t48 = Ifges(6,6) * t107 + (-Ifges(6,2) * t148 + t178) * t108;
t163 = t47 / 0.2e1 + t48 / 0.2e1;
t49 = Ifges(7,5) * t107 + (Ifges(7,1) * t151 - t177) * t108;
t50 = Ifges(6,5) * t107 + (Ifges(6,1) * t151 - t179) * t108;
t162 = t49 / 0.2e1 + t50 / 0.2e1;
t161 = Ifges(3,5) * t170 + Ifges(3,6) * t169 + Ifges(3,3) * t147;
t160 = t117 / 0.2e1 + t118 / 0.2e1;
t119 = Ifges(7,2) * t151 + t177;
t120 = Ifges(6,2) * t151 + t179;
t159 = t119 / 0.2e1 + t120 / 0.2e1;
t121 = Ifges(7,1) * t148 + t176;
t122 = Ifges(6,1) * t148 + t178;
t158 = t121 / 0.2e1 + t122 / 0.2e1;
t68 = -t97 * mrSges(4,1) + t98 * mrSges(4,2);
t30 = t63 * mrSges(5,1) + t64 * mrSges(5,2);
t20 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t3 = -t148 * t8 + t151 * t19;
t33 = -t148 * t82 + t151 * t71;
t109 = -t146 * mrSges(4,1) + t144 * mrSges(4,2);
t77 = t107 * mrSges(5,1) + t108 * mrSges(5,2);
t135 = t148 * mrSges(7,2);
t114 = -t151 * mrSges(7,1) + t135;
t15 = -t149 * t37 + t186 * t32;
t156 = mrSges(6,1) * t148 + mrSges(6,2) * t151;
t155 = -t53 * t144 + t54 * t146;
t7 = pkin(4) * t169 - t15;
t132 = -pkin(5) * t151 - pkin(4);
t116 = t182 * t151;
t113 = t182 * t148;
t112 = Ifges(4,1) * t144 + Ifges(4,4) * t146;
t111 = Ifges(4,4) * t144 + Ifges(4,2) * t146;
t84 = -mrSges(4,1) * t169 - mrSges(4,3) * t98;
t83 = mrSges(4,2) * t169 + mrSges(4,3) * t97;
t79 = Ifges(5,1) * t108 - Ifges(5,4) * t107;
t78 = Ifges(5,4) * t108 - Ifges(5,2) * t107;
t75 = mrSges(6,1) * t107 - mrSges(6,3) * t171;
t74 = mrSges(7,1) * t107 - mrSges(7,3) * t171;
t73 = -mrSges(6,2) * t107 - mrSges(6,3) * t172;
t72 = -mrSges(7,2) * t107 - mrSges(7,3) * t172;
t70 = t156 * t108;
t57 = Ifges(4,1) * t98 + Ifges(4,4) * t97 - Ifges(4,5) * t169;
t56 = Ifges(4,4) * t98 + Ifges(4,2) * t97 - Ifges(4,6) * t169;
t55 = pkin(5) * t172 + t80;
t52 = -mrSges(5,1) * t169 - mrSges(5,3) * t64;
t51 = mrSges(5,2) * t169 - mrSges(5,3) * t63;
t46 = -Ifges(6,6) * t172 + t173;
t45 = -Ifges(7,6) * t172 + t174;
t29 = -qJ(6) * t172 + t34;
t28 = Ifges(5,1) * t64 - Ifges(5,4) * t63 - Ifges(5,5) * t169;
t27 = Ifges(5,4) * t64 - Ifges(5,2) * t63 - Ifges(5,6) * t169;
t26 = mrSges(6,1) * t63 - mrSges(6,3) * t44;
t25 = mrSges(7,1) * t63 - mrSges(7,3) * t44;
t24 = -mrSges(6,2) * t63 + mrSges(6,3) * t43;
t23 = -mrSges(7,2) * t63 + mrSges(7,3) * t43;
t22 = pkin(5) * t107 - qJ(6) * t171 + t33;
t21 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t5 = -t43 * pkin(5) + t7;
t2 = qJ(6) * t43 + t4;
t1 = pkin(5) * t63 - qJ(6) * t44 + t3;
t6 = [t97 * t56 + t98 * t57 + 0.2e1 * t54 * t83 + 0.2e1 * t53 * t84 + 0.2e1 * t91 * t68 + 0.2e1 * t67 * t30 + t64 * t28 + 0.2e1 * t16 * t51 + 0.2e1 * t15 * t52 + 0.2e1 * t1 * t25 + 0.2e1 * t3 * t26 + 0.2e1 * t5 * t20 + 0.2e1 * t7 * t21 + 0.2e1 * t2 * t23 + 0.2e1 * t4 * t24 + (-t27 + t193) * t63 + (t13 + t14) * t44 + (t11 + t12) * t43 + m(3) * (pkin(1) ^ 2 * t145 ^ 2 + t100 ^ 2 + t99 ^ 2) + m(4) * (t53 ^ 2 + t54 ^ 2 + t91 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t67 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) + ((-0.2e1 * t99 * mrSges(3,3) + Ifges(3,5) * t147 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t150) * t145) * t150 + (0.2e1 * t100 * mrSges(3,3) - Ifges(4,5) * t98 + Ifges(3,6) * t147 - Ifges(4,6) * t97 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t150 + (Ifges(5,3) + Ifges(3,2) + Ifges(4,3)) * t152) * t145 + t181) * t152) * t145 + (t161 - 0.2e1 * t175 + 0.2e1 * t184) * t147 + Ifges(2,3); t184 - t175 + t144 * t57 / 0.2e1 + t146 * t56 / 0.2e1 + t131 * t30 + t91 * t109 + t97 * t111 / 0.2e1 + t98 * t112 / 0.2e1 + t82 * t51 - pkin(2) * t68 + t5 * t69 + t7 * t70 + t2 * t72 + t4 * t73 + t1 * t74 + t3 * t75 + t67 * t77 + t64 * t79 / 0.2e1 + t55 * t20 + t29 * t23 + t33 * t26 + t34 * t24 + t22 * t25 + (t21 - t52) * t80 + m(5) * (t131 * t67 - t15 * t80 + t16 * t82) + m(6) * (t3 * t33 + t34 * t4 + t7 * t80) + m(7) * (t1 * t22 + t2 * t29 + t5 * t55) + (-t78 / 0.2e1 + t45 / 0.2e1 + t46 / 0.2e1) * t63 + (-t144 * t84 + t146 * t83) * qJ(3) + t161 - (Ifges(4,5) * t144 + Ifges(4,6) * t146 + t168) * t169 / 0.2e1 + t162 * t44 + t163 * t43 + (t28 / 0.2e1 - t15 * mrSges(5,3) + t164 * t151 + t165 * t148) * t108 + t155 * mrSges(4,3) + m(4) * (-pkin(2) * t91 + qJ(3) * t155) + (t9 / 0.2e1 + t10 / 0.2e1 - t27 / 0.2e1 - t16 * mrSges(5,3)) * t107; -0.2e1 * pkin(2) * t109 + t146 * t111 + t144 * t112 + 0.2e1 * t131 * t77 + 0.2e1 * t22 * t74 + 0.2e1 * t29 * t72 + 0.2e1 * t33 * t75 + 0.2e1 * t34 * t73 + 0.2e1 * t55 * t69 + t70 * t190 + Ifges(3,3) + 0.2e1 * t167 * qJ(3) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t82 + t45 + t46 - t78) * t107 + m(4) * (t167 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t131 ^ 2 + t82 ^ 2 + t191) + m(6) * (t33 ^ 2 + t34 ^ 2 + t191) + m(7) * (t22 ^ 2 + t29 ^ 2 + t55 ^ 2) + (mrSges(5,3) * t190 + t79 + (t49 + t50) * t151 + (-t47 - t48) * t148) * t108; (t25 + t26) * t151 + (t23 + t24) * t148 + m(7) * (t1 * t151 + t148 * t2) + m(6) * (t148 * t4 + t151 * t3) + m(5) * t67 + m(4) * t91 + t30 + t68; -m(4) * pkin(2) + (t74 + t75) * t151 + (t72 + t73) * t148 + m(6) * (t148 * t34 + t151 * t33) + m(7) * (t148 * t29 + t151 * t22) + m(5) * t131 + t77 + t109; m(4) + m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t194; -Ifges(5,3) * t169 + t15 * mrSges(5,1) - t16 * mrSges(5,2) - pkin(4) * t21 + t113 * t25 + t5 * t114 - t116 * t23 + t132 * t20 + t160 * t63 + t158 * t44 + t159 * t43 + m(7) * (t1 * t113 - t116 * t2 + t132 * t5) + (t4 * mrSges(6,3) + t2 * mrSges(7,3) + (m(6) * t4 + t24) * pkin(10) - t165) * t151 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t26) * pkin(10) + t164) * t148 - t181 + t192 * t7; -t116 * t72 + t132 * t69 + t113 * t74 + t55 * t114 - t82 * mrSges(5,2) - pkin(4) * t70 + m(7) * (t113 * t22 - t116 * t29 + t132 * t55) + t160 * t107 + (-mrSges(5,1) + t192) * t80 + (t34 * mrSges(6,3) + t29 * mrSges(7,3) + t158 * t108 + (m(6) * t34 + t73) * pkin(10) + t163) * t151 + (-t33 * mrSges(6,3) - t22 * mrSges(7,3) - t159 * t108 + (-m(6) * t33 - t75) * pkin(10) + t162) * t148 + t168; m(7) * (t113 * t151 - t116 * t148); -0.2e1 * pkin(4) * t115 + 0.2e1 * t132 * t114 + Ifges(5,3) + pkin(10) * mrSges(6,3) * t194 + m(7) * (t113 ^ 2 + t116 ^ 2 + t132 ^ 2) + m(6) * (t166 * pkin(10) ^ 2 + pkin(4) ^ 2) + (t116 * t195 + t119 + t120) * t151 + (t113 * t195 + t121 + t122) * t148; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t193; mrSges(6,1) * t33 + mrSges(7,1) * t22 - mrSges(6,2) * t34 - mrSges(7,2) * t29 + (-Ifges(6,6) - Ifges(7,6)) * t172 + (m(7) * t22 + t74) * pkin(5) + t173 + t174; -t180 - t135 + (mrSges(6,1) + mrSges(7,1) + t188) * t151; mrSges(7,1) * t113 + mrSges(7,2) * t116 - t156 * pkin(10) + (m(7) * t113 - t148 * mrSges(7,3)) * pkin(5) + t118 + t117; Ifges(6,3) + Ifges(7,3) + (0.2e1 * mrSges(7,1) + t188) * pkin(5); m(7) * t5 + t20; m(7) * t55 + t69; 0; m(7) * t132 + t114; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
