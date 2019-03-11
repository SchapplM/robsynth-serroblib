% Calculate joint inertia matrix for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:18
% EndTime: 2019-03-09 14:29:21
% DurationCPUTime: 1.27s
% Computational Cost: add. (2142->297), mult. (3818->418), div. (0->0), fcn. (3865->8), ass. (0->122)
t170 = pkin(3) + pkin(7);
t112 = sin(qJ(2));
t116 = cos(qJ(2));
t169 = t112 ^ 2 + t116 ^ 2;
t109 = sin(qJ(6));
t113 = cos(qJ(6));
t110 = sin(qJ(5));
t111 = sin(qJ(4));
t114 = cos(qJ(5));
t115 = cos(qJ(4));
t127 = t110 * t111 - t114 * t115;
t74 = -t110 * t115 - t114 * t111;
t138 = t109 * t74 - t113 * t127;
t40 = -t109 * t127 - t113 * t74;
t168 = t138 ^ 2 + t40 ^ 2;
t117 = -pkin(2) - pkin(8);
t139 = -qJ(3) * t112 - pkin(1);
t66 = t117 * t116 + t139;
t87 = t170 * t112;
t78 = t115 * t87;
t30 = pkin(4) * t112 + t78 + (pkin(9) * t116 - t66) * t111;
t147 = t115 * t116;
t43 = t111 * t87 + t115 * t66;
t33 = -pkin(9) * t147 + t43;
t13 = -t110 * t33 + t114 * t30;
t59 = t74 * t116;
t4 = pkin(5) * t112 - pkin(10) * t59 + t13;
t14 = t110 * t30 + t114 * t33;
t58 = t127 * t116;
t5 = pkin(10) * t58 + t14;
t2 = -t109 * t5 + t113 * t4;
t3 = t109 * t4 + t113 * t5;
t167 = t138 * t2 + t3 * t40;
t152 = -pkin(9) + t117;
t81 = t152 * t111;
t82 = t152 * t115;
t48 = -t110 * t81 + t114 * t82;
t28 = pkin(10) * t127 + t48;
t49 = t110 * t82 + t114 * t81;
t29 = pkin(10) * t74 + t49;
t10 = t109 * t28 + t113 * t29;
t9 = -t109 * t29 + t113 * t28;
t166 = t10 * t40 + t138 * t9;
t155 = pkin(4) * t110;
t93 = pkin(4) * t114 + pkin(5);
t61 = -t109 * t155 + t113 * t93;
t62 = t109 * t93 + t113 * t155;
t165 = t138 * t61 + t40 * t62;
t164 = t109 * t40 + t113 * t138;
t158 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + Ifges(6,3) * t112;
t157 = -m(4) * pkin(2) + mrSges(4,2);
t156 = -t111 / 0.2e1;
t154 = t62 * mrSges(7,2);
t153 = Ifges(6,3) + Ifges(7,3);
t151 = Ifges(5,4) * t111;
t150 = Ifges(5,4) * t115;
t149 = t109 * mrSges(7,2);
t148 = t115 * mrSges(5,1);
t90 = t111 * pkin(4) + qJ(3);
t146 = t169 * pkin(7) ^ 2;
t88 = t170 * t116;
t145 = -t111 ^ 2 - t115 ^ 2;
t144 = t127 ^ 2 + t74 ^ 2;
t143 = pkin(5) * t149;
t26 = -t109 * t59 + t113 * t58;
t27 = t109 * t58 + t113 * t59;
t142 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t112;
t63 = pkin(4) * t147 + t88;
t141 = mrSges(7,1) * t138 - t40 * mrSges(7,2);
t140 = m(5) * t145;
t137 = t145 * mrSges(5,3);
t57 = t61 * mrSges(7,1);
t136 = Ifges(7,3) + t57 - t154;
t134 = -t127 * t13 - t14 * t74;
t133 = -t127 * t48 - t49 * t74;
t36 = Ifges(7,6) * t40;
t37 = Ifges(7,5) * t138;
t132 = t9 * mrSges(7,1) - t10 * mrSges(7,2) - t36 + t37;
t131 = -t111 * mrSges(5,2) + t148;
t130 = -Ifges(5,5) * t111 - Ifges(5,6) * t115;
t129 = t110 * t74 + t114 * t127;
t42 = -t111 * t66 + t78;
t128 = t43 * t111 + t42 * t115;
t126 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t142;
t125 = (mrSges(6,1) * t114 - mrSges(6,2) * t110) * pkin(4);
t124 = -mrSges(6,1) * t127 + t74 * mrSges(6,2) + t141;
t68 = Ifges(6,6) * t74;
t69 = Ifges(6,5) * t127;
t123 = t48 * mrSges(6,1) - t49 * mrSges(6,2) + t132 + t68 - t69;
t122 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t126 + t158;
t118 = qJ(3) ^ 2;
t98 = Ifges(5,5) * t115;
t97 = Ifges(5,3) * t112;
t94 = t113 * pkin(5) * mrSges(7,1);
t86 = Ifges(5,1) * t115 - t151;
t85 = -Ifges(5,2) * t111 + t150;
t84 = mrSges(5,1) * t111 + mrSges(5,2) * t115;
t83 = -pkin(2) * t116 + t139;
t80 = -mrSges(5,2) * t112 - mrSges(5,3) * t147;
t79 = mrSges(5,3) * t111 * t116 + mrSges(5,1) * t112;
t65 = t131 * t116;
t56 = Ifges(5,5) * t112 + (-Ifges(5,1) * t111 - t150) * t116;
t55 = t112 * Ifges(5,6) + (-Ifges(5,2) * t115 - t151) * t116;
t52 = -pkin(5) * t74 + t90;
t51 = mrSges(6,1) * t112 - mrSges(6,3) * t59;
t50 = -mrSges(6,2) * t112 + mrSges(6,3) * t58;
t47 = -Ifges(6,1) * t127 + Ifges(6,4) * t74;
t46 = -Ifges(6,4) * t127 + Ifges(6,2) * t74;
t45 = -mrSges(6,1) * t74 - mrSges(6,2) * t127;
t34 = -pkin(5) * t58 + t63;
t32 = -mrSges(6,1) * t58 + mrSges(6,2) * t59;
t25 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t112;
t24 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t112;
t19 = mrSges(7,1) * t112 - mrSges(7,3) * t27;
t18 = -mrSges(7,2) * t112 + mrSges(7,3) * t26;
t17 = Ifges(7,1) * t138 - Ifges(7,4) * t40;
t16 = Ifges(7,4) * t138 - Ifges(7,2) * t40;
t15 = mrSges(7,1) * t40 + mrSges(7,2) * t138;
t11 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t112;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t112;
t1 = [0.2e1 * t34 * t11 + 0.2e1 * t13 * t51 + 0.2e1 * t14 * t50 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 + t58 * t24 + t59 * t25 + t26 * t7 + t27 * t8 + 0.2e1 * t63 * t32 + 0.2e1 * t42 * t79 + 0.2e1 * t43 * t80 + 0.2e1 * t88 * t65 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t83 * mrSges(4,2) - t111 * t56 - t115 * t55 + (Ifges(4,3) + Ifges(3,2)) * t116) * t116 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t83 * mrSges(4,3) + t97 + (Ifges(4,2) + Ifges(3,1)) * t112 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t130) * t116 + t142 + t158) * t112 + m(4) * (t83 ^ 2 + t146) + m(3) * (pkin(1) ^ 2 + t146) + m(5) * (t42 ^ 2 + t43 ^ 2 + t88 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t63 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t169; t138 * t8 / 0.2e1 - t127 * t25 / 0.2e1 + (Ifges(5,6) * t156 + t98 / 0.2e1 - t69 / 0.2e1 + t68 / 0.2e1 + t37 / 0.2e1 - t36 / 0.2e1 + Ifges(3,5) - Ifges(4,4) - pkin(2) * mrSges(4,1)) * t112 - t40 * t7 / 0.2e1 + m(6) * (t13 * t48 + t14 * t49 + t63 * t90) + m(7) * (t10 * t3 + t2 * t9 + t34 * t52) + t88 * t84 + t90 * t32 + t74 * t24 / 0.2e1 + t63 * t45 + qJ(3) * t65 + t49 * t50 + t48 * t51 + t52 * t11 + t58 * t46 / 0.2e1 + t59 * t47 / 0.2e1 + t26 * t16 / 0.2e1 + t27 * t17 / 0.2e1 + t34 * t15 + (-t55 / 0.2e1 + t117 * t80 - t43 * mrSges(5,3)) * t111 - t167 * mrSges(7,3) + (t56 / 0.2e1 - t42 * mrSges(5,3) + t117 * t79) * t115 + t10 * t18 + t9 * t19 + m(5) * (qJ(3) * t88 + t128 * t117) - t134 * mrSges(6,3) + (Ifges(3,6) - Ifges(4,5) - t115 * t85 / 0.2e1 + t86 * t156 + qJ(3) * mrSges(4,1)) * t116 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t116 + (-mrSges(3,1) + t157) * t112) * pkin(7); -0.2e1 * pkin(2) * mrSges(4,2) - t111 * t85 + t115 * t86 + 0.2e1 * t52 * t15 - t40 * t16 + t138 * t17 + 0.2e1 * t90 * t45 + t74 * t46 - t127 * t47 + Ifges(4,1) + Ifges(3,3) + m(7) * (t10 ^ 2 + t52 ^ 2 + t9 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2 + t90 ^ 2) + m(5) * (-t145 * t117 ^ 2 + t118) + m(4) * (pkin(2) ^ 2 + t118) + 0.2e1 * (t84 + mrSges(4,3)) * qJ(3) - 0.2e1 * t166 * mrSges(7,3) - 0.2e1 * t133 * mrSges(6,3) + 0.2e1 * t117 * t137; t111 * t80 + t115 * t79 + t40 * t18 + t138 * t19 - t74 * t50 - t127 * t51 + (m(4) * pkin(7) + mrSges(4,1)) * t112 + m(7) * t167 + m(6) * t134 + m(5) * t128; m(6) * t133 + m(7) * t166 - t144 * mrSges(6,3) - t168 * mrSges(7,3) - t117 * t140 + t137 + t157; m(6) * t144 + m(7) * t168 + m(4) - t140; (t114 * t51 + t110 * t50 + m(6) * (t110 * t14 + t114 * t13)) * pkin(4) + m(7) * (t2 * t61 + t3 * t62) + t130 * t116 + t97 + t61 * t19 + t62 * t18 + t42 * mrSges(5,1) - t43 * mrSges(5,2) + t122; m(7) * (t10 * t62 + t61 * t9) + t98 + t117 * t148 + (-mrSges(5,2) * t117 - Ifges(5,6)) * t111 - t165 * mrSges(7,3) + (m(6) * (t110 * t49 + t114 * t48) + t129 * mrSges(6,3)) * pkin(4) + t123; -m(6) * t129 * pkin(4) + m(7) * t165 + t124 + t131; -0.2e1 * t154 + Ifges(5,3) + 0.2e1 * t57 + 0.2e1 * t125 + m(7) * (t61 ^ 2 + t62 ^ 2) + m(6) * (t110 ^ 2 + t114 ^ 2) * pkin(4) ^ 2 + t153; (m(7) * (t109 * t3 + t113 * t2) + t109 * t18 + t113 * t19) * pkin(5) + t122; (m(7) * (t10 * t109 + t113 * t9) - t164 * mrSges(7,3)) * pkin(5) + t123; m(7) * t164 * pkin(5) + t124; Ifges(6,3) + t94 + t125 + (m(7) * (t109 * t62 + t113 * t61) - t149) * pkin(5) + t136; -0.2e1 * t143 + 0.2e1 * t94 + m(7) * (t109 ^ 2 + t113 ^ 2) * pkin(5) ^ 2 + t153; t126; t132; t141; t136; Ifges(7,3) + t94 - t143; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
