% Calculate joint inertia matrix for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:53:22
% EndTime: 2018-11-23 17:53:23
% DurationCPUTime: 1.09s
% Computational Cost: add. (2062->280), mult. (3642->381), div. (0->0), fcn. (3835->8), ass. (0->113)
t108 = sin(qJ(6));
t109 = sin(qJ(5));
t112 = cos(qJ(6));
t113 = cos(qJ(5));
t128 = t108 * t113 + t112 * t109;
t80 = -t108 * t109 + t112 * t113;
t176 = t128 ^ 2 + t80 ^ 2;
t115 = cos(qJ(2));
t100 = -pkin(2) * t115 - pkin(1);
t110 = sin(qJ(3));
t111 = sin(qJ(2));
t114 = cos(qJ(3));
t81 = t110 * t115 + t111 * t114;
t127 = -qJ(4) * t81 + t100;
t165 = pkin(3) + pkin(9);
t78 = t110 * t111 - t114 * t115;
t26 = t165 * t78 + t127;
t164 = -pkin(8) - pkin(7);
t89 = t164 * t111;
t90 = t164 * t115;
t56 = -t110 * t90 - t114 * t89;
t32 = pkin(4) * t81 + t56;
t29 = t113 * t32;
t5 = pkin(5) * t81 + t29 + (-pkin(10) * t78 - t26) * t109;
t10 = t109 * t32 + t113 * t26;
t150 = t113 * t78;
t6 = pkin(10) * t150 + t10;
t3 = -t108 * t6 + t112 * t5;
t4 = t108 * t5 + t112 * t6;
t175 = t4 * t128 + t3 * t80;
t86 = mrSges(6,1) * t109 + mrSges(6,2) * t113;
t174 = mrSges(5,3) + t86;
t147 = t109 ^ 2 + t113 ^ 2;
t139 = t147 * mrSges(6,3);
t101 = Ifges(6,5) * t113;
t154 = t108 * t128;
t173 = -pkin(5) * mrSges(7,3) * t154 + t101;
t153 = t109 * t78;
t172 = Ifges(6,5) * t153 + Ifges(6,6) * t150 + Ifges(6,3) * t81;
t93 = pkin(2) * t110 + qJ(4);
t170 = t93 ^ 2;
t42 = pkin(3) * t78 + t127;
t169 = -0.2e1 * t42;
t46 = mrSges(7,1) * t128 + mrSges(7,2) * t80;
t168 = 0.2e1 * t46;
t167 = 0.2e1 * t100;
t99 = -pkin(2) * t114 - pkin(3);
t92 = -pkin(9) + t99;
t161 = -pkin(10) + t92;
t160 = mrSges(5,1) + mrSges(4,3);
t159 = -pkin(10) - t165;
t158 = Ifges(7,5) * t80 - Ifges(7,6) * t128;
t157 = Ifges(6,4) * t109;
t156 = Ifges(6,4) * t113;
t155 = qJ(4) * t93;
t152 = t112 * t80;
t151 = t113 * mrSges(6,1);
t146 = t111 ^ 2 + t115 ^ 2;
t145 = 0.2e1 * mrSges(7,3);
t37 = t80 * t78;
t38 = t128 * t78;
t144 = Ifges(7,5) * t38 + Ifges(7,6) * t37 + Ifges(7,3) * t81;
t143 = mrSges(7,3) * t152;
t58 = t110 * t89 - t114 * t90;
t142 = t56 ^ 2 + t58 ^ 2;
t141 = t80 * mrSges(7,1) - mrSges(7,2) * t128;
t140 = m(6) * t147;
t138 = t147 * t165;
t137 = 0.2e1 * t174;
t64 = t161 * t109;
t65 = t161 * t113;
t39 = -t108 * t64 + t112 * t65;
t40 = t108 * t65 + t112 * t64;
t136 = t128 * t40 + t80 * t39;
t83 = t159 * t109;
t84 = t159 * t113;
t49 = -t108 * t83 + t112 * t84;
t50 = t108 * t84 + t112 * t83;
t135 = t128 * t50 + t80 * t49;
t9 = -t109 * t26 + t29;
t134 = t10 * t109 + t9 * t113;
t133 = t109 * mrSges(6,2) - t151;
t43 = mrSges(6,1) * t81 - mrSges(6,3) * t153;
t44 = -mrSges(6,2) * t81 + mrSges(6,3) * t150;
t132 = t109 * t44 + t113 * t43;
t131 = t39 * mrSges(7,1) - t40 * mrSges(7,2) + t158;
t130 = t49 * mrSges(7,1) - t50 * mrSges(7,2) + t158;
t129 = -0.2e1 * t139;
t126 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t144;
t125 = (mrSges(4,1) * t114 - mrSges(4,2) * t110) * pkin(2);
t124 = (mrSges(7,1) * t112 - mrSges(7,2) * t108) * pkin(5);
t47 = Ifges(7,4) * t80 - Ifges(7,2) * t128;
t48 = Ifges(7,1) * t80 - Ifges(7,4) * t128;
t87 = -Ifges(6,2) * t109 + t156;
t88 = Ifges(6,1) * t113 - t157;
t123 = -t109 * t87 + t113 * t88 - t128 * t47 + t80 * t48 + Ifges(5,1) + Ifges(4,3);
t122 = -t176 * mrSges(7,3) + mrSges(5,2) - t139;
t11 = Ifges(7,4) * t38 + Ifges(7,2) * t37 + Ifges(7,6) * t81;
t12 = Ifges(7,1) * t38 + Ifges(7,4) * t37 + Ifges(7,5) * t81;
t17 = (-pkin(5) * t113 - pkin(4)) * t78 + t58;
t24 = Ifges(6,6) * t81 + (Ifges(6,2) * t113 + t157) * t78;
t25 = Ifges(6,5) * t81 + (Ifges(6,1) * t109 + t156) * t78;
t33 = -pkin(4) * t78 + t58;
t121 = t17 * t46 + t37 * t47 / 0.2e1 + t38 * t48 / 0.2e1 - t109 * t24 / 0.2e1 + t113 * t25 / 0.2e1 + t33 * t86 + t88 * t153 / 0.2e1 + t87 * t150 / 0.2e1 - t128 * t11 / 0.2e1 + t80 * t12 / 0.2e1 + (-mrSges(4,2) + mrSges(5,3)) * t58 - t134 * mrSges(6,3) + (Ifges(4,5) - Ifges(5,4)) * t81 + (Ifges(5,5) - Ifges(4,6)) * t78 + (mrSges(5,2) - mrSges(4,1)) * t56 + (-Ifges(6,6) * t109 + t101 + t158) * t81 / 0.2e1 - t175 * mrSges(7,3);
t117 = qJ(4) ^ 2;
t102 = t109 * pkin(5);
t95 = qJ(4) + t102;
t85 = t102 + t93;
t41 = t133 * t78;
t19 = mrSges(7,1) * t81 - mrSges(7,3) * t38;
t18 = -mrSges(7,2) * t81 + mrSges(7,3) * t37;
t14 = -mrSges(7,1) * t37 + mrSges(7,2) * t38;
t1 = [-0.2e1 * pkin(1) * (-t115 * mrSges(3,1) + t111 * mrSges(3,2)) + t111 * (Ifges(3,1) * t111 + Ifges(3,4) * t115) + t115 * (Ifges(3,4) * t111 + Ifges(3,2) * t115) + 0.2e1 * t33 * t41 + 0.2e1 * t9 * t43 + 0.2e1 * t10 * t44 + t37 * t11 + t38 * t12 + 0.2e1 * t17 * t14 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + Ifges(2,3) + 0.2e1 * t146 * pkin(7) * mrSges(3,3) + (mrSges(4,2) * t167 + mrSges(5,3) * t169 + (Ifges(4,1) + Ifges(5,2)) * t81 + 0.2e1 * t160 * t56 + t144 + t172) * t81 + (mrSges(4,1) * t167 + mrSges(5,2) * t169 + t109 * t25 + t113 * t24 + (Ifges(4,2) + Ifges(5,3)) * t78 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t81 - 0.2e1 * t160 * t58) * t78 + m(3) * (t146 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t100 ^ 2 + t142) + m(5) * (t42 ^ 2 + t142) + m(6) * (t10 ^ 2 + t33 ^ 2 + t9 ^ 2) + m(7) * (t17 ^ 2 + t3 ^ 2 + t4 ^ 2); (-t93 * t78 + t99 * t81) * mrSges(5,1) + m(6) * (t134 * t92 + t33 * t93) + Ifges(3,5) * t111 + Ifges(3,6) * t115 + t93 * t41 + t85 * t14 + t121 + t39 * t19 + t40 * t18 + t132 * t92 + m(5) * (t56 * t99 + t58 * t93) + m(7) * (t17 * t85 + t3 * t39 + t4 * t40) + (-t111 * mrSges(3,1) - t115 * mrSges(3,2)) * pkin(7) + (m(4) * (t110 * t58 - t114 * t56) + (-t110 * t78 - t114 * t81) * mrSges(4,3)) * pkin(2); 0.2e1 * t99 * mrSges(5,2) + t85 * t168 + Ifges(3,3) + t93 * t137 + 0.2e1 * t125 - t136 * t145 + t92 * t129 + m(7) * (t39 ^ 2 + t40 ^ 2 + t85 ^ 2) + m(6) * (t147 * t92 ^ 2 + t170) + m(5) * (t99 ^ 2 + t170) + m(4) * (t110 ^ 2 + t114 ^ 2) * pkin(2) ^ 2 + t123; -t132 * t165 + m(7) * (t17 * t95 + t3 * t49 + t4 * t50) + m(5) * (-pkin(3) * t56 + qJ(4) * t58) + (-pkin(3) * t81 - qJ(4) * t78) * mrSges(5,1) + m(6) * (qJ(4) * t33 - t134 * t165) + t95 * t14 + t121 + qJ(4) * t41 + t49 * t19 + t50 * t18; (t95 + t85) * t46 + t125 + (t99 - pkin(3)) * mrSges(5,2) + m(6) * (-t92 * t138 + t155) + m(7) * (t39 * t49 + t40 * t50 + t85 * t95) + m(5) * (-pkin(3) * t99 + t155) + ((-t39 - t49) * t80 - (t40 + t50) * t128) * mrSges(7,3) + t123 + (t165 - t92) * t139 + t174 * (t93 + qJ(4)); -0.2e1 * pkin(3) * mrSges(5,2) + t95 * t168 + qJ(4) * t137 - t135 * t145 - t165 * t129 + m(7) * (t49 ^ 2 + t50 ^ 2 + t95 ^ 2) + m(6) * (t147 * t165 ^ 2 + t117) + m(5) * (pkin(3) ^ 2 + t117) + t123; m(5) * t56 + m(6) * t134 + m(7) * t175 + t81 * mrSges(5,1) + t128 * t18 + t80 * t19 + t132; m(5) * t99 + m(7) * t136 + t92 * t140 + t122; -m(5) * pkin(3) - m(6) * t138 + m(7) * t135 + t122; m(7) * t176 + m(5) + t140; t9 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (t108 * t4 + t112 * t3) + t108 * t18 + t112 * t19) * pkin(5) + t126 + t172; t92 * t151 + (-mrSges(6,2) * t92 - Ifges(6,6)) * t109 + (m(7) * (t108 * t40 + t112 * t39) - t143) * pkin(5) + t131 + t173; -t165 * t151 + (mrSges(6,2) * t165 - Ifges(6,6)) * t109 + (m(7) * (t108 * t50 + t112 * t49) - t143) * pkin(5) + t130 + t173; m(7) * (t152 + t154) * pkin(5) - t133 + t141; Ifges(6,3) + Ifges(7,3) + m(7) * (t108 ^ 2 + t112 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t124; t126; t131; t130; t141; Ifges(7,3) + t124; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
