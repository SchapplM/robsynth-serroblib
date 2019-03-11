% Calculate joint inertia matrix for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:38
% EndTime: 2019-03-09 16:03:42
% DurationCPUTime: 1.41s
% Computational Cost: add. (1610->349), mult. (3450->459), div. (0->0), fcn. (3374->8), ass. (0->120)
t163 = Ifges(6,5) - Ifges(5,6);
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t162 = t121 ^ 2 + t124 ^ 2;
t161 = pkin(9) - qJ(5);
t160 = -m(5) * pkin(3) - mrSges(5,1);
t75 = t161 * t124;
t159 = t75 ^ 2;
t158 = -2 * mrSges(6,3);
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t117 = sin(pkin(6));
t125 = cos(qJ(2));
t143 = t117 * t125;
t118 = cos(pkin(6));
t122 = sin(qJ(2));
t144 = t117 * t122;
t61 = -t118 * t124 + t121 * t144;
t33 = -t120 * t61 + t123 * t143;
t157 = t33 / 0.2e1;
t34 = t120 * t143 + t123 * t61;
t156 = t34 / 0.2e1;
t155 = -pkin(4) - pkin(10);
t154 = t120 / 0.2e1;
t153 = -t123 / 0.2e1;
t139 = t120 ^ 2 + t123 ^ 2;
t67 = m(7) * t139;
t152 = m(6) + t67;
t151 = pkin(1) * t118;
t63 = -pkin(8) * t144 + t125 * t151;
t150 = t63 * mrSges(3,1);
t64 = pkin(8) * t143 + t122 * t151;
t149 = t64 * mrSges(3,2);
t148 = -mrSges(6,3) + mrSges(5,2);
t119 = qJ(4) + pkin(5);
t46 = pkin(9) * t118 + t64;
t47 = (-pkin(2) * t125 - pkin(9) * t122 - pkin(1)) * t117;
t19 = t121 * t47 + t124 * t46;
t62 = t118 * t121 + t124 * t144;
t39 = mrSges(5,1) * t143 + t62 * mrSges(5,2);
t147 = Ifges(7,4) * t120;
t146 = Ifges(7,4) * t123;
t145 = qJ(5) * t61;
t142 = t120 * t124;
t141 = t123 * t124;
t140 = t162 * pkin(9) ^ 2;
t78 = -Ifges(7,5) * t120 - Ifges(7,6) * t123;
t138 = t78 / 0.2e1 + Ifges(6,6);
t137 = Ifges(5,2) + Ifges(4,3) + Ifges(6,3);
t70 = -t124 * pkin(3) - t121 * qJ(4) - pkin(2);
t5 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t62;
t136 = Ifges(3,5) * t144 + Ifges(3,6) * t143 + Ifges(3,3) * t118;
t45 = -t118 * pkin(2) - t63;
t135 = Ifges(6,6) * t143;
t28 = t62 * mrSges(6,1) + t61 * mrSges(6,2);
t18 = -t121 * t46 + t124 * t47;
t134 = t139 * mrSges(7,3);
t77 = t121 * mrSges(6,1) - t124 * mrSges(6,2);
t66 = t124 * pkin(4) - t70;
t15 = pkin(3) * t143 - t18;
t40 = -mrSges(6,2) * t143 - t62 * mrSges(6,3);
t13 = t61 * pkin(3) - t62 * qJ(4) + t45;
t3 = pkin(5) * t62 + t155 * t61 - t13;
t9 = pkin(4) * t143 - qJ(5) * t62 + t15;
t4 = pkin(10) * t143 + t9;
t1 = -t120 * t4 + t123 * t3;
t2 = t120 * t3 + t123 * t4;
t132 = -t120 * t1 + t123 * t2;
t131 = -mrSges(7,1) * t120 - t123 * mrSges(7,2);
t41 = pkin(5) * t121 + pkin(10) * t124 + t66;
t72 = t161 * t121;
t26 = -t120 * t72 + t123 * t41;
t27 = t120 * t41 + t123 * t72;
t130 = -t120 * t26 + t123 * t27;
t48 = -Ifges(7,5) * t141 + Ifges(7,6) * t142 + Ifges(7,3) * t121;
t14 = -qJ(4) * t143 + t19;
t129 = (-Ifges(5,4) - Ifges(4,5)) * t62 + (Ifges(4,6) + t163) * t61;
t127 = qJ(4) ^ 2;
t126 = -pkin(3) - pkin(4);
t112 = -pkin(3) + t155;
t106 = Ifges(5,4) * t121;
t105 = Ifges(4,5) * t121;
t104 = Ifges(4,6) * t124;
t86 = Ifges(4,1) * t121 + Ifges(4,4) * t124;
t85 = Ifges(5,1) * t121 - Ifges(5,5) * t124;
t84 = -Ifges(6,1) * t124 - Ifges(6,4) * t121;
t83 = -Ifges(7,1) * t120 - t146;
t82 = Ifges(4,4) * t121 + Ifges(4,2) * t124;
t81 = -Ifges(6,4) * t124 - Ifges(6,2) * t121;
t80 = -Ifges(7,2) * t123 - t147;
t79 = Ifges(5,5) * t121 - Ifges(5,3) * t124;
t74 = -mrSges(4,1) * t124 + mrSges(4,2) * t121;
t73 = -mrSges(5,1) * t124 - mrSges(5,3) * t121;
t71 = mrSges(7,1) * t123 - mrSges(7,2) * t120;
t69 = mrSges(7,1) * t121 + mrSges(7,3) * t141;
t68 = -mrSges(7,2) * t121 + mrSges(7,3) * t142;
t65 = t131 * t124;
t50 = Ifges(7,5) * t121 + (-Ifges(7,1) * t123 + t147) * t124;
t49 = Ifges(7,6) * t121 + (Ifges(7,2) * t120 - t146) * t124;
t38 = -mrSges(4,1) * t143 - mrSges(4,3) * t62;
t37 = mrSges(4,2) * t143 - mrSges(4,3) * t61;
t36 = mrSges(6,1) * t143 - mrSges(6,3) * t61;
t35 = -mrSges(5,2) * t61 - mrSges(5,3) * t143;
t30 = mrSges(4,1) * t61 + mrSges(4,2) * t62;
t29 = mrSges(5,1) * t61 - mrSges(5,3) * t62;
t25 = Ifges(4,1) * t62 - Ifges(4,4) * t61 - Ifges(4,5) * t143;
t24 = Ifges(5,1) * t62 - Ifges(5,4) * t143 + Ifges(5,5) * t61;
t23 = Ifges(6,1) * t61 - Ifges(6,4) * t62 + Ifges(6,5) * t143;
t22 = Ifges(4,4) * t62 - Ifges(4,2) * t61 - Ifges(4,6) * t143;
t21 = Ifges(6,4) * t61 - Ifges(6,2) * t62 + t135;
t20 = Ifges(5,5) * t62 - Ifges(5,6) * t143 + Ifges(5,3) * t61;
t17 = mrSges(7,1) * t62 - mrSges(7,3) * t34;
t16 = -mrSges(7,2) * t62 + mrSges(7,3) * t33;
t12 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t11 = -t14 - t145;
t10 = -pkin(4) * t61 - t13;
t8 = -t119 * t143 + t145 + t19;
t7 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t62;
t6 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t62;
t31 = [0.2e1 * t1 * t17 + 0.2e1 * t10 * t28 + 0.2e1 * t11 * t36 + 0.2e1 * t8 * t12 + 0.2e1 * t13 * t29 + 0.2e1 * t14 * t35 + 0.2e1 * t15 * t39 + 0.2e1 * t2 * t16 + 0.2e1 * t18 * t38 + 0.2e1 * t19 * t37 + 0.2e1 * t45 * t30 + t33 * t6 + t34 * t7 + 0.2e1 * t9 * t40 + Ifges(2,3) + (t20 + t23 - t22) * t61 + (t136 - 0.2e1 * t149 + 0.2e1 * t150) * t118 + (t24 + t25 + t5 - t21) * t62 + ((-0.2e1 * t63 * mrSges(3,3) + Ifges(3,5) * t118 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t122) * t117) * t122 + (0.2e1 * t64 * mrSges(3,3) + Ifges(3,6) * t118 - Ifges(6,6) * t62 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t122 + (Ifges(3,2) + t137) * t125) * t117 + t129) * t125) * t117 + m(3) * (pkin(1) ^ 2 * t117 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(4) * (t18 ^ 2 + t19 ^ 2 + t45 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t9 ^ 2); (t86 / 0.2e1 - t81 / 0.2e1 + t85 / 0.2e1 + t48 / 0.2e1) * t62 + (t79 / 0.2e1 - t82 / 0.2e1 + t84 / 0.2e1) * t61 + m(5) * (t13 * t70 + (t15 * t121 + t14 * t124) * pkin(9)) + m(4) * (-pkin(2) * t45 + (-t18 * t121 + t19 * t124) * pkin(9)) + (t6 * t154 + t7 * t153 + t19 * mrSges(4,3) + t11 * mrSges(6,3) + t14 * mrSges(5,2) - t20 / 0.2e1 + t22 / 0.2e1 - t23 / 0.2e1 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t143 + (t35 + t37) * pkin(9)) * t124 + (-t135 / 0.2e1 + t15 * mrSges(5,2) - t18 * mrSges(4,3) - t9 * mrSges(6,3) - t21 / 0.2e1 + t24 / 0.2e1 + t25 / 0.2e1 + t5 / 0.2e1 + (-t38 + t39) * pkin(9)) * t121 + t136 + t2 * t68 + t1 * t69 + t70 * t29 + t72 * t40 + t13 * t73 + t45 * t74 + t10 * t77 + t8 * t65 + t66 * t28 - pkin(2) * t30 + t26 * t17 + t27 * t16 + (t12 - t36) * t75 + t50 * t156 + t49 * t157 + (-t105 / 0.2e1 - t104 / 0.2e1 - t106 / 0.2e1) * t143 + m(7) * (t1 * t26 + t2 * t27 + t75 * t8) + m(6) * (t10 * t66 - t11 * t75 + t72 * t9) - t149 + t150; -0.2e1 * pkin(2) * t74 + 0.2e1 * t26 * t69 + 0.2e1 * t27 * t68 + 0.2e1 * t75 * t65 + 0.2e1 * t66 * t77 + 0.2e1 * t70 * t73 + Ifges(3,3) + (t72 * t158 + t48 - t81 + t85 + t86) * t121 + m(5) * (t70 ^ 2 + t140) + m(4) * (pkin(2) ^ 2 + t140) + m(7) * (t26 ^ 2 + t27 ^ 2 + t159) + m(6) * (t66 ^ 2 + t72 ^ 2 + t159) + (t120 * t49 - t123 * t50 + t75 * t158 - t79 + t82 - t84) * t124 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(9) * t162; (t112 * t16 - t2 * mrSges(7,3) - t6 / 0.2e1) * t123 + (-t112 * t17 + t1 * mrSges(7,3) - t7 / 0.2e1) * t120 + (t35 - t36) * qJ(4) + m(6) * (-qJ(4) * t11 + t126 * t9) + m(5) * (-pkin(3) * t15 + qJ(4) * t14) + t138 * t62 - t137 * t143 - t129 + m(7) * (t132 * t112 + t119 * t8) + t126 * t40 + t119 * t12 + t8 * t71 + t80 * t157 + t83 * t156 - pkin(3) * t39 + t18 * mrSges(4,1) - t19 * mrSges(4,2) + t9 * mrSges(6,2) - t11 * mrSges(6,1) + t14 * mrSges(5,3) - t15 * mrSges(5,1); t72 * mrSges(6,2) + t119 * t65 + t104 + t105 + t106 + (t71 + mrSges(6,1)) * t75 + (t112 * t68 - t27 * mrSges(7,3) - t49 / 0.2e1) * t123 + (-t112 * t69 + t26 * mrSges(7,3) - t50 / 0.2e1) * t120 + m(6) * (qJ(4) * t75 + t126 * t72) + m(7) * (t130 * t112 + t119 * t75) + (-pkin(3) * mrSges(5,2) - t126 * mrSges(6,3) + t138) * t121 + (t148 * qJ(4) + t83 * t153 + t80 * t154 + t163) * t124 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t124 + (-mrSges(4,1) + t160) * t121) * pkin(9); 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * t126 * mrSges(6,2) + 0.2e1 * t119 * t71 - t120 * t83 - t123 * t80 + m(7) * (t139 * t112 ^ 2 + t119 ^ 2) + m(5) * (pkin(3) ^ 2 + t127) + m(6) * (t126 ^ 2 + t127) + t137 + 0.2e1 * (mrSges(5,3) + mrSges(6,1)) * qJ(4) - 0.2e1 * t112 * t134; m(5) * t15 + m(6) * t9 + m(7) * t132 - t120 * t17 + t123 * t16 + t39 + t40; -t120 * t69 + t123 * t68 + m(7) * t130 + m(6) * t72 + (m(5) * pkin(9) + t148) * t121; m(6) * t126 + t112 * t67 + mrSges(6,2) - t134 + t160; m(5) + t152; t120 * t16 + t123 * t17 + m(7) * (t1 * t123 + t120 * t2) + m(6) * t10 + t28; t120 * t68 + t123 * t69 + m(7) * (t120 * t27 + t123 * t26) + m(6) * t66 + t77; 0; 0; t152; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t26 - mrSges(7,2) * t27 + t48; t131 * t112 + t78; t131; t71; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t31(1) t31(2) t31(4) t31(7) t31(11) t31(16); t31(2) t31(3) t31(5) t31(8) t31(12) t31(17); t31(4) t31(5) t31(6) t31(9) t31(13) t31(18); t31(7) t31(8) t31(9) t31(10) t31(14) t31(19); t31(11) t31(12) t31(13) t31(14) t31(15) t31(20); t31(16) t31(17) t31(18) t31(19) t31(20) t31(21);];
Mq  = res;
