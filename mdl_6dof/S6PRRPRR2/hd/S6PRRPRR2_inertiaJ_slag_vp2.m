% Calculate joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:15:03
% EndTime: 2018-11-23 15:15:04
% DurationCPUTime: 0.86s
% Computational Cost: add. (1576->258), mult. (3280->387), div. (0->0), fcn. (3665->12), ass. (0->95)
t82 = sin(pkin(12));
t84 = cos(pkin(12));
t88 = sin(qJ(3));
t92 = cos(qJ(3));
t59 = t82 * t92 + t84 * t88;
t91 = cos(qJ(5));
t114 = t59 * t91;
t58 = t82 * t88 - t84 * t92;
t124 = Ifges(6,5) * t114 + Ifges(6,3) * t58;
t83 = sin(pkin(6));
t89 = sin(qJ(2));
t113 = t83 * t89;
t85 = cos(pkin(6));
t48 = -t113 * t88 + t85 * t92;
t49 = t113 * t92 + t85 * t88;
t24 = -t84 * t48 + t49 * t82;
t23 = t24 ^ 2;
t111 = -qJ(4) - pkin(8);
t103 = t111 * t88;
t66 = t111 * t92;
t43 = -t84 * t103 - t66 * t82;
t123 = t43 ^ 2;
t81 = t92 ^ 2;
t122 = 0.2e1 * t43;
t121 = m(7) * pkin(5);
t72 = pkin(3) * t82 + pkin(9);
t119 = pkin(10) + t72;
t87 = sin(qJ(5));
t118 = Ifges(6,4) * t87;
t117 = Ifges(6,4) * t91;
t116 = t24 * t43;
t115 = t59 * t87;
t93 = cos(qJ(2));
t112 = t83 * t93;
t74 = -pkin(3) * t92 - pkin(2);
t35 = pkin(4) * t58 - pkin(9) * t59 + t74;
t45 = t103 * t82 - t84 * t66;
t13 = t87 * t35 + t91 * t45;
t86 = sin(qJ(6));
t90 = cos(qJ(6));
t61 = -t86 * t87 + t90 * t91;
t62 = t86 * t91 + t87 * t90;
t110 = Ifges(7,5) * t62 + Ifges(7,6) * t61;
t64 = -t91 * mrSges(6,1) + t87 * mrSges(6,2);
t109 = t64 - mrSges(5,1);
t108 = Ifges(6,5) * t87 + Ifges(6,6) * t91;
t107 = t87 ^ 2 + t91 ^ 2;
t106 = t88 ^ 2 + t81;
t27 = t62 * t59;
t28 = t61 * t59;
t105 = Ifges(7,5) * t28 - Ifges(7,6) * t27 + Ifges(7,3) * t58;
t26 = t48 * t82 + t49 * t84;
t16 = -t112 * t91 - t26 * t87;
t17 = -t112 * t87 + t26 * t91;
t5 = t16 * t90 - t17 * t86;
t6 = t16 * t86 + t17 * t90;
t104 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t73 = -pkin(3) * t84 - pkin(4);
t38 = t58 * mrSges(5,1) + t59 * mrSges(5,2);
t40 = -t61 * mrSges(7,1) + t62 * mrSges(7,2);
t12 = t91 * t35 - t45 * t87;
t102 = mrSges(6,1) * t87 + mrSges(6,2) * t91;
t101 = -t16 * t87 + t17 * t91;
t100 = -t48 * t88 + t49 * t92;
t53 = t119 * t87;
t54 = t119 * t91;
t33 = -t53 * t90 - t54 * t86;
t34 = -t53 * t86 + t54 * t90;
t99 = t33 * mrSges(7,1) - t34 * mrSges(7,2) + t110;
t10 = -pkin(10) * t115 + t13;
t7 = pkin(5) * t58 - pkin(10) * t114 + t12;
t2 = -t10 * t86 + t7 * t90;
t3 = t10 * t90 + t7 * t86;
t98 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t105;
t97 = (mrSges(7,1) * t90 - mrSges(7,2) * t86) * pkin(5);
t77 = t83 ^ 2;
t70 = t77 * t93 ^ 2;
t68 = Ifges(6,1) * t87 + t117;
t67 = Ifges(6,2) * t91 + t118;
t65 = -mrSges(4,1) * t92 + mrSges(4,2) * t88;
t63 = -pkin(5) * t91 + t73;
t42 = Ifges(7,1) * t62 + Ifges(7,4) * t61;
t41 = Ifges(7,4) * t62 + Ifges(7,2) * t61;
t37 = mrSges(6,1) * t58 - mrSges(6,3) * t114;
t36 = -mrSges(6,2) * t58 - mrSges(6,3) * t115;
t32 = t102 * t59;
t20 = pkin(5) * t115 + t43;
t19 = Ifges(6,5) * t58 + (Ifges(6,1) * t91 - t118) * t59;
t18 = Ifges(6,6) * t58 + (-Ifges(6,2) * t87 + t117) * t59;
t15 = mrSges(7,1) * t58 - mrSges(7,3) * t28;
t14 = -mrSges(7,2) * t58 - mrSges(7,3) * t27;
t11 = mrSges(7,1) * t27 + mrSges(7,2) * t28;
t9 = Ifges(7,1) * t28 - Ifges(7,4) * t27 + Ifges(7,5) * t58;
t8 = Ifges(7,4) * t28 - Ifges(7,2) * t27 + Ifges(7,6) * t58;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t23) + m(6) * (t16 ^ 2 + t17 ^ 2 + t23) + m(5) * (t26 ^ 2 + t23 + t70) + m(4) * (t48 ^ 2 + t49 ^ 2 + t70) + m(3) * (t77 * t89 ^ 2 + t85 ^ 2 + t70); -t26 * t58 * mrSges(5,3) + t6 * t14 + t5 * t15 + t16 * t37 + t17 * t36 + t100 * mrSges(4,3) + (t59 * mrSges(5,3) + t11 + t32) * t24 + (-t89 * mrSges(3,2) + (mrSges(3,1) - t38 - t65) * t93) * t83 + m(7) * (t2 * t5 + t20 * t24 + t3 * t6) + m(6) * (t12 * t16 + t13 * t17 + t116) + m(5) * (-t112 * t74 + t26 * t45 + t116) + m(4) * (pkin(2) * t112 + pkin(8) * t100); Ifges(4,2) * t81 - 0.2e1 * pkin(2) * t65 + 0.2e1 * t20 * t11 + 0.2e1 * t12 * t37 + 0.2e1 * t13 * t36 + 0.2e1 * t3 * t14 + 0.2e1 * t2 * t15 - t27 * t8 + t28 * t9 + t32 * t122 + 0.2e1 * t74 * t38 + Ifges(3,3) + (Ifges(4,1) * t88 + 0.2e1 * Ifges(4,4) * t92) * t88 + 0.2e1 * t106 * pkin(8) * mrSges(4,3) + (mrSges(5,3) * t122 + Ifges(5,1) * t59 - t18 * t87 + t19 * t91) * t59 + (-0.2e1 * t45 * mrSges(5,3) + Ifges(5,2) * t58 + (-Ifges(6,6) * t87 - (2 * Ifges(5,4))) * t59 + t105 + t124) * t58 + m(4) * (pkin(8) ^ 2 * t106 + pkin(2) ^ 2) + m(5) * (t45 ^ 2 + t74 ^ 2 + t123) + m(6) * (t12 ^ 2 + t13 ^ 2 + t123) + m(7) * (t2 ^ 2 + t20 ^ 2 + t3 ^ 2); t48 * mrSges(4,1) - t49 * mrSges(4,2) - t26 * mrSges(5,2) + (-t5 * t62 + t6 * t61) * mrSges(7,3) + t101 * mrSges(6,3) + (t40 + t109) * t24 + m(7) * (t24 * t63 + t33 * t5 + t34 * t6) + m(6) * (t101 * t72 + t24 * t73) + m(5) * (-t24 * t84 + t26 * t82) * pkin(3); Ifges(4,6) * t92 + Ifges(4,5) * t88 + t62 * t9 / 0.2e1 + t63 * t11 + t73 * t32 - Ifges(5,6) * t58 + t61 * t8 / 0.2e1 - t45 * mrSges(5,2) + t33 * t15 + t34 * t14 + t20 * t40 - t27 * t41 / 0.2e1 + t28 * t42 / 0.2e1 + t109 * t43 + (-t88 * mrSges(4,1) - t92 * mrSges(4,2)) * pkin(8) + (-t2 * t62 + t3 * t61) * mrSges(7,3) + (t18 / 0.2e1 + t72 * t36 + t13 * mrSges(6,3)) * t91 + (t19 / 0.2e1 - t72 * t37 - t12 * mrSges(6,3)) * t87 + m(6) * (t43 * t73 + (-t12 * t87 + t13 * t91) * t72) + m(7) * (t2 * t33 + t20 * t63 + t3 * t34) + (Ifges(5,5) + t91 * t68 / 0.2e1 - t87 * t67 / 0.2e1) * t59 + (m(5) * (-t43 * t84 + t45 * t82) + (-t82 * t58 - t84 * t59) * mrSges(5,3)) * pkin(3) + (t108 + t110) * t58 / 0.2e1; 0.2e1 * t63 * t40 + t61 * t41 + t62 * t42 + 0.2e1 * t73 * t64 + t91 * t67 + t87 * t68 + Ifges(4,3) + Ifges(5,3) + m(7) * (t33 ^ 2 + t34 ^ 2 + t63 ^ 2) + m(6) * (t107 * t72 ^ 2 + t73 ^ 2) + m(5) * (t82 ^ 2 + t84 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t84 - mrSges(5,2) * t82) * pkin(3) + 0.2e1 * (-t33 * t62 + t34 * t61) * mrSges(7,3) + 0.2e1 * t107 * t72 * mrSges(6,3); -m(5) * t112 + m(7) * (t5 * t61 + t6 * t62) + m(6) * (t16 * t91 + t17 * t87); t62 * t14 + t61 * t15 + t87 * t36 + t91 * t37 + m(7) * (t2 * t61 + t3 * t62) + m(6) * (t12 * t91 + t13 * t87) + m(5) * t74 + t38; m(7) * (t33 * t61 + t34 * t62); m(5) + m(6) * t107 + m(7) * (t61 ^ 2 + t62 ^ 2); t16 * mrSges(6,1) - t17 * mrSges(6,2) + (t5 * t90 + t6 * t86) * t121 + t104; -Ifges(6,6) * t115 + t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t2 * t90 + t3 * t86) + t86 * t14 + t90 * t15) * pkin(5) + t98 + t124; -t102 * t72 + (m(7) * (t33 * t90 + t34 * t86) + (t61 * t86 - t62 * t90) * mrSges(7,3)) * pkin(5) + t99 + t108; (t61 * t90 + t62 * t86) * t121 - t64 - t40; Ifges(6,3) + Ifges(7,3) + m(7) * (t86 ^ 2 + t90 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t97; t104; t98; t99; -t40; Ifges(7,3) + t97; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
