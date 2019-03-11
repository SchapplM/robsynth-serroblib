% Calculate joint inertia matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:31
% EndTime: 2019-03-09 03:04:33
% DurationCPUTime: 0.72s
% Computational Cost: add. (1009->205), mult. (1862->281), div. (0->0), fcn. (1812->8), ass. (0->85)
t113 = Ifges(7,2) + Ifges(6,3);
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t87 = t67 ^ 2 + t69 ^ 2;
t70 = cos(qJ(3));
t112 = t70 ^ 2;
t111 = m(6) + m(7);
t105 = m(5) * pkin(3);
t110 = (mrSges(7,2) + mrSges(6,3)) * t87;
t64 = sin(pkin(9));
t53 = t64 * pkin(1) + pkin(7);
t85 = qJ(4) + t53;
t34 = t85 * t70;
t63 = sin(pkin(10));
t65 = cos(pkin(10));
t68 = sin(qJ(3));
t83 = t85 * t68;
t14 = t63 * t34 + t65 * t83;
t109 = t14 ^ 2;
t36 = t63 * t68 - t65 * t70;
t108 = t36 ^ 2;
t107 = 0.2e1 * t14;
t66 = cos(pkin(9));
t55 = -t66 * pkin(1) - pkin(2);
t39 = -t70 * pkin(3) + t55;
t106 = 0.2e1 * t39;
t104 = t63 * pkin(3);
t103 = t65 * pkin(3);
t102 = Ifges(6,4) * t67;
t101 = Ifges(6,4) * t69;
t100 = Ifges(7,5) * t67;
t99 = Ifges(7,5) * t69;
t98 = Ifges(6,6) * t36;
t97 = Ifges(7,6) * t36;
t96 = t14 * t36;
t38 = t63 * t70 + t65 * t68;
t94 = t38 * t67;
t93 = t38 * t69;
t13 = t36 * pkin(4) - t38 * pkin(8) + t39;
t16 = t65 * t34 - t63 * t83;
t4 = t67 * t13 + t69 * t16;
t19 = -t36 * mrSges(6,2) - mrSges(6,3) * t94;
t22 = -mrSges(7,2) * t94 + t36 * mrSges(7,3);
t92 = t19 + t22;
t20 = t36 * mrSges(6,1) - mrSges(6,3) * t93;
t21 = -t36 * mrSges(7,1) + mrSges(7,2) * t93;
t91 = t20 - t21;
t52 = pkin(8) + t104;
t90 = t87 * t38 * t52;
t41 = -t69 * mrSges(6,1) + t67 * mrSges(6,2);
t89 = t41 - mrSges(5,1);
t88 = t87 * t52 ^ 2;
t86 = t68 ^ 2 + t112;
t54 = -pkin(4) - t103;
t81 = Ifges(7,6) * t94 + (Ifges(7,4) + Ifges(6,5)) * t93 + t113 * t36;
t1 = t36 * qJ(6) + t4;
t3 = t69 * t13 - t67 * t16;
t2 = -t36 * pkin(5) - t3;
t80 = t1 * t69 + t2 * t67;
t79 = -t3 * t67 + t4 * t69;
t78 = -t70 * mrSges(4,1) + t68 * mrSges(4,2);
t77 = t67 * mrSges(6,1) + t69 * mrSges(6,2);
t40 = -t69 * mrSges(7,1) - t67 * mrSges(7,3);
t76 = t67 * mrSges(7,1) - t69 * mrSges(7,3);
t75 = t69 * pkin(5) + t67 * qJ(6);
t74 = pkin(5) * t67 - qJ(6) * t69;
t73 = -m(7) * t74 - t76 - t77;
t58 = Ifges(7,4) * t67;
t57 = Ifges(6,5) * t67;
t56 = Ifges(6,6) * t69;
t45 = Ifges(6,1) * t67 + t101;
t44 = Ifges(7,1) * t67 - t99;
t43 = Ifges(6,2) * t69 + t102;
t42 = -Ifges(7,3) * t69 + t100;
t35 = t38 ^ 2;
t33 = t54 - t75;
t30 = t38 * mrSges(5,2);
t18 = t77 * t38;
t17 = t76 * t38;
t9 = Ifges(6,5) * t36 + (Ifges(6,1) * t69 - t102) * t38;
t8 = Ifges(7,4) * t36 + (Ifges(7,1) * t69 + t100) * t38;
t7 = t98 + (-Ifges(6,2) * t67 + t101) * t38;
t6 = t97 + (Ifges(7,3) * t67 + t99) * t38;
t5 = t74 * t38 + t14;
t10 = [Ifges(4,2) * t112 + 0.2e1 * t55 * t78 + t30 * t106 + 0.2e1 * t5 * t17 + t18 * t107 + 0.2e1 * t4 * t19 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + 0.2e1 * t1 * t22 + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t106 - 0.2e1 * t16 * mrSges(5,3) + Ifges(5,2) * t36 + t81) * t36 + (mrSges(5,3) * t107 + Ifges(5,1) * t38 - 0.2e1 * Ifges(5,4) * t36 + (t8 + t9) * t69 + (t6 - t7 - t98) * t67) * t38 + m(4) * (t86 * t53 ^ 2 + t55 ^ 2) + m(5) * (t16 ^ 2 + t39 ^ 2 + t109) + m(6) * (t3 ^ 2 + t4 ^ 2 + t109) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (t64 ^ 2 + t66 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t68 + 0.2e1 * Ifges(4,4) * t70) * t68 + 0.2e1 * (t66 * mrSges(3,1) - t64 * mrSges(3,2)) * pkin(1) + 0.2e1 * t86 * t53 * mrSges(4,3); (t17 + t18) * t36 + (-t91 * t67 + t92 * t69) * t38 + m(7) * (t5 * t36 + t80 * t38) + m(5) * (t16 * t38 + t96) + m(6) * (t79 * t38 + t96); m(3) + m(5) * (t35 + t108) + m(4) * t86 + (t87 * t35 + t108) * t111; -t16 * mrSges(5,2) + Ifges(4,5) * t68 + Ifges(4,6) * t70 + t33 * t17 + t54 * t18 + t5 * t40 + (-t68 * mrSges(4,1) - t70 * mrSges(4,2)) * t53 + t89 * t14 + (t58 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1 - Ifges(5,6) - mrSges(5,3) * t104) * t36 + (t4 * mrSges(6,3) + t1 * mrSges(7,2) - t97 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t92 * t52) * t69 + (-t3 * mrSges(6,3) + t2 * mrSges(7,2) + t8 / 0.2e1 + t9 / 0.2e1 - t91 * t52) * t67 + m(6) * (t54 * t14 + t79 * t52) + m(7) * (t33 * t5 + t80 * t52) + (-t14 * t65 + t16 * t63) * t105 + (-mrSges(5,3) * t103 + Ifges(5,5) + (t44 / 0.2e1 + t45 / 0.2e1) * t69 + (t42 / 0.2e1 - t43 / 0.2e1) * t67) * t38; -t30 + (t40 + t89) * t36 + m(7) * (t33 * t36 + t90) + m(6) * (t54 * t36 + t90) + (-t36 * t65 + t38 * t63) * t105 - t78 + t38 * t110; 0.2e1 * t33 * t40 + 0.2e1 * t54 * t41 + Ifges(4,3) + Ifges(5,3) + (-t42 + t43) * t69 + (t44 + t45) * t67 + m(6) * (t54 ^ 2 + t88) + m(7) * (t33 ^ 2 + t88) + 0.2e1 * t52 * t110 + (0.2e1 * t65 * mrSges(5,1) - 0.2e1 * t63 * mrSges(5,2) + (t63 ^ 2 + t65 ^ 2) * t105) * pkin(3); t36 * mrSges(5,1) + t30 + t91 * t69 + t92 * t67 + m(7) * (t67 * t1 - t69 * t2) + m(6) * (t69 * t3 + t67 * t4) + m(5) * t39; 0; 0; t87 * t111 + m(5); -Ifges(6,6) * t94 - pkin(5) * t21 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t22 + t1 * mrSges(7,3) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t3 * mrSges(6,1) + t81; t73 * t38; -t74 * mrSges(7,2) - Ifges(7,6) * t69 + t73 * t52 + t56 + t57 + t58; m(7) * t75 - t40 - t41; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t113; m(7) * t2 + t21; m(7) * t94; (m(7) * t52 + mrSges(7,2)) * t67; -m(7) * t69; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
