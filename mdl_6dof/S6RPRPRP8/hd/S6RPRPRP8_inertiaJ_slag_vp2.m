% Calculate joint inertia matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:27
% EndTime: 2019-03-09 03:24:29
% DurationCPUTime: 0.76s
% Computational Cost: add. (986->211), mult. (1734->280), div. (0->0), fcn. (1680->6), ass. (0->84)
t113 = Ifges(7,2) + Ifges(6,3);
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t87 = t65 ^ 2 + t67 ^ 2;
t68 = cos(qJ(3));
t112 = t68 ^ 2;
t111 = 2 * mrSges(5,1);
t110 = m(6) + m(7);
t105 = m(5) * pkin(3);
t109 = (mrSges(7,2) + mrSges(6,3)) * t87;
t66 = sin(qJ(3));
t69 = -pkin(1) - pkin(7);
t85 = -qJ(4) + t69;
t39 = t85 * t66;
t63 = sin(pkin(9));
t64 = cos(pkin(9));
t81 = t85 * t68;
t20 = t63 * t39 - t64 * t81;
t108 = t20 ^ 2;
t36 = t63 * t66 - t64 * t68;
t35 = t36 ^ 2;
t107 = -2 * mrSges(5,3);
t104 = pkin(3) * t63;
t103 = pkin(3) * t64;
t102 = Ifges(6,4) * t65;
t101 = Ifges(6,4) * t67;
t100 = Ifges(7,5) * t65;
t99 = Ifges(7,5) * t67;
t37 = t63 * t68 + t64 * t66;
t98 = Ifges(6,6) * t37;
t97 = Ifges(7,6) * t37;
t96 = t36 * t20;
t95 = t36 * t65;
t94 = t36 * t67;
t15 = -t37 * mrSges(6,2) + mrSges(6,3) * t95;
t18 = mrSges(7,2) * t95 + t37 * mrSges(7,3);
t92 = t15 + t18;
t16 = t37 * mrSges(6,1) + mrSges(6,3) * t94;
t17 = -t37 * mrSges(7,1) - mrSges(7,2) * t94;
t91 = t16 - t17;
t53 = t66 * pkin(3) + qJ(2);
t14 = pkin(4) * t37 + pkin(8) * t36 + t53;
t22 = t64 * t39 + t63 * t81;
t4 = t65 * t14 + t67 * t22;
t51 = pkin(8) + t104;
t90 = t87 * t37 * t51;
t41 = -t67 * mrSges(6,1) + t65 * mrSges(6,2);
t89 = t41 - mrSges(5,1);
t88 = t87 * t51 ^ 2;
t86 = t66 ^ 2 + t112;
t52 = -pkin(4) - t103;
t83 = m(4) * t86;
t82 = t86 * mrSges(4,3);
t79 = -Ifges(7,6) * t95 + (-Ifges(7,4) - Ifges(6,5)) * t94 + t113 * t37;
t1 = qJ(6) * t37 + t4;
t3 = t14 * t67 - t22 * t65;
t2 = -pkin(5) * t37 - t3;
t78 = t1 * t67 + t2 * t65;
t77 = -t3 * t65 + t4 * t67;
t76 = -t65 * mrSges(6,1) - t67 * mrSges(6,2);
t40 = -t67 * mrSges(7,1) - t65 * mrSges(7,3);
t75 = -t65 * mrSges(7,1) + t67 * mrSges(7,3);
t74 = pkin(5) * t67 + qJ(6) * t65;
t73 = -pkin(5) * t65 + qJ(6) * t67;
t72 = m(7) * t73 + t75 + t76;
t70 = qJ(2) ^ 2;
t56 = Ifges(7,4) * t65;
t55 = Ifges(6,5) * t65;
t54 = Ifges(6,6) * t67;
t45 = Ifges(6,1) * t65 + t101;
t44 = Ifges(7,1) * t65 - t99;
t43 = Ifges(6,2) * t67 + t102;
t42 = -Ifges(7,3) * t67 + t100;
t34 = t37 ^ 2;
t32 = t52 - t74;
t29 = t36 * mrSges(5,2);
t13 = t76 * t36;
t12 = t75 * t36;
t9 = Ifges(6,5) * t37 + (-Ifges(6,1) * t67 + t102) * t36;
t8 = Ifges(7,4) * t37 + (-Ifges(7,1) * t67 - t100) * t36;
t7 = t98 + (Ifges(6,2) * t65 - t101) * t36;
t6 = t97 + (-Ifges(7,3) * t65 - t99) * t36;
t5 = t73 * t36 + t20;
t10 = [Ifges(4,1) * t112 - 0.2e1 * t53 * t29 + 0.2e1 * t5 * t12 + 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t20 * t13 - (2 * pkin(1) * mrSges(3,2)) + Ifges(3,1) + Ifges(2,3) - 0.2e1 * t69 * t82 + (Ifges(5,2) * t37 + t22 * t107 + t53 * t111 + t79) * t37 + (t20 * t107 + Ifges(5,1) * t36 + 0.2e1 * Ifges(5,4) * t37 + (-t8 - t9) * t67 + (-t6 + t7 + t98) * t65) * t36 + m(4) * (t86 * t69 ^ 2 + t70) + m(3) * ((pkin(1) ^ 2) + t70) + m(5) * (t22 ^ 2 + t53 ^ 2 + t108) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t108) + (-0.2e1 * Ifges(4,4) * t68 + Ifges(4,2) * t66) * t66 + 0.2e1 * (mrSges(4,1) * t66 + mrSges(4,2) * t68 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t35 * mrSges(5,3) + mrSges(3,2) + (t12 + t13) * t36 - t82 + (-mrSges(5,3) * t37 - t91 * t65 + t92 * t67) * t37 + m(6) * (t77 * t37 + t96) + m(7) * (t36 * t5 + t78 * t37) + m(5) * (t22 * t37 + t96) + t69 * t83; m(3) + m(5) * (t34 + t35) + t83 + (t87 * t34 + t35) * t110; -t22 * mrSges(5,2) + t32 * t12 + t52 * t13 + t5 * t40 + (mrSges(4,1) * t69 + Ifges(4,5)) * t68 + (-mrSges(4,2) * t69 - Ifges(4,6)) * t66 + t89 * t20 + (-mrSges(5,3) * t104 + t56 / 0.2e1 + t55 / 0.2e1 + t54 / 0.2e1 - Ifges(5,6)) * t37 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) - t97 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t92 * t51) * t67 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t8 / 0.2e1 + t9 / 0.2e1 - t91 * t51) * t65 + m(7) * (t32 * t5 + t78 * t51) + m(6) * (t52 * t20 + t77 * t51) + (-t20 * t64 + t22 * t63) * t105 + (mrSges(5,3) * t103 - Ifges(5,5) + (-t44 / 0.2e1 - t45 / 0.2e1) * t67 + (-t42 / 0.2e1 + t43 / 0.2e1) * t65) * t36; t68 * mrSges(4,1) - t66 * mrSges(4,2) + (t40 + t89) * t36 + (-mrSges(5,2) + t109) * t37 + m(6) * (t36 * t52 + t90) + m(7) * (t32 * t36 + t90) + (-t36 * t64 + t37 * t63) * t105; 0.2e1 * t32 * t40 + 0.2e1 * t52 * t41 + Ifges(4,3) + Ifges(5,3) + (-t42 + t43) * t67 + (t45 + t44) * t65 + m(7) * (t32 ^ 2 + t88) + m(6) * (t52 ^ 2 + t88) + 0.2e1 * t51 * t109 + (t64 * t111 - 0.2e1 * mrSges(5,2) * t63 + (t63 ^ 2 + t64 ^ 2) * t105) * pkin(3); t37 * mrSges(5,1) - t29 + t91 * t67 + t92 * t65 + m(6) * (t3 * t67 + t4 * t65) + m(7) * (t1 * t65 - t2 * t67) + m(5) * t53; 0; 0; t87 * t110 + m(5); Ifges(6,6) * t95 - pkin(5) * t17 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t18 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t79; t72 * t37; t73 * mrSges(7,2) - Ifges(7,6) * t67 + t72 * t51 + t54 + t55 + t56; m(7) * t74 - t40 - t41; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t113; m(7) * t2 + t17; m(7) * t65 * t37; (m(7) * t51 + mrSges(7,2)) * t65; -m(7) * t67; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
