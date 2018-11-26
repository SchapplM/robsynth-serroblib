% Calculate joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2018-11-23 15:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:48:05
% EndTime: 2018-11-23 15:48:05
% DurationCPUTime: 0.73s
% Computational Cost: add. (1467->217), mult. (2752->319), div. (0->0), fcn. (2949->10), ass. (0->85)
t73 = cos(pkin(11));
t77 = sin(qJ(4));
t103 = cos(qJ(4));
t71 = sin(pkin(11));
t90 = t103 * t71;
t50 = t77 * t73 + t90;
t75 = sin(qJ(6));
t76 = sin(qJ(5));
t78 = cos(qJ(6));
t79 = cos(qJ(5));
t52 = t75 * t79 + t76 * t78;
t23 = t52 * t50;
t51 = -t75 * t76 + t78 * t79;
t24 = t51 * t50;
t10 = t23 * mrSges(7,1) + t24 * mrSges(7,2);
t86 = mrSges(6,1) * t76 + mrSges(6,2) * t79;
t28 = t86 * t50;
t113 = t10 + t28;
t97 = t71 * t77;
t48 = -t103 * t73 + t97;
t98 = t50 * t79;
t112 = Ifges(6,5) * t98 + Ifges(6,3) * t48;
t72 = sin(pkin(10));
t60 = pkin(1) * t72 + qJ(3);
t104 = pkin(7) + t60;
t39 = t104 * t73;
t25 = t104 * t90 + t39 * t77;
t111 = t25 ^ 2;
t46 = t48 ^ 2;
t68 = t73 ^ 2;
t110 = 0.2e1 * t25;
t74 = cos(pkin(10));
t62 = -pkin(1) * t74 - pkin(2);
t53 = -pkin(3) * t73 + t62;
t109 = 0.2e1 * t53;
t108 = m(7) * pkin(5);
t106 = -pkin(9) - pkin(8);
t105 = pkin(4) * t48;
t102 = Ifges(6,4) * t76;
t101 = Ifges(6,4) * t79;
t100 = t25 * t48;
t99 = t50 * t76;
t22 = -pkin(8) * t50 + t105 + t53;
t27 = t103 * t39 - t104 * t97;
t9 = t76 * t22 + t79 * t27;
t96 = Ifges(7,5) * t52 + Ifges(7,6) * t51;
t54 = -t79 * mrSges(6,1) + t76 * mrSges(6,2);
t95 = t54 - mrSges(5,1);
t94 = Ifges(6,5) * t76 + Ifges(6,6) * t79;
t93 = t71 ^ 2 + t68;
t92 = t76 ^ 2 + t79 ^ 2;
t91 = Ifges(7,5) * t24 - Ifges(7,6) * t23 + Ifges(7,3) * t48;
t89 = t92 * t50;
t88 = -t73 * mrSges(4,1) + t71 * mrSges(4,2);
t31 = -t51 * mrSges(7,1) + t52 * mrSges(7,2);
t8 = t79 * t22 - t27 * t76;
t87 = -t76 * t8 + t79 * t9;
t57 = t106 * t76;
t58 = t106 * t79;
t35 = t57 * t78 + t58 * t75;
t36 = t57 * t75 - t58 * t78;
t85 = t35 * mrSges(7,1) - t36 * mrSges(7,2) + t96;
t4 = pkin(5) * t48 - pkin(9) * t98 + t8;
t5 = -pkin(9) * t99 + t9;
t2 = t4 * t78 - t5 * t75;
t3 = t4 * t75 + t5 * t78;
t84 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t91;
t83 = (mrSges(7,1) * t78 - mrSges(7,2) * t75) * pkin(5);
t63 = -pkin(5) * t79 - pkin(4);
t56 = Ifges(6,1) * t76 + t101;
t55 = Ifges(6,2) * t79 + t102;
t47 = t50 ^ 2;
t40 = t50 * mrSges(5,2);
t33 = Ifges(7,1) * t52 + Ifges(7,4) * t51;
t32 = Ifges(7,4) * t52 + Ifges(7,2) * t51;
t30 = mrSges(6,1) * t48 - mrSges(6,3) * t98;
t29 = -mrSges(6,2) * t48 - mrSges(6,3) * t99;
t15 = Ifges(6,5) * t48 + (Ifges(6,1) * t79 - t102) * t50;
t14 = Ifges(6,6) * t48 + (-Ifges(6,2) * t76 + t101) * t50;
t13 = pkin(5) * t99 + t25;
t12 = mrSges(7,1) * t48 - mrSges(7,3) * t24;
t11 = -mrSges(7,2) * t48 - mrSges(7,3) * t23;
t7 = Ifges(7,1) * t24 - Ifges(7,4) * t23 + Ifges(7,5) * t48;
t6 = Ifges(7,4) * t24 - Ifges(7,2) * t23 + Ifges(7,6) * t48;
t1 = [Ifges(4,2) * t68 + 0.2e1 * t62 * t88 + t40 * t109 + t28 * t110 + 0.2e1 * t9 * t29 + 0.2e1 * t8 * t30 - t23 * t6 + t24 * t7 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t13 * t10 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t71 + 0.2e1 * Ifges(4,4) * t73) * t71 + (mrSges(5,1) * t109 - 0.2e1 * mrSges(5,3) * t27 + Ifges(5,2) * t48 + t112 + t91) * t48 + (mrSges(5,3) * t110 + Ifges(5,1) * t50 - t76 * t14 + t79 * t15 + (-Ifges(6,6) * t76 - (2 * Ifges(5,4))) * t48) * t50 + m(5) * (t27 ^ 2 + t53 ^ 2 + t111) + m(4) * (t93 * t60 ^ 2 + t62 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2 + t111) + m(7) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (t72 ^ 2 + t74 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t74 - mrSges(3,2) * t72) * pkin(1) + 0.2e1 * t93 * t60 * mrSges(4,3); t24 * t11 - t23 * t12 + (t79 * t29 - t76 * t30) * t50 + t113 * t48 + m(7) * (t13 * t48 - t2 * t23 + t24 * t3) + m(6) * (t87 * t50 + t100) + m(5) * (t27 * t50 + t100); m(3) + m(7) * (t23 ^ 2 + t24 ^ 2 + t46) + m(6) * (t92 * t47 + t46) + m(5) * (t47 + t46) + m(4) * t93; t48 * mrSges(5,1) + t52 * t11 + t51 * t12 + t76 * t29 + t79 * t30 + t40 + m(7) * (t2 * t51 + t3 * t52) + m(6) * (t76 * t9 + t79 * t8) + m(5) * t53 + m(4) * t62 + t88; m(7) * (-t23 * t51 + t24 * t52); m(4) + m(5) + m(6) * t92 + m(7) * (t51 ^ 2 + t52 ^ 2); t63 * t10 - Ifges(5,6) * t48 + t51 * t6 / 0.2e1 + t52 * t7 / 0.2e1 - t27 * mrSges(5,2) - pkin(4) * t28 + t13 * t31 - t23 * t32 / 0.2e1 + t24 * t33 / 0.2e1 + t35 * t12 + t36 * t11 + t95 * t25 + (-t2 * t52 + t3 * t51) * mrSges(7,3) + (pkin(8) * t29 + t9 * mrSges(6,3) + t14 / 0.2e1) * t79 + (-pkin(8) * t30 - t8 * mrSges(6,3) + t15 / 0.2e1) * t76 + m(6) * (-pkin(4) * t25 + t87 * pkin(8)) + m(7) * (t13 * t63 + t2 * t35 + t3 * t36) + (-t76 * t55 / 0.2e1 + t79 * t56 / 0.2e1 + Ifges(5,5)) * t50 + (t94 + t96) * t48 / 0.2e1; -t40 + (t23 * t52 + t24 * t51) * mrSges(7,3) + mrSges(6,3) * t89 + (t31 + t95) * t48 + m(7) * (-t23 * t35 + t24 * t36 + t48 * t63) + m(6) * (pkin(8) * t89 - t105); m(7) * (t35 * t51 + t36 * t52); -0.2e1 * pkin(4) * t54 + 0.2e1 * t63 * t31 + t51 * t32 + t52 * t33 + t79 * t55 + t76 * t56 + Ifges(5,3) + m(7) * (t35 ^ 2 + t36 ^ 2 + t63 ^ 2) + m(6) * (t92 * pkin(8) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t35 * t52 + t36 * t51) * mrSges(7,3) + 0.2e1 * t92 * pkin(8) * mrSges(6,3); -Ifges(6,6) * t99 + t8 * mrSges(6,1) - t9 * mrSges(6,2) + (t78 * t12 + t75 * t11 + m(7) * (t2 * t78 + t3 * t75)) * pkin(5) + t84 + t112; (-t23 * t78 + t24 * t75) * t108 - t113; (t51 * t78 + t52 * t75) * t108 - t54 - t31; -t86 * pkin(8) + (m(7) * (t35 * t78 + t36 * t75) + (t51 * t75 - t52 * t78) * mrSges(7,3)) * pkin(5) + t85 + t94; Ifges(6,3) + Ifges(7,3) + m(7) * (t75 ^ 2 + t78 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t83; t84; -t10; -t31; t85; Ifges(7,3) + t83; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
