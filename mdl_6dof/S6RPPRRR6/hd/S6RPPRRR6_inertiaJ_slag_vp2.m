% Calculate joint inertia matrix for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2018-11-23 15:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:50:23
% EndTime: 2018-11-23 15:50:24
% DurationCPUTime: 0.60s
% Computational Cost: add. (772->212), mult. (1409->297), div. (0->0), fcn. (1243->6), ass. (0->81)
t62 = (qJ(2) - pkin(7));
t102 = -2 * t62;
t65 = sin(qJ(6));
t66 = sin(qJ(5));
t68 = cos(qJ(6));
t69 = cos(qJ(5));
t35 = t65 * t66 - t68 * t69;
t36 = t65 * t69 + t68 * t66;
t101 = t35 * mrSges(7,1) + t36 * mrSges(7,2);
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t88 = t69 * t70;
t100 = Ifges(6,5) * t88 + Ifges(6,3) * t67;
t63 = (pkin(1) + qJ(3));
t99 = t63 ^ 2;
t98 = 2 * t63;
t97 = m(7) * pkin(5);
t96 = t69 / 0.2e1;
t95 = -pkin(9) - pkin(8);
t94 = Ifges(6,4) * t66;
t93 = Ifges(6,4) * t69;
t92 = Ifges(6,6) * t67;
t90 = t62 * t67;
t89 = t66 * t70;
t41 = -t69 * mrSges(6,1) + t66 * mrSges(6,2);
t87 = mrSges(5,1) - t41;
t40 = t67 * pkin(4) - t70 * pkin(8) + t63;
t18 = t66 * t40 + t69 * t90;
t86 = t66 ^ 2 + t69 ^ 2;
t59 = t67 ^ 2;
t61 = t70 ^ 2;
t85 = t61 + t59;
t25 = t36 * t70;
t27 = t35 * t70;
t84 = -Ifges(7,5) * t27 - Ifges(7,6) * t25 + Ifges(7,3) * t67;
t83 = t86 * mrSges(6,3);
t82 = t85 * mrSges(5,3);
t24 = t36 * t67;
t26 = t35 * t67;
t81 = -t24 * mrSges(7,1) + t26 * mrSges(7,2);
t79 = t66 * mrSges(6,1) + t69 * mrSges(6,2);
t30 = t69 * t40;
t17 = -t66 * t90 + t30;
t78 = -t17 * t66 + t18 * t69;
t38 = -t67 * mrSges(6,2) - mrSges(6,3) * t89;
t39 = t67 * mrSges(6,1) - mrSges(6,3) * t88;
t77 = t69 * t38 - t66 * t39;
t44 = t95 * t66;
t45 = t95 * t69;
t13 = t68 * t44 + t65 * t45;
t14 = t65 * t44 - t68 * t45;
t32 = Ifges(7,6) * t35;
t33 = Ifges(7,5) * t36;
t76 = t13 * mrSges(7,1) - t14 * mrSges(7,2) - t32 + t33;
t7 = -pkin(9) * t88 + t30 + (-t62 * t66 + pkin(5)) * t67;
t8 = -pkin(9) * t89 + t18;
t2 = -t65 * t8 + t68 * t7;
t3 = t65 * t7 + t68 * t8;
t75 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t84;
t74 = (t68 * mrSges(7,1) - t65 * mrSges(7,2)) * pkin(5);
t71 = (qJ(2) ^ 2);
t57 = t62 ^ 2;
t56 = Ifges(6,5) * t66;
t55 = Ifges(6,6) * t69;
t51 = -t69 * pkin(5) - pkin(4);
t49 = t61 * t62;
t48 = t61 * t57;
t43 = Ifges(6,1) * t66 + t93;
t42 = Ifges(6,2) * t69 + t94;
t34 = (pkin(5) * t66 - t62) * t70;
t28 = t79 * t70;
t23 = Ifges(6,5) * t67 + (Ifges(6,1) * t69 - t94) * t70;
t22 = t92 + (-Ifges(6,2) * t66 + t93) * t70;
t16 = t67 * mrSges(7,1) + t27 * mrSges(7,3);
t15 = -t67 * mrSges(7,2) - t25 * mrSges(7,3);
t11 = Ifges(7,1) * t36 - Ifges(7,4) * t35;
t10 = Ifges(7,4) * t36 - Ifges(7,2) * t35;
t6 = t25 * mrSges(7,1) - t27 * mrSges(7,2);
t5 = -Ifges(7,1) * t27 - Ifges(7,4) * t25 + Ifges(7,5) * t67;
t4 = -Ifges(7,4) * t27 - Ifges(7,2) * t25 + Ifges(7,6) * t67;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(4,3) * t98) + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + 0.2e1 * t17 * t39 + 0.2e1 * t18 * t38 - t25 * t4 - t27 * t5 + 0.2e1 * t34 * t6 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + (mrSges(5,1) * t98 + Ifges(5,2) * t67 + t100 + t84) * t67 + ((mrSges(5,2) * t98) + Ifges(5,1) * t70 - 0.2e1 * Ifges(5,4) * t67 + t69 * t23 + t28 * t102 + (-t22 - t92) * t66) * t70 + m(4) * (t71 + t99) + (m(3) * (pkin(1) ^ 2 + t71)) + m(5) * (t59 * t57 + t48 + t99) + m(6) * (t17 ^ 2 + t18 ^ 2 + t48) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + (2 * (mrSges(4,2) + mrSges(3,3)) * qJ(2)) + t82 * t102; -(m(3) * pkin(1)) - t67 * mrSges(5,1) - t70 * mrSges(5,2) - t36 * t15 + t35 * t16 - t66 * t38 - t69 * t39 + mrSges(3,2) - mrSges(4,3) + m(7) * (t35 * t2 - t36 * t3) + m(6) * (-t69 * t17 - t66 * t18) + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t98; m(3) + m(4) + m(5) + m(6) * t86 + m(7) * (t35 ^ 2 + t36 ^ 2); m(4) * qJ(2) - t26 * t15 - t24 * t16 + mrSges(4,2) + (-t28 - t6) * t70 + t77 * t67 - t82 + m(7) * (-t24 * t2 - t26 * t3 - t70 * t34) + m(6) * (t78 * t67 + t49) + m(5) * (t59 * t62 + t49); m(7) * (-t24 * t35 + t26 * t36); m(4) + m(5) * t85 + m(6) * (t86 * t59 + t61) + m(7) * (t24 ^ 2 + t26 ^ 2 + t61); t66 * t23 / 0.2e1 + t22 * t96 + t51 * t6 - t35 * t4 / 0.2e1 + t36 * t5 / 0.2e1 - t25 * t10 / 0.2e1 - t27 * t11 / 0.2e1 - pkin(4) * t28 + t34 * t101 + t14 * t15 + t13 * t16 + m(7) * (t13 * t2 + t14 * t3 + t51 * t34) + (t56 / 0.2e1 + t55 / 0.2e1 + t33 / 0.2e1 - t32 / 0.2e1 - Ifges(5,6) - (t62 * mrSges(5,2))) * t67 + (-t2 * t36 - t3 * t35) * mrSges(7,3) + t78 * mrSges(6,3) + (m(6) * t78 + t77) * pkin(8) + (Ifges(5,5) + t43 * t96 - t66 * t42 / 0.2e1 + (m(6) * pkin(4) + t87) * t62) * t70; m(7) * (t13 * t35 - t14 * t36); (t24 * t36 + t26 * t35) * mrSges(7,3) + (-t101 + t87) * t70 + (-mrSges(5,2) + t83) * t67 + m(6) * (t86 * t67 * pkin(8) + pkin(4) * t70) + m(7) * (-t13 * t24 - t14 * t26 - t51 * t70); -0.2e1 * pkin(4) * t41 - t35 * t10 + t36 * t11 + t69 * t42 + t66 * t43 + 0.2e1 * t51 * t101 + Ifges(5,3) + m(7) * (t13 ^ 2 + t14 ^ 2 + t51 ^ 2) + m(6) * (t86 * pkin(8) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t13 * t36 - t14 * t35) * mrSges(7,3) + 0.2e1 * pkin(8) * t83; -Ifges(6,6) * t89 + t17 * mrSges(6,1) - t18 * mrSges(6,2) + (m(7) * (t2 * t68 + t3 * t65) + t65 * t15 + t68 * t16) * pkin(5) + t75 + t100; (t35 * t68 - t36 * t65) * t97 + t101 + t41; -t79 * t67 + (-t24 * t68 - t26 * t65) * t97 + t81; t55 + t56 - t79 * pkin(8) + (m(7) * (t13 * t68 + t14 * t65) + (-t65 * t35 - t68 * t36) * mrSges(7,3)) * pkin(5) + t76; Ifges(6,3) + Ifges(7,3) + m(7) * (t65 ^ 2 + t68 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t74; t75; t101; t81; t76; Ifges(7,3) + t74; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
