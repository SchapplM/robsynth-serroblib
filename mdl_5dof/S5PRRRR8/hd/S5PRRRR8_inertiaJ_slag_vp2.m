% Calculate joint inertia matrix for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:52
% EndTime: 2019-12-05 17:14:55
% DurationCPUTime: 0.54s
% Computational Cost: add. (648->157), mult. (1418->241), div. (0->0), fcn. (1429->10), ass. (0->71)
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t41 = -t65 * mrSges(6,1) + t61 * mrSges(6,2);
t104 = t41 - mrSges(5,1);
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t66 = cos(qJ(4));
t67 = cos(qJ(3));
t39 = t62 * t67 + t66 * t63;
t77 = mrSges(6,1) * t61 + mrSges(6,2) * t65;
t17 = t77 * t39;
t98 = -pkin(8) - pkin(7);
t45 = t98 * t67;
t81 = t98 * t63;
t25 = -t62 * t45 - t66 * t81;
t103 = m(6) * t25 + t17;
t38 = t62 * t63 - t66 * t67;
t89 = t61 * mrSges(6,3);
t19 = -t38 * mrSges(6,2) - t39 * t89;
t51 = -t67 * pkin(3) - pkin(2);
t18 = t38 * pkin(4) - t39 * pkin(9) + t51;
t27 = -t66 * t45 + t62 * t81;
t2 = t65 * t18 - t61 * t27;
t92 = t39 * t65;
t20 = t38 * mrSges(6,1) - mrSges(6,3) * t92;
t3 = t61 * t18 + t65 * t27;
t102 = m(6) * (-t2 * t61 + t3 * t65) + t65 * t19 - t61 * t20;
t60 = cos(pkin(5));
t59 = sin(pkin(5));
t64 = sin(qJ(2));
t91 = t59 * t64;
t31 = t60 * t67 - t63 * t91;
t32 = t60 * t63 + t67 * t91;
t13 = -t66 * t31 + t62 * t32;
t101 = t13 ^ 2;
t100 = t25 ^ 2;
t58 = t67 ^ 2;
t99 = 0.2e1 * t25;
t97 = mrSges(6,3) * t65;
t96 = Ifges(6,4) * t61;
t95 = Ifges(6,4) * t65;
t94 = t25 * t13;
t93 = t39 * t61;
t68 = cos(qJ(2));
t90 = t59 * t68;
t86 = Ifges(6,5) * t92 + Ifges(6,3) * t38;
t85 = Ifges(6,5) * t61 + Ifges(6,6) * t65;
t84 = t61 ^ 2 + t65 ^ 2;
t83 = t63 ^ 2 + t58;
t43 = Ifges(6,2) * t65 + t96;
t44 = Ifges(6,1) * t61 + t95;
t82 = t65 * t43 + t61 * t44 + Ifges(5,3);
t49 = t62 * pkin(3) + pkin(9);
t80 = t84 * t49;
t15 = t62 * t31 + t66 * t32;
t7 = -t61 * t15 - t65 * t90;
t8 = t65 * t15 - t61 * t90;
t78 = -t61 * t7 + t65 * t8;
t76 = -t31 * t63 + t32 * t67;
t75 = 0.2e1 * t84 * mrSges(6,3);
t74 = (t66 * mrSges(5,1) - t62 * mrSges(5,2)) * pkin(3);
t73 = -t15 * mrSges(5,2) + t104 * t13 - t7 * t89 + t8 * t97;
t10 = Ifges(6,6) * t38 + (-Ifges(6,2) * t61 + t95) * t39;
t11 = Ifges(6,5) * t38 + (Ifges(6,1) * t65 - t96) * t39;
t72 = -t27 * mrSges(5,2) + t61 * t11 / 0.2e1 - t2 * t89 + t3 * t97 - t43 * t93 / 0.2e1 + t44 * t92 / 0.2e1 + Ifges(5,5) * t39 + t65 * t10 / 0.2e1 + (t85 / 0.2e1 - Ifges(5,6)) * t38 + t104 * t25;
t54 = t59 ^ 2;
t50 = -t66 * pkin(3) - pkin(4);
t47 = t54 * t68 ^ 2;
t42 = -t67 * mrSges(4,1) + t63 * mrSges(4,2);
t21 = t38 * mrSges(5,1) + t39 * mrSges(5,2);
t1 = [m(2) + m(6) * (t7 ^ 2 + t8 ^ 2 + t101) + m(5) * (t15 ^ 2 + t101 + t47) + m(4) * (t31 ^ 2 + t32 ^ 2 + t47) + m(3) * (t54 * t64 ^ 2 + t60 ^ 2 + t47); t13 * t17 + t8 * t19 + t7 * t20 + (t13 * t39 - t15 * t38) * mrSges(5,3) + t76 * mrSges(4,3) + (-t64 * mrSges(3,2) + (mrSges(3,1) - t21 - t42) * t68) * t59 + m(6) * (t2 * t7 + t3 * t8 + t94) + m(5) * (t27 * t15 - t51 * t90 + t94) + m(4) * (pkin(2) * t90 + t76 * pkin(7)); Ifges(4,2) * t58 - 0.2e1 * pkin(2) * t42 + t17 * t99 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t51 * t21 + Ifges(3,3) + (Ifges(4,1) * t63 + 0.2e1 * Ifges(4,4) * t67) * t63 + 0.2e1 * t83 * pkin(7) * mrSges(4,3) + (-0.2e1 * t27 * mrSges(5,3) + Ifges(5,2) * t38 + t86) * t38 + m(6) * (t2 ^ 2 + t3 ^ 2 + t100) + m(5) * (t27 ^ 2 + t51 ^ 2 + t100) + m(4) * (t83 * pkin(7) ^ 2 + pkin(2) ^ 2) + (mrSges(5,3) * t99 + Ifges(5,1) * t39 - t61 * t10 + t65 * t11 + (-Ifges(6,6) * t61 - (2 * Ifges(5,4))) * t38) * t39; t31 * mrSges(4,1) - t32 * mrSges(4,2) + m(6) * (t50 * t13 + t49 * t78) + m(5) * (-t13 * t66 + t15 * t62) * pkin(3) + t73; t72 + (m(5) * (-t25 * t66 + t27 * t62) + (-t62 * t38 - t66 * t39) * mrSges(5,3)) * pkin(3) + (-t63 * mrSges(4,1) - t67 * mrSges(4,2)) * pkin(7) + Ifges(4,6) * t67 + Ifges(4,5) * t63 + t103 * t50 + t102 * t49; 0.2e1 * t50 * t41 + Ifges(4,3) + 0.2e1 * t74 + t49 * t75 + m(6) * (t49 ^ 2 * t84 + t50 ^ 2) + m(5) * (t62 ^ 2 + t66 ^ 2) * pkin(3) ^ 2 + t82; m(6) * (-pkin(4) * t13 + pkin(9) * t78) + t73; -t103 * pkin(4) + t102 * pkin(9) + t72; m(6) * (-pkin(4) * t50 + pkin(9) * t80) + (-pkin(4) + t50) * t41 + t74 + (pkin(9) * t84 + t80) * mrSges(6,3) + t82; -0.2e1 * pkin(4) * t41 + m(6) * (pkin(9) ^ 2 * t84 + pkin(4) ^ 2) + pkin(9) * t75 + t82; t7 * mrSges(6,1) - t8 * mrSges(6,2); t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t93 + t86; -t49 * t77 + t85; -t77 * pkin(9) + t85; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
