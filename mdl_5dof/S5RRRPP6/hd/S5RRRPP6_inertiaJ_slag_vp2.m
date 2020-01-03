% Calculate joint inertia matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:00
% EndTime: 2019-12-31 21:00:02
% DurationCPUTime: 0.72s
% Computational Cost: add. (801->218), mult. (1645->306), div. (0->0), fcn. (1514->6), ass. (0->79)
t103 = 2 * pkin(6);
t85 = cos(qJ(2));
t102 = pkin(6) * t85;
t83 = sin(qJ(2));
t75 = t83 * pkin(6);
t82 = sin(qJ(3));
t101 = Ifges(4,4) * t82;
t84 = cos(qJ(3));
t100 = Ifges(4,4) * t84;
t99 = t82 * t83;
t98 = t83 * t84;
t96 = -qJ(4) - pkin(7);
t59 = -t85 * pkin(2) - t83 * pkin(7) - pkin(1);
t54 = t84 * t59;
t94 = qJ(4) * t83;
t14 = -t84 * t94 + t54 + (-pkin(6) * t82 - pkin(3)) * t85;
t32 = t84 * t102 + t82 * t59;
t22 = -t82 * t94 + t32;
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t4 = t80 * t14 + t81 * t22;
t51 = t80 * t82 - t81 * t84;
t43 = t51 * t83;
t30 = t85 * mrSges(6,1) - t43 * mrSges(6,2);
t58 = pkin(3) * t99 + t75;
t95 = t82 ^ 2 + t84 ^ 2;
t61 = t96 * t84;
t91 = t96 * t82;
t24 = -t80 * t61 - t81 * t91;
t26 = -t81 * t61 + t80 * t91;
t93 = t24 ^ 2 + t26 ^ 2;
t92 = Ifges(6,2) + Ifges(4,3) + Ifges(5,3);
t70 = -t84 * pkin(3) - pkin(2);
t52 = t80 * t84 + t81 * t82;
t42 = t52 * t83;
t11 = t42 * mrSges(5,1) - t43 * mrSges(5,2);
t17 = t51 * mrSges(5,1) + t52 * mrSges(5,2);
t10 = t42 * mrSges(6,1) + t43 * mrSges(6,3);
t16 = t51 * mrSges(6,1) - t52 * mrSges(6,3);
t90 = -Ifges(4,5) * t98 + (Ifges(6,4) + Ifges(5,5)) * t43 + (Ifges(5,6) - Ifges(6,6)) * t42;
t89 = t82 * mrSges(4,1) + t84 * mrSges(4,2);
t3 = t81 * t14 - t80 * t22;
t87 = pkin(6) ^ 2;
t79 = t85 ^ 2;
t77 = t83 ^ 2;
t74 = t77 * t87;
t73 = Ifges(4,5) * t82;
t72 = Ifges(4,6) * t84;
t69 = -t81 * pkin(3) - pkin(4);
t67 = t80 * pkin(3) + qJ(5);
t63 = Ifges(4,1) * t82 + t100;
t62 = Ifges(4,2) * t84 + t101;
t60 = -t84 * mrSges(4,1) + t82 * mrSges(4,2);
t57 = -t85 * mrSges(4,1) - mrSges(4,3) * t98;
t56 = t85 * mrSges(4,2) - mrSges(4,3) * t99;
t50 = t89 * t83;
t49 = Ifges(6,4) * t52;
t48 = Ifges(5,5) * t52;
t47 = Ifges(5,6) * t51;
t46 = Ifges(6,6) * t51;
t41 = -Ifges(4,5) * t85 + (Ifges(4,1) * t84 - t101) * t83;
t40 = -Ifges(4,6) * t85 + (-Ifges(4,2) * t82 + t100) * t83;
t31 = -t82 * t102 + t54;
t29 = -t85 * mrSges(5,1) + t43 * mrSges(5,3);
t28 = t85 * mrSges(5,2) - t42 * mrSges(5,3);
t27 = -t42 * mrSges(6,2) - t85 * mrSges(6,3);
t21 = Ifges(5,1) * t52 - Ifges(5,4) * t51;
t20 = Ifges(6,1) * t52 + Ifges(6,5) * t51;
t19 = Ifges(5,4) * t52 - Ifges(5,2) * t51;
t18 = Ifges(6,5) * t52 + Ifges(6,3) * t51;
t12 = t51 * pkin(4) - t52 * qJ(5) + t70;
t9 = -Ifges(5,1) * t43 - Ifges(5,4) * t42 - Ifges(5,5) * t85;
t8 = -Ifges(6,1) * t43 - Ifges(6,4) * t85 + Ifges(6,5) * t42;
t7 = -Ifges(5,4) * t43 - Ifges(5,2) * t42 - Ifges(5,6) * t85;
t6 = -Ifges(6,5) * t43 - Ifges(6,6) * t85 + Ifges(6,3) * t42;
t5 = t42 * pkin(4) + t43 * qJ(5) + t58;
t2 = t85 * pkin(4) - t3;
t1 = -t85 * qJ(5) + t4;
t13 = [0.2e1 * t1 * t27 + 0.2e1 * t5 * t10 + 0.2e1 * t58 * t11 + 0.2e1 * t2 * t30 + 0.2e1 * t4 * t28 + 0.2e1 * t3 * t29 + 0.2e1 * t31 * t57 + 0.2e1 * t32 * t56 + Ifges(2,3) - (t8 + t9) * t43 + (t6 - t7) * t42 + (t77 + t79) * mrSges(3,3) * t103 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t83 + t103 * t50 - t82 * t40 + t84 * t41) * t83 + m(3) * (pkin(1) ^ 2 + t79 * t87 + t74) + m(4) * (t31 ^ 2 + t32 ^ 2 + t74) + m(5) * (t3 ^ 2 + t4 ^ 2 + t58 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t92) * t85 + (Ifges(4,6) * t82 + (2 * Ifges(3,4))) * t83 + t90) * t85; -pkin(2) * t50 + t12 * t10 + t70 * t11 + t5 * t16 + t58 * t17 + (-t73 / 0.2e1 - t72 / 0.2e1 - t48 / 0.2e1 + t47 / 0.2e1 - t49 / 0.2e1 - t46 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t85 - (t20 / 0.2e1 + t21 / 0.2e1) * t43 + (t18 / 0.2e1 - t19 / 0.2e1) * t42 + (t27 + t28) * t26 + (-t29 + t30) * t24 + (t40 / 0.2e1 + pkin(7) * t56 + t32 * mrSges(4,3)) * t84 + (-pkin(7) * t57 - t31 * mrSges(4,3) + t41 / 0.2e1) * t82 + (Ifges(3,5) - t82 * t62 / 0.2e1 + t84 * t63 / 0.2e1 + (t60 - mrSges(3,1)) * pkin(6)) * t83 + m(4) * (-pkin(2) * t75 + (-t31 * t82 + t32 * t84) * pkin(7)) + m(5) * (-t24 * t3 + t26 * t4 + t70 * t58) + m(6) * (t26 * t1 + t12 * t5 + t24 * t2) + (t2 * mrSges(6,2) - t3 * mrSges(5,3) + t8 / 0.2e1 + t9 / 0.2e1) * t52 + (-t4 * mrSges(5,3) - t1 * mrSges(6,2) + t6 / 0.2e1 - t7 / 0.2e1) * t51; -0.2e1 * pkin(2) * t60 + 0.2e1 * t12 * t16 + 0.2e1 * t70 * t17 + t84 * t62 + t82 * t63 + Ifges(3,3) + m(6) * (t12 ^ 2 + t93) + m(5) * (t70 ^ 2 + t93) + m(4) * (t95 * pkin(7) ^ 2 + pkin(2) ^ 2) + (t20 + t21) * t52 + (t18 - t19) * t51 + 0.2e1 * t95 * pkin(7) * mrSges(4,3) + 0.2e1 * (t24 * t52 - t26 * t51) * (mrSges(6,2) + mrSges(5,3)); -Ifges(4,6) * t99 - t2 * mrSges(6,1) + t3 * mrSges(5,1) - t4 * mrSges(5,2) + t69 * t30 + t1 * mrSges(6,3) + t67 * t27 + m(6) * (t67 * t1 + t69 * t2) + t31 * mrSges(4,1) - t32 * mrSges(4,2) - t92 * t85 + (t81 * t29 + t80 * t28 + m(5) * (t3 * t81 + t4 * t80)) * pkin(3) - t90; m(6) * (t69 * t24 + t67 * t26) - t26 * mrSges(5,2) - t24 * mrSges(5,1) - t24 * mrSges(6,1) + t26 * mrSges(6,3) + t46 + t48 + t49 - t47 + t73 + t72 - t89 * pkin(7) + (-t67 * t51 + t69 * t52) * mrSges(6,2) + (m(5) * (-t24 * t81 + t26 * t80) + (-t80 * t51 - t81 * t52) * mrSges(5,3)) * pkin(3); -0.2e1 * t69 * mrSges(6,1) + 0.2e1 * t67 * mrSges(6,3) + m(6) * (t67 ^ 2 + t69 ^ 2) + t92 + (0.2e1 * t81 * mrSges(5,1) - 0.2e1 * t80 * mrSges(5,2) + m(5) * (t80 ^ 2 + t81 ^ 2) * pkin(3)) * pkin(3); m(5) * t58 + m(6) * t5 + t10 + t11; m(5) * t70 + m(6) * t12 + t16 + t17; 0; m(5) + m(6); m(6) * t2 + t30; m(6) * t24 + t52 * mrSges(6,2); m(6) * t69 - mrSges(6,1); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
