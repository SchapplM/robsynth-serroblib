% Calculate joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:10
% DurationCPUTime: 0.51s
% Computational Cost: add. (698->157), mult. (1278->213), div. (0->0), fcn. (1171->6), ass. (0->64)
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t82 = t52 ^ 2 + t55 ^ 2;
t74 = t82 * mrSges(6,3);
t35 = t52 * mrSges(6,1) + t55 * mrSges(6,2);
t94 = mrSges(5,3) + t35;
t53 = sin(qJ(3));
t41 = t53 * pkin(2) + qJ(4);
t93 = t41 ^ 2;
t54 = sin(qJ(2));
t56 = cos(qJ(3));
t57 = cos(qJ(2));
t32 = t53 * t54 - t56 * t57;
t33 = t53 * t57 + t56 * t54;
t45 = -t57 * pkin(2) - pkin(1);
t64 = -t33 * qJ(4) + t45;
t12 = t32 * pkin(3) + t64;
t92 = -0.2e1 * t12;
t91 = 0.2e1 * t45;
t90 = pkin(3) + pkin(8);
t89 = -pkin(7) - pkin(6);
t88 = Ifges(6,4) * t52;
t87 = Ifges(6,4) * t55;
t86 = t32 * t52;
t85 = t32 * t55;
t84 = t55 * mrSges(6,1);
t83 = mrSges(5,1) + mrSges(4,3);
t81 = t54 ^ 2 + t57 ^ 2;
t80 = qJ(4) * t41;
t78 = Ifges(6,5) * t86 + Ifges(6,6) * t85 + Ifges(6,3) * t33;
t77 = t89 * t54;
t38 = t89 * t57;
t19 = -t53 * t38 - t56 * t77;
t21 = -t56 * t38 + t53 * t77;
t76 = t19 ^ 2 + t21 ^ 2;
t44 = -t56 * pkin(2) - pkin(3);
t75 = m(6) * t82;
t73 = t82 * t90;
t46 = Ifges(6,5) * t55;
t72 = -Ifges(6,6) * t52 + t46;
t71 = 0.2e1 * t94;
t8 = t90 * t32 + t64;
t9 = t33 * pkin(4) + t19;
t1 = -t52 * t8 + t55 * t9;
t2 = t52 * t9 + t55 * t8;
t70 = t55 * t1 + t52 * t2;
t69 = mrSges(5,2) - t74;
t68 = -t52 * mrSges(6,2) + t84;
t13 = t33 * mrSges(6,1) - mrSges(6,3) * t86;
t14 = -t33 * mrSges(6,2) + mrSges(6,3) * t85;
t67 = t55 * t13 + t52 * t14;
t66 = -0.2e1 * t74;
t36 = -Ifges(6,2) * t52 + t87;
t37 = Ifges(6,1) * t55 - t88;
t65 = -t52 * t36 + t55 * t37 + Ifges(5,1) + Ifges(4,3);
t63 = (t56 * mrSges(4,1) - t53 * mrSges(4,2)) * pkin(2);
t10 = -t32 * pkin(4) + t21;
t6 = Ifges(6,6) * t33 + (Ifges(6,2) * t55 + t88) * t32;
t7 = Ifges(6,5) * t33 + (Ifges(6,1) * t52 + t87) * t32;
t62 = t37 * t86 / 0.2e1 + t36 * t85 / 0.2e1 - t52 * t6 / 0.2e1 + t55 * t7 / 0.2e1 + t10 * t35 + (-mrSges(4,2) + mrSges(5,3)) * t21 - t70 * mrSges(6,3) + (Ifges(5,5) - Ifges(4,6)) * t32 + (mrSges(5,2) - mrSges(4,1)) * t19 + (t72 / 0.2e1 + Ifges(4,5) - Ifges(5,4)) * t33;
t59 = qJ(4) ^ 2;
t40 = -pkin(8) + t44;
t11 = t68 * t32;
t3 = [t54 * (Ifges(3,1) * t54 + Ifges(3,4) * t57) - 0.2e1 * pkin(1) * (-t57 * mrSges(3,1) + t54 * mrSges(3,2)) + t57 * (Ifges(3,4) * t54 + Ifges(3,2) * t57) - 0.2e1 * t10 * t11 + 0.2e1 * t1 * t13 + 0.2e1 * t2 * t14 + Ifges(2,3) + 0.2e1 * t81 * pkin(6) * mrSges(3,3) + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(5) * (t12 ^ 2 + t76) + m(4) * (t45 ^ 2 + t76) + m(3) * (t81 * pkin(6) ^ 2 + pkin(1) ^ 2) + (mrSges(4,2) * t91 + mrSges(5,3) * t92 + (Ifges(5,2) + Ifges(4,1)) * t33 + 0.2e1 * t83 * t19 + t78) * t33 + (mrSges(4,1) * t91 + mrSges(5,2) * t92 + t52 * t7 + t55 * t6 + (Ifges(5,3) + Ifges(4,2)) * t32 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t33 - 0.2e1 * t83 * t21) * t32; t62 + (m(4) * (-t19 * t56 + t21 * t53) + (-t53 * t32 - t56 * t33) * mrSges(4,3)) * pkin(2) + m(6) * (t41 * t10 + t40 * t70) + t67 * t40 + m(5) * (t44 * t19 + t41 * t21) + (-t54 * mrSges(3,1) - t57 * mrSges(3,2)) * pkin(6) + (-t41 * t32 + t44 * t33) * mrSges(5,1) + Ifges(3,5) * t54 + Ifges(3,6) * t57 - t41 * t11; 0.2e1 * t44 * mrSges(5,2) + Ifges(3,3) + t41 * t71 + 0.2e1 * t63 + t40 * t66 + m(6) * (t82 * t40 ^ 2 + t93) + m(5) * (t44 ^ 2 + t93) + m(4) * (t53 ^ 2 + t56 ^ 2) * pkin(2) ^ 2 + t65; t62 + m(6) * (qJ(4) * t10 - t70 * t90) - t67 * t90 + m(5) * (-pkin(3) * t19 + qJ(4) * t21) + (-pkin(3) * t33 - qJ(4) * t32) * mrSges(5,1) - qJ(4) * t11; t63 + (-pkin(3) + t44) * mrSges(5,2) + m(6) * (-t40 * t73 + t80) + m(5) * (-pkin(3) * t44 + t80) + t65 + (-t40 + t90) * t74 + t94 * (t41 + qJ(4)); -0.2e1 * pkin(3) * mrSges(5,2) + qJ(4) * t71 - t90 * t66 + m(6) * (t82 * t90 ^ 2 + t59) + m(5) * (pkin(3) ^ 2 + t59) + t65; m(5) * t19 + m(6) * t70 + t33 * mrSges(5,1) + t67; m(5) * t44 + t40 * t75 + t69; -m(5) * pkin(3) - m(6) * t73 + t69; m(5) + t75; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t78; t40 * t68 + t72; -t90 * t84 + t46 + (mrSges(6,2) * t90 - Ifges(6,6)) * t52; t68; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
