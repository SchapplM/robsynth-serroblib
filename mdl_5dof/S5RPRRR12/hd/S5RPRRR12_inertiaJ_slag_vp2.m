% Calculate joint inertia matrix for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:30
% EndTime: 2019-12-31 19:12:31
% DurationCPUTime: 0.51s
% Computational Cost: add. (622->139), mult. (1120->203), div. (0->0), fcn. (1021->6), ass. (0->63)
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t77 = t53 ^ 2 + t56 ^ 2;
t99 = mrSges(6,3) * t77;
t37 = -t56 * mrSges(6,1) + t53 * mrSges(6,2);
t98 = -mrSges(5,1) + t37;
t58 = cos(qJ(3));
t97 = t58 ^ 2;
t54 = sin(qJ(4));
t55 = sin(qJ(3));
t57 = cos(qJ(4));
t33 = t54 * t55 - t57 * t58;
t34 = t54 * t58 + t57 * t55;
t82 = t53 * mrSges(6,3);
t11 = -t34 * mrSges(6,2) + t33 * t82;
t83 = t33 * t56;
t12 = t34 * mrSges(6,1) + mrSges(6,3) * t83;
t96 = t56 * t11 - t53 * t12;
t59 = -pkin(1) - pkin(6);
t89 = -pkin(7) + t59;
t36 = t89 * t55;
t73 = t89 * t58;
t14 = t54 * t36 - t57 * t73;
t68 = -mrSges(6,1) * t53 - mrSges(6,2) * t56;
t9 = t68 * t33;
t95 = m(6) * t14 + t9;
t41 = t55 * pkin(3) + qJ(2);
t10 = t34 * pkin(4) + t33 * pkin(8) + t41;
t16 = t57 * t36 + t54 * t73;
t2 = t56 * t10 - t53 * t16;
t3 = t53 * t10 + t56 * t16;
t90 = t3 * t56;
t69 = -t2 * t53 + t90;
t94 = m(6) * t69 + t96;
t93 = t14 ^ 2;
t30 = t33 ^ 2;
t92 = -2 * mrSges(5,3);
t87 = Ifges(6,4) * t53;
t86 = Ifges(6,4) * t56;
t85 = t33 * t14;
t84 = t33 * t53;
t79 = -Ifges(6,5) * t83 + Ifges(6,3) * t34;
t78 = Ifges(6,5) * t53 + Ifges(6,6) * t56;
t76 = t55 ^ 2 + t97;
t38 = Ifges(6,2) * t56 + t87;
t39 = Ifges(6,1) * t53 + t86;
t75 = t56 * t38 + t53 * t39 + Ifges(5,3);
t74 = m(4) * t76;
t72 = t77 * pkin(8);
t43 = t54 * pkin(3) + pkin(8);
t71 = t77 * t43;
t70 = t76 * mrSges(4,3);
t67 = t33 * t57 - t34 * t54;
t66 = 0.2e1 * t99;
t65 = t98 * t33 + (-mrSges(5,2) + t99) * t34;
t64 = (t57 * mrSges(5,1) - t54 * mrSges(5,2)) * pkin(3);
t6 = Ifges(6,6) * t34 + (Ifges(6,2) * t53 - t86) * t33;
t7 = Ifges(6,5) * t34 + (-Ifges(6,1) * t56 + t87) * t33;
t63 = -t16 * mrSges(5,2) + mrSges(6,3) * t90 - t2 * t82 + t53 * t7 / 0.2e1 + t38 * t84 / 0.2e1 - t39 * t83 / 0.2e1 - Ifges(5,5) * t33 + t56 * t6 / 0.2e1 + (t78 / 0.2e1 - Ifges(5,6)) * t34 + t98 * t14;
t60 = qJ(2) ^ 2;
t44 = -t57 * pkin(3) - pkin(4);
t29 = t34 ^ 2;
t1 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t14 * t9 + Ifges(4,1) * t97 - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t59 * t70 + (0.2e1 * t41 * mrSges(5,1) + Ifges(5,2) * t34 + t16 * t92 + t79) * t34 + (-0.2e1 * t41 * mrSges(5,2) + t14 * t92 + Ifges(5,1) * t33 + t53 * t6 - t56 * t7 + (Ifges(6,6) * t53 + (2 * Ifges(5,4))) * t34) * t33 + m(6) * (t2 ^ 2 + t3 ^ 2 + t93) + m(5) * (t16 ^ 2 + t41 ^ 2 + t93) + m(4) * (t76 * t59 ^ 2 + t60) + m(3) * ((pkin(1) ^ 2) + t60) + (-0.2e1 * Ifges(4,4) * t58 + Ifges(4,2) * t55) * t55 + 0.2e1 * (t55 * mrSges(4,1) + t58 * mrSges(4,2) + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t30 * mrSges(5,3) + t33 * t9 + mrSges(3,2) - t70 + (-mrSges(5,3) * t34 + t96) * t34 + m(6) * (t69 * t34 + t85) + m(5) * (t34 * t16 + t85) + t59 * t74; m(3) + t74 + m(5) * (t29 + t30) + m(6) * (t77 * t29 + t30); t63 + (m(5) * (-t14 * t57 + t16 * t54) + t67 * mrSges(5,3)) * pkin(3) + (t59 * mrSges(4,1) + Ifges(4,5)) * t58 + (-t59 * mrSges(4,2) - Ifges(4,6)) * t55 + t95 * t44 + t94 * t43; t58 * mrSges(4,1) - t55 * mrSges(4,2) + m(6) * (t44 * t33 + t34 * t71) - m(5) * t67 * pkin(3) + t65; 0.2e1 * t44 * t37 + Ifges(4,3) + 0.2e1 * t64 + t43 * t66 + m(6) * (t77 * t43 ^ 2 + t44 ^ 2) + m(5) * (t54 ^ 2 + t57 ^ 2) * pkin(3) ^ 2 + t75; -t95 * pkin(4) + t94 * pkin(8) + t63; m(6) * (-pkin(4) * t33 + t34 * t72) + t65; m(6) * (-pkin(4) * t44 + pkin(8) * t71) + (t44 - pkin(4)) * t37 + t64 + (t71 + t72) * mrSges(6,3) + t75; -0.2e1 * pkin(4) * t37 + m(6) * (t77 * pkin(8) ^ 2 + pkin(4) ^ 2) + pkin(8) * t66 + t75; t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,6) * t84 + t79; t68 * t34; t68 * t43 + t78; t68 * pkin(8) + t78; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
