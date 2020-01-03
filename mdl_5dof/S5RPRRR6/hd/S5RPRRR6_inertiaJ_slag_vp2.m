% Calculate joint inertia matrix for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:58
% EndTime: 2019-12-31 19:00:59
% DurationCPUTime: 0.54s
% Computational Cost: add. (629->136), mult. (1199->205), div. (0->0), fcn. (1105->8), ass. (0->65)
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t76 = t55 ^ 2 + t58 ^ 2;
t99 = mrSges(6,3) * t76;
t60 = cos(qJ(3));
t98 = t60 ^ 2;
t56 = sin(qJ(4));
t57 = sin(qJ(3));
t59 = cos(qJ(4));
t33 = t56 * t57 - t59 * t60;
t35 = t56 * t60 + t59 * t57;
t81 = t55 * mrSges(6,3);
t15 = -t33 * mrSges(6,2) - t35 * t81;
t82 = t35 * t58;
t16 = t33 * mrSges(6,1) - mrSges(6,3) * t82;
t97 = t58 * t15 - t55 * t16;
t53 = sin(pkin(9));
t42 = t53 * pkin(1) + pkin(6);
t88 = pkin(7) + t42;
t25 = t88 * t60;
t73 = t88 * t57;
t11 = t56 * t25 + t59 * t73;
t68 = mrSges(6,1) * t55 + mrSges(6,2) * t58;
t14 = t68 * t35;
t96 = m(6) * t11 + t14;
t13 = t59 * t25 - t56 * t73;
t54 = cos(pkin(9));
t43 = -t54 * pkin(1) - pkin(2);
t36 = -t60 * pkin(3) + t43;
t90 = pkin(4) * t33;
t9 = -t35 * pkin(8) + t36 + t90;
t2 = -t55 * t13 + t58 * t9;
t3 = t58 * t13 + t55 * t9;
t89 = t3 * t58;
t70 = -t2 * t55 + t89;
t95 = m(6) * t70 + t97;
t94 = t11 ^ 2;
t93 = t33 ^ 2;
t92 = 0.2e1 * t11;
t91 = 0.2e1 * t36;
t86 = Ifges(6,4) * t55;
t85 = Ifges(6,4) * t58;
t84 = t11 * t33;
t83 = t35 * t55;
t78 = Ifges(6,5) * t82 + Ifges(6,3) * t33;
t77 = Ifges(6,5) * t55 + Ifges(6,6) * t58;
t75 = t57 ^ 2 + t98;
t38 = Ifges(6,2) * t58 + t86;
t39 = Ifges(6,1) * t55 + t85;
t74 = t58 * t38 + t55 * t39 + Ifges(5,3);
t72 = t76 * pkin(8);
t45 = t56 * pkin(3) + pkin(8);
t71 = t76 * t45;
t69 = -t60 * mrSges(4,1) + t57 * mrSges(4,2);
t67 = 0.2e1 * t99;
t26 = t33 * mrSges(5,1);
t37 = -t58 * mrSges(6,1) + t55 * mrSges(6,2);
t66 = t33 * t37 - t26 + (-mrSges(5,2) + t99) * t35;
t65 = (t59 * mrSges(5,1) - t56 * mrSges(5,2)) * pkin(3);
t6 = Ifges(6,6) * t33 + (-Ifges(6,2) * t55 + t85) * t35;
t7 = Ifges(6,5) * t33 + (Ifges(6,1) * t58 - t86) * t35;
t64 = -t13 * mrSges(5,2) + mrSges(6,3) * t89 - t2 * t81 + t55 * t7 / 0.2e1 - t38 * t83 / 0.2e1 + t39 * t82 / 0.2e1 + Ifges(5,5) * t35 + t58 * t6 / 0.2e1 + (t77 / 0.2e1 - Ifges(5,6)) * t33 + (-mrSges(5,1) + t37) * t11;
t46 = -t59 * pkin(3) - pkin(4);
t30 = t35 ^ 2;
t1 = [Ifges(2,3) + Ifges(3,3) + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + t14 * t92 + t26 * t91 + 0.2e1 * t43 * t69 + Ifges(4,2) * t98 + (mrSges(5,2) * t91 + mrSges(5,3) * t92 + Ifges(5,1) * t35 - t55 * t6 + t58 * t7) * t35 + (-0.2e1 * t13 * mrSges(5,3) + Ifges(5,2) * t33 + (-Ifges(6,6) * t55 - (2 * Ifges(5,4))) * t35 + t78) * t33 + m(6) * (t2 ^ 2 + t3 ^ 2 + t94) + m(5) * (t13 ^ 2 + t36 ^ 2 + t94) + m(4) * (t42 ^ 2 * t75 + t43 ^ 2) + m(3) * (t53 ^ 2 + t54 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t57 + 0.2e1 * Ifges(4,4) * t60) * t57 + 0.2e1 * (t54 * mrSges(3,1) - t53 * mrSges(3,2)) * pkin(1) + 0.2e1 * t75 * t42 * mrSges(4,3); t33 * t14 + t97 * t35 + m(6) * (t35 * t70 + t84) + m(5) * (t13 * t35 + t84); m(3) + m(6) * (t30 * t76 + t93) + m(5) * (t30 + t93) + m(4) * t75; t64 + (m(5) * (-t11 * t59 + t13 * t56) + (-t56 * t33 - t59 * t35) * mrSges(5,3)) * pkin(3) + (-t57 * mrSges(4,1) - t60 * mrSges(4,2)) * t42 + Ifges(4,6) * t60 + Ifges(4,5) * t57 + t96 * t46 + t95 * t45; m(6) * (t46 * t33 + t35 * t71) + m(5) * (-t33 * t59 + t35 * t56) * pkin(3) + t66 - t69; 0.2e1 * t46 * t37 + Ifges(4,3) + 0.2e1 * t65 + t45 * t67 + m(6) * (t45 ^ 2 * t76 + t46 ^ 2) + m(5) * (t56 ^ 2 + t59 ^ 2) * pkin(3) ^ 2 + t74; -t96 * pkin(4) + t95 * pkin(8) + t64; m(6) * (t35 * t72 - t90) + t66; m(6) * (-pkin(4) * t46 + pkin(8) * t71) + (t46 - pkin(4)) * t37 + t65 + (t71 + t72) * mrSges(6,3) + t74; -0.2e1 * pkin(4) * t37 + m(6) * (pkin(8) ^ 2 * t76 + pkin(4) ^ 2) + pkin(8) * t67 + t74; t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t83 + t78; -t14; -t45 * t68 + t77; -pkin(8) * t68 + t77; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
