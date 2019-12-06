% Calculate joint inertia matrix for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:49:47
% EndTime: 2019-12-05 18:49:49
% DurationCPUTime: 0.41s
% Computational Cost: add. (674->134), mult. (1327->208), div. (0->0), fcn. (1405->8), ass. (0->71)
t57 = sin(qJ(5));
t55 = t57 ^ 2;
t89 = mrSges(6,3) * t55;
t49 = pkin(6) * t89;
t61 = cos(qJ(5));
t56 = t61 ^ 2;
t88 = mrSges(6,3) * t56;
t50 = pkin(6) * t88;
t81 = t61 * mrSges(6,1);
t39 = t57 * mrSges(6,2) - t81;
t91 = pkin(4) * t39;
t95 = t49 + t50 - t91;
t62 = cos(qJ(4));
t90 = t62 * pkin(3);
t46 = -pkin(4) - t90;
t30 = t46 * t39;
t58 = sin(qJ(4));
t45 = t58 * pkin(3) + pkin(6);
t37 = t45 * t89;
t38 = t45 * t88;
t51 = mrSges(5,1) * t90;
t94 = t30 + t37 + t38 + t51;
t59 = sin(qJ(3));
t60 = sin(qJ(2));
t63 = cos(qJ(3));
t64 = cos(qJ(2));
t35 = t59 * t60 - t63 * t64;
t48 = t64 * pkin(2) + pkin(1);
t22 = -t35 * pkin(3) + t48;
t93 = 0.2e1 * t22;
t92 = pkin(2) * t59;
t87 = Ifges(6,4) * t57;
t86 = Ifges(6,4) * t61;
t36 = -t59 * t64 - t63 * t60;
t18 = t58 * t35 + t62 * t36;
t85 = t18 * t57;
t84 = t18 * t61;
t47 = t63 * pkin(2) + pkin(3);
t28 = t58 * t47 + t62 * t92;
t83 = t28 * mrSges(5,2);
t82 = t58 * mrSges(5,2);
t17 = -t62 * t35 + t58 * t36;
t80 = Ifges(6,5) * t84 + Ifges(6,3) * t17;
t79 = Ifges(6,5) * t57 + Ifges(6,6) * t61;
t78 = t55 + t56;
t77 = pkin(3) * t82;
t40 = Ifges(6,2) * t61 + t87;
t41 = Ifges(6,1) * t57 + t86;
t76 = t61 * t40 + t57 * t41 + Ifges(5,3);
t75 = t78 * t45;
t74 = Ifges(4,3) + t76;
t8 = -t17 * mrSges(6,2) - mrSges(6,3) * t85;
t9 = t17 * mrSges(6,1) - mrSges(6,3) * t84;
t73 = -t57 * t9 + t61 * t8;
t72 = mrSges(6,1) * t57 + mrSges(6,2) * t61;
t27 = t62 * t47 - t58 * t92;
t3 = Ifges(6,6) * t17 + (-Ifges(6,2) * t57 + t86) * t18;
t4 = Ifges(6,5) * t17 + (Ifges(6,1) * t61 - t87) * t18;
t71 = t57 * t4 / 0.2e1 - t40 * t85 / 0.2e1 + t41 * t84 / 0.2e1 + Ifges(5,5) * t18 + t61 * t3 / 0.2e1 + (t79 / 0.2e1 - Ifges(5,6)) * t17;
t70 = (t63 * mrSges(4,1) - t59 * mrSges(4,2)) * pkin(2);
t69 = Ifges(4,5) * t36 + Ifges(4,6) * t35 + t71;
t25 = -pkin(4) - t27;
t19 = t25 * t39;
t26 = pkin(6) + t28;
t20 = t26 * t89;
t21 = t26 * t88;
t23 = t27 * mrSges(5,1);
t68 = t19 + t20 + t21 + t23 + t76 - t83;
t7 = t72 * t18;
t6 = t17 * pkin(4) - t18 * pkin(6) + t22;
t1 = [Ifges(2,3) + m(5) * t22 ^ 2 + m(4) * t48 ^ 2 + m(3) * pkin(1) ^ 2 + 0.2e1 * t48 * (-t35 * mrSges(4,1) + t36 * mrSges(4,2)) + t35 * (Ifges(4,4) * t36 + Ifges(4,2) * t35) + t36 * (Ifges(4,1) * t36 + Ifges(4,4) * t35) + 0.2e1 * pkin(1) * (t64 * mrSges(3,1) - t60 * mrSges(3,2)) - t60 * (-Ifges(3,1) * t60 - Ifges(3,4) * t64) - t64 * (-Ifges(3,4) * t60 - Ifges(3,2) * t64) + (mrSges(5,1) * t93 + Ifges(5,2) * t17 + t80) * t17 + (mrSges(5,2) * t93 + Ifges(5,1) * t18 - t57 * t3 + t61 * t4 + (-Ifges(6,6) * t57 - (2 * Ifges(5,4))) * t17) * t18 + (m(6) * t78 * t6 + 0.2e1 * t57 * t8 + 0.2e1 * t61 * t9) * t6; t69 + t73 * t26 + (-t28 * t17 - t27 * t18) * mrSges(5,3) + (t35 * t59 - t36 * t63) * pkin(2) * mrSges(4,3) - Ifges(3,5) * t60 - Ifges(3,6) * t64 + t25 * t7; -0.2e1 * t83 + Ifges(3,3) + 0.2e1 * t19 + 0.2e1 * t20 + 0.2e1 * t21 + 0.2e1 * t23 + 0.2e1 * t70 + m(6) * (t78 * t26 ^ 2 + t25 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t59 ^ 2 + t63 ^ 2) * pkin(2) ^ 2 + t74; t46 * t7 + t73 * t45 + (-t17 * t58 - t18 * t62) * pkin(3) * mrSges(5,3) + t69; t68 + m(6) * (t46 * t25 + t26 * t75) + (m(5) * (t27 * t62 + t28 * t58) - t82) * pkin(3) + t70 + Ifges(4,3) + t94; -0.2e1 * t77 + 0.2e1 * t30 + 0.2e1 * t37 + 0.2e1 * t38 + 0.2e1 * t51 + m(6) * (t78 * t45 ^ 2 + t46 ^ 2) + m(5) * (t58 ^ 2 + t62 ^ 2) * pkin(3) ^ 2 + t74; -pkin(4) * t7 + t73 * pkin(6) + t71; m(6) * (t78 * t26 * pkin(6) - pkin(4) * t25) + t68 + t95; -t77 + m(6) * (-pkin(4) * t46 + pkin(6) * t75) + t76 + t94 + t95; 0.2e1 * t50 + 0.2e1 * t49 - 0.2e1 * t91 + m(6) * (t78 * pkin(6) ^ 2 + pkin(4) ^ 2) + t76; t6 * t81 + (-mrSges(6,2) * t6 - Ifges(6,6) * t18) * t57 + t80; -t72 * t26 + t79; -t72 * t45 + t79; -t72 * pkin(6) + t79; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
