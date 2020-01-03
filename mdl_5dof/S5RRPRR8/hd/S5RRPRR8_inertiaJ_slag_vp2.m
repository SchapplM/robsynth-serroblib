% Calculate joint inertia matrix for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:54
% EndTime: 2019-12-31 20:16:56
% DurationCPUTime: 0.51s
% Computational Cost: add. (1092->156), mult. (2075->232), div. (0->0), fcn. (2264->8), ass. (0->66)
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t39 = -t58 * t62 + t59 * t64;
t40 = t58 * t64 + t59 * t62;
t61 = sin(qJ(4));
t92 = cos(qJ(4));
t27 = -t92 * t39 + t61 * t40;
t28 = t61 * t39 + t92 * t40;
t60 = sin(qJ(5));
t85 = t60 * mrSges(6,3);
t15 = -t27 * mrSges(6,2) - t28 * t85;
t63 = cos(qJ(5));
t88 = t28 * t63;
t16 = t27 * mrSges(6,1) - mrSges(6,3) * t88;
t99 = t63 * t15 - t60 * t16;
t82 = -qJ(3) - pkin(6);
t75 = t82 * t62;
t76 = t82 * t64;
t30 = t58 * t75 - t59 * t76;
t21 = t39 * pkin(7) + t30;
t29 = t58 * t76 + t59 * t75;
t69 = -t40 * pkin(7) + t29;
t11 = t61 * t21 - t92 * t69;
t98 = t11 ^ 2;
t97 = 0.2e1 * t11;
t51 = -t64 * pkin(2) - pkin(1);
t31 = -t39 * pkin(3) + t51;
t96 = 0.2e1 * t31;
t95 = 0.2e1 * t39;
t94 = pkin(2) * t58;
t10 = t27 * pkin(4) - t28 * pkin(8) + t31;
t13 = t92 * t21 + t61 * t69;
t3 = t60 * t10 + t63 * t13;
t93 = t3 * t63;
t91 = Ifges(6,4) * t60;
t90 = Ifges(6,4) * t63;
t89 = t28 * t60;
t50 = t59 * pkin(2) + pkin(3);
t35 = t92 * t50 - t61 * t94;
t87 = t35 * mrSges(5,1);
t36 = t61 * t50 + t92 * t94;
t86 = t36 * mrSges(5,2);
t81 = Ifges(6,5) * t88 + Ifges(6,3) * t27;
t80 = Ifges(6,5) * t60 + Ifges(6,6) * t63;
t79 = t60 ^ 2 + t63 ^ 2;
t78 = t62 ^ 2 + t64 ^ 2;
t46 = Ifges(6,2) * t63 + t91;
t47 = Ifges(6,1) * t60 + t90;
t77 = t63 * t46 + t60 * t47 + Ifges(5,3);
t34 = pkin(8) + t36;
t74 = t79 * t34;
t73 = -t39 * mrSges(4,1) + t40 * mrSges(4,2);
t2 = t63 * t10 - t60 * t13;
t72 = -t2 * t60 + t93;
t71 = mrSges(6,1) * t60 + mrSges(6,2) * t63;
t70 = 0.2e1 * t79 * mrSges(6,3);
t45 = -t63 * mrSges(6,1) + t60 * mrSges(6,2);
t6 = Ifges(6,6) * t27 + (-Ifges(6,2) * t60 + t90) * t28;
t7 = Ifges(6,5) * t27 + (Ifges(6,1) * t63 - t91) * t28;
t68 = -t13 * mrSges(5,2) + mrSges(6,3) * t93 - t2 * t85 + t60 * t7 / 0.2e1 - t46 * t89 / 0.2e1 + t47 * t88 / 0.2e1 + Ifges(5,5) * t28 + t63 * t6 / 0.2e1 + (t80 / 0.2e1 - Ifges(5,6)) * t27 + (t45 - mrSges(5,1)) * t11;
t33 = -pkin(4) - t35;
t23 = t28 * mrSges(5,2);
t14 = t71 * t28;
t1 = [t62 * (Ifges(3,1) * t62 + Ifges(3,4) * t64) + t64 * (Ifges(3,4) * t62 + Ifges(3,2) * t64) - 0.2e1 * pkin(1) * (-t64 * mrSges(3,1) + t62 * mrSges(3,2)) + 0.2e1 * t51 * t73 + Ifges(4,2) * t39 ^ 2 + t23 * t96 + t14 * t97 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + t30 * mrSges(4,3) * t95 + (mrSges(5,1) * t96 - 0.2e1 * t13 * mrSges(5,3) + Ifges(5,2) * t27 + t81) * t27 + 0.2e1 * t78 * pkin(6) * mrSges(3,3) + (-0.2e1 * t29 * mrSges(4,3) + Ifges(4,1) * t40 + Ifges(4,4) * t95) * t40 + (mrSges(5,3) * t97 + Ifges(5,1) * t28 - t60 * t6 + t63 * t7 + (-Ifges(6,6) * t60 - (2 * Ifges(5,4))) * t27) * t28 + m(3) * (t78 * pkin(6) ^ 2 + pkin(1) ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t98) + m(5) * (t13 ^ 2 + t31 ^ 2 + t98) + m(4) * (t29 ^ 2 + t30 ^ 2 + t51 ^ 2); t68 + (m(4) * (t29 * t59 + t30 * t58) + (t58 * t39 - t59 * t40) * mrSges(4,3)) * pkin(2) + m(6) * (t33 * t11 + t34 * t72) + t99 * t34 + m(5) * (-t35 * t11 + t36 * t13) + (-t62 * mrSges(3,1) - t64 * mrSges(3,2)) * pkin(6) + (-t36 * t27 - t35 * t28) * mrSges(5,3) + Ifges(3,6) * t64 + Ifges(3,5) * t62 + t33 * t14 + Ifges(4,6) * t39 + Ifges(4,5) * t40 + t29 * mrSges(4,1) - t30 * mrSges(4,2); 0.2e1 * t87 - 0.2e1 * t86 + 0.2e1 * t33 * t45 + Ifges(3,3) + Ifges(4,3) + t34 * t70 + m(6) * (t79 * t34 ^ 2 + t33 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2) + t77 + (0.2e1 * t59 * mrSges(4,1) - 0.2e1 * t58 * mrSges(4,2) + m(4) * (t58 ^ 2 + t59 ^ 2) * pkin(2)) * pkin(2); t27 * mrSges(5,1) + t60 * t15 + t63 * t16 + t23 + m(6) * (t63 * t2 + t60 * t3) + m(5) * t31 + m(4) * t51 + t73; 0; m(6) * t79 + m(4) + m(5); t68 + (-m(6) * t11 - t14) * pkin(4) + (m(6) * t72 + t99) * pkin(8); m(6) * (-pkin(4) * t33 + pkin(8) * t74) + t87 - t86 + (-pkin(4) + t33) * t45 + (t79 * pkin(8) + t74) * mrSges(6,3) + t77; 0; -0.2e1 * pkin(4) * t45 + m(6) * (t79 * pkin(8) ^ 2 + pkin(4) ^ 2) + pkin(8) * t70 + t77; t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t89 + t81; -t34 * t71 + t80; -t45; -pkin(8) * t71 + t80; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
