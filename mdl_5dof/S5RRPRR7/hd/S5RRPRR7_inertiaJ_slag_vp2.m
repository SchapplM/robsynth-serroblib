% Calculate joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (478->117), mult. (762->153), div. (0->0), fcn. (599->6), ass. (0->58)
t40 = sin(qJ(5));
t41 = sin(qJ(4));
t43 = cos(qJ(5));
t44 = cos(qJ(4));
t21 = -t40 * t41 + t43 * t44;
t17 = t21 ^ 2;
t82 = t44 ^ 2;
t19 = -t40 * t44 - t43 * t41;
t81 = t19 ^ 2 + t17;
t69 = t41 ^ 2 + t82;
t63 = t69 * mrSges(5,3);
t80 = t41 * mrSges(5,1) + t44 * mrSges(5,2) + mrSges(4,3);
t73 = t19 * t40;
t79 = pkin(4) * mrSges(6,3) * t73 + Ifges(5,5) * t44;
t7 = -t19 * mrSges(6,1) + t21 * mrSges(6,2);
t78 = 0.2e1 * t7;
t42 = sin(qJ(2));
t29 = t42 * pkin(1) + qJ(3);
t76 = t29 ^ 2;
t45 = cos(qJ(2));
t34 = -t45 * pkin(1) - pkin(2);
t28 = -pkin(7) + t34;
t75 = -pkin(8) + t28;
t46 = -pkin(2) - pkin(7);
t74 = -pkin(8) + t46;
t72 = t21 * t43;
t71 = t44 * mrSges(5,1);
t70 = Ifges(6,5) * t21 + Ifges(6,6) * t19;
t68 = qJ(3) * t29;
t66 = 0.2e1 * mrSges(6,3);
t65 = mrSges(6,3) * t72;
t64 = m(5) * t69;
t62 = t69 * t46;
t61 = t21 * mrSges(6,1) + t19 * mrSges(6,2);
t60 = 0.2e1 * t80;
t12 = t75 * t41;
t13 = t75 * t44;
t4 = -t40 * t12 + t43 * t13;
t5 = t43 * t12 + t40 * t13;
t59 = t19 * t5 - t21 * t4;
t22 = t74 * t41;
t23 = t74 * t44;
t8 = -t40 * t22 + t43 * t23;
t9 = t43 * t22 + t40 * t23;
t58 = t19 * t9 - t21 * t8;
t57 = -t41 * mrSges(5,2) + t71;
t56 = t4 * mrSges(6,1) - t5 * mrSges(6,2) + t70;
t55 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + t70;
t54 = -0.2e1 * t63;
t53 = (t45 * mrSges(3,1) - t42 * mrSges(3,2)) * pkin(1);
t52 = (t43 * mrSges(6,1) - t40 * mrSges(6,2)) * pkin(4);
t51 = Ifges(5,1) * t82 + Ifges(6,1) * t17 + Ifges(4,1) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t44 + Ifges(5,2) * t41) * t41 + (0.2e1 * Ifges(6,4) * t21 + Ifges(6,2) * t19) * t19;
t50 = -t81 * mrSges(6,3) + mrSges(4,2) - t63;
t47 = qJ(3) ^ 2;
t36 = t41 * pkin(4);
t31 = qJ(3) + t36;
t24 = t29 + t36;
t1 = [0.2e1 * t34 * mrSges(4,2) + t24 * t78 + Ifges(2,3) + t29 * t60 + 0.2e1 * t53 + t59 * t66 + t28 * t54 + m(6) * (t24 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(5) * (t69 * t28 ^ 2 + t76) + m(4) * (t34 ^ 2 + t76) + m(3) * (t42 ^ 2 + t45 ^ 2) * pkin(1) ^ 2 + t51; (t31 + t24) * t7 + t53 + (t34 - pkin(2)) * mrSges(4,2) + m(6) * (t31 * t24 + t8 * t4 + t9 * t5) + m(4) * (-pkin(2) * t34 + t68) + m(5) * (t28 * t62 + t68) + ((-t4 - t8) * t21 + (t5 + t9) * t19) * mrSges(6,3) + t51 + (-t28 - t46) * t63 + t80 * (t29 + qJ(3)); -0.2e1 * pkin(2) * mrSges(4,2) + t31 * t78 + qJ(3) * t60 + t58 * t66 + t46 * t54 + m(6) * (t31 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t69 * t46 ^ 2 + t47) + m(4) * (pkin(2) ^ 2 + t47) + t51; m(4) * t34 - m(6) * t59 + t28 * t64 + t50; -m(4) * pkin(2) + m(5) * t62 - m(6) * t58 + t50; m(6) * t81 + m(4) + t64; -Ifges(5,6) * t41 + t57 * t28 + (-t65 + m(6) * (t4 * t43 + t40 * t5)) * pkin(4) + t56 + t79; t46 * t71 + (-t46 * mrSges(5,2) - Ifges(5,6)) * t41 + (-t65 + m(6) * (t40 * t9 + t43 * t8)) * pkin(4) + t55 + t79; m(6) * (t72 - t73) * pkin(4) + t57 + t61; Ifges(5,3) + Ifges(6,3) + m(6) * (t40 ^ 2 + t43 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t52; t56; t55; t61; Ifges(6,3) + t52; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
