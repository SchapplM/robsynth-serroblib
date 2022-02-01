% Calculate joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:26
% EndTime: 2022-01-20 10:05:26
% DurationCPUTime: 0.38s
% Computational Cost: add. (454->116), mult. (847->158), div. (0->0), fcn. (654->8), ass. (0->61)
t42 = sin(pkin(9));
t48 = cos(qJ(5));
t66 = t42 * t48;
t78 = Ifges(6,5) * t66;
t38 = t42 ^ 2;
t44 = cos(pkin(9));
t39 = t44 ^ 2;
t61 = t38 + t39;
t77 = t61 * mrSges(5,3);
t49 = cos(qJ(2));
t36 = t49 * pkin(1) + pkin(2);
t43 = sin(pkin(8));
t45 = cos(pkin(8));
t47 = sin(qJ(2));
t72 = pkin(1) * t47;
t16 = t43 * t36 + t45 * t72;
t13 = qJ(4) + t16;
t34 = t43 * pkin(2) + qJ(4);
t76 = t13 + t34;
t46 = sin(qJ(5));
t67 = t42 * t46;
t23 = t44 * mrSges(6,2) - mrSges(6,3) * t67;
t75 = 0.2e1 * t23;
t24 = -t44 * mrSges(6,1) - mrSges(6,3) * t66;
t74 = 0.2e1 * t24;
t25 = -t44 * mrSges(5,1) + t42 * mrSges(5,2);
t73 = 0.2e1 * t25;
t71 = t45 * pkin(2);
t70 = t13 * t34;
t15 = t45 * t36 - t43 * t72;
t69 = t15 * mrSges(4,1);
t68 = t16 * mrSges(4,2);
t65 = t44 * t46;
t64 = t44 * t48;
t63 = t46 * (-Ifges(6,6) * t44 + (Ifges(6,4) * t48 - Ifges(6,2) * t46) * t42);
t62 = t46 * t24;
t60 = t46 ^ 2 + t48 ^ 2;
t19 = -mrSges(6,1) * t67 - mrSges(6,2) * t66;
t59 = t44 * t19 + t23 * t66;
t58 = t45 * mrSges(4,1) - t43 * mrSges(4,2);
t57 = t46 * t23 + t48 * t24 + t25;
t56 = -t44 * pkin(4) - t42 * pkin(7) - pkin(3);
t55 = (t49 * mrSges(3,1) - t47 * mrSges(3,2)) * pkin(1);
t9 = -Ifges(6,6) * t67 - Ifges(6,3) * t44 + t78;
t54 = Ifges(3,3) + Ifges(4,3) + ((Ifges(6,1) * t48 - Ifges(6,4) * t46) * t66 + Ifges(5,1) * t42) * t42 + (0.2e1 * Ifges(5,4) * t42 + Ifges(5,2) * t44 - t78 - t9) * t44;
t53 = -0.2e1 * t42 * t19 + 0.2e1 * t77;
t52 = -t42 * t63 + t54;
t35 = -pkin(3) - t71;
t33 = t34 ^ 2;
t26 = t38 * t33;
t20 = t56 - t71;
t14 = -pkin(3) - t15;
t11 = t13 ^ 2;
t8 = t38 * t11;
t6 = -t15 + t56;
t5 = t38 * t70;
t4 = t46 * t20 + t34 * t64;
t3 = t48 * t20 - t34 * t65;
t2 = t13 * t64 + t46 * t6;
t1 = -t13 * t65 + t48 * t6;
t7 = [t52 + t53 * t13 + 0.2e1 * t55 + m(3) * (t47 ^ 2 + t49 ^ 2) * pkin(1) ^ 2 + m(6) * (t1 ^ 2 + t2 ^ 2 + t8) + m(5) * (t39 * t11 + t14 ^ 2 + t8) + m(4) * (t15 ^ 2 + t16 ^ 2) + Ifges(2,3) + 0.2e1 * t69 - 0.2e1 * t68 + t2 * t75 + t1 * t74 + t14 * t73; t69 - t68 + (t14 + t35) * t25 + (t3 + t1) * t24 + (t4 + t2) * t23 + t55 + m(6) * (t3 * t1 + t4 * t2 + t5) + m(5) * (t35 * t14 + t39 * t70 + t5) + (-t76 * t19 - t63) * t42 + (m(4) * (t15 * t45 + t16 * t43) + t58) * pkin(2) + t54 + t76 * t77; t4 * t75 + t3 * t74 + t35 * t73 + t53 * t34 + m(6) * (t3 ^ 2 + t4 ^ 2 + t26) + m(5) * (t39 * t33 + t35 ^ 2 + t26) + t52 + (0.2e1 * t58 + m(4) * (t43 ^ 2 + t45 ^ 2) * pkin(2)) * pkin(2); (m(6) * (-t1 * t46 - t13 * t44 + t2 * t48) - t62) * t42 + t59; (m(6) * (-t3 * t46 - t34 * t44 + t4 * t48) - t62) * t42 + t59; m(4) + m(5) * t61 + m(6) * (t60 * t38 + t39); m(6) * (t48 * t1 + t46 * t2) + m(5) * t14 + t57; m(6) * (t48 * t3 + t46 * t4) + m(5) * t35 + t57; 0; m(6) * t60 + m(5); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t9; t3 * mrSges(6,1) - t4 * mrSges(6,2) + t9; t19; t48 * mrSges(6,1) - t46 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
