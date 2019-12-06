% Calculate joint inertia matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:48
% EndTime: 2019-12-05 18:37:49
% DurationCPUTime: 0.45s
% Computational Cost: add. (1138->147), mult. (2138->219), div. (0->0), fcn. (2340->8), ass. (0->58)
t56 = sin(qJ(3));
t57 = sin(qJ(2));
t59 = cos(qJ(3));
t60 = cos(qJ(2));
t40 = -t56 * t57 + t59 * t60;
t41 = t56 * t60 + t57 * t59;
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t26 = t40 * t54 - t41 * t53;
t27 = t40 * t53 + t41 * t54;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t12 = t26 * t58 - t27 * t55;
t80 = 0.2e1 * t12;
t79 = 0.2e1 * t26;
t78 = -pkin(7) - pkin(6);
t77 = pkin(2) * t56;
t76 = pkin(3) * t53;
t49 = pkin(2) * t59 + pkin(3);
t34 = t54 * t49 - t53 * t77;
t33 = pkin(4) + t34;
t36 = t49 * t53 + t54 * t77;
t20 = t33 * t55 + t36 * t58;
t75 = t20 * mrSges(6,2);
t74 = t34 * mrSges(5,1);
t73 = t36 * mrSges(5,2);
t48 = pkin(3) * t54 + pkin(4);
t37 = t48 * t55 + t58 * t76;
t72 = t37 * mrSges(6,2);
t46 = t78 * t57;
t47 = t78 * t60;
t29 = t59 * t46 + t47 * t56;
t21 = -qJ(4) * t41 + t29;
t30 = t56 * t46 - t59 * t47;
t22 = qJ(4) * t40 + t30;
t8 = t53 * t21 + t54 * t22;
t71 = t57 ^ 2 + t60 ^ 2;
t70 = Ifges(4,3) + Ifges(5,3) + Ifges(6,3);
t50 = -pkin(2) * t60 - pkin(1);
t13 = t26 * t55 + t27 * t58;
t69 = -t12 * mrSges(6,1) + t13 * mrSges(6,2);
t68 = -t26 * mrSges(5,1) + t27 * mrSges(5,2);
t7 = t54 * t21 - t22 * t53;
t67 = t54 * mrSges(5,1) - t53 * mrSges(5,2);
t4 = -pkin(8) * t27 + t7;
t5 = pkin(8) * t26 + t8;
t2 = t4 * t58 - t5 * t55;
t3 = t4 * t55 + t5 * t58;
t66 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t12;
t31 = -pkin(3) * t40 + t50;
t65 = (mrSges(4,1) * t59 - mrSges(4,2) * t56) * pkin(2);
t64 = t29 * mrSges(4,1) + t7 * mrSges(5,1) - t30 * mrSges(4,2) - t8 * mrSges(5,2) + Ifges(4,5) * t41 + Ifges(5,5) * t27 + Ifges(4,6) * t40 + Ifges(5,6) * t26 + t66;
t35 = t48 * t58 - t55 * t76;
t32 = t35 * mrSges(6,1);
t19 = t33 * t58 - t36 * t55;
t15 = t19 * mrSges(6,1);
t14 = -pkin(4) * t26 + t31;
t1 = [t60 * (Ifges(3,4) * t57 + Ifges(3,2) * t60) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t60 + mrSges(3,2) * t57) + t57 * (Ifges(3,1) * t57 + Ifges(3,4) * t60) + 0.2e1 * t50 * (-mrSges(4,1) * t40 + mrSges(4,2) * t41) + t40 * (Ifges(4,4) * t41 + Ifges(4,2) * t40) + t41 * (Ifges(4,1) * t41 + Ifges(4,4) * t40) + Ifges(5,2) * t26 ^ 2 + 0.2e1 * t31 * t68 + Ifges(6,2) * t12 ^ 2 + 0.2e1 * t14 * t69 + Ifges(2,3) + t3 * mrSges(6,3) * t80 + t8 * mrSges(5,3) * t79 + (-0.2e1 * t7 * mrSges(5,3) + Ifges(5,1) * t27 + Ifges(5,4) * t79) * t27 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t13 + Ifges(6,4) * t80) * t13 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t31 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2 + t50 ^ 2) + m(3) * (t71 * pkin(6) ^ 2 + pkin(1) ^ 2) + 0.2e1 * (-t29 * t41 + t30 * t40) * mrSges(4,3) + 0.2e1 * t71 * pkin(6) * mrSges(3,3); t64 + (m(4) * (t29 * t59 + t30 * t56) + (t56 * t40 - t59 * t41) * mrSges(4,3)) * pkin(2) + m(6) * (t19 * t2 + t20 * t3) + m(5) * (t34 * t7 + t36 * t8) + (-mrSges(3,1) * t57 - mrSges(3,2) * t60) * pkin(6) + (t12 * t20 - t13 * t19) * mrSges(6,3) + (t26 * t36 - t27 * t34) * mrSges(5,3) + Ifges(3,6) * t60 + Ifges(3,5) * t57; 0.2e1 * t74 - 0.2e1 * t73 - 0.2e1 * t75 + Ifges(3,3) + 0.2e1 * t15 + 0.2e1 * t65 + m(5) * (t34 ^ 2 + t36 ^ 2) + m(6) * (t19 ^ 2 + t20 ^ 2) + m(4) * (t56 ^ 2 + t59 ^ 2) * pkin(2) ^ 2 + t70; t64 + (m(5) * (t53 * t8 + t54 * t7) + (t26 * t53 - t27 * t54) * mrSges(5,3)) * pkin(3) + m(6) * (t2 * t35 + t3 * t37) + (t12 * t37 - t13 * t35) * mrSges(6,3); t15 + m(6) * (t19 * t35 + t20 * t37) + t32 - t73 + t74 + t65 + (-t20 - t37) * mrSges(6,2) + (m(5) * (t34 * t54 + t36 * t53) + t67) * pkin(3) + t70; -0.2e1 * t72 + 0.2e1 * t32 + m(6) * (t35 ^ 2 + t37 ^ 2) + t70 + (0.2e1 * t67 + m(5) * (t53 ^ 2 + t54 ^ 2) * pkin(3)) * pkin(3); m(5) * t31 + m(6) * t14 + t68 + t69; 0; 0; m(5) + m(6); t66; Ifges(6,3) + t15 - t75; Ifges(6,3) + t32 - t72; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
