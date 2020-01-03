% Calculate joint inertia matrix for
% S5RPRRR5
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:48
% EndTime: 2020-01-03 11:53:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (458->107), mult. (854->151), div. (0->0), fcn. (728->8), ass. (0->48)
t41 = sin(qJ(5));
t42 = sin(qJ(4));
t44 = cos(qJ(5));
t45 = cos(qJ(4));
t25 = t41 * t45 + t44 * t42;
t68 = t25 ^ 2;
t67 = t45 ^ 2;
t24 = -t41 * t42 + t44 * t45;
t66 = t41 * pkin(4) * t24 * mrSges(6,3) + Ifges(5,5) * t42 + Ifges(5,6) * t45;
t6 = -t24 * mrSges(6,1) + t25 * mrSges(6,2);
t65 = 0.2e1 * t6;
t39 = sin(pkin(9));
t64 = pkin(1) * t39;
t63 = t45 * pkin(4);
t40 = cos(pkin(9));
t32 = t40 * pkin(1) + pkin(2);
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t17 = t46 * t32 - t43 * t64;
t62 = t17 * mrSges(4,1);
t18 = t43 * t32 + t46 * t64;
t61 = t18 * mrSges(4,2);
t60 = Ifges(6,5) * t25 + Ifges(6,6) * t24;
t59 = t42 ^ 2 + t67;
t58 = 0.2e1 * mrSges(6,3);
t57 = t44 * t25 * mrSges(6,3);
t15 = pkin(7) + t18;
t56 = t59 * t15;
t14 = -pkin(3) - t17;
t55 = Ifges(6,1) * t68 + Ifges(5,2) * t67 + Ifges(4,3) + (Ifges(5,1) * t42 + 0.2e1 * Ifges(5,4) * t45) * t42 + (0.2e1 * Ifges(6,4) * t25 + Ifges(6,2) * t24) * t24;
t26 = -t45 * mrSges(5,1) + t42 * mrSges(5,2);
t54 = -t42 * mrSges(5,1) - t45 * mrSges(5,2);
t8 = (-pkin(8) - t15) * t42;
t36 = t45 * pkin(8);
t9 = t45 * t15 + t36;
t2 = -t41 * t9 + t44 * t8;
t3 = t41 * t8 + t44 * t9;
t53 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t60;
t52 = 0.2e1 * t59 * mrSges(5,3);
t28 = (-pkin(8) - pkin(7)) * t42;
t29 = t45 * pkin(7) + t36;
t10 = t44 * t28 - t41 * t29;
t11 = t41 * t28 + t44 * t29;
t51 = t10 * mrSges(6,1) - t11 * mrSges(6,2) + t60;
t50 = (t44 * mrSges(6,1) - t41 * mrSges(6,2)) * pkin(4);
t33 = -pkin(3) - t63;
t12 = t14 - t63;
t1 = [0.2e1 * t62 - 0.2e1 * t61 + t12 * t65 + 0.2e1 * t14 * t26 + Ifges(2,3) + Ifges(3,3) + (-t2 * t25 + t3 * t24) * t58 + t15 * t52 + m(6) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t59 * t15 ^ 2 + t14 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + t55 + (0.2e1 * t40 * mrSges(3,1) - 0.2e1 * t39 * mrSges(3,2) + m(3) * (t39 ^ 2 + t40 ^ 2) * pkin(1)) * pkin(1); m(6) * (t2 * t24 + t3 * t25); m(3) + m(4) + m(5) * t59 + m(6) * (t24 ^ 2 + t68); t62 - t61 + (t12 + t33) * t6 + (t14 - pkin(3)) * t26 + m(6) * (t10 * t2 + t11 * t3 + t33 * t12) + m(5) * (-pkin(3) * t14 + pkin(7) * t56) + ((-t10 - t2) * t25 + (t11 + t3) * t24) * mrSges(6,3) + (t59 * pkin(7) + t56) * mrSges(5,3) + t55; m(6) * (t10 * t24 + t11 * t25); -0.2e1 * pkin(3) * t26 + t33 * t65 + (-t10 * t25 + t11 * t24) * t58 + pkin(7) * t52 + m(6) * (t10 ^ 2 + t11 ^ 2 + t33 ^ 2) + m(5) * (t59 * pkin(7) ^ 2 + pkin(3) ^ 2) + t55; t54 * t15 + (-t57 + m(6) * (t2 * t44 + t3 * t41)) * pkin(4) + t53 + t66; m(6) * (t24 * t44 + t25 * t41) * pkin(4) - t26 - t6; t54 * pkin(7) + (-t57 + m(6) * (t10 * t44 + t11 * t41)) * pkin(4) + t51 + t66; Ifges(5,3) + Ifges(6,3) + m(6) * (t41 ^ 2 + t44 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t50; t53; -t6; t51; Ifges(6,3) + t50; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
