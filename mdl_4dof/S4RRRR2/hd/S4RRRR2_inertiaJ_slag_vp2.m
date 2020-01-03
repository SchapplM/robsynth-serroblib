% Calculate joint inertia matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:09
% DurationCPUTime: 0.26s
% Computational Cost: add. (300->88), mult. (582->127), div. (0->0), fcn. (474->6), ass. (0->40)
t36 = cos(qJ(3));
t56 = t36 ^ 2;
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t35 = cos(qJ(4));
t17 = -t32 * t33 + t35 * t36;
t55 = t32 * pkin(3) * t17 * mrSges(5,3) + Ifges(4,5) * t33 + Ifges(4,6) * t36;
t18 = t32 * t36 + t35 * t33;
t6 = -t17 * mrSges(5,1) + t18 * mrSges(5,2);
t54 = 0.2e1 * t6;
t37 = cos(qJ(2));
t53 = t37 * pkin(1);
t52 = Ifges(5,5) * t18 + Ifges(5,6) * t17;
t51 = t33 ^ 2 + t56;
t50 = 0.2e1 * mrSges(5,3);
t49 = t35 * t18 * mrSges(5,3);
t26 = -t36 * pkin(3) - pkin(2);
t34 = sin(qJ(2));
t24 = t34 * pkin(1) + pkin(6);
t48 = t51 * t24;
t47 = Ifges(5,1) * t18 ^ 2 + Ifges(4,2) * t56 + Ifges(3,3) + (Ifges(4,1) * t33 + 0.2e1 * Ifges(4,4) * t36) * t33 + (0.2e1 * Ifges(5,4) * t18 + Ifges(5,2) * t17) * t17;
t46 = -t33 * mrSges(4,1) - t36 * mrSges(4,2);
t13 = (-pkin(7) - t24) * t33;
t29 = t36 * pkin(7);
t14 = t36 * t24 + t29;
t4 = t35 * t13 - t32 * t14;
t5 = t32 * t13 + t35 * t14;
t45 = t4 * mrSges(5,1) - t5 * mrSges(5,2) + t52;
t21 = (-pkin(7) - pkin(6)) * t33;
t22 = t36 * pkin(6) + t29;
t8 = t35 * t21 - t32 * t22;
t9 = t32 * t21 + t35 * t22;
t44 = t8 * mrSges(5,1) - t9 * mrSges(5,2) + t52;
t43 = 0.2e1 * t51 * mrSges(4,3);
t42 = (t37 * mrSges(3,1) - t34 * mrSges(3,2)) * pkin(1);
t41 = (t35 * mrSges(5,1) - t32 * mrSges(5,2)) * pkin(3);
t25 = -pkin(2) - t53;
t20 = -t36 * mrSges(4,1) + t33 * mrSges(4,2);
t19 = t26 - t53;
t1 = [t19 * t54 + 0.2e1 * t25 * t20 + Ifges(2,3) + 0.2e1 * t42 + (t5 * t17 - t4 * t18) * t50 + t24 * t43 + m(5) * (t19 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(4) * (t51 * t24 ^ 2 + t25 ^ 2) + m(3) * (t34 ^ 2 + t37 ^ 2) * pkin(1) ^ 2 + t47; (t19 + t26) * t6 + (t25 - pkin(2)) * t20 + t42 + m(5) * (t26 * t19 + t8 * t4 + t9 * t5) + m(4) * (-pkin(2) * t25 + pkin(6) * t48) + ((-t4 - t8) * t18 + (t5 + t9) * t17) * mrSges(5,3) + (t51 * pkin(6) + t48) * mrSges(4,3) + t47; -0.2e1 * pkin(2) * t20 + t26 * t54 + (t9 * t17 - t8 * t18) * t50 + pkin(6) * t43 + m(5) * (t26 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(4) * (t51 * pkin(6) ^ 2 + pkin(2) ^ 2) + t47; t46 * t24 + (-t49 + m(5) * (t32 * t5 + t35 * t4)) * pkin(3) + t45 + t55; t46 * pkin(6) + (-t49 + m(5) * (t32 * t9 + t35 * t8)) * pkin(3) + t44 + t55; Ifges(4,3) + Ifges(5,3) + m(5) * (t32 ^ 2 + t35 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t41; t45; t44; Ifges(5,3) + t41; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
