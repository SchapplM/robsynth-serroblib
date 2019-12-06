% Calculate joint inertia matrix for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:48
% EndTime: 2019-12-05 16:50:50
% DurationCPUTime: 0.36s
% Computational Cost: add. (321->114), mult. (712->155), div. (0->0), fcn. (600->6), ass. (0->39)
t58 = -mrSges(5,1) - mrSges(6,1);
t57 = mrSges(5,2) - mrSges(6,3);
t54 = mrSges(6,2) + mrSges(5,3);
t38 = cos(qJ(3));
t32 = t38 ^ 2;
t56 = 2 * mrSges(6,3);
t55 = -pkin(7) - pkin(6);
t53 = Ifges(6,2) + Ifges(5,3);
t35 = sin(qJ(3));
t52 = t35 ^ 2 + t32;
t25 = t55 * t38;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t50 = t55 * t35;
t10 = -t37 * t25 + t34 * t50;
t8 = -t34 * t25 - t37 * t50;
t51 = t10 ^ 2 + t8 ^ 2;
t29 = -t38 * pkin(3) - pkin(2);
t22 = t34 * t38 + t37 * t35;
t36 = sin(qJ(2));
t14 = t22 * t36;
t21 = t34 * t35 - t37 * t38;
t16 = t21 * t36;
t49 = -t10 * t16 + t8 * t14;
t48 = t52 * mrSges(4,3);
t46 = -t35 * mrSges(4,1) - t38 * mrSges(4,2);
t44 = (t37 * mrSges(5,1) - t34 * mrSges(5,2)) * pkin(3);
t43 = t58 * t14 + t57 * t16;
t42 = t58 * t8 + (Ifges(6,4) + Ifges(5,5)) * t22 + (-Ifges(5,6) + Ifges(6,6)) * t21 - t57 * t10;
t39 = cos(qJ(2));
t33 = t39 ^ 2;
t31 = t36 ^ 2;
t28 = -t37 * pkin(3) - pkin(4);
t26 = t34 * pkin(3) + qJ(5);
t24 = -t38 * mrSges(4,1) + t35 * mrSges(4,2);
t4 = t21 * mrSges(5,1) + t22 * mrSges(5,2);
t3 = t21 * mrSges(6,1) - t22 * mrSges(6,3);
t2 = t21 * pkin(4) - t22 * qJ(5) + t29;
t1 = [m(2) + m(3) * (t31 + t33) + m(4) * (t52 * t31 + t33) + (m(5) + m(6)) * (t14 ^ 2 + t16 ^ 2 + t33); (-mrSges(3,2) + t48) * t36 + (mrSges(3,1) - t24 - t3 - t4) * t39 + m(4) * (t52 * t36 * pkin(6) + t39 * pkin(2)) + m(5) * (-t29 * t39 + t49) + m(6) * (-t2 * t39 + t49) + t54 * (t14 * t22 + t16 * t21); Ifges(4,2) * t32 - 0.2e1 * pkin(2) * t24 + 0.2e1 * t2 * t3 + 0.2e1 * t29 * t4 + Ifges(3,3) + (Ifges(4,1) * t35 + 0.2e1 * Ifges(4,4) * t38) * t35 + 0.2e1 * pkin(6) * t48 + m(6) * (t2 ^ 2 + t51) + m(5) * (t29 ^ 2 + t51) + m(4) * (t52 * pkin(6) ^ 2 + pkin(2) ^ 2) + ((Ifges(6,1) + Ifges(5,1)) * t22 + 0.2e1 * t54 * t8) * t22 + ((Ifges(6,3) + Ifges(5,2)) * t21 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t22 - 0.2e1 * t54 * t10) * t21; t46 * t36 + m(6) * (t28 * t14 - t26 * t16) + m(5) * (-t14 * t37 - t16 * t34) * pkin(3) + t43; m(6) * (t26 * t10 + t28 * t8) + Ifges(4,5) * t35 + Ifges(4,6) * t38 + t46 * pkin(6) + (-t26 * t21 + t28 * t22) * mrSges(6,2) + (m(5) * (t10 * t34 - t37 * t8) + (-t34 * t21 - t37 * t22) * mrSges(5,3)) * pkin(3) + t42; -0.2e1 * t28 * mrSges(6,1) + t26 * t56 + Ifges(4,3) + 0.2e1 * t44 + m(6) * (t26 ^ 2 + t28 ^ 2) + m(5) * (t34 ^ 2 + t37 ^ 2) * pkin(3) ^ 2 + t53; m(6) * (-pkin(4) * t14 - t16 * qJ(5)) + t43; m(6) * (-pkin(4) * t8 + qJ(5) * t10) + (-pkin(4) * t22 - qJ(5) * t21) * mrSges(6,2) + t42; m(6) * (-pkin(4) * t28 + qJ(5) * t26) + t44 + (t26 + qJ(5)) * mrSges(6,3) + (-t28 + pkin(4)) * mrSges(6,1) + t53; 0.2e1 * pkin(4) * mrSges(6,1) + qJ(5) * t56 + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t53; m(6) * t14; m(6) * t8 + t22 * mrSges(6,2); m(6) * t28 - mrSges(6,1); -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
