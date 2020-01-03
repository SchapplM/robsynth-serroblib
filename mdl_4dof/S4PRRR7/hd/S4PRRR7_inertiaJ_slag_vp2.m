% Calculate joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:04
% EndTime: 2019-12-31 16:36:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (183->95), mult. (455->147), div. (0->0), fcn. (383->8), ass. (0->51)
t56 = 2 * pkin(6);
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t14 = -t34 * mrSges(5,1) + t31 * mrSges(5,2);
t55 = -m(5) * pkin(3) - mrSges(4,1) + t14;
t30 = cos(pkin(4));
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t29 = sin(pkin(4));
t33 = sin(qJ(2));
t47 = t29 * t33;
t7 = -t30 * t35 + t32 * t47;
t54 = t7 ^ 2;
t53 = t34 / 0.2e1;
t52 = pkin(6) * t35;
t51 = t7 * t32;
t50 = Ifges(5,4) * t31;
t49 = Ifges(5,4) * t34;
t48 = Ifges(5,6) * t35;
t36 = cos(qJ(2));
t46 = t29 * t36;
t45 = t31 * t32;
t44 = t32 * t34;
t43 = t31 ^ 2 + t34 ^ 2;
t13 = -t35 * pkin(3) - t32 * pkin(7) - pkin(2);
t3 = t34 * t13 - t31 * t52;
t4 = t31 * t13 + t34 * t52;
t41 = -t3 * t31 + t4 * t34;
t9 = t30 * t32 + t35 * t47;
t40 = t9 * t35 + t51;
t39 = mrSges(5,1) * t31 + mrSges(5,2) * t34;
t38 = pkin(6) ^ 2;
t28 = t35 ^ 2;
t26 = t32 ^ 2;
t24 = t29 ^ 2;
t23 = t26 * t38;
t22 = Ifges(5,5) * t31;
t21 = Ifges(5,6) * t34;
t19 = t24 * t36 ^ 2;
t18 = Ifges(5,5) * t44;
t17 = Ifges(5,1) * t31 + t49;
t16 = Ifges(5,2) * t34 + t50;
t15 = -t35 * mrSges(4,1) + t32 * mrSges(4,2);
t12 = -t35 * mrSges(5,1) - mrSges(5,3) * t44;
t11 = t35 * mrSges(5,2) - mrSges(5,3) * t45;
t10 = t39 * t32;
t6 = -Ifges(5,5) * t35 + (Ifges(5,1) * t34 - t50) * t32;
t5 = -t48 + (-Ifges(5,2) * t31 + t49) * t32;
t2 = -t31 * t46 + t9 * t34;
t1 = -t9 * t31 - t34 * t46;
t8 = [m(2) + m(3) * (t24 * t33 ^ 2 + t30 ^ 2 + t19) + m(4) * (t9 ^ 2 + t19 + t54) + m(5) * (t1 ^ 2 + t2 ^ 2 + t54); t1 * t12 + t7 * t10 + t2 * t11 + t40 * mrSges(4,3) + (-t33 * mrSges(3,2) + (mrSges(3,1) - t15) * t36) * t29 + m(4) * (pkin(2) * t46 + pkin(6) * t40) + m(5) * (pkin(6) * t51 + t3 * t1 + t4 * t2); -0.2e1 * pkin(2) * t15 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t12 + Ifges(3,3) + (t26 + t28) * mrSges(4,3) * t56 + m(5) * (t3 ^ 2 + t4 ^ 2 + t23) + m(4) * (pkin(2) ^ 2 + t28 * t38 + t23) + (-t18 + (Ifges(5,3) + Ifges(4,2)) * t35) * t35 + (Ifges(4,1) * t32 + 0.2e1 * Ifges(4,4) * t35 + t10 * t56 + t34 * t6 + (-t5 + t48) * t31) * t32; -t9 * mrSges(4,2) + (m(5) * pkin(7) + mrSges(5,3)) * (-t1 * t31 + t2 * t34) + t55 * t7; -pkin(3) * t10 + t31 * t6 / 0.2e1 + t5 * t53 + (-pkin(6) * mrSges(4,2) - t22 / 0.2e1 - t21 / 0.2e1 + Ifges(4,6)) * t35 + t41 * mrSges(5,3) + (m(5) * t41 + t34 * t11 - t31 * t12) * pkin(7) + (t17 * t53 - t31 * t16 / 0.2e1 + Ifges(4,5) + t55 * pkin(6)) * t32; Ifges(4,3) + t31 * t17 + t34 * t16 + m(5) * (t43 * pkin(7) ^ 2 + pkin(3) ^ 2) - 0.2e1 * pkin(3) * t14 + 0.2e1 * t43 * pkin(7) * mrSges(5,3); t1 * mrSges(5,1) - t2 * mrSges(5,2); t3 * mrSges(5,1) - t4 * mrSges(5,2) - Ifges(5,6) * t45 - Ifges(5,3) * t35 + t18; -pkin(7) * t39 + t21 + t22; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
