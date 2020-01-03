% Calculate joint inertia matrix for
% S4RRRR1
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:09
% DurationCPUTime: 0.19s
% Computational Cost: add. (215->65), mult. (422->97), div. (0->0), fcn. (252->6), ass. (0->45)
t32 = cos(qJ(4));
t28 = t32 ^ 2;
t29 = sin(qJ(4));
t27 = t29 ^ 2;
t49 = mrSges(5,3) * t27;
t22 = pkin(7) * t49;
t48 = mrSges(5,3) * t28;
t23 = pkin(7) * t48;
t15 = -t32 * mrSges(5,1) + t29 * mrSges(5,2);
t52 = pkin(3) * t15;
t55 = t22 + t23 - t52;
t33 = cos(qJ(3));
t51 = t33 * pkin(2);
t20 = -pkin(3) - t51;
t10 = t20 * t15;
t30 = sin(qJ(3));
t19 = t30 * pkin(2) + pkin(7);
t13 = t19 * t49;
t14 = t19 * t48;
t24 = mrSges(4,1) * t51;
t54 = t10 + t13 + t14 + t24;
t31 = sin(qJ(2));
t53 = pkin(1) * t31;
t34 = cos(qJ(2));
t21 = t34 * pkin(1) + pkin(2);
t9 = t30 * t21 + t33 * t53;
t50 = t9 * mrSges(4,2);
t47 = t30 * mrSges(4,2);
t46 = Ifges(5,5) * t29 + Ifges(5,6) * t32;
t45 = t27 + t28;
t44 = pkin(2) * t47;
t43 = Ifges(5,2) * t28 + Ifges(4,3) + (Ifges(5,1) * t29 + 0.2e1 * Ifges(5,4) * t32) * t29;
t42 = t45 * t19;
t41 = Ifges(3,3) + t43;
t40 = -mrSges(5,1) * t29 - mrSges(5,2) * t32;
t8 = t33 * t21 - t30 * t53;
t39 = (t34 * mrSges(3,1) - t31 * mrSges(3,2)) * pkin(1);
t6 = -pkin(3) - t8;
t1 = t6 * t15;
t7 = pkin(7) + t9;
t2 = t7 * t49;
t3 = t7 * t48;
t4 = t8 * mrSges(4,1);
t38 = t1 + t2 + t3 + t4 + t43 - t50;
t5 = [-0.2e1 * t50 + Ifges(2,3) + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 + 0.2e1 * t4 + 0.2e1 * t39 + m(5) * (t45 * t7 ^ 2 + t6 ^ 2) + m(4) * (t8 ^ 2 + t9 ^ 2) + m(3) * (t31 ^ 2 + t34 ^ 2) * pkin(1) ^ 2 + t41; t38 + m(5) * (t20 * t6 + t7 * t42) + (m(4) * (t30 * t9 + t33 * t8) - t47) * pkin(2) + t39 + Ifges(3,3) + t54; -0.2e1 * t44 + 0.2e1 * t10 + 0.2e1 * t13 + 0.2e1 * t14 + 0.2e1 * t24 + m(5) * (t45 * t19 ^ 2 + t20 ^ 2) + m(4) * (t30 ^ 2 + t33 ^ 2) * pkin(2) ^ 2 + t41; m(5) * (t45 * t7 * pkin(7) - pkin(3) * t6) + t38 + t55; -t44 + m(5) * (-pkin(3) * t20 + pkin(7) * t42) + t43 + t54 + t55; 0.2e1 * t22 + 0.2e1 * t23 - 0.2e1 * t52 + m(5) * (pkin(7) ^ 2 * t45 + pkin(3) ^ 2) + t43; t40 * t7 + t46; t19 * t40 + t46; pkin(7) * t40 + t46; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;
