% Calculate joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:16
% EndTime: 2020-01-03 11:38:18
% DurationCPUTime: 0.38s
% Computational Cost: add. (485->101), mult. (892->155), div. (0->0), fcn. (882->8), ass. (0->42)
t39 = cos(qJ(3));
t54 = t39 ^ 2;
t53 = m(5) * pkin(3);
t32 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(3));
t23 = -t32 * t37 + t34 * t39;
t52 = t23 ^ 2;
t51 = 0.2e1 * t23;
t50 = pkin(3) * t32;
t28 = pkin(3) * t34 + pkin(4);
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t18 = t28 * t38 - t36 * t50;
t49 = t18 * mrSges(6,1);
t19 = t28 * t36 + t38 * t50;
t48 = t19 * mrSges(6,2);
t33 = sin(pkin(8));
t27 = pkin(1) * t33 + pkin(6);
t46 = qJ(4) + t27;
t21 = t46 * t37;
t22 = t46 * t39;
t8 = -t32 * t21 + t34 * t22;
t47 = t37 ^ 2 + t54;
t35 = cos(pkin(8));
t29 = -pkin(1) * t35 - pkin(2);
t24 = t32 * t39 + t34 * t37;
t12 = t23 * t38 - t24 * t36;
t13 = t23 * t36 + t24 * t38;
t4 = -t12 * mrSges(6,1) + t13 * mrSges(6,2);
t45 = -t23 * mrSges(5,1) + t24 * mrSges(5,2);
t7 = -t34 * t21 - t22 * t32;
t44 = -t39 * mrSges(4,1) + t37 * mrSges(4,2);
t5 = -pkin(7) * t24 + t7;
t6 = pkin(7) * t23 + t8;
t2 = -t36 * t6 + t38 * t5;
t3 = t36 * t5 + t38 * t6;
t43 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t12;
t25 = -pkin(3) * t39 + t29;
t42 = -t4 - t45;
t14 = -pkin(4) * t23 + t25;
t1 = [Ifges(2,3) + Ifges(3,3) + t8 * mrSges(5,3) * t51 + 0.2e1 * t14 * t4 + 0.2e1 * t29 * t44 + 0.2e1 * t25 * t45 + Ifges(5,2) * t52 + Ifges(4,2) * t54 + (-0.2e1 * mrSges(6,3) * t2 + Ifges(6,1) * t13) * t13 + (-0.2e1 * mrSges(5,3) * t7 + Ifges(5,1) * t24 + Ifges(5,4) * t51) * t24 + (0.2e1 * mrSges(6,3) * t3 + 0.2e1 * Ifges(6,4) * t13 + Ifges(6,2) * t12) * t12 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t25 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t27 ^ 2 * t47 + t29 ^ 2) + m(3) * (t33 ^ 2 + t35 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t37 + 0.2e1 * Ifges(4,4) * t39) * t37 + 0.2e1 * (t35 * mrSges(3,1) - t33 * mrSges(3,2)) * pkin(1) + 0.2e1 * t47 * t27 * mrSges(4,3); m(5) * (t23 * t7 + t24 * t8) + m(6) * (t12 * t2 + t13 * t3); m(3) + m(4) * t47 + m(5) * (t24 ^ 2 + t52) + m(6) * (t12 ^ 2 + t13 ^ 2); m(6) * (t18 * t2 + t19 * t3) - t8 * mrSges(5,2) + t7 * mrSges(5,1) + Ifges(5,5) * t24 + Ifges(5,6) * t23 + Ifges(4,5) * t37 + Ifges(4,6) * t39 + (-mrSges(4,1) * t37 - mrSges(4,2) * t39) * t27 + (t12 * t19 - t13 * t18) * mrSges(6,3) + (m(5) * (t32 * t8 + t34 * t7) + (t23 * t32 - t24 * t34) * mrSges(5,3)) * pkin(3) + t43; m(6) * (t12 * t18 + t13 * t19) + (t23 * t34 + t24 * t32) * t53 + t42 - t44; 0.2e1 * t49 - 0.2e1 * t48 + Ifges(4,3) + Ifges(5,3) + Ifges(6,3) + m(6) * (t18 ^ 2 + t19 ^ 2) + (0.2e1 * t34 * mrSges(5,1) - 0.2e1 * t32 * mrSges(5,2) + (t32 ^ 2 + t34 ^ 2) * t53) * pkin(3); m(5) * t25 + m(6) * t14 - t42; 0; 0; m(5) + m(6); t43; -t4; Ifges(6,3) - t48 + t49; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
