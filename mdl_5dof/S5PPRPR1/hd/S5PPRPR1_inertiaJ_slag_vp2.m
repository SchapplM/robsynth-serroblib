% Calculate joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:54
% EndTime: 2019-12-05 15:00:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (174->62), mult. (405->99), div. (0->0), fcn. (386->8), ass. (0->29)
t21 = sin(pkin(9));
t23 = cos(pkin(9));
t14 = -t23 * mrSges(5,1) + t21 * mrSges(5,2);
t25 = sin(qJ(5));
t27 = cos(qJ(5));
t11 = t27 * t21 + t25 * t23;
t8 = -t25 * t21 + t27 * t23;
t3 = -t8 * mrSges(6,1) + t11 * mrSges(6,2);
t33 = t14 + t3;
t22 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(3));
t32 = cos(qJ(3));
t9 = t26 * t22 - t32 * t24;
t6 = t9 ^ 2;
t20 = t23 ^ 2;
t31 = pkin(6) + qJ(4);
t30 = t21 ^ 2 + t20;
t29 = qJ(4) * t30;
t17 = -t23 * pkin(4) - pkin(3);
t15 = t31 * t23;
t13 = t31 * t21;
t12 = t32 * t22 + t26 * t24;
t7 = t12 ^ 2;
t5 = -t25 * t13 + t27 * t15;
t4 = -t27 * t13 - t25 * t15;
t2 = t8 * t12;
t1 = t11 * t12;
t10 = [m(2) + m(3) * (t22 ^ 2 + t24 ^ 2) + m(4) * (t7 + t6) + m(5) * (t30 * t7 + t6) + m(6) * (t1 ^ 2 + t2 ^ 2 + t6); m(6) * (-t8 * t1 + t11 * t2); m(3) + m(4) + m(5) * t30 + m(6) * (t11 ^ 2 + t8 ^ 2); (t1 * t11 + t2 * t8) * mrSges(6,3) + (-mrSges(4,1) + t33) * t9 + (t30 * mrSges(5,3) - mrSges(4,2)) * t12 + m(5) * (-pkin(3) * t9 + t12 * t29) + m(6) * (-t4 * t1 + t17 * t9 + t5 * t2); m(6) * (t5 * t11 + t4 * t8); Ifges(5,2) * t20 - 0.2e1 * pkin(3) * t14 + 0.2e1 * t17 * t3 + Ifges(4,3) + (Ifges(5,1) * t21 + 0.2e1 * Ifges(5,4) * t23) * t21 + (0.2e1 * t5 * mrSges(6,3) + Ifges(6,2) * t8) * t8 + 0.2e1 * mrSges(5,3) * t29 + m(6) * (t17 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(5) * (t30 * qJ(4) ^ 2 + pkin(3) ^ 2) + (-0.2e1 * t4 * mrSges(6,3) + Ifges(6,1) * t11 + 0.2e1 * Ifges(6,4) * t8) * t11; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t9; 0; -m(5) * pkin(3) + m(6) * t17 + t33; m(5) + m(6); -t1 * mrSges(6,1) - t2 * mrSges(6,2); -t3; t4 * mrSges(6,1) - t5 * mrSges(6,2) + Ifges(6,5) * t11 + Ifges(6,6) * t8; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
