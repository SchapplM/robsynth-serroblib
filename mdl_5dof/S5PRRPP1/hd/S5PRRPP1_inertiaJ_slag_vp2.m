% Calculate joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:12
% EndTime: 2019-12-05 16:06:13
% DurationCPUTime: 0.30s
% Computational Cost: add. (231->84), mult. (480->118), div. (0->0), fcn. (397->4), ass. (0->27)
t24 = cos(qJ(3));
t39 = t24 ^ 2;
t38 = m(5) * pkin(3);
t37 = m(5) + m(6);
t18 = -t24 * pkin(3) - pkin(2);
t36 = 0.2e1 * t18;
t35 = mrSges(6,2) + mrSges(5,3);
t34 = -qJ(4) - pkin(6);
t23 = sin(qJ(3));
t33 = t23 ^ 2 + t39;
t14 = t34 * t24;
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t29 = t34 * t23;
t4 = -t21 * t14 - t22 * t29;
t6 = -t22 * t14 + t21 * t29;
t32 = t4 ^ 2 + t6 ^ 2;
t28 = -t24 * mrSges(4,1) + t23 * mrSges(4,2);
t10 = t21 * t23 - t22 * t24;
t12 = t21 * t24 + t22 * t23;
t7 = t10 * mrSges(6,1);
t8 = t12 * mrSges(5,2);
t27 = -t10 * mrSges(5,1) + t12 * mrSges(6,3) - t7 - t8;
t17 = -t22 * pkin(3) - pkin(4);
t15 = t21 * pkin(3) + qJ(5);
t2 = t10 * pkin(4) - t12 * qJ(5) + t18;
t1 = [m(4) * t33 + m(2) + m(3) + t37 * (t10 ^ 2 + t12 ^ 2); t37 * (t4 * t10 + t6 * t12); Ifges(3,3) + Ifges(4,2) * t39 + 0.2e1 * t2 * t7 - 0.2e1 * pkin(2) * t28 + t8 * t36 + 0.2e1 * t33 * pkin(6) * mrSges(4,3) + m(6) * (t2 ^ 2 + t32) + m(5) * (t18 ^ 2 + t32) + m(4) * (t33 * pkin(6) ^ 2 + pkin(2) ^ 2) + (-0.2e1 * t2 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t12 + 0.2e1 * t35 * t4) * t12 + (mrSges(5,1) * t36 + (Ifges(6,3) + Ifges(5,2)) * t10 - 0.2e1 * t35 * t6 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t12) * t10 + (Ifges(4,1) * t23 + 0.2e1 * Ifges(4,4) * t24) * t23; m(6) * (t17 * t10 + t15 * t12) + (-t10 * t22 + t12 * t21) * t38 + t27 - t28; Ifges(4,5) * t23 + Ifges(4,6) * t24 - t4 * mrSges(6,1) + t6 * mrSges(6,3) + m(6) * (t15 * t6 + t17 * t4) - t6 * mrSges(5,2) - t4 * mrSges(5,1) + (-t23 * mrSges(4,1) - t24 * mrSges(4,2)) * pkin(6) + (t17 * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t12 + (-t15 * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t10 + (m(5) * (t21 * t6 - t22 * t4) + (-t21 * t10 - t22 * t12) * mrSges(5,3)) * pkin(3); -0.2e1 * t17 * mrSges(6,1) + 0.2e1 * t15 * mrSges(6,3) + Ifges(6,2) + Ifges(4,3) + Ifges(5,3) + m(6) * (t15 ^ 2 + t17 ^ 2) + (0.2e1 * t22 * mrSges(5,1) - 0.2e1 * t21 * mrSges(5,2) + (t21 ^ 2 + t22 ^ 2) * t38) * pkin(3); 0; m(5) * t18 + m(6) * t2 - t27; 0; t37; m(6) * t10; m(6) * t4 + t12 * mrSges(6,2); m(6) * t17 - mrSges(6,1); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
