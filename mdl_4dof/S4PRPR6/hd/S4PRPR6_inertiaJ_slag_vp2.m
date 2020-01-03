% Calculate joint inertia matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (110->51), mult. (255->79), div. (0->0), fcn. (212->6), ass. (0->26)
t18 = sin(pkin(7));
t19 = cos(pkin(7));
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t7 = -t20 * t18 + t22 * t19;
t8 = t22 * t18 + t20 * t19;
t1 = -t7 * mrSges(5,1) + t8 * mrSges(5,2);
t10 = -t19 * mrSges(4,1) + t18 * mrSges(4,2);
t29 = -t1 - t10;
t15 = t19 ^ 2;
t28 = m(4) + m(5);
t27 = pkin(5) + qJ(3);
t26 = t18 ^ 2 + t15;
t25 = qJ(3) * t26;
t23 = cos(qJ(2));
t21 = sin(qJ(2));
t17 = t23 ^ 2;
t16 = t21 ^ 2;
t12 = -t19 * pkin(3) - pkin(2);
t11 = t27 * t19;
t9 = t27 * t18;
t5 = t7 * t21;
t4 = t8 * t21;
t3 = t22 * t11 - t20 * t9;
t2 = -t20 * t11 - t22 * t9;
t6 = [m(2) + m(3) * (t16 + t17) + m(4) * (t26 * t16 + t17) + m(5) * (t4 ^ 2 + t5 ^ 2 + t17); (t4 * t8 + t5 * t7) * mrSges(5,3) + (mrSges(3,1) + t29) * t23 + (t26 * mrSges(4,3) - mrSges(3,2)) * t21 + m(4) * (t23 * pkin(2) + t21 * t25) + m(5) * (-t12 * t23 - t2 * t4 + t3 * t5); Ifges(4,2) * t15 - 0.2e1 * pkin(2) * t10 + 0.2e1 * t12 * t1 + Ifges(3,3) + (Ifges(4,1) * t18 + 0.2e1 * Ifges(4,4) * t19) * t18 + (-0.2e1 * t2 * mrSges(5,3) + Ifges(5,1) * t8) * t8 + 0.2e1 * mrSges(4,3) * t25 + m(5) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t26 * qJ(3) ^ 2 + pkin(2) ^ 2) + (0.2e1 * t3 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t8 + Ifges(5,2) * t7) * t7; -t28 * t23; -m(4) * pkin(2) + m(5) * t12 - t29; t28; -t4 * mrSges(5,1) - t5 * mrSges(5,2); t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t8 + Ifges(5,6) * t7; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
