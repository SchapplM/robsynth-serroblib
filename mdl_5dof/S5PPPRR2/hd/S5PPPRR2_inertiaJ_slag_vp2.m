% Calculate joint inertia matrix for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:27
% DurationCPUTime: 0.20s
% Computational Cost: add. (126->64), mult. (386->101), div. (0->0), fcn. (336->8), ass. (0->37)
t41 = m(6) * pkin(6) + mrSges(6,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t8 = -t25 * mrSges(6,1) + t23 * mrSges(6,2);
t40 = m(6) * pkin(4) + mrSges(5,1) - t8;
t22 = cos(pkin(8));
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t20 = sin(pkin(8));
t21 = cos(pkin(9));
t33 = t20 * t21;
t3 = t22 * t26 + t24 * t33;
t39 = t3 ^ 2;
t17 = t25 ^ 2;
t38 = t24 * t3;
t37 = t26 * t3;
t19 = sin(pkin(9));
t36 = t19 * t20;
t34 = t19 * t26;
t32 = t23 ^ 2 + t17;
t5 = -t22 * t24 + t26 * t33;
t1 = -t5 * t23 + t25 * t36;
t2 = t23 * t36 + t5 * t25;
t30 = -t1 * t23 + t2 * t25;
t6 = -t25 * t21 - t23 * t34;
t7 = -t23 * t21 + t25 * t34;
t29 = -t6 * t23 + t7 * t25;
t28 = -mrSges(6,1) * t23 - mrSges(6,2) * t25;
t18 = t26 ^ 2;
t16 = t24 ^ 2;
t14 = t22 ^ 2;
t13 = t21 ^ 2;
t12 = t20 ^ 2;
t11 = t19 ^ 2;
t10 = t16 * t11;
t9 = t11 * t12;
t4 = [m(2) + m(3) * (t12 + t14) + m(4) * (t13 * t12 + t14 + t9) + m(5) * (t5 ^ 2 + t39 + t9) + m(6) * (t1 ^ 2 + t2 ^ 2 + t39); m(5) * (t26 * t5 - t33 + t38) * t19 + (t6 * t1 + t38 * t19 + t7 * t2) * m(6); m(3) + m(4) * (t11 + t13) + m(5) * (t18 * t11 + t10 + t13) + m(6) * (t6 ^ 2 + t7 ^ 2 + t10); -m(4) * t22 + m(5) * (t24 * t5 - t37) + m(6) * (t30 * t24 - t37); m(6) * (t29 - t34) * t24; m(4) + m(5) * (t16 + t18) + m(6) * (t32 * t16 + t18); -t5 * mrSges(5,2) - t40 * t3 + t41 * t30; -t40 * t19 * t24 - mrSges(5,2) * t34 + t41 * t29; t40 * t26 + (t32 * t41 - mrSges(5,2)) * t24; Ifges(5,3) + Ifges(6,2) * t17 - 0.2e1 * pkin(4) * t8 + m(6) * (t32 * pkin(6) ^ 2 + pkin(4) ^ 2) + (Ifges(6,1) * t23 + 0.2e1 * Ifges(6,4) * t25) * t23 + 0.2e1 * t32 * pkin(6) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2); t6 * mrSges(6,1) - t7 * mrSges(6,2); t28 * t24; Ifges(6,5) * t23 + Ifges(6,6) * t25 + t28 * pkin(6); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
