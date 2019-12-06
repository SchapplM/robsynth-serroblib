% Calculate joint inertia matrix for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:50
% EndTime: 2019-12-05 14:57:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (115->50), mult. (326->75), div. (0->0), fcn. (300->8), ass. (0->31)
t24 = cos(qJ(5));
t17 = t24 ^ 2;
t22 = sin(qJ(5));
t30 = t22 ^ 2 + t17;
t37 = m(6) * t30;
t10 = -t24 * mrSges(6,1) + t22 * mrSges(6,2);
t36 = m(6) * pkin(4) + mrSges(5,1) - t10;
t19 = sin(pkin(8));
t18 = sin(pkin(9));
t20 = cos(pkin(9));
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t9 = t25 * t18 + t23 * t20;
t3 = t9 * t19;
t35 = t3 ^ 2;
t7 = t23 * t18 - t25 * t20;
t34 = t7 ^ 2;
t33 = t7 * t3;
t32 = -m(4) - m(5);
t31 = t18 ^ 2 + t20 ^ 2;
t29 = mrSges(6,3) * t30;
t21 = cos(pkin(8));
t5 = t7 * t19;
t1 = -t21 * t24 + t22 * t5;
t2 = -t21 * t22 - t24 * t5;
t28 = -t1 * t22 + t2 * t24;
t27 = -mrSges(6,1) * t22 - mrSges(6,2) * t24;
t15 = t21 ^ 2;
t13 = t19 ^ 2;
t6 = t9 ^ 2;
t4 = [m(2) + m(3) * (t13 + t15) + m(4) * (t31 * t13 + t15) + m(5) * (t5 ^ 2 + t15 + t35) + m(6) * (t1 ^ 2 + t2 ^ 2 + t35); m(5) * (-t9 * t5 + t33) + m(6) * (t28 * t9 + t33); m(3) + m(4) * t31 + m(5) * (t6 + t34) + m(6) * (t30 * t6 + t34); m(6) * (t24 * t1 + t22 * t2) + t32 * t21; 0; -t32 + t37; t5 * mrSges(5,2) + (m(6) * pkin(6) + mrSges(6,3)) * t28 - t36 * t3; (pkin(6) * t37 - mrSges(5,2) + t29) * t9 - t36 * t7; 0; Ifges(5,3) + Ifges(6,2) * t17 - 0.2e1 * pkin(4) * t10 + m(6) * (t30 * pkin(6) ^ 2 + pkin(4) ^ 2) + (Ifges(6,1) * t22 + 0.2e1 * Ifges(6,4) * t24) * t22 + 0.2e1 * pkin(6) * t29; t1 * mrSges(6,1) - t2 * mrSges(6,2); t27 * t9; -t10; Ifges(6,5) * t22 + Ifges(6,6) * t24 + t27 * pkin(6); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
