% Calculate joint inertia matrix for
% S5PPRPR3
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:48
% EndTime: 2019-12-05 15:04:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (144->65), mult. (389->102), div. (0->0), fcn. (347->8), ass. (0->32)
t22 = sin(pkin(8));
t21 = sin(pkin(9));
t23 = cos(pkin(9));
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t9 = t21 * t28 + t23 * t26;
t3 = t9 * t22;
t39 = t3 ^ 2;
t7 = t21 * t26 - t23 * t28;
t38 = t7 ^ 2;
t27 = cos(qJ(5));
t19 = t27 ^ 2;
t37 = m(5) * pkin(3);
t36 = t7 * t3;
t25 = sin(qJ(5));
t10 = -t27 * mrSges(6,1) + t25 * mrSges(6,2);
t35 = -mrSges(5,1) + t10;
t34 = t25 ^ 2 + t19;
t33 = t26 ^ 2 + t28 ^ 2;
t32 = t34 * mrSges(6,3);
t24 = cos(pkin(8));
t5 = t7 * t22;
t1 = -t24 * t27 + t25 * t5;
t2 = -t24 * t25 - t27 * t5;
t31 = -t1 * t25 + t2 * t27;
t30 = -mrSges(6,1) * t25 - mrSges(6,2) * t27;
t16 = t24 ^ 2;
t15 = t22 ^ 2;
t14 = -t23 * pkin(3) - pkin(4);
t13 = t21 * pkin(3) + pkin(6);
t6 = t9 ^ 2;
t4 = [m(2) + m(3) * (t15 + t16) + m(4) * (t33 * t15 + t16) + m(5) * (t5 ^ 2 + t16 + t39) + m(6) * (t1 ^ 2 + t2 ^ 2 + t39); m(5) * (-t9 * t5 + t36) + m(6) * (t31 * t9 + t36); m(3) + m(4) * t33 + m(5) * (t6 + t38) + m(6) * (t34 * t6 + t38); t5 * mrSges(5,2) + t35 * t3 + (-t26 * mrSges(4,1) - t28 * mrSges(4,2)) * t22 + t31 * mrSges(6,3) + m(6) * (t31 * t13 + t14 * t3) + (-t21 * t5 - t23 * t3) * t37; t28 * mrSges(4,1) - t26 * mrSges(4,2) + t35 * t7 + (-mrSges(5,2) + t32) * t9 + m(6) * (t34 * t9 * t13 + t14 * t7) + (t21 * t9 - t23 * t7) * t37; Ifges(6,2) * t19 + 0.2e1 * t14 * t10 + Ifges(4,3) + Ifges(5,3) + (Ifges(6,1) * t25 + 0.2e1 * Ifges(6,4) * t27) * t25 + m(6) * (t34 * t13 ^ 2 + t14 ^ 2) + m(5) * (t21 ^ 2 + t23 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (t23 * mrSges(5,1) - t21 * mrSges(5,2)) * pkin(3) + 0.2e1 * t13 * t32; -m(5) * t24 + m(6) * (t27 * t1 + t25 * t2); 0; 0; m(6) * t34 + m(5); t1 * mrSges(6,1) - t2 * mrSges(6,2); t30 * t9; Ifges(6,5) * t25 + Ifges(6,6) * t27 + t30 * t13; -t10; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
