% Calculate joint inertia matrix for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:38
% EndTime: 2019-12-05 15:12:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (181->58), mult. (401->84), div. (0->0), fcn. (361->8), ass. (0->29)
t27 = cos(qJ(5));
t21 = t27 ^ 2;
t24 = sin(qJ(5));
t39 = t24 ^ 2 + t21;
t43 = mrSges(6,3) * t39;
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t10 = -t26 * t22 + t29 * t23;
t11 = t29 * t22 + t26 * t23;
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t6 = -t28 * t10 + t25 * t11;
t42 = t6 ^ 2;
t40 = Ifges(6,5) * t24 + Ifges(6,6) * t27;
t38 = Ifges(6,2) * t21 + Ifges(5,3) + (Ifges(6,1) * t24 + 0.2e1 * Ifges(6,4) * t27) * t24;
t37 = t39 * pkin(7);
t16 = t25 * pkin(3) + pkin(7);
t36 = t39 * t16;
t35 = -mrSges(6,1) * t24 - mrSges(6,2) * t27;
t34 = 0.2e1 * t43;
t14 = -t27 * mrSges(6,1) + t24 * mrSges(6,2);
t8 = t25 * t10 + t28 * t11;
t33 = (-mrSges(5,1) + t14) * t6 + (-mrSges(5,2) + t43) * t8;
t32 = (t28 * mrSges(5,1) - t25 * mrSges(5,2)) * pkin(3);
t17 = -t28 * pkin(3) - pkin(4);
t5 = t8 ^ 2;
t1 = [m(2) + m(6) * (t39 * t5 + t42) + m(5) * (t5 + t42) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t22 ^ 2 + t23 ^ 2); 0; m(6) * t39 + m(3) + m(4) + m(5); t10 * mrSges(4,1) - t11 * mrSges(4,2) + m(6) * (t17 * t6 + t8 * t36) + m(5) * (t25 * t8 - t28 * t6) * pkin(3) + t33; 0; 0.2e1 * t17 * t14 + Ifges(4,3) + 0.2e1 * t32 + t16 * t34 + m(6) * (t39 * t16 ^ 2 + t17 ^ 2) + m(5) * (t25 ^ 2 + t28 ^ 2) * pkin(3) ^ 2 + t38; m(6) * (-pkin(4) * t6 + t8 * t37) + t33; 0; m(6) * (-pkin(4) * t17 + pkin(7) * t36) + (t17 - pkin(4)) * t14 + t32 + (t36 + t37) * mrSges(6,3) + t38; -0.2e1 * pkin(4) * t14 + m(6) * (t39 * pkin(7) ^ 2 + pkin(4) ^ 2) + pkin(7) * t34 + t38; t35 * t8; -t14; t35 * t16 + t40; t35 * pkin(7) + t40; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
