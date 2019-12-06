% Calculate joint inertia matrix for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:19
% EndTime: 2019-12-05 15:26:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (120->52), mult. (253->72), div. (0->0), fcn. (174->6), ass. (0->26)
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t3 = t17 * t20 - t18 * t22;
t1 = t3 ^ 2;
t5 = t17 * t22 + t18 * t20;
t2 = t5 ^ 2;
t11 = t17 * pkin(2) + qJ(4);
t32 = t11 ^ 2;
t21 = cos(qJ(5));
t16 = t21 ^ 2;
t19 = sin(qJ(5));
t28 = t19 ^ 2 + t16;
t7 = m(6) * t28;
t31 = m(5) + t7;
t30 = t11 * t5;
t8 = t19 * mrSges(6,1) + t21 * mrSges(6,2);
t29 = t8 + mrSges(5,3);
t14 = -t18 * pkin(2) - pkin(3);
t27 = t28 * mrSges(6,3);
t10 = -pkin(6) + t14;
t26 = t28 * t10;
t25 = mrSges(5,2) - t27;
t24 = t21 * mrSges(6,1) - t19 * mrSges(6,2);
t4 = [m(2) + m(3) * (t20 ^ 2 + t22 ^ 2) + m(6) * (t28 * t1 + t2) + (m(4) + m(5)) * (t1 + t2); t22 * mrSges(3,1) - t20 * mrSges(3,2) + (-mrSges(4,2) + t29) * t5 + (-mrSges(4,1) + t25) * t3 + m(5) * (t14 * t3 + t30) + m(6) * (t3 * t26 + t30) + m(4) * (t17 * t5 - t18 * t3) * pkin(2); Ifges(6,1) * t16 + 0.2e1 * t14 * mrSges(5,2) + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + (-0.2e1 * Ifges(6,4) * t21 + Ifges(6,2) * t19) * t19 + m(6) * (t28 * t10 ^ 2 + t32) + m(5) * (t14 ^ 2 + t32) + m(4) * (t17 ^ 2 + t18 ^ 2) * pkin(2) ^ 2 + 0.2e1 * t29 * t11 + 0.2e1 * (t18 * mrSges(4,1) - t17 * mrSges(4,2)) * pkin(2) - 0.2e1 * t10 * t27; 0; 0; m(4) + t31; 0.2e1 * (m(5) / 0.2e1 + t7 / 0.2e1) * t3; m(5) * t14 + m(6) * t26 + t25; 0; t31; t24 * t3; Ifges(6,5) * t21 - Ifges(6,6) * t19 + t24 * t10; -t8; t24; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
