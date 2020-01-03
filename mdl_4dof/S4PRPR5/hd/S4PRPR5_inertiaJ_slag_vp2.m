% Calculate joint inertia matrix for
% S4PRPR5
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:22:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (70->39), mult. (166->61), div. (0->0), fcn. (119->6), ass. (0->17)
t12 = sin(pkin(7));
t13 = cos(pkin(7));
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t2 = t12 * t15 - t13 * t17;
t21 = t2 ^ 2;
t16 = cos(qJ(4));
t11 = t16 ^ 2;
t14 = sin(qJ(4));
t20 = t14 ^ 2 + t11;
t19 = t20 * mrSges(5,3);
t9 = -t13 * pkin(2) - pkin(3);
t8 = t12 * pkin(2) + pkin(5);
t5 = -t16 * mrSges(5,1) + t14 * mrSges(5,2);
t4 = t12 * t17 + t13 * t15;
t1 = t4 ^ 2;
t3 = [m(2) + m(3) * (t15 ^ 2 + t17 ^ 2) + m(4) * (t1 + t21) + m(5) * (t20 * t1 + t21); t17 * mrSges(3,1) - t15 * mrSges(3,2) + (-mrSges(4,1) + t5) * t2 + (-mrSges(4,2) + t19) * t4 + m(5) * (t20 * t8 * t4 + t9 * t2) + m(4) * (t12 * t4 - t13 * t2) * pkin(2); Ifges(5,2) * t11 + 0.2e1 * t9 * t5 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t14 + 0.2e1 * Ifges(5,4) * t16) * t14 + m(5) * (t20 * t8 ^ 2 + t9 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t13 * mrSges(4,1) - t12 * mrSges(4,2)) * pkin(2) + 0.2e1 * t8 * t19; 0; 0; m(5) * t20 + m(4); (-mrSges(5,1) * t14 - mrSges(5,2) * t16) * t4; (-mrSges(5,2) * t8 + Ifges(5,6)) * t16 + (-mrSges(5,1) * t8 + Ifges(5,5)) * t14; -t5; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
