% Calculate joint inertia matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:52
% DurationCPUTime: 0.15s
% Computational Cost: add. (77->53), mult. (184->68), div. (0->0), fcn. (102->4), ass. (0->18)
t22 = -2 * mrSges(5,3);
t21 = m(5) * pkin(3);
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t20 = t12 ^ 2 + t14 ^ 2;
t19 = -qJ(4) - pkin(5);
t18 = mrSges(5,1) + t21;
t17 = t20 * mrSges(4,3);
t3 = -t14 * mrSges(5,1) + t12 * mrSges(5,2);
t15 = cos(qJ(2));
t13 = sin(qJ(2));
t11 = t15 ^ 2;
t9 = t13 ^ 2;
t6 = -t14 * pkin(3) - pkin(2);
t5 = t19 * t14;
t4 = -t14 * mrSges(4,1) + t12 * mrSges(4,2);
t2 = t19 * t12;
t1 = [m(2) + m(3) * (t9 + t11) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t20 * t9 + t11); (mrSges(3,1) - t3 - t4) * t15 + (t20 * mrSges(5,3) - mrSges(3,2) + t17) * t13 + m(4) * (t20 * t13 * pkin(5) + t15 * pkin(2)) + m(5) * (-t6 * t15 + (-t12 * t2 - t14 * t5) * t13); -0.2e1 * pkin(2) * t4 + 0.2e1 * t6 * t3 + Ifges(3,3) + m(5) * (t2 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t20 * pkin(5) ^ 2 + pkin(2) ^ 2) + (t5 * t22 + (Ifges(5,2) + Ifges(4,2)) * t14) * t14 + 0.2e1 * pkin(5) * t17 + (t2 * t22 + 0.2e1 * (Ifges(4,4) + Ifges(5,4)) * t14 + (Ifges(5,1) + Ifges(4,1)) * t12) * t12; ((-mrSges(4,2) - mrSges(5,2)) * t14 + (-mrSges(4,1) - t18) * t12) * t13; t5 * mrSges(5,2) + t18 * t2 + (-mrSges(4,2) * pkin(5) + Ifges(4,6) + Ifges(5,6)) * t14 + (-mrSges(4,1) * pkin(5) - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t12; Ifges(4,3) + Ifges(5,3) + (0.2e1 * mrSges(5,1) + t21) * pkin(3); -m(5) * t15; m(5) * t6 + t3; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
