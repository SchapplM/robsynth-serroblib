% Calculate joint inertia matrix for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (69->35), mult. (150->49), div. (0->0), fcn. (72->4), ass. (0->16)
t13 = cos(qJ(4));
t24 = t13 ^ 2;
t11 = sin(qJ(4));
t23 = Ifges(5,5) * t11 + Ifges(5,6) * t13;
t22 = t11 ^ 2 + t24;
t21 = Ifges(5,2) * t24 + Ifges(4,3) + (Ifges(5,1) * t11 + 0.2e1 * Ifges(5,4) * t13) * t11;
t12 = sin(qJ(3));
t5 = t12 * pkin(2) + pkin(6);
t20 = t22 * t5;
t19 = 0.2e1 * t22 * mrSges(5,3);
t18 = -mrSges(5,1) * t11 - mrSges(5,2) * t13;
t14 = cos(qJ(3));
t17 = (t14 * mrSges(4,1) - t12 * mrSges(4,2)) * pkin(2);
t6 = -t14 * pkin(2) - pkin(3);
t3 = -t13 * mrSges(5,1) + t11 * mrSges(5,2);
t1 = [m(5) * t22 + m(2) + m(3) + m(4); 0; 0.2e1 * t6 * t3 + Ifges(3,3) + 0.2e1 * t17 + t5 * t19 + m(5) * (t22 * t5 ^ 2 + t6 ^ 2) + m(4) * (t12 ^ 2 + t14 ^ 2) * pkin(2) ^ 2 + t21; 0; m(5) * (-pkin(3) * t6 + pkin(6) * t20) + (-pkin(3) + t6) * t3 + t17 + (t22 * pkin(6) + t20) * mrSges(5,3) + t21; -0.2e1 * pkin(3) * t3 + m(5) * (t22 * pkin(6) ^ 2 + pkin(3) ^ 2) + pkin(6) * t19 + t21; -t3; t18 * t5 + t23; t18 * pkin(6) + t23; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
