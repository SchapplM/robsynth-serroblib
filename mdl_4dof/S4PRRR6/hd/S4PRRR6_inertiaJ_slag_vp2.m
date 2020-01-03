% Calculate joint inertia matrix for
% S4PRRR6
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:41
% DurationCPUTime: 0.15s
% Computational Cost: add. (143->64), mult. (345->102), div. (0->0), fcn. (283->6), ass. (0->28)
t24 = cos(qJ(3));
t18 = t24 ^ 2;
t34 = -pkin(6) - pkin(5);
t21 = sin(qJ(3));
t33 = t21 ^ 2 + t18;
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t11 = t20 * t24 + t23 * t21;
t22 = sin(qJ(2));
t6 = t11 * t22;
t10 = -t20 * t21 + t23 * t24;
t7 = t10 * t22;
t32 = -t6 * mrSges(5,1) - t7 * mrSges(5,2);
t31 = t33 * mrSges(4,3);
t13 = t34 * t21;
t14 = t34 * t24;
t3 = t23 * t13 + t20 * t14;
t4 = t20 * t13 - t23 * t14;
t30 = t3 * mrSges(5,1) - t4 * mrSges(5,2) + Ifges(5,5) * t11 + Ifges(5,6) * t10;
t29 = -t21 * mrSges(4,1) - t24 * mrSges(4,2);
t28 = (t23 * mrSges(5,1) - t20 * mrSges(5,2)) * pkin(3);
t25 = cos(qJ(2));
t19 = t25 ^ 2;
t17 = t22 ^ 2;
t15 = -t24 * pkin(3) - pkin(2);
t12 = -t24 * mrSges(4,1) + t21 * mrSges(4,2);
t1 = -t10 * mrSges(5,1) + t11 * mrSges(5,2);
t2 = [m(2) + m(3) * (t17 + t19) + m(4) * (t33 * t17 + t19) + m(5) * (t6 ^ 2 + t7 ^ 2 + t19); (t7 * t10 + t6 * t11) * mrSges(5,3) + (mrSges(3,1) - t1 - t12) * t25 + (-mrSges(3,2) + t31) * t22 + m(4) * (t33 * t22 * pkin(5) + t25 * pkin(2)) + m(5) * (-t15 * t25 - t3 * t6 + t4 * t7); Ifges(4,2) * t18 - 0.2e1 * pkin(2) * t12 + 0.2e1 * t15 * t1 + Ifges(3,3) + (Ifges(4,1) * t21 + 0.2e1 * Ifges(4,4) * t24) * t21 + (-0.2e1 * t3 * mrSges(5,3) + Ifges(5,1) * t11) * t11 + 0.2e1 * pkin(5) * t31 + m(5) * (t15 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t33 * pkin(5) ^ 2 + pkin(2) ^ 2) + (0.2e1 * t4 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t11 + Ifges(5,2) * t10) * t10; t29 * t22 + m(5) * (t20 * t7 - t23 * t6) * pkin(3) + t32; Ifges(4,5) * t21 + Ifges(4,6) * t24 + t29 * pkin(5) + (m(5) * (t20 * t4 + t23 * t3) + (t20 * t10 - t23 * t11) * mrSges(5,3)) * pkin(3) + t30; Ifges(4,3) + Ifges(5,3) + m(5) * (t20 ^ 2 + t23 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t28; t32; t30; Ifges(5,3) + t28; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
