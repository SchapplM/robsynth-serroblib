% Calculate joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:24
% DurationCPUTime: 0.24s
% Computational Cost: add. (309->81), mult. (604->127), div. (0->0), fcn. (594->6), ass. (0->34)
t28 = sin(pkin(7));
t29 = cos(pkin(7));
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t17 = -t28 * t31 + t29 * t33;
t18 = t28 * t33 + t29 * t31;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t9 = t32 * t17 - t30 * t18;
t45 = 0.2e1 * t9;
t44 = 0.2e1 * t17;
t43 = pkin(2) * t28;
t24 = t29 * pkin(2) + pkin(3);
t14 = t32 * t24 - t30 * t43;
t42 = t14 * mrSges(5,1);
t15 = t30 * t24 + t32 * t43;
t41 = t15 * mrSges(5,2);
t40 = -qJ(3) - pkin(5);
t22 = t40 * t31;
t23 = t40 * t33;
t12 = t28 * t22 - t29 * t23;
t39 = t31 ^ 2 + t33 ^ 2;
t10 = t30 * t17 + t32 * t18;
t38 = -t9 * mrSges(5,1) + t10 * mrSges(5,2);
t25 = -t33 * pkin(2) - pkin(1);
t37 = -t17 * mrSges(4,1) + t18 * mrSges(4,2);
t11 = t29 * t22 + t28 * t23;
t4 = -t18 * pkin(6) + t11;
t5 = t17 * pkin(6) + t12;
t2 = -t30 * t5 + t32 * t4;
t3 = t30 * t4 + t32 * t5;
t36 = t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t10 + Ifges(5,6) * t9;
t13 = -t17 * pkin(3) + t25;
t1 = [Ifges(2,3) + t3 * mrSges(5,3) * t45 + t12 * mrSges(4,3) * t44 + 0.2e1 * t13 * t38 + Ifges(5,2) * t9 ^ 2 + Ifges(4,2) * t17 ^ 2 + t31 * (Ifges(3,1) * t31 + Ifges(3,4) * t33) + t33 * (Ifges(3,4) * t31 + Ifges(3,2) * t33) - 0.2e1 * pkin(1) * (-t33 * mrSges(3,1) + t31 * mrSges(3,2)) + 0.2e1 * t25 * t37 + 0.2e1 * t39 * pkin(5) * mrSges(3,3) + (-0.2e1 * t11 * mrSges(4,3) + Ifges(4,1) * t18 + Ifges(4,4) * t44) * t18 + (-0.2e1 * t2 * mrSges(5,3) + Ifges(5,1) * t10 + Ifges(5,4) * t45) * t10 + m(5) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) + m(3) * (t39 * pkin(5) ^ 2 + pkin(1) ^ 2); m(5) * (t14 * t2 + t15 * t3) + t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,6) * t17 + Ifges(3,5) * t31 + Ifges(3,6) * t33 + Ifges(4,5) * t18 + (-t31 * mrSges(3,1) - t33 * mrSges(3,2)) * pkin(5) + (-t14 * t10 + t15 * t9) * mrSges(5,3) + (m(4) * (t11 * t29 + t12 * t28) + (t28 * t17 - t29 * t18) * mrSges(4,3)) * pkin(2) + t36; 0.2e1 * t42 - 0.2e1 * t41 + Ifges(3,3) + Ifges(4,3) + Ifges(5,3) + m(5) * (t14 ^ 2 + t15 ^ 2) + (0.2e1 * t29 * mrSges(4,1) - 0.2e1 * t28 * mrSges(4,2) + m(4) * (t28 ^ 2 + t29 ^ 2) * pkin(2)) * pkin(2); m(4) * t25 + m(5) * t13 + t37 + t38; 0; m(4) + m(5); t36; Ifges(5,3) - t41 + t42; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
