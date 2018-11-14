% Calculate joint inertia matrix for
% S4RRPR1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:18
% EndTime: 2018-11-14 13:53:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (142->48), mult. (266->62), div. (0->0), fcn. (186->6), ass. (0->27)
t16 = sin(qJ(2));
t29 = pkin(1) * t16;
t13 = sin(pkin(7));
t28 = t13 * pkin(2);
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t18 = cos(qJ(2));
t12 = t18 * pkin(1) + pkin(2);
t14 = cos(pkin(7));
t6 = t14 * t12 - t13 * t29;
t5 = pkin(3) + t6;
t8 = t13 * t12 + t14 * t29;
t3 = t15 * t5 + t17 * t8;
t27 = t3 * mrSges(5,2);
t26 = t6 * mrSges(4,1);
t25 = t8 * mrSges(4,2);
t11 = t14 * pkin(2) + pkin(3);
t9 = t15 * t11 + t17 * t28;
t24 = t9 * mrSges(5,2);
t23 = Ifges(3,3) + Ifges(4,3) + Ifges(5,3);
t22 = t14 * mrSges(4,1) - t13 * mrSges(4,2);
t21 = (t18 * mrSges(3,1) - t16 * mrSges(3,2)) * pkin(1);
t7 = t17 * t11 - t15 * t28;
t4 = t7 * mrSges(5,1);
t2 = -t15 * t8 + t17 * t5;
t1 = t2 * mrSges(5,1);
t10 = [0.2e1 * t26 - 0.2e1 * t25 - 0.2e1 * t27 + Ifges(2,3) + 0.2e1 * t1 + 0.2e1 * t21 + m(5) * (t2 ^ 2 + t3 ^ 2) + m(4) * (t6 ^ 2 + t8 ^ 2) + m(3) * (t16 ^ 2 + t18 ^ 2) * pkin(1) ^ 2 + t23; m(5) * (t7 * t2 + t9 * t3) + t4 + t1 - t25 + t26 + t21 + (-t9 - t3) * mrSges(5,2) + (m(4) * (t13 * t8 + t14 * t6) + t22) * pkin(2) + t23; -0.2e1 * t24 + 0.2e1 * t4 + m(5) * (t7 ^ 2 + t9 ^ 2) + t23 + (0.2e1 * t22 + m(4) * (t13 ^ 2 + t14 ^ 2) * pkin(2)) * pkin(2); 0; 0; m(4) + m(5); Ifges(5,3) + t1 - t27; Ifges(5,3) + t4 - t24; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1) t10(2) t10(4) t10(7); t10(2) t10(3) t10(5) t10(8); t10(4) t10(5) t10(6) t10(9); t10(7) t10(8) t10(9) t10(10);];
Mq  = res;
