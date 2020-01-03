% Calculate joint inertia matrix for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (272->66), mult. (527->103), div. (0->0), fcn. (534->6), ass. (0->31)
t26 = sin(pkin(7));
t27 = cos(pkin(7));
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t18 = -t29 * t26 + t31 * t27;
t19 = t31 * t26 + t29 * t27;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t9 = t30 * t18 - t28 * t19;
t42 = 0.2e1 * t9;
t25 = t27 ^ 2;
t41 = 0.2e1 * t18;
t40 = pkin(5) + qJ(2);
t20 = t40 * t26;
t21 = t40 * t27;
t12 = -t29 * t20 + t31 * t21;
t39 = t26 ^ 2 + t25;
t10 = t28 * t18 + t30 * t19;
t38 = -t9 * mrSges(5,1) + t10 * mrSges(5,2);
t22 = -t27 * pkin(2) - pkin(1);
t37 = -t27 * mrSges(3,1) + t26 * mrSges(3,2);
t36 = -t18 * mrSges(4,1) + t19 * mrSges(4,2);
t11 = -t31 * t20 - t29 * t21;
t4 = -t19 * pkin(6) + t11;
t5 = t18 * pkin(6) + t12;
t2 = -t28 * t5 + t30 * t4;
t3 = t28 * t4 + t30 * t5;
t35 = t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t10 + Ifges(5,6) * t9;
t34 = (t30 * mrSges(5,1) - t28 * mrSges(5,2)) * pkin(3);
t13 = -t18 * pkin(3) + t22;
t1 = [Ifges(2,3) + t3 * mrSges(5,3) * t42 + t12 * mrSges(4,3) * t41 + 0.2e1 * t13 * t38 + Ifges(5,2) * t9 ^ 2 + Ifges(4,2) * t18 ^ 2 + 0.2e1 * t22 * t36 + Ifges(3,2) * t25 - 0.2e1 * pkin(1) * t37 + (Ifges(3,1) * t26 + 0.2e1 * Ifges(3,4) * t27) * t26 + 0.2e1 * t39 * qJ(2) * mrSges(3,3) + (-0.2e1 * t11 * mrSges(4,3) + Ifges(4,1) * t19 + Ifges(4,4) * t41) * t19 + (-0.2e1 * t2 * mrSges(5,3) + Ifges(5,1) * t10 + Ifges(5,4) * t42) * t10 + m(5) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t11 ^ 2 + t12 ^ 2 + t22 ^ 2) + m(3) * (t39 * qJ(2) ^ 2 + pkin(1) ^ 2); -m(3) * pkin(1) + m(4) * t22 + m(5) * t13 + t36 + t37 + t38; m(3) + m(4) + m(5); t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,5) * t19 + Ifges(4,6) * t18 + (m(5) * (t2 * t30 + t28 * t3) + (-t30 * t10 + t28 * t9) * mrSges(5,3)) * pkin(3) + t35; 0; Ifges(4,3) + Ifges(5,3) + m(5) * (t28 ^ 2 + t30 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t34; t35; 0; Ifges(5,3) + t34; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
