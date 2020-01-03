% Calculate joint inertia matrix for
% S4RPRR3
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:06
% EndTime: 2019-12-31 16:49:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (157->57), mult. (310->94), div. (0->0), fcn. (248->6), ass. (0->23)
t23 = cos(qJ(3));
t31 = t23 ^ 2;
t18 = sin(pkin(7));
t14 = t18 * pkin(1) + pkin(5);
t30 = pkin(6) + t14;
t21 = sin(qJ(3));
t29 = t21 ^ 2 + t31;
t19 = cos(pkin(7));
t15 = -t19 * pkin(1) - pkin(2);
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t10 = -t20 * t21 + t22 * t23;
t11 = t20 * t23 + t22 * t21;
t4 = -t10 * mrSges(5,1) + t11 * mrSges(5,2);
t5 = t30 * t21;
t6 = t30 * t23;
t2 = -t20 * t6 - t22 * t5;
t3 = -t20 * t5 + t22 * t6;
t28 = t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t11 + Ifges(5,6) * t10;
t27 = -t23 * mrSges(4,1) + t21 * mrSges(4,2);
t26 = (t22 * mrSges(5,1) - t20 * mrSges(5,2)) * pkin(3);
t12 = -t23 * pkin(3) + t15;
t1 = [Ifges(2,3) + Ifges(3,3) + 0.2e1 * t12 * t4 + 0.2e1 * t15 * t27 + Ifges(4,2) * t31 + (-0.2e1 * t2 * mrSges(5,3) + Ifges(5,1) * t11) * t11 + (0.2e1 * t3 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t11 + Ifges(5,2) * t10) * t10 + m(5) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t29 * t14 ^ 2 + t15 ^ 2) + m(3) * (t18 ^ 2 + t19 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t21 + 0.2e1 * Ifges(4,4) * t23) * t21 + 0.2e1 * (t19 * mrSges(3,1) - t18 * mrSges(3,2)) * pkin(1) + 0.2e1 * t29 * t14 * mrSges(4,3); m(5) * (t2 * t10 + t3 * t11); m(3) + m(4) * t29 + m(5) * (t10 ^ 2 + t11 ^ 2); Ifges(4,5) * t21 + Ifges(4,6) * t23 + (-t21 * mrSges(4,1) - t23 * mrSges(4,2)) * t14 + (m(5) * (t2 * t22 + t20 * t3) + (t20 * t10 - t22 * t11) * mrSges(5,3)) * pkin(3) + t28; m(5) * (t10 * t22 + t11 * t20) * pkin(3) - t27 - t4; Ifges(4,3) + Ifges(5,3) + m(5) * (t20 ^ 2 + t22 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t26; t28; -t4; Ifges(5,3) + t26; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
