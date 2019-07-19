% Calculate joint inertia matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:31
% EndTime: 2019-07-18 18:16:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (139->57), mult. (194->69), div. (0->0), fcn. (99->4), ass. (0->22)
t23 = 2 * mrSges(4,3);
t11 = cos(qJ(4));
t12 = cos(qJ(2));
t7 = -t12 * pkin(1) - pkin(2);
t5 = -pkin(3) + t7;
t10 = sin(qJ(2));
t6 = t10 * pkin(1) + qJ(3);
t9 = sin(qJ(4));
t1 = t11 * t5 - t9 * t6;
t22 = t1 * mrSges(5,1);
t2 = t11 * t6 + t9 * t5;
t21 = t2 * mrSges(5,2);
t13 = -pkin(2) - pkin(3);
t3 = -t9 * qJ(3) + t11 * t13;
t20 = t3 * mrSges(5,1);
t4 = t11 * qJ(3) + t9 * t13;
t19 = t4 * mrSges(5,2);
t18 = Ifges(4,2) + Ifges(3,3) + Ifges(5,3);
t17 = t11 * mrSges(5,1) - t9 * mrSges(5,2);
t16 = -mrSges(4,1) - t17;
t15 = (t12 * mrSges(3,1) - t10 * mrSges(3,2)) * pkin(1);
t8 = [-0.2e1 * t7 * mrSges(4,1) - 0.2e1 * t22 + 0.2e1 * t21 + t6 * t23 + Ifges(2,3) + 0.2e1 * t15 + m(5) * (t1 ^ 2 + t2 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2) + m(3) * (t10 ^ 2 + t12 ^ 2) * pkin(1) ^ 2 + t18; t15 + (qJ(3) + t6) * mrSges(4,3) + (t4 + t2) * mrSges(5,2) + (-t1 - t3) * mrSges(5,1) + (-t7 + pkin(2)) * mrSges(4,1) + m(5) * (t3 * t1 + t4 * t2) + m(4) * (-pkin(2) * t7 + qJ(3) * t6) + t18; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t20 + 0.2e1 * t19 + qJ(3) * t23 + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t18; m(5) * (t11 * t1 + t9 * t2) + m(4) * t7 + t16; -m(4) * pkin(2) + m(5) * (t11 * t3 + t9 * t4) + t16; m(4) + m(5) * (t11 ^ 2 + t9 ^ 2); -Ifges(5,3) - t21 + t22; -Ifges(5,3) - t19 + t20; t17; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq  = res;
