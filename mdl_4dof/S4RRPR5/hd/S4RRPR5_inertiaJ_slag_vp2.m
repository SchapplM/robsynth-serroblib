% Calculate joint inertia matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (129->54), mult. (220->66), div. (0->0), fcn. (91->4), ass. (0->25)
t15 = cos(qJ(4));
t35 = t15 ^ 2;
t13 = sin(qJ(4));
t34 = t13 * mrSges(5,1) + t15 * mrSges(5,2) + mrSges(4,3);
t28 = t13 ^ 2 + t35;
t26 = t28 * mrSges(5,3);
t14 = sin(qJ(2));
t6 = t14 * pkin(1) + qJ(3);
t33 = t6 ^ 2;
t31 = qJ(3) * t6;
t30 = t15 * mrSges(5,1);
t16 = cos(qJ(2));
t8 = -t16 * pkin(1) - pkin(2);
t27 = m(5) * t28;
t17 = -pkin(2) - pkin(6);
t25 = t28 * t17;
t24 = 0.2e1 * t34;
t23 = mrSges(4,2) - t26;
t22 = Ifges(5,1) * t35 + Ifges(4,1) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t15 + Ifges(5,2) * t13) * t13;
t21 = -0.2e1 * t26;
t20 = (t16 * mrSges(3,1) - t14 * mrSges(3,2)) * pkin(1);
t18 = qJ(3) ^ 2;
t9 = Ifges(5,5) * t15;
t5 = -pkin(6) + t8;
t1 = [0.2e1 * t8 * mrSges(4,2) + Ifges(2,3) + t6 * t24 + 0.2e1 * t20 + t5 * t21 + m(5) * (t28 * t5 ^ 2 + t33) + m(4) * (t8 ^ 2 + t33) + m(3) * (t14 ^ 2 + t16 ^ 2) * pkin(1) ^ 2 + t22; t20 + (-pkin(2) + t8) * mrSges(4,2) + m(5) * (t5 * t25 + t31) + m(4) * (-pkin(2) * t8 + t31) + t22 + (-t17 - t5) * t26 + t34 * (qJ(3) + t6); -0.2e1 * pkin(2) * mrSges(4,2) + qJ(3) * t24 + t17 * t21 + m(5) * (t28 * t17 ^ 2 + t18) + m(4) * (pkin(2) ^ 2 + t18) + t22; m(4) * t8 + t5 * t27 + t23; -m(4) * pkin(2) + m(5) * t25 + t23; m(4) + t27; t5 * t30 + t9 + (-mrSges(5,2) * t5 - Ifges(5,6)) * t13; t17 * t30 + t9 + (-mrSges(5,2) * t17 - Ifges(5,6)) * t13; -t13 * mrSges(5,2) + t30; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
