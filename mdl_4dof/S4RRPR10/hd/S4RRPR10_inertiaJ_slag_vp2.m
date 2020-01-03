% Calculate joint inertia matrix for
% S4RRPR10
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:11
% DurationCPUTime: 0.30s
% Computational Cost: add. (194->95), mult. (375->129), div. (0->0), fcn. (244->4), ass. (0->40)
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t48 = t25 ^ 2 + t27 ^ 2;
t47 = -m(4) * pkin(2) + mrSges(4,2);
t24 = sin(qJ(4));
t46 = -t24 / 0.2e1;
t28 = -pkin(2) - pkin(6);
t45 = pkin(3) + pkin(5);
t26 = cos(qJ(4));
t43 = mrSges(5,3) * t27;
t8 = -t25 * mrSges(5,2) - t26 * t43;
t44 = t24 * t8;
t42 = Ifges(5,4) * t24;
t41 = Ifges(5,4) * t26;
t40 = t26 * t28;
t39 = t48 * pkin(5) ^ 2;
t38 = t24 ^ 2 + t26 ^ 2;
t37 = m(5) * t38;
t36 = t38 * mrSges(5,3);
t35 = -t25 * qJ(3) - pkin(1);
t13 = t45 * t25;
t6 = t28 * t27 + t35;
t1 = t26 * t13 - t24 * t6;
t2 = t24 * t13 + t26 * t6;
t33 = t26 * t1 + t24 * t2;
t32 = t26 * mrSges(5,1) - t24 * mrSges(5,2);
t31 = -Ifges(5,5) * t24 - Ifges(5,6) * t26;
t29 = qJ(3) ^ 2;
t16 = Ifges(5,5) * t26;
t15 = Ifges(5,3) * t25;
t14 = t45 * t27;
t12 = Ifges(5,1) * t26 - t42;
t11 = -Ifges(5,2) * t24 + t41;
t10 = t24 * mrSges(5,1) + t26 * mrSges(5,2);
t9 = -t27 * pkin(2) + t35;
t7 = t25 * mrSges(5,1) + t24 * t43;
t5 = t32 * t27;
t4 = Ifges(5,5) * t25 + (-Ifges(5,1) * t24 - t41) * t27;
t3 = Ifges(5,6) * t25 + (-Ifges(5,2) * t26 - t42) * t27;
t17 = [0.2e1 * t1 * t7 + 0.2e1 * t14 * t5 + 0.2e1 * t2 * t8 + Ifges(2,3) + m(5) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(4) * (t9 ^ 2 + t39) + m(3) * (pkin(1) ^ 2 + t39) + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t9 * mrSges(4,3) + t15 + (Ifges(4,2) + Ifges(3,1)) * t25) * t25 + (-t24 * t4 - t26 * t3 + 0.2e1 * t9 * mrSges(4,2) + 0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t27 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t31) * t25) * t27 + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(5) * t48; t26 * t4 / 0.2e1 + t3 * t46 + qJ(3) * t5 + m(5) * (qJ(3) * t14 + t33 * t28) + t14 * t10 + t28 * t44 + t7 * t40 - t33 * mrSges(5,3) + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,6) * t46 + t16 / 0.2e1 - pkin(2) * mrSges(4,1)) * t25 + (-Ifges(4,5) + Ifges(3,6) + qJ(3) * mrSges(4,1) + t12 * t46 - t26 * t11 / 0.2e1) * t27 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t27 + (-mrSges(3,1) + t47) * t25) * pkin(5); -0.2e1 * pkin(2) * mrSges(4,2) - t24 * t11 + t26 * t12 + Ifges(4,1) + Ifges(3,3) + m(5) * (t38 * t28 ^ 2 + t29) + m(4) * (pkin(2) ^ 2 + t29) + 0.2e1 * (t10 + mrSges(4,3)) * qJ(3) - 0.2e1 * t28 * t36; m(5) * t33 + t44 + t26 * t7 + (m(4) * pkin(5) + mrSges(4,1)) * t25; t28 * t37 - t36 + t47; m(4) + t37; t1 * mrSges(5,1) - t2 * mrSges(5,2) + t31 * t27 + t15; mrSges(5,1) * t40 + t16 + (-mrSges(5,2) * t28 - Ifges(5,6)) * t24; t32; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1), t17(2), t17(4), t17(7); t17(2), t17(3), t17(5), t17(8); t17(4), t17(5), t17(6), t17(9); t17(7), t17(8), t17(9), t17(10);];
Mq = res;
