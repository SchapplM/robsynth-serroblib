% Calculate joint inertia matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:11
% DurationCPUTime: 0.25s
% Computational Cost: add. (220->78), mult. (429->109), div. (0->0), fcn. (364->4), ass. (0->28)
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t41 = (t27 * mrSges(4,1) + (-mrSges(4,2) - mrSges(5,2)) * t25) * pkin(2);
t40 = 2 * mrSges(5,1);
t39 = m(5) * pkin(3);
t38 = -pkin(6) - pkin(5);
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t13 = -t25 * t26 + t27 * t28;
t37 = t25 * t13;
t35 = Ifges(4,3) + Ifges(5,3);
t18 = t38 * t26;
t19 = t38 * t28;
t6 = t25 * t18 - t27 * t19;
t34 = t26 ^ 2 + t28 ^ 2;
t21 = -t28 * pkin(2) - pkin(1);
t5 = t27 * t18 + t25 * t19;
t14 = t25 * t28 + t27 * t26;
t2 = -t14 * qJ(4) + t5;
t32 = m(5) * t2 - t14 * mrSges(5,3);
t3 = t13 * qJ(4) + t6;
t31 = t5 * mrSges(4,1) + t2 * mrSges(5,1) - t6 * mrSges(4,2) - t3 * mrSges(5,2) + (Ifges(4,5) + Ifges(5,5)) * t14 + (Ifges(4,6) + Ifges(5,6)) * t13;
t30 = pkin(2) ^ 2;
t22 = t25 ^ 2 * t30;
t20 = t27 * pkin(2) + pkin(3);
t8 = t14 * mrSges(5,2);
t7 = -t13 * pkin(3) + t21;
t1 = [Ifges(2,3) + 0.2e1 * t7 * t8 - 0.2e1 * pkin(1) * (-t28 * mrSges(3,1) + t26 * mrSges(3,2)) + t26 * (Ifges(3,1) * t26 + Ifges(3,4) * t28) + t28 * (Ifges(3,4) * t26 + Ifges(3,2) * t28) + m(5) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + m(4) * (t21 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(3) * (pkin(5) ^ 2 * t34 + pkin(1) ^ 2) + (0.2e1 * t21 * mrSges(4,2) - 0.2e1 * t5 * mrSges(4,3) - 0.2e1 * t2 * mrSges(5,3) + (Ifges(4,1) + Ifges(5,1)) * t14) * t14 + 0.2e1 * t34 * pkin(5) * mrSges(3,3) + (-0.2e1 * t21 * mrSges(4,1) - 0.2e1 * t7 * mrSges(5,1) + 0.2e1 * t6 * mrSges(4,3) + 0.2e1 * t3 * mrSges(5,3) + 0.2e1 * (Ifges(4,4) + Ifges(5,4)) * t14 + (Ifges(5,2) + Ifges(4,2)) * t13) * t13; Ifges(3,5) * t26 + Ifges(3,6) * t28 + t32 * t20 + (-t26 * mrSges(3,1) - t28 * mrSges(3,2)) * pkin(5) + (mrSges(5,3) * t37 + (-t27 * t14 + t37) * mrSges(4,3) + m(5) * t25 * t3 + m(4) * (t25 * t6 + t27 * t5)) * pkin(2) + t31; t20 * t40 + Ifges(3,3) + m(5) * (t20 ^ 2 + t22) + m(4) * (t27 ^ 2 * t30 + t22) + 0.2e1 * t41 + t35; pkin(3) * t32 + t31; t20 * t39 + (pkin(3) + t20) * mrSges(5,1) + t41 + t35; (t40 + t39) * pkin(3) + t35; m(5) * t7 - t13 * mrSges(5,1) + t8; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
