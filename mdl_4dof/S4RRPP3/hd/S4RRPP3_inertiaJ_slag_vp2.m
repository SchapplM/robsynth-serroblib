% Calculate joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (181->75), mult. (366->105), div. (0->0), fcn. (297->4), ass. (0->22)
t21 = cos(qJ(2));
t15 = -t21 * pkin(2) - pkin(1);
t29 = 0.2e1 * t15;
t28 = mrSges(5,2) + mrSges(4,3);
t27 = -qJ(3) - pkin(5);
t20 = sin(qJ(2));
t26 = t20 ^ 2 + t21 ^ 2;
t11 = t27 * t21;
t18 = sin(pkin(6));
t19 = cos(pkin(6));
t24 = t27 * t20;
t3 = -t18 * t11 - t19 * t24;
t5 = -t19 * t11 + t18 * t24;
t25 = t3 ^ 2 + t5 ^ 2;
t14 = -t19 * pkin(2) - pkin(3);
t12 = t18 * pkin(2) + qJ(4);
t9 = t18 * t21 + t19 * t20;
t8 = t18 * t20 - t19 * t21;
t7 = t9 * mrSges(4,2);
t6 = t8 * mrSges(5,1);
t1 = t8 * pkin(3) - t9 * qJ(4) + t15;
t2 = [Ifges(2,3) + 0.2e1 * t1 * t6 + t20 * (Ifges(3,1) * t20 + Ifges(3,4) * t21) + t21 * (Ifges(3,4) * t20 + Ifges(3,2) * t21) - 0.2e1 * pkin(1) * (-t21 * mrSges(3,1) + t20 * mrSges(3,2)) + t7 * t29 + 0.2e1 * t26 * pkin(5) * mrSges(3,3) + m(5) * (t1 ^ 2 + t25) + m(4) * (t15 ^ 2 + t25) + m(3) * (t26 * pkin(5) ^ 2 + pkin(1) ^ 2) + (-0.2e1 * t1 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t9 + 0.2e1 * t28 * t3) * t9 + (mrSges(4,1) * t29 + (Ifges(5,3) + Ifges(4,2)) * t8 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t9 - 0.2e1 * t28 * t5) * t8; m(5) * (t12 * t5 + t14 * t3) - t5 * mrSges(4,2) - t3 * mrSges(4,1) - t3 * mrSges(5,1) + t5 * mrSges(5,3) + Ifges(3,5) * t20 + Ifges(3,6) * t21 + (-t20 * mrSges(3,1) - t21 * mrSges(3,2)) * pkin(5) + (t14 * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t9 + (-t12 * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t8 + (m(4) * (t18 * t5 - t19 * t3) + (-t18 * t8 - t19 * t9) * mrSges(4,3)) * pkin(2); -0.2e1 * t14 * mrSges(5,1) + 0.2e1 * t12 * mrSges(5,3) + Ifges(5,2) + Ifges(3,3) + Ifges(4,3) + m(5) * (t12 ^ 2 + t14 ^ 2) + (0.2e1 * t19 * mrSges(4,1) - 0.2e1 * t18 * mrSges(4,2) + m(4) * (t18 ^ 2 + t19 ^ 2) * pkin(2)) * pkin(2); m(4) * t15 + m(5) * t1 + t8 * mrSges(4,1) - t9 * mrSges(5,3) + t6 + t7; 0; m(4) + m(5); m(5) * t3 + t9 * mrSges(5,2); m(5) * t14 - mrSges(5,1); 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
