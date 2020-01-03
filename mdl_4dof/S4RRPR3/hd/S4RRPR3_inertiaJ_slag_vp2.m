% Calculate joint inertia matrix for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (142->52), mult. (278->73), div. (0->0), fcn. (168->6), ass. (0->28)
t23 = cos(qJ(4));
t39 = t23 ^ 2;
t21 = sin(qJ(4));
t8 = -t23 * mrSges(5,1) + t21 * mrSges(5,2);
t38 = 0.2e1 * t8;
t22 = sin(qJ(2));
t37 = pkin(1) * t22;
t24 = cos(qJ(2));
t14 = t24 * pkin(1) + pkin(2);
t19 = sin(pkin(7));
t20 = cos(pkin(7));
t4 = t20 * t14 - t19 * t37;
t36 = t4 * mrSges(4,1);
t5 = t19 * t14 + t20 * t37;
t35 = t5 * mrSges(4,2);
t34 = Ifges(5,5) * t21 + Ifges(5,6) * t23;
t33 = t21 ^ 2 + t39;
t12 = t19 * pkin(2) + pkin(6);
t32 = t33 * t12;
t31 = Ifges(5,2) * t39 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t21 + 0.2e1 * Ifges(5,4) * t23) * t21;
t30 = t20 * mrSges(4,1) - t19 * mrSges(4,2);
t29 = -mrSges(5,1) * t21 - mrSges(5,2) * t23;
t28 = 0.2e1 * t33 * mrSges(5,3);
t27 = (t24 * mrSges(3,1) - t22 * mrSges(3,2)) * pkin(1);
t13 = -t20 * pkin(2) - pkin(3);
t3 = pkin(6) + t5;
t2 = -pkin(3) - t4;
t1 = [0.2e1 * t36 - 0.2e1 * t35 + t2 * t38 + Ifges(2,3) + 0.2e1 * t27 + t3 * t28 + m(5) * (t33 * t3 ^ 2 + t2 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(3) * (t22 ^ 2 + t24 ^ 2) * pkin(1) ^ 2 + t31; m(5) * (t13 * t2 + t3 * t32) - t35 + t36 + (t2 + t13) * t8 + t27 + (m(4) * (t19 * t5 + t20 * t4) + t30) * pkin(2) + (t33 * t3 + t32) * mrSges(5,3) + t31; t13 * t38 + t12 * t28 + m(5) * (t33 * t12 ^ 2 + t13 ^ 2) + t31 + (0.2e1 * t30 + m(4) * (t19 ^ 2 + t20 ^ 2) * pkin(2)) * pkin(2); 0; 0; m(5) * t33 + m(4); t29 * t3 + t34; t29 * t12 + t34; -t8; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
