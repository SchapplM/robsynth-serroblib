% Calculate joint inertia matrix for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (248->75), mult. (474->104), div. (0->0), fcn. (386->6), ass. (0->35)
t30 = cos(pkin(7));
t47 = t30 ^ 2;
t29 = sin(pkin(7));
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t15 = -t31 * t29 + t33 * t30;
t16 = t33 * t29 + t31 * t30;
t5 = -t15 * mrSges(5,1) + t16 * mrSges(5,2);
t46 = 0.2e1 * t5;
t34 = cos(qJ(2));
t45 = t34 * pkin(1);
t44 = Ifges(5,5) * t16 + Ifges(5,6) * t15;
t43 = t29 ^ 2 + t47;
t42 = 2 * mrSges(5,3);
t23 = -t30 * pkin(3) - pkin(2);
t19 = -t30 * mrSges(4,1) + t29 * mrSges(4,2);
t41 = t43 * qJ(3);
t40 = Ifges(5,1) * t16 ^ 2 + Ifges(4,2) * t47 + Ifges(3,3) + (Ifges(4,1) * t29 + 0.2e1 * Ifges(4,4) * t30) * t29 + (0.2e1 * Ifges(5,4) * t16 + Ifges(5,2) * t15) * t15;
t39 = 0.2e1 * t43 * mrSges(4,3);
t32 = sin(qJ(2));
t38 = (t34 * mrSges(3,1) - t32 * mrSges(3,2)) * pkin(1);
t37 = t19 + t5;
t26 = t30 * pkin(6);
t24 = -pkin(2) - t45;
t22 = t32 * pkin(1) + qJ(3);
t20 = t30 * qJ(3) + t26;
t18 = (-pkin(6) - qJ(3)) * t29;
t17 = t23 - t45;
t9 = t30 * t22 + t26;
t8 = (-pkin(6) - t22) * t29;
t7 = t31 * t18 + t33 * t20;
t6 = t33 * t18 - t31 * t20;
t4 = t31 * t8 + t33 * t9;
t3 = -t31 * t9 + t33 * t8;
t1 = [t17 * t46 + 0.2e1 * t24 * t19 + Ifges(2,3) + 0.2e1 * t38 + (t4 * t15 - t3 * t16) * t42 + t22 * t39 + m(5) * (t17 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t43 * t22 ^ 2 + t24 ^ 2) + m(3) * (t32 ^ 2 + t34 ^ 2) * pkin(1) ^ 2 + t40; (t17 + t23) * t5 + (t24 - pkin(2)) * t19 + t38 + m(5) * (t23 * t17 + t6 * t3 + t7 * t4) + m(4) * (-pkin(2) * t24 + t22 * t41) + ((-t3 - t6) * t16 + (t4 + t7) * t15) * mrSges(5,3) + (t43 * t22 + t41) * mrSges(4,3) + t40; -0.2e1 * pkin(2) * t19 + t23 * t46 + (t7 * t15 - t6 * t16) * t42 + qJ(3) * t39 + m(5) * (t23 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(4) * (t43 * qJ(3) ^ 2 + pkin(2) ^ 2) + t40; m(4) * t24 + m(5) * t17 + t37; -m(4) * pkin(2) + m(5) * t23 + t37; m(4) + m(5); t3 * mrSges(5,1) - t4 * mrSges(5,2) + t44; t6 * mrSges(5,1) - t7 * mrSges(5,2) + t44; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
