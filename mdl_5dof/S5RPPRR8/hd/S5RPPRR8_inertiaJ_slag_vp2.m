% Calculate joint inertia matrix for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:05
% DurationCPUTime: 0.19s
% Computational Cost: add. (304->78), mult. (459->104), div. (0->0), fcn. (342->6), ass. (0->31)
t27 = cos(qJ(5));
t42 = t27 ^ 2;
t25 = sin(qJ(5));
t36 = t25 ^ 2 + t42;
t34 = t36 * mrSges(6,3);
t15 = -t27 * mrSges(6,1) + t25 * mrSges(6,2);
t41 = mrSges(5,1) - t15;
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t9 = t26 * t23 - t28 * t24;
t40 = t9 ^ 2;
t39 = -0.2e1 * t15;
t29 = -pkin(1) - pkin(2);
t13 = -t23 * qJ(2) + t24 * t29;
t12 = -pkin(3) + t13;
t14 = t24 * qJ(2) + t23 * t29;
t5 = t26 * t12 + t28 * t14;
t4 = t28 * t12 - t26 * t14;
t38 = t4 * mrSges(5,1);
t37 = t5 * mrSges(5,2);
t3 = -pkin(7) + t5;
t35 = t36 * t3;
t11 = t28 * t23 + t26 * t24;
t33 = t36 * t11;
t32 = -mrSges(6,1) * t25 - mrSges(6,2) * t27;
t31 = Ifges(6,2) * t42 + Ifges(5,3) + (Ifges(6,1) * t25 + 0.2e1 * Ifges(6,4) * t27) * t25;
t8 = t11 ^ 2;
t2 = pkin(4) - t4;
t1 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t13 * mrSges(4,1) - 0.2e1 * t38 + 0.2e1 * t14 * mrSges(4,2) + 0.2e1 * t37 + 0.2e1 * qJ(2) * mrSges(3,3) + t2 * t39 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) - 0.2e1 * t3 * t34 + m(6) * (t36 * t3 ^ 2 + t2 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(4) * (t13 ^ 2 + t14 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + t31; -m(3) * pkin(1) - t24 * mrSges(4,1) + t23 * mrSges(4,2) - mrSges(3,1) + t41 * t9 + (mrSges(5,2) - t34) * t11 + m(6) * (t9 * t2 + t3 * t33) + m(5) * (t11 * t5 - t9 * t4) + m(4) * (t24 * t13 + t23 * t14); m(3) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(5) * (t8 + t40) + m(6) * (t36 * t8 + t40); 0; 0; m(6) * t36 + m(4) + m(5); m(6) * (-pkin(4) * t2 + pkin(7) * t35) - t37 + t38 + (t2 + pkin(4)) * t15 + (-t36 * pkin(7) + t35) * mrSges(6,3) - t31; -t11 * mrSges(5,2) + (m(6) * pkin(7) + mrSges(6,3)) * t33 + (-m(6) * pkin(4) - t41) * t9; 0; pkin(4) * t39 + m(6) * (t36 * pkin(7) ^ 2 + pkin(4) ^ 2) + 0.2e1 * pkin(7) * t34 + t31; (-mrSges(6,2) * t3 - Ifges(6,6)) * t27 + (-mrSges(6,1) * t3 - Ifges(6,5)) * t25; t32 * t11; -t15; Ifges(6,5) * t25 + Ifges(6,6) * t27 + t32 * pkin(7); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
