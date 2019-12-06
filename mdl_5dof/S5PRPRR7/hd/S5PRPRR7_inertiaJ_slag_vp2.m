% Calculate joint inertia matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:45
% EndTime: 2019-12-05 15:59:46
% DurationCPUTime: 0.26s
% Computational Cost: add. (238->87), mult. (483->127), div. (0->0), fcn. (388->6), ass. (0->39)
t32 = cos(qJ(4));
t25 = t32 ^ 2;
t29 = sin(qJ(4));
t48 = t29 ^ 2 + t25;
t45 = m(5) * t48;
t53 = -m(4) - t45;
t52 = -t48 * mrSges(5,3) + mrSges(4,2);
t51 = m(6) * pkin(4);
t34 = -pkin(2) - pkin(6);
t50 = -pkin(7) + t34;
t49 = mrSges(5,1) * t29 + mrSges(5,2) * t32 + mrSges(4,3);
t28 = sin(qJ(5));
t31 = cos(qJ(5));
t12 = -t28 * t32 - t31 * t29;
t38 = t28 * t29 - t31 * t32;
t47 = t12 ^ 2 + t38 ^ 2;
t33 = cos(qJ(2));
t6 = t38 * t33;
t7 = t12 * t33;
t46 = t6 * mrSges(6,1) - t7 * mrSges(6,2);
t44 = -mrSges(6,1) * t38 + t12 * mrSges(6,2);
t43 = t48 * t34;
t42 = t12 * t7 + t38 * t6;
t15 = t50 * t29;
t16 = t50 * t32;
t3 = -t15 * t28 + t16 * t31;
t4 = t15 * t31 + t16 * t28;
t41 = t3 * mrSges(6,1) - t4 * mrSges(6,2) - Ifges(6,5) * t38 + Ifges(6,6) * t12;
t40 = -t32 * mrSges(5,1) + t29 * mrSges(5,2);
t39 = t12 * t28 + t31 * t38;
t37 = (mrSges(6,1) * t31 - mrSges(6,2) * t28) * pkin(4);
t35 = qJ(3) ^ 2;
t30 = sin(qJ(2));
t26 = t33 ^ 2;
t24 = t30 ^ 2;
t22 = t30 * qJ(3);
t19 = pkin(4) * t29 + qJ(3);
t2 = -mrSges(6,1) * t12 - mrSges(6,2) * t38;
t1 = [m(2) + m(5) * (t48 * t26 + t24) + m(6) * (t6 ^ 2 + t7 ^ 2 + t24) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t24 + t26); t42 * mrSges(6,3) + (mrSges(3,1) - t52) * t33 + (-mrSges(3,2) + t2 + t49) * t30 + m(4) * (pkin(2) * t33 + t22) + m(5) * (-t33 * t43 + t22) + m(6) * (t19 * t30 + t3 * t6 + t4 * t7); Ifges(5,1) * t25 - 0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t19 * t2 + Ifges(4,1) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t32 + Ifges(5,2) * t29) * t29 - (-0.2e1 * mrSges(6,3) * t3 - Ifges(6,1) * t38) * t38 + (0.2e1 * mrSges(6,3) * t4 - 0.2e1 * Ifges(6,4) * t38 + Ifges(6,2) * t12) * t12 + m(6) * (t19 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t48 * t34 ^ 2 + t35) + m(4) * (pkin(2) ^ 2 + t35) - 0.2e1 * mrSges(5,3) * t43 + 0.2e1 * t49 * qJ(3); -m(6) * t42 + t53 * t33; -m(4) * pkin(2) - t47 * mrSges(6,3) + m(6) * (-t12 * t4 - t3 * t38) + t34 * t45 + t52; m(6) * t47 - t53; t40 * t33 + (t28 * t7 + t31 * t6) * t51 + t46; (mrSges(5,1) * t34 + Ifges(5,5)) * t32 + (-mrSges(5,2) * t34 - Ifges(5,6)) * t29 + (m(6) * (t28 * t4 + t3 * t31) + t39 * mrSges(6,3)) * pkin(4) + t41; -t39 * t51 - t40 + t44; Ifges(5,3) + Ifges(6,3) + m(6) * (t28 ^ 2 + t31 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t37; t46; t41; t44; Ifges(6,3) + t37; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
