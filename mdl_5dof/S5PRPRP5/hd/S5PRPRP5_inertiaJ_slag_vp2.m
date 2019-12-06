% Calculate joint inertia matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:36
% DurationCPUTime: 0.26s
% Computational Cost: add. (240->79), mult. (533->104), div. (0->0), fcn. (468->6), ass. (0->35)
t49 = m(5) + m(6);
t43 = mrSges(6,2) + mrSges(5,3);
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t18 = -t28 * mrSges(4,1) + t27 * mrSges(4,2);
t29 = sin(qJ(4));
t44 = cos(qJ(4));
t17 = t44 * t27 + t29 * t28;
t33 = -t29 * t27 + t44 * t28;
t3 = -mrSges(6,1) * t33 - t17 * mrSges(6,3);
t4 = -mrSges(5,1) * t33 + t17 * mrSges(5,2);
t48 = -t18 - t3 - t4;
t47 = -m(6) * pkin(4) - mrSges(6,1);
t46 = mrSges(5,1) - t47;
t45 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t24 = t28 ^ 2;
t42 = pkin(6) + qJ(3);
t41 = t27 ^ 2 + t24;
t19 = t42 * t28;
t36 = t42 * t27;
t6 = t29 * t19 + t44 * t36;
t8 = t44 * t19 - t29 * t36;
t40 = t6 ^ 2 + t8 ^ 2;
t39 = m(4) + t49;
t21 = -t28 * pkin(3) - pkin(2);
t30 = sin(qJ(2));
t10 = t17 * t30;
t12 = t33 * t30;
t38 = t6 * t10 + t8 * t12;
t35 = qJ(3) * t41;
t31 = cos(qJ(2));
t26 = t31 ^ 2;
t25 = t30 ^ 2;
t2 = -pkin(4) * t33 - t17 * qJ(5) + t21;
t1 = [m(2) + m(3) * (t25 + t26) + m(4) * (t41 * t25 + t26) + t49 * (t10 ^ 2 + t12 ^ 2 + t26); (t41 * mrSges(4,3) - mrSges(3,2)) * t30 + (mrSges(3,1) + t48) * t31 + m(4) * (t31 * pkin(2) + t30 * t35) + m(5) * (-t21 * t31 + t38) + m(6) * (-t2 * t31 + t38) + t43 * (t10 * t17 + t12 * t33); Ifges(4,2) * t24 - 0.2e1 * pkin(2) * t18 + 0.2e1 * t2 * t3 + 0.2e1 * t21 * t4 + Ifges(3,3) + (Ifges(4,1) * t27 + 0.2e1 * Ifges(4,4) * t28) * t27 + 0.2e1 * mrSges(4,3) * t35 + m(5) * (t21 ^ 2 + t40) + m(6) * (t2 ^ 2 + t40) + m(4) * (t41 * qJ(3) ^ 2 + pkin(2) ^ 2) + ((Ifges(6,1) + Ifges(5,1)) * t17 + 0.2e1 * t43 * t6) * t17 - (-(Ifges(6,3) + Ifges(5,2)) * t33 - 0.2e1 * t43 * t8 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t17) * t33; -t39 * t31; -m(4) * pkin(2) + m(5) * t21 + m(6) * t2 - t48; t39; -t46 * t10 + t45 * t12; (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t17 - (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t33 + t45 * t8 - t46 * t6; 0; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t10; m(6) * t6 + t17 * mrSges(6,2); 0; t47; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
