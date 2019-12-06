% Calculate joint inertia matrix for
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:32:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (154->71), mult. (334->94), div. (0->0), fcn. (245->6), ass. (0->28)
t36 = 2 * mrSges(6,3);
t35 = m(5) + m(6);
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t5 = t19 * t22 - t20 * t24;
t34 = t5 ^ 2;
t33 = m(6) * pkin(4);
t21 = sin(qJ(4));
t32 = t21 * mrSges(5,2);
t23 = cos(qJ(4));
t31 = t21 ^ 2 + t23 ^ 2;
t14 = t19 * pkin(2) + pkin(6);
t30 = qJ(5) + t14;
t28 = mrSges(6,1) + t33;
t15 = -t20 * pkin(2) - pkin(3);
t27 = t31 * mrSges(5,3);
t16 = t21 * mrSges(6,2);
t9 = -t23 * mrSges(6,1) + t16;
t26 = mrSges(5,1) + t28;
t10 = -t23 * mrSges(5,1) + t32;
t8 = -t23 * pkin(4) + t15;
t7 = t19 * t24 + t20 * t22;
t4 = t7 ^ 2;
t3 = t30 * t23;
t2 = t30 * t21;
t1 = [m(2) + m(3) * (t22 ^ 2 + t24 ^ 2) + m(4) * (t4 + t34) + (t31 * t4 + t34) * t35; t24 * mrSges(3,1) - t22 * mrSges(3,2) + (-mrSges(4,1) + t10 + t9) * t5 + (t31 * mrSges(6,3) - mrSges(4,2) + t27) * t7 + m(5) * (t31 * t7 * t14 + t15 * t5) + m(6) * (t8 * t5 + (t2 * t21 + t23 * t3) * t7) + m(4) * (t19 * t7 - t20 * t5) * pkin(2); 0.2e1 * t15 * t10 + 0.2e1 * t8 * t9 + Ifges(3,3) + Ifges(4,3) + m(6) * (t2 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(5) * (t31 * t14 ^ 2 + t15 ^ 2) + m(4) * (t19 ^ 2 + t20 ^ 2) * pkin(2) ^ 2 + (t3 * t36 + (Ifges(5,2) + Ifges(6,2)) * t23) * t23 + 0.2e1 * (t20 * mrSges(4,1) - t19 * mrSges(4,2)) * pkin(2) + 0.2e1 * t14 * t27 + (t2 * t36 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t23 + (Ifges(5,1) + Ifges(6,1)) * t21) * t21; 0; m(6) * (-t23 * t2 + t21 * t3); t31 * t35 + m(4); ((-mrSges(5,2) - mrSges(6,2)) * t23 - t26 * t21) * t7; -t3 * mrSges(6,2) - t28 * t2 + (-mrSges(5,2) * t14 + Ifges(5,6) + Ifges(6,6)) * t23 + (-mrSges(5,1) * t14 - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t21; t26 * t23 - t16 - t32; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t33) * pkin(4); m(6) * t5; m(6) * t8 + t9; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
