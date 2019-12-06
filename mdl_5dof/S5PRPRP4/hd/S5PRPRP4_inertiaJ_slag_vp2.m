% Calculate joint inertia matrix for
% S5PRPRP4
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:05
% DurationCPUTime: 0.19s
% Computational Cost: add. (151->67), mult. (338->84), div. (0->0), fcn. (237->6), ass. (0->27)
t24 = cos(qJ(4));
t19 = t24 ^ 2;
t22 = sin(qJ(4));
t31 = t22 ^ 2 + t19;
t40 = m(5) + m(6);
t39 = m(4) * pkin(2);
t10 = -t24 * mrSges(5,1) + t22 * mrSges(5,2);
t9 = -t24 * mrSges(6,1) - t22 * mrSges(6,3);
t38 = t10 + t9;
t37 = -m(6) * pkin(4) - mrSges(6,1);
t36 = (mrSges(6,2) + mrSges(5,3)) * t31;
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t6 = t20 * t23 - t21 * t25;
t35 = t6 ^ 2;
t16 = t20 * pkin(2) + pkin(6);
t8 = t20 * t25 + t21 * t23;
t34 = t31 * t16 * t8;
t32 = t31 * t16 ^ 2;
t17 = -t21 * pkin(2) - pkin(3);
t28 = t24 * pkin(4) + t22 * qJ(5);
t27 = (m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t24 + (-mrSges(5,1) + t37) * t22;
t5 = t8 ^ 2;
t4 = t17 - t28;
t1 = [m(2) + m(3) * (t23 ^ 2 + t25 ^ 2) + m(4) * (t5 + t35) + (t31 * t5 + t35) * t40; t25 * mrSges(3,1) - t23 * mrSges(3,2) + (-mrSges(4,1) + t38) * t6 + (-mrSges(4,2) + t36) * t8 + m(5) * (t17 * t6 + t34) + m(6) * (t4 * t6 + t34) + (t20 * t8 - t21 * t6) * t39; 0.2e1 * t17 * t10 + 0.2e1 * t4 * t9 + Ifges(3,3) + Ifges(4,3) + (Ifges(6,3) + Ifges(5,2)) * t19 + m(6) * (t4 ^ 2 + t32) + m(5) * (t17 ^ 2 + t32) + ((Ifges(6,1) + Ifges(5,1)) * t22 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t24) * t22 + 0.2e1 * t16 * t36 + (0.2e1 * t21 * mrSges(4,1) - 0.2e1 * t20 * mrSges(4,2) + (t20 ^ 2 + t21 ^ 2) * t39) * pkin(2); 0; 0; t31 * t40 + m(4); t27 * t8; (qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t24 + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t22 + t27 * t16; m(6) * t28 - t38; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t22 * t8; (m(6) * t16 + mrSges(6,2)) * t22; -m(6) * t24; t37; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
