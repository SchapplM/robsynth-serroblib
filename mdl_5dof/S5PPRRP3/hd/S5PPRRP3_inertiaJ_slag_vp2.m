% Calculate joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:35
% DurationCPUTime: 0.24s
% Computational Cost: add. (135->68), mult. (383->92), div. (0->0), fcn. (262->6), ass. (0->33)
t26 = cos(qJ(4));
t20 = t26 ^ 2;
t24 = sin(qJ(4));
t36 = t24 ^ 2 + t20;
t23 = cos(pkin(8));
t22 = sin(pkin(8));
t27 = cos(qJ(3));
t39 = t22 * t27;
t5 = t23 * t26 + t24 * t39;
t7 = -t23 * t24 + t26 * t39;
t50 = t24 * t5 + t7 * t26;
t49 = m(6) + m(5);
t48 = mrSges(6,2) + mrSges(5,3);
t47 = -m(6) * pkin(4) - mrSges(6,1);
t46 = t48 * t36;
t45 = -mrSges(5,1) + t47;
t44 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t25 = sin(qJ(3));
t40 = t22 * t25;
t38 = t36 * pkin(6) * t25;
t37 = t36 * pkin(6) ^ 2;
t10 = -t26 * mrSges(5,1) + t24 * mrSges(5,2);
t9 = -t26 * mrSges(6,1) - t24 * mrSges(6,3);
t35 = mrSges(4,1) - t10 - t9;
t33 = t50 * pkin(6);
t29 = t45 * t24 + t44 * t26;
t21 = t27 ^ 2;
t19 = t25 ^ 2;
t17 = t23 ^ 2;
t16 = t22 ^ 2;
t13 = t19 * t16;
t8 = -t26 * pkin(4) - t24 * qJ(5) - pkin(3);
t1 = [m(2) + m(3) * (t16 + t17) + m(4) * (t21 * t16 + t13 + t17) + t49 * (t5 ^ 2 + t7 ^ 2 + t13); t49 * (-t39 + t50) * t25; m(3) + m(4) * (t19 + t21) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t19 * t36 + t21); (-t27 * mrSges(4,2) - t25 * t35) * t22 + m(5) * (-pkin(3) * t40 + t33) + m(6) * (t8 * t40 + t33) + t48 * t50; t35 * t27 + m(6) * (-t8 * t27 + t38) + m(5) * (pkin(3) * t27 + t38) + (-mrSges(4,2) + t46) * t25; -0.2e1 * pkin(3) * t10 + 0.2e1 * t8 * t9 + Ifges(4,3) + m(6) * (t8 ^ 2 + t37) + m(5) * (pkin(3) ^ 2 + t37) + (Ifges(5,2) + Ifges(6,3)) * t20 + ((Ifges(6,1) + Ifges(5,1)) * t24 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t26) * t24 + 0.2e1 * pkin(6) * t46; t44 * t7 + t45 * t5; t29 * t25; (qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t26 + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t24 + t29 * pkin(6); Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t5; m(6) * t24 * t25; (m(6) * pkin(6) + mrSges(6,2)) * t24; t47; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
