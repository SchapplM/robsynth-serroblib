% Calculate joint inertia matrix for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:06:40
% DurationCPUTime: 0.17s
% Computational Cost: add. (118->61), mult. (277->77), div. (0->0), fcn. (203->6), ass. (0->26)
t33 = -2 * mrSges(6,3);
t32 = m(5) + m(6);
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t19 = sin(qJ(3));
t29 = cos(qJ(3));
t3 = t19 * t16 - t29 * t17;
t31 = t3 ^ 2;
t30 = m(6) * pkin(4);
t18 = sin(qJ(4));
t28 = t18 * mrSges(5,2);
t27 = -qJ(5) - pkin(6);
t20 = cos(qJ(4));
t26 = t18 ^ 2 + t20 ^ 2;
t24 = mrSges(6,1) + t30;
t23 = t26 * mrSges(5,3);
t13 = t18 * mrSges(6,2);
t7 = -t20 * mrSges(6,1) + t13;
t22 = mrSges(5,1) + t24;
t12 = -t20 * pkin(4) - pkin(3);
t9 = t27 * t20;
t8 = -t20 * mrSges(5,1) + t28;
t6 = t27 * t18;
t5 = t29 * t16 + t19 * t17;
t2 = t5 ^ 2;
t1 = [m(2) + m(3) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t2 + t31) + (t26 * t2 + t31) * t32; 0; t26 * t32 + m(3) + m(4); (-mrSges(4,1) + t7 + t8) * t3 + (t26 * mrSges(6,3) - mrSges(4,2) + t23) * t5 + m(5) * (t26 * t5 * pkin(6) - pkin(3) * t3) + m(6) * (t12 * t3 + (-t18 * t6 - t20 * t9) * t5); m(6) * (-t9 * t18 + t6 * t20); -0.2e1 * pkin(3) * t8 + 0.2e1 * t12 * t7 + Ifges(4,3) + m(5) * (t26 * pkin(6) ^ 2 + pkin(3) ^ 2) + m(6) * (t12 ^ 2 + t6 ^ 2 + t9 ^ 2) + (t9 * t33 + (Ifges(6,2) + Ifges(5,2)) * t20) * t20 + 0.2e1 * pkin(6) * t23 + (t6 * t33 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t20 + (Ifges(6,1) + Ifges(5,1)) * t18) * t18; ((-mrSges(5,2) - mrSges(6,2)) * t20 - t22 * t18) * t5; t22 * t20 - t13 - t28; t9 * mrSges(6,2) + t24 * t6 + (-mrSges(5,2) * pkin(6) + Ifges(5,6) + Ifges(6,6)) * t20 + (-mrSges(5,1) * pkin(6) - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t18; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t30) * pkin(4); m(6) * t3; 0; m(6) * t12 + t7; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
