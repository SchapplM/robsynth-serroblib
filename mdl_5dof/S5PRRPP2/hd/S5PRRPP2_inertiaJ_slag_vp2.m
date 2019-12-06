% Calculate joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:44
% EndTime: 2019-12-05 16:08:45
% DurationCPUTime: 0.34s
% Computational Cost: add. (273->101), mult. (615->137), div. (0->0), fcn. (526->6), ass. (0->35)
t48 = m(5) * pkin(3);
t46 = m(5) + m(6);
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t15 = t28 * t30 - t29 * t32;
t16 = t28 * t32 + t29 * t30;
t3 = t15 * mrSges(6,1) - t16 * mrSges(6,3);
t4 = t15 * mrSges(5,1) + t16 * mrSges(5,2);
t47 = -t3 - t4;
t45 = mrSges(6,2) + mrSges(5,3);
t26 = t32 ^ 2;
t44 = -qJ(4) - pkin(6);
t43 = t30 ^ 2 + t26;
t19 = t44 * t32;
t39 = t44 * t30;
t6 = -t28 * t19 - t29 * t39;
t8 = -t29 * t19 + t28 * t39;
t42 = t6 ^ 2 + t8 ^ 2;
t23 = -t32 * pkin(3) - pkin(2);
t31 = sin(qJ(2));
t10 = t16 * t31;
t12 = t15 * t31;
t41 = t6 * t10 - t8 * t12;
t38 = t43 * mrSges(4,3);
t37 = -t30 * mrSges(4,1) - t32 * mrSges(4,2);
t33 = cos(qJ(2));
t27 = t33 ^ 2;
t25 = t31 ^ 2;
t22 = -t29 * pkin(3) - pkin(4);
t20 = t28 * pkin(3) + qJ(5);
t18 = -t32 * mrSges(4,1) + t30 * mrSges(4,2);
t2 = t15 * pkin(4) - t16 * qJ(5) + t23;
t1 = [m(2) + m(3) * (t25 + t27) + m(4) * (t43 * t25 + t27) + t46 * (t10 ^ 2 + t12 ^ 2 + t27); (-mrSges(3,2) + t38) * t31 + (mrSges(3,1) - t18 + t47) * t33 + m(4) * (t43 * t31 * pkin(6) + t33 * pkin(2)) + m(5) * (-t23 * t33 + t41) + m(6) * (-t2 * t33 + t41) + t45 * (t10 * t16 + t12 * t15); Ifges(4,2) * t26 - 0.2e1 * pkin(2) * t18 + 0.2e1 * t2 * t3 + 0.2e1 * t23 * t4 + Ifges(3,3) + (Ifges(4,1) * t30 + 0.2e1 * Ifges(4,4) * t32) * t30 + 0.2e1 * pkin(6) * t38 + m(6) * (t2 ^ 2 + t42) + m(5) * (t23 ^ 2 + t42) + m(4) * (t43 * pkin(6) ^ 2 + pkin(2) ^ 2) + ((Ifges(6,1) + Ifges(5,1)) * t16 + 0.2e1 * t45 * t6) * t16 + ((Ifges(6,3) + Ifges(5,2)) * t15 - 0.2e1 * t45 * t8 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t16) * t15; t37 * t31 - (-mrSges(5,2) + mrSges(6,3)) * t12 + (-mrSges(5,1) - mrSges(6,1)) * t10 + m(6) * (t22 * t10 - t20 * t12) + (-t10 * t29 - t12 * t28) * t48; m(6) * (t20 * t8 + t22 * t6) - t8 * mrSges(5,2) - t6 * mrSges(5,1) - t6 * mrSges(6,1) + t8 * mrSges(6,3) + Ifges(4,5) * t30 + Ifges(4,6) * t32 + t37 * pkin(6) + (t22 * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t16 + (-t20 * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t15 + (m(5) * (t28 * t8 - t29 * t6) + (-t28 * t15 - t29 * t16) * mrSges(5,3)) * pkin(3); -0.2e1 * t22 * mrSges(6,1) + 0.2e1 * t20 * mrSges(6,3) + Ifges(6,2) + Ifges(4,3) + Ifges(5,3) + m(6) * (t20 ^ 2 + t22 ^ 2) + (0.2e1 * t29 * mrSges(5,1) - 0.2e1 * t28 * mrSges(5,2) + (t28 ^ 2 + t29 ^ 2) * t48) * pkin(3); -t46 * t33; m(5) * t23 + m(6) * t2 - t47; 0; t46; m(6) * t10; m(6) * t6 + t16 * mrSges(6,2); m(6) * t22 - mrSges(6,1); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
