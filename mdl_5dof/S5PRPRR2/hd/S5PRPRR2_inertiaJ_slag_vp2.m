% Calculate joint inertia matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:48
% EndTime: 2019-12-05 15:44:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (241->70), mult. (504->99), div. (0->0), fcn. (452->8), ass. (0->35)
t32 = cos(qJ(5));
t26 = t32 ^ 2;
t29 = sin(qJ(5));
t43 = t29 ^ 2 + t26;
t51 = mrSges(6,3) * t43;
t50 = m(4) * pkin(2);
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t15 = -t27 * t31 + t28 * t34;
t16 = t27 * t34 + t28 * t31;
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t6 = -t33 * t15 + t16 * t30;
t49 = t6 ^ 2;
t47 = pkin(2) * t27;
t22 = pkin(2) * t28 + pkin(3);
t12 = t22 * t33 - t30 * t47;
t46 = t12 * mrSges(5,1);
t13 = t30 * t22 + t33 * t47;
t45 = t13 * mrSges(5,2);
t44 = Ifges(6,5) * t29 + Ifges(6,6) * t32;
t42 = Ifges(6,2) * t26 + Ifges(5,3) + (Ifges(6,1) * t29 + 0.2e1 * Ifges(6,4) * t32) * t29;
t41 = t43 * pkin(7);
t11 = pkin(7) + t13;
t40 = t43 * t11;
t39 = -mrSges(6,1) * t29 - mrSges(6,2) * t32;
t38 = 0.2e1 * t51;
t19 = -mrSges(6,1) * t32 + mrSges(6,2) * t29;
t8 = t15 * t30 + t16 * t33;
t37 = (-mrSges(5,1) + t19) * t6 + (-mrSges(5,2) + t51) * t8;
t10 = -pkin(4) - t12;
t5 = t8 ^ 2;
t1 = [m(2) + m(6) * (t43 * t5 + t49) + m(5) * (t5 + t49) + m(4) * (t15 ^ 2 + t16 ^ 2) + m(3) * (t31 ^ 2 + t34 ^ 2); t34 * mrSges(3,1) + t15 * mrSges(4,1) - t31 * mrSges(3,2) - t16 * mrSges(4,2) + m(6) * (t10 * t6 + t8 * t40) + m(5) * (-t12 * t6 + t13 * t8) + (t15 * t28 + t16 * t27) * t50 + t37; 0.2e1 * t46 - 0.2e1 * t45 + 0.2e1 * t10 * t19 + Ifges(3,3) + Ifges(4,3) + t11 * t38 + m(6) * (t43 * t11 ^ 2 + t10 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + t42 + (0.2e1 * mrSges(4,1) * t28 - 0.2e1 * mrSges(4,2) * t27 + (t27 ^ 2 + t28 ^ 2) * t50) * pkin(2); 0; 0; m(6) * t43 + m(4) + m(5); m(6) * (-pkin(4) * t6 + t8 * t41) + t37; m(6) * (-pkin(4) * t10 + pkin(7) * t40) - t45 + t46 + (-pkin(4) + t10) * t19 + (t40 + t41) * mrSges(6,3) + t42; 0; -0.2e1 * pkin(4) * t19 + m(6) * (t43 * pkin(7) ^ 2 + pkin(4) ^ 2) + pkin(7) * t38 + t42; t39 * t8; t39 * t11 + t44; -t19; t39 * pkin(7) + t44; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
