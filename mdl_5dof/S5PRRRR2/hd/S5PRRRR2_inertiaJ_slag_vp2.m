% Calculate joint inertia matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:12
% EndTime: 2019-07-18 13:30:12
% DurationCPUTime: 0.19s
% Computational Cost: add. (198->67), mult. (410->86), div. (0->0), fcn. (245->6), ass. (0->46)
t31 = cos(qJ(5));
t27 = t31 ^ 2;
t28 = sin(qJ(5));
t26 = t28 ^ 2;
t46 = t26 + t27;
t43 = m(6) * t46;
t50 = mrSges(6,3) * t26;
t20 = pkin(6) * t50;
t49 = mrSges(6,3) * t27;
t21 = pkin(6) * t49;
t57 = t20 + t21;
t29 = sin(qJ(4));
t53 = pkin(3) * t29;
t19 = pkin(6) + t53;
t12 = t19 * t50;
t13 = t19 * t49;
t32 = cos(qJ(4));
t22 = t32 * pkin(3) * mrSges(5,1);
t56 = t12 + t13 + t22;
t33 = cos(qJ(3));
t44 = pkin(2) * t33 + pkin(3);
t30 = sin(qJ(3));
t54 = pkin(2) * t30;
t7 = t29 * t54 - t32 * t44;
t55 = t7 ^ 2;
t52 = t32 * t7;
t9 = t29 * t44 + t32 * t54;
t51 = t9 * mrSges(5,2);
t14 = -mrSges(6,1) * t31 + mrSges(6,2) * t28;
t48 = t32 * t14;
t47 = Ifges(6,5) * t28 + Ifges(6,6) * t31;
t45 = Ifges(6,2) * t27 + Ifges(5,3) + (Ifges(6,1) * t28 + 0.2e1 * Ifges(6,4) * t31) * t28;
t42 = Ifges(4,3) + t45;
t41 = pkin(6) * t43;
t40 = -mrSges(6,1) * t28 - mrSges(6,2) * t31;
t1 = t7 * t14;
t6 = pkin(6) + t9;
t2 = t6 * t50;
t3 = t6 * t49;
t4 = t7 * mrSges(5,1);
t39 = t1 + t2 + t3 - t4 + t45;
t38 = (mrSges(4,1) * t33 - mrSges(4,2) * t30) * pkin(2);
t37 = (-t29 * mrSges(5,2) - t48) * pkin(3);
t35 = pkin(3) ^ 2;
t25 = t32 ^ 2 * t35;
t5 = [m(2) + m(3) + m(4) + m(5) + t43; 0; -0.2e1 * t51 + Ifges(3,3) + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 + 0.2e1 * t38 + m(6) * (t46 * t6 ^ 2 + t55) + m(5) * (t9 ^ 2 + t55) + m(4) * (t30 ^ 2 + t33 ^ 2) * pkin(2) ^ 2 + t42; 0; Ifges(4,3) + t38 + (-t9 - t53) * mrSges(5,2) + t6 * t19 * t43 + t39 + (-t48 - m(6) * t52 + m(5) * (t29 * t9 - t52)) * pkin(3) + t56; 0.2e1 * t12 + 0.2e1 * t13 + 0.2e1 * t22 + 0.2e1 * t37 + m(6) * (t46 * t19 ^ 2 + t25) + m(5) * (t29 ^ 2 * t35 + t25) + t42; 0; t6 * t41 + t39 - t51 + t57; t19 * t41 + t37 + t45 + t56 + t57; pkin(6) ^ 2 * t43 + 0.2e1 * t20 + 0.2e1 * t21 + t45; -t14; t40 * t6 + t47; t40 * t19 + t47; t40 * pkin(6) + t47; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq  = res;
