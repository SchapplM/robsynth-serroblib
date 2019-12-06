% Calculate joint inertia matrix for
% S5RPPRR2
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:36
% DurationCPUTime: 0.35s
% Computational Cost: add. (451->90), mult. (756->128), div. (0->0), fcn. (774->6), ass. (0->41)
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t56 = sin(qJ(4));
t57 = cos(qJ(4));
t23 = -t57 * t37 - t56 * t38;
t25 = -t56 * t37 + t57 * t38;
t40 = sin(qJ(5));
t41 = cos(qJ(5));
t11 = t41 * t23 - t40 * t25;
t66 = t11 ^ 2;
t46 = t40 * t23 + t41 * t25;
t65 = t46 ^ 2 + t66;
t64 = t40 * t11 - t41 * t46;
t61 = t23 ^ 2;
t35 = t38 ^ 2;
t60 = 0.2e1 * t11;
t39 = -pkin(1) - qJ(3);
t58 = -pkin(6) + t39;
t26 = t58 * t37;
t27 = t58 * t38;
t15 = t57 * t26 + t56 * t27;
t55 = t23 * t15;
t54 = t37 * mrSges(4,1) + t38 * mrSges(4,2);
t53 = t37 ^ 2 + t35;
t28 = t37 * pkin(3) + qJ(2);
t52 = t25 ^ 2 + t61;
t51 = m(4) * t53;
t50 = -t11 * mrSges(6,1) + mrSges(6,2) * t46;
t49 = mrSges(6,1) * t46 + mrSges(6,2) * t11;
t48 = t53 * mrSges(4,3);
t47 = -t23 * mrSges(5,1) + t25 * mrSges(5,2);
t14 = -t56 * t26 + t57 * t27;
t4 = -t25 * pkin(7) + t14;
t5 = t23 * pkin(7) + t15;
t2 = t41 * t4 - t40 * t5;
t3 = t40 * t4 + t41 * t5;
t45 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t46 + Ifges(6,6) * t11;
t44 = (t41 * mrSges(6,1) - t40 * mrSges(6,2)) * pkin(4);
t42 = qJ(2) ^ 2;
t16 = -t23 * pkin(4) + t28;
t1 = [Ifges(3,1) + Ifges(2,3) + Ifges(6,2) * t66 + 0.2e1 * t16 * t50 + 0.2e1 * t28 * t47 + Ifges(5,2) * t61 - (2 * pkin(1) * mrSges(3,2)) + Ifges(4,1) * t35 + 0.2e1 * mrSges(5,3) * t55 + t3 * mrSges(6,3) * t60 + (-0.2e1 * Ifges(4,4) * t38 + Ifges(4,2) * t37) * t37 - 0.2e1 * t39 * t48 + (-0.2e1 * t14 * mrSges(5,3) + Ifges(5,1) * t25 + 0.2e1 * Ifges(5,4) * t23) * t25 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t46 + Ifges(6,4) * t60) * t46 + m(6) * (t16 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2 + t28 ^ 2) + m(4) * (t53 * t39 ^ 2 + t42) + m(3) * ((pkin(1) ^ 2) + t42) + 0.2e1 * (mrSges(3,3) + t54) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) - t65 * mrSges(6,3) - t52 * mrSges(5,3) - t48 + m(6) * (-t3 * t11 + t2 * t46) + m(5) * (t25 * t14 - t55) + t39 * t51; m(5) * t52 + m(6) * t65 + m(3) + t51; m(4) * qJ(2) + m(5) * t28 + m(6) * t16 + t47 + t50 + t54; 0; m(4) + m(5) + m(6); t14 * mrSges(5,1) - t15 * mrSges(5,2) + Ifges(5,5) * t25 + Ifges(5,6) * t23 + (m(6) * (t2 * t41 + t3 * t40) + t64 * mrSges(6,3)) * pkin(4) + t45; -m(6) * t64 * pkin(4) + t25 * mrSges(5,1) + t23 * mrSges(5,2) + t49; 0; Ifges(5,3) + Ifges(6,3) + m(6) * (t40 ^ 2 + t41 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t44; t45; t49; 0; Ifges(6,3) + t44; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
