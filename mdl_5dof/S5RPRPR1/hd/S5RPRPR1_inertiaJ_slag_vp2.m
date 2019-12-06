% Calculate joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:04
% EndTime: 2019-12-05 17:47:05
% DurationCPUTime: 0.35s
% Computational Cost: add. (504->107), mult. (853->152), div. (0->0), fcn. (852->6), ass. (0->47)
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t25 = -t38 * t43 - t39 * t41;
t27 = -t38 * t41 + t39 * t43;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t11 = t42 * t25 - t40 * t27;
t71 = t11 ^ 2;
t49 = t40 * t25 + t42 * t27;
t70 = t49 ^ 2 + t71;
t32 = t39 * pkin(3) + pkin(4);
t61 = pkin(3) * t38;
t17 = t42 * t32 - t40 * t61;
t18 = t40 * t32 + t42 * t61;
t69 = t18 * t11 - t17 * t49;
t68 = t43 ^ 2;
t66 = m(5) * pkin(3);
t64 = t25 ^ 2;
t63 = 0.2e1 * t11;
t60 = t17 * mrSges(6,1);
t59 = t18 * mrSges(6,2);
t44 = -pkin(1) - pkin(6);
t56 = -qJ(4) + t44;
t28 = t56 * t41;
t29 = t56 * t43;
t15 = t39 * t28 + t38 * t29;
t58 = t25 * t15;
t57 = t41 ^ 2 + t68;
t33 = t41 * pkin(3) + qJ(2);
t55 = t27 ^ 2 + t64;
t54 = m(4) * t57;
t53 = -t11 * mrSges(6,1) + mrSges(6,2) * t49;
t52 = mrSges(6,1) * t49 + mrSges(6,2) * t11;
t51 = t57 * mrSges(4,3);
t50 = -t25 * mrSges(5,1) + t27 * mrSges(5,2);
t14 = -t38 * t28 + t39 * t29;
t4 = -t27 * pkin(7) + t14;
t5 = t25 * pkin(7) + t15;
t2 = t42 * t4 - t40 * t5;
t3 = t40 * t4 + t42 * t5;
t48 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t49 + Ifges(6,6) * t11;
t47 = t25 * t38 - t27 * t39;
t45 = qJ(2) ^ 2;
t16 = -t25 * pkin(4) + t33;
t1 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * mrSges(5,3) * t58 + t3 * mrSges(6,3) * t63 + 0.2e1 * t16 * t53 + Ifges(6,2) * t71 + 0.2e1 * t33 * t50 + Ifges(5,2) * t64 + Ifges(4,1) * t68 - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t44 * t51 + (-0.2e1 * t14 * mrSges(5,3) + Ifges(5,1) * t27 + 0.2e1 * Ifges(5,4) * t25) * t27 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t49 + Ifges(6,4) * t63) * t49 + m(6) * (t16 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2 + t33 ^ 2) + m(4) * (t57 * t44 ^ 2 + t45) + m(3) * ((pkin(1) ^ 2) + t45) + (-0.2e1 * Ifges(4,4) * t43 + Ifges(4,2) * t41) * t41 + 0.2e1 * (t41 * mrSges(4,1) + t43 * mrSges(4,2) + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) - t70 * mrSges(6,3) - t55 * mrSges(5,3) - t51 + m(6) * (-t3 * t11 + t2 * t49) + m(5) * (t27 * t14 - t58) + t44 * t54; m(5) * t55 + m(6) * t70 + m(3) + t54; m(6) * (t17 * t2 + t18 * t3) - t15 * mrSges(5,2) + t14 * mrSges(5,1) + Ifges(5,5) * t27 + Ifges(5,6) * t25 + (t44 * mrSges(4,1) + Ifges(4,5)) * t43 + (-t44 * mrSges(4,2) - Ifges(4,6)) * t41 + t69 * mrSges(6,3) + (m(5) * (t14 * t39 + t15 * t38) + t47 * mrSges(5,3)) * pkin(3) + t48; -m(6) * t69 + t43 * mrSges(4,1) + t27 * mrSges(5,1) - t41 * mrSges(4,2) + t25 * mrSges(5,2) - t47 * t66 + t52; 0.2e1 * t60 - 0.2e1 * t59 + Ifges(4,3) + Ifges(5,3) + Ifges(6,3) + m(6) * (t17 ^ 2 + t18 ^ 2) + (0.2e1 * t39 * mrSges(5,1) - 0.2e1 * t38 * mrSges(5,2) + (t38 ^ 2 + t39 ^ 2) * t66) * pkin(3); m(5) * t33 + m(6) * t16 + t50 + t53; 0; 0; m(5) + m(6); t48; t52; Ifges(6,3) - t59 + t60; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
