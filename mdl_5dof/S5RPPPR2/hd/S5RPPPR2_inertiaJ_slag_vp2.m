% Calculate joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:58
% EndTime: 2019-12-05 17:31:01
% DurationCPUTime: 0.47s
% Computational Cost: add. (561->134), mult. (1218->210), div. (0->0), fcn. (1202->8), ass. (0->58)
t56 = cos(pkin(7));
t71 = 0.2e1 * t56;
t70 = 2 * qJ(2);
t51 = sin(pkin(9));
t54 = cos(pkin(9));
t53 = sin(pkin(7));
t55 = cos(pkin(8));
t64 = t55 * t53;
t28 = -t56 * t51 + t54 * t64;
t57 = sin(qJ(5));
t52 = sin(pkin(8));
t58 = cos(qJ(5));
t65 = t52 * t58;
t15 = -t28 * t57 + t53 * t65;
t66 = t52 * t57;
t16 = t28 * t58 + t53 * t66;
t69 = Ifges(6,5) * t16 + Ifges(6,6) * t15;
t47 = t55 ^ 2;
t67 = t52 * t53;
t20 = mrSges(5,1) * t67 - t28 * mrSges(5,3);
t6 = -t15 * mrSges(6,1) + t16 * mrSges(6,2);
t68 = t20 - t6;
t34 = -t56 * pkin(2) - t53 * qJ(3) - pkin(1);
t61 = qJ(2) * t56;
t23 = t52 * t34 + t55 * t61;
t17 = -t56 * qJ(4) + t23;
t24 = (pkin(3) * t52 - qJ(4) * t55 + qJ(2)) * t53;
t8 = t54 * t17 + t51 * t24;
t63 = mrSges(4,1) * t67 + mrSges(4,2) * t64;
t62 = t57 ^ 2 + t58 ^ 2;
t27 = t51 * t64 + t56 * t54;
t60 = Ifges(6,3) * t27 + t69;
t11 = t27 * mrSges(5,1) + t28 * mrSges(5,2);
t22 = t55 * t34 - t52 * t61;
t18 = t56 * pkin(3) - t22;
t7 = -t51 * t17 + t54 * t24;
t59 = qJ(2) ^ 2;
t48 = t56 ^ 2;
t46 = t54 ^ 2;
t45 = t53 ^ 2;
t44 = t52 ^ 2;
t43 = t51 ^ 2;
t41 = t53 * mrSges(3,2);
t40 = t45 * t59;
t39 = t43 * t44;
t33 = -t56 * mrSges(4,1) - mrSges(4,3) * t64;
t32 = t56 * mrSges(4,2) - mrSges(4,3) * t67;
t30 = t54 * t65 - t57 * t55;
t29 = -t54 * t66 - t58 * t55;
t19 = -mrSges(5,2) * t67 - t27 * mrSges(5,3);
t10 = t27 * mrSges(6,1) - t16 * mrSges(6,3);
t9 = -t27 * mrSges(6,2) + t15 * mrSges(6,3);
t5 = pkin(6) * t67 + t8;
t4 = -pkin(4) * t67 - t7;
t3 = t27 * pkin(4) - t28 * pkin(6) + t18;
t2 = t57 * t3 + t58 * t5;
t1 = t58 * t3 - t57 * t5;
t12 = [Ifges(5,1) * t28 ^ 2 + Ifges(6,1) * t16 ^ 2 - 0.2e1 * pkin(1) * t41 + 0.2e1 * t1 * t10 + 0.2e1 * t18 * t11 + 0.2e1 * t8 * t19 + 0.2e1 * t2 * t9 + 0.2e1 * t7 * t20 + 0.2e1 * t22 * t33 + 0.2e1 * t23 * t32 + 0.2e1 * t4 * t6 + Ifges(2,3) + (0.2e1 * Ifges(6,4) * t16 + Ifges(6,2) * t15) * t15 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3)) * t56) * t56 + (-0.2e1 * Ifges(5,4) * t28 + Ifges(5,2) * t27 + t60 + t69) * t27 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(5) * (t18 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2 + t40) + m(3) * (pkin(1) ^ 2 + t48 * t59 + t40) + (t45 + t48) * mrSges(3,3) * t70 + (t63 * t70 + (0.2e1 * Ifges(5,5) * t28 + Ifges(4,6) * t71 - 0.2e1 * Ifges(5,6) * t27) * t52 + (-Ifges(4,5) * t55 + Ifges(3,4)) * t71 + (Ifges(4,1) * t47 + Ifges(3,1) + (-0.2e1 * Ifges(4,4) * t55 + (Ifges(5,3) + Ifges(4,2)) * t52) * t52) * t53) * t53; -m(3) * pkin(1) - t56 * mrSges(3,1) + t29 * t10 + t30 * t9 + t41 + (-t11 + t33) * t55 + (t54 * t19 - t68 * t51 + t32) * t52 + m(6) * (t51 * t52 * t4 + t29 * t1 + t30 * t2) + m(5) * (-t55 * t18 + (-t51 * t7 + t54 * t8) * t52) + m(4) * (t55 * t22 + t52 * t23); m(3) + m(4) * (t44 + t47) + m(5) * (t46 * t44 + t39 + t47) + m(6) * (t29 ^ 2 + t30 ^ 2 + t39); m(4) * t53 * qJ(2) + t68 * t54 + (-t57 * t10 + t58 * t9 + t19) * t51 + m(6) * (-t54 * t4 + (-t1 * t57 + t2 * t58) * t51) + m(5) * (t51 * t8 + t54 * t7) + t63; m(6) * (-t29 * t57 + t30 * t58 - t52 * t54) * t51; m(4) + m(5) * (t43 + t46) + m(6) * (t62 * t43 + t46); t58 * t10 + t57 * t9 + m(6) * (t58 * t1 + t57 * t2) + m(5) * t18 + t11; -m(5) * t55 + m(6) * (t58 * t29 + t57 * t30); 0; m(6) * t62 + m(5); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t60; t29 * mrSges(6,1) - t30 * mrSges(6,2); (-mrSges(6,1) * t57 - mrSges(6,2) * t58) * t51; t58 * mrSges(6,1) - t57 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;
