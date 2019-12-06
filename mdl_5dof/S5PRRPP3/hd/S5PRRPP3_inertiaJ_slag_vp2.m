% Calculate joint inertia matrix for
% S5PRRPP3
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:02
% DurationCPUTime: 0.47s
% Computational Cost: add. (308->153), mult. (724->201), div. (0->0), fcn. (548->6), ass. (0->70)
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t88 = t49 ^ 2 + t50 ^ 2;
t54 = cos(qJ(2));
t52 = sin(qJ(2));
t53 = cos(qJ(3));
t68 = t52 * t53;
t12 = t49 * t68 + t54 * t50;
t14 = -t54 * t49 + t50 * t68;
t87 = t12 * t49 + t14 * t50;
t86 = 2 * pkin(6);
t83 = m(5) + m(6);
t85 = mrSges(5,3) + mrSges(6,2);
t47 = t53 ^ 2;
t84 = m(5) * pkin(3);
t22 = -t50 * pkin(4) - t49 * qJ(5) - pkin(3);
t82 = m(6) * t22;
t81 = pkin(6) * t52;
t80 = pkin(6) * t53;
t79 = Ifges(5,4) * t49;
t78 = Ifges(5,4) * t50;
t77 = Ifges(6,5) * t49;
t76 = Ifges(6,5) * t50;
t73 = t49 * mrSges(6,3);
t51 = sin(qJ(3));
t72 = t49 * t51;
t23 = -t53 * pkin(3) - t51 * qJ(4) - pkin(2);
t71 = t50 * t23;
t70 = t50 * t51;
t69 = t51 * t52;
t18 = t53 * mrSges(5,2) - mrSges(5,3) * t72;
t21 = -mrSges(6,2) * t72 - t53 * mrSges(6,3);
t67 = t18 + t21;
t19 = -t53 * mrSges(5,1) - mrSges(5,3) * t70;
t20 = t53 * mrSges(6,1) + mrSges(6,2) * t70;
t66 = -t19 + t20;
t40 = t49 * mrSges(5,2);
t25 = -t50 * mrSges(5,1) + t40;
t65 = t25 - mrSges(4,1);
t16 = mrSges(5,1) * t72 + mrSges(5,2) * t70;
t4 = t49 * t23 + t50 * t80;
t64 = t88 * qJ(4) ^ 2;
t45 = t51 ^ 2;
t63 = (t45 + t47) * mrSges(4,3);
t60 = t87 * qJ(4);
t1 = -t53 * qJ(5) + t4;
t2 = -t71 + (pkin(6) * t49 + pkin(4)) * t53;
t59 = t1 * t50 + t2 * t49;
t3 = -t49 * t80 + t71;
t58 = -t3 * t49 + t4 * t50;
t56 = pkin(6) ^ 2;
t48 = t54 ^ 2;
t46 = t52 ^ 2;
t42 = t45 * t56;
t37 = t45 * t46;
t36 = t45 * t81;
t31 = mrSges(6,1) * t72;
t30 = -t53 * mrSges(4,1) + t51 * mrSges(4,2);
t29 = Ifges(5,1) * t49 + t78;
t28 = Ifges(6,1) * t49 - t76;
t27 = Ifges(5,2) * t50 + t79;
t26 = -Ifges(6,3) * t50 + t77;
t24 = -t50 * mrSges(6,1) - t73;
t15 = -mrSges(6,3) * t70 + t31;
t10 = (pkin(4) * t49 - qJ(5) * t50 + pkin(6)) * t51;
t9 = -Ifges(5,5) * t53 + (Ifges(5,1) * t50 - t79) * t51;
t8 = -Ifges(6,4) * t53 + (Ifges(6,1) * t50 + t77) * t51;
t7 = -Ifges(5,6) * t53 + (-Ifges(5,2) * t49 + t78) * t51;
t6 = -Ifges(6,6) * t53 + (Ifges(6,3) * t49 + t76) * t51;
t5 = [m(2) + m(3) * (t46 + t48) + m(4) * (t47 * t46 + t37 + t48) + t83 * (t12 ^ 2 + t14 ^ 2 + t37); (mrSges(3,1) - t30) * t54 + t67 * t14 + t66 * t12 + (-mrSges(3,2) + (t15 + t16) * t51 + t63) * t52 + m(4) * (t54 * pkin(2) + t47 * t81 + t36) + m(5) * (-t3 * t12 + t4 * t14 + t36) + m(6) * (t1 * t14 + t10 * t69 + t2 * t12); -0.2e1 * pkin(2) * t30 + 0.2e1 * t1 * t21 + 0.2e1 * t10 * t15 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + Ifges(3,3) + t63 * t86 + m(5) * (t3 ^ 2 + t4 ^ 2 + t42) + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(4) * (pkin(2) ^ 2 + t47 * t56 + t42) + (Ifges(6,2) + Ifges(5,3) + Ifges(4,2)) * t47 + (Ifges(4,1) * t51 + t16 * t86 + (t8 + t9) * t50 + (t6 - t7) * t49 + ((2 * Ifges(4,4)) + (-Ifges(6,4) - Ifges(5,5)) * t50 + (Ifges(5,6) - Ifges(6,6)) * t49) * t53) * t51; (-t53 * mrSges(4,2) + (t24 + t65) * t51) * t52 + m(5) * (-pkin(3) * t69 + t60) + m(6) * (t22 * t69 + t60) + t85 * t87; -pkin(3) * t16 + t22 * t15 + (-t6 / 0.2e1 + t7 / 0.2e1) * t50 + (t8 / 0.2e1 + t9 / 0.2e1) * t49 + (-pkin(6) * mrSges(4,2) + Ifges(4,6) + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t50 + (-Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t49) * t53 + (t24 + t82) * t10 + t58 * mrSges(5,3) + t59 * mrSges(6,2) + (m(5) * t58 + m(6) * t59 + t66 * t49 + t67 * t50) * qJ(4) + (Ifges(4,5) + (t28 / 0.2e1 + t29 / 0.2e1) * t50 + (t26 / 0.2e1 - t27 / 0.2e1) * t49 + (t65 - t84) * pkin(6)) * t51; -0.2e1 * pkin(3) * t25 + 0.2e1 * t22 * t24 + Ifges(4,3) + (-t26 + t27) * t50 + (t28 + t29) * t49 + m(5) * (pkin(3) ^ 2 + t64) + m(6) * (t22 ^ 2 + t64) + 0.2e1 * t85 * qJ(4) * t88; t83 * t69; m(6) * t10 + t31 + (m(5) * pkin(6) - t50 * mrSges(6,3)) * t51 + t16; t40 - t84 - t73 + t82 + (-mrSges(5,1) - mrSges(6,1)) * t50; t83; m(6) * t12; m(6) * t2 + t20; (m(6) * qJ(4) + mrSges(6,2)) * t49; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
