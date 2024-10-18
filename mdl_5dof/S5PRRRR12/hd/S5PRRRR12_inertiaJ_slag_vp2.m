% Calculate joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:16
% EndTime: 2024-09-28 18:07:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (649->135), mult. (1613->202), div. (0->0), fcn. (1618->12), ass. (0->70)
t52 = cos(pkin(6));
t58 = cos(qJ(5));
t75 = t52 * t58;
t50 = sin(pkin(6));
t54 = sin(qJ(5));
t78 = t50 * t54;
t32 = pkin(4) * t75 - pkin(10) * t78;
t34 = t52 * mrSges(6,1) - mrSges(6,3) * t78;
t16 = t32 * t34;
t76 = t52 * t54;
t77 = t50 * t58;
t33 = pkin(4) * t76 + pkin(10) * t77;
t35 = -t52 * mrSges(6,2) + mrSges(6,3) * t77;
t17 = t33 * t35;
t86 = t16 + t17;
t46 = t50 * pkin(10);
t55 = sin(qJ(4));
t36 = t55 * pkin(3) + t46;
t59 = cos(qJ(4));
t82 = t59 * pkin(3);
t42 = pkin(4) + t82;
t20 = -t54 * t36 + t42 * t75;
t14 = t20 * t34;
t21 = t58 * t36 + t42 * t76;
t15 = t21 * t35;
t44 = mrSges(5,1) * t82;
t85 = t14 + t15 + t44;
t56 = sin(qJ(3));
t84 = pkin(2) * t56;
t51 = sin(pkin(5));
t57 = sin(qJ(2));
t60 = cos(qJ(3));
t61 = cos(qJ(2));
t25 = (-t56 * t57 + t60 * t61) * t51;
t26 = (t56 * t61 + t57 * t60) * t51;
t12 = t59 * t25 - t55 * t26;
t53 = cos(pkin(5));
t8 = -t50 * t12 + t53 * t52;
t83 = t50 * t8;
t43 = t60 * pkin(2) + pkin(3);
t30 = t55 * t43 + t59 * t84;
t81 = t30 * mrSges(5,2);
t47 = t50 ^ 2;
t80 = t42 * t47;
t31 = (-mrSges(6,1) * t58 + mrSges(6,2) * t54) * t50;
t79 = t50 * t31;
t74 = t55 * mrSges(5,2);
t73 = -0.2e1 * t79;
t72 = pkin(3) * t74;
t71 = Ifges(6,5) * t78 + Ifges(6,6) * t77 + Ifges(6,3) * t52;
t70 = Ifges(5,3) + (Ifges(6,5) * t52 + (Ifges(6,1) * t54 + Ifges(6,4) * t58) * t50) * t78 + (Ifges(6,6) * t52 + (Ifges(6,4) * t54 + Ifges(6,2) * t58) * t50) * t77 + t52 * t71;
t29 = t59 * t43 - t55 * t84;
t69 = Ifges(4,3) + t70;
t68 = t12 * t52 + t50 * t53;
t13 = t55 * t25 + t59 * t26;
t3 = -t54 * t13 + t68 * t58;
t4 = t58 * t13 + t68 * t54;
t67 = t12 * mrSges(5,1) - t13 * mrSges(5,2) + t3 * t34 + t8 * t31 + t4 * t35;
t66 = (t60 * mrSges(4,1) - t56 * mrSges(4,2)) * pkin(2);
t27 = t29 * mrSges(5,1);
t23 = t46 + t30;
t28 = pkin(4) + t29;
t9 = -t54 * t23 + t28 * t75;
t6 = t9 * t34;
t10 = t58 * t23 + t28 * t76;
t7 = t10 * t35;
t65 = t27 + t6 + t7 + t70 - t81;
t64 = t25 * mrSges(4,1) - t26 * mrSges(4,2) + t67;
t49 = t53 ^ 2;
t1 = [m(2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t49) + m(4) * (t25 ^ 2 + t26 ^ 2 + t49) + m(3) * (t49 + (t57 ^ 2 + t61 ^ 2) * t51 ^ 2); (t61 * mrSges(3,1) - t57 * mrSges(3,2)) * t51 + m(6) * (t10 * t4 - t28 * t83 + t9 * t3) + m(5) * (t29 * t12 + t30 * t13) + m(4) * (t25 * t60 + t26 * t56) * pkin(2) + t64; t28 * t73 - 0.2e1 * t81 + Ifges(3,3) + 0.2e1 * t27 + 0.2e1 * t6 + 0.2e1 * t7 + 0.2e1 * t66 + m(6) * (t47 * t28 ^ 2 + t10 ^ 2 + t9 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(4) * (t56 ^ 2 + t60 ^ 2) * pkin(2) ^ 2 + t69; m(6) * (t20 * t3 + t21 * t4 - t42 * t83) + m(5) * (t12 * t59 + t13 * t55) * pkin(3) + t64; Ifges(4,3) + t65 + (m(5) * (t29 * t59 + t30 * t55) - t74) * pkin(3) + m(6) * (t21 * t10 + t20 * t9 + t28 * t80) + t66 + (-t28 - t42) * t79 + t85; -0.2e1 * t72 + t42 * t73 + 0.2e1 * t14 + 0.2e1 * t15 + 0.2e1 * t44 + m(6) * (t47 * t42 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(5) * (t55 ^ 2 + t59 ^ 2) * pkin(3) ^ 2 + t69; m(6) * (-pkin(4) * t83 + t32 * t3 + t33 * t4) + t67; m(6) * (t47 * pkin(4) * t28 + t33 * t10 + t32 * t9) + (-pkin(4) - t28) * t79 + t65 + t86; m(6) * (pkin(4) * t80 + t32 * t20 + t33 * t21) - t72 + (-pkin(4) - t42) * t79 + t70 + t85 + t86; 0.2e1 * t17 + 0.2e1 * t16 + pkin(4) * t73 + m(6) * (t47 * pkin(4) ^ 2 + t32 ^ 2 + t33 ^ 2) + t70; t3 * mrSges(6,1) - t4 * mrSges(6,2); t9 * mrSges(6,1) - t10 * mrSges(6,2) + t71; t20 * mrSges(6,1) - t21 * mrSges(6,2) + t71; t32 * mrSges(6,1) - t33 * mrSges(6,2) + t71; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
