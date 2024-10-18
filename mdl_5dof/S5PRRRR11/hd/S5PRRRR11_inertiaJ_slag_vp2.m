% Calculate joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:12
% EndTime: 2024-09-27 21:45:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (750->147), mult. (1930->216), div. (0->0), fcn. (1952->8), ass. (0->66)
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t85 = t56 * Ifges(4,5) + t59 * Ifges(4,6);
t52 = sin(pkin(5));
t34 = (t59 * mrSges(4,1) - t56 * mrSges(4,2)) * t52;
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t29 = (-t55 * t56 + t58 * t59) * t52;
t30 = (t55 * t59 + t56 * t58) * t52;
t17 = -t29 * mrSges(5,1) + t30 * mrSges(5,2);
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t15 = t57 * t29 - t54 * t30;
t16 = t54 * t29 + t57 * t30;
t4 = t15 * mrSges(6,1) - t16 * mrSges(6,2);
t63 = -t17 + t4;
t84 = -t34 - t63;
t83 = t59 ^ 2;
t82 = Ifges(6,5) * t16 + Ifges(6,6) * t15;
t53 = cos(pkin(5));
t81 = Ifges(4,3) * t53 + t85 * t52;
t80 = t16 ^ 2;
t79 = t30 ^ 2;
t78 = pkin(2) * t53;
t77 = pkin(3) * t55;
t28 = Ifges(5,5) * t30;
t27 = Ifges(5,6) * t29;
t45 = t58 * pkin(3) + pkin(4);
t33 = t54 * t45 + t57 * t77;
t76 = t33 * mrSges(6,2);
t75 = t52 * t56;
t74 = t52 * t59;
t73 = t54 * mrSges(6,2);
t70 = Ifges(5,3) + Ifges(6,3);
t44 = t59 * t78;
t21 = t53 * pkin(3) + t44 + (-pkin(7) - pkin(8)) * t75;
t36 = pkin(7) * t74 + t56 * t78;
t25 = pkin(8) * t74 + t36;
t9 = t55 * t21 + t58 * t25;
t69 = pkin(4) * t73;
t68 = Ifges(6,3) * t53 + t82;
t8 = t58 * t21 - t55 * t25;
t32 = t57 * t45 - t54 * t77;
t31 = t32 * mrSges(6,1);
t67 = Ifges(6,3) + t31 - t76;
t39 = (-pkin(3) * t59 - pkin(2)) * t52;
t5 = t53 * pkin(4) - t30 * pkin(9) + t8;
t6 = t29 * pkin(9) + t9;
t2 = t57 * t5 - t54 * t6;
t3 = t54 * t5 + t57 * t6;
t65 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t68;
t64 = (t58 * mrSges(5,1) - t55 * mrSges(5,2)) * pkin(3);
t48 = Ifges(5,3) * t53;
t62 = t8 * mrSges(5,1) - t9 * mrSges(5,2) + t27 + t28 + t48 + t65;
t51 = t53 ^ 2;
t50 = t52 ^ 2;
t46 = t57 * pkin(4) * mrSges(6,1);
t38 = -t53 * mrSges(4,2) + mrSges(4,3) * t74;
t37 = t53 * mrSges(4,1) - mrSges(4,3) * t75;
t35 = -pkin(7) * t75 + t44;
t23 = t53 * mrSges(5,1) - t30 * mrSges(5,3);
t22 = -t53 * mrSges(5,2) + t29 * mrSges(5,3);
t18 = -t29 * pkin(4) + t39;
t11 = t53 * mrSges(6,1) - t16 * mrSges(6,3);
t10 = -t53 * mrSges(6,2) + t15 * mrSges(6,3);
t1 = [m(3) + m(2) + m(6) * (t15 ^ 2 + t51 + t80) + m(5) * (t29 ^ 2 + t51 + t79) + m(4) * (t51 + (t56 ^ 2 + t83) * t50); t16 * t10 + t15 * t11 + t30 * t22 + t29 * t23 + t84 * t53 + m(6) * (t2 * t15 + t3 * t16 + t18 * t53) + m(5) * (t8 * t29 + t9 * t30 + t39 * t53) + (t56 * t38 + t59 * t37 + m(4) * (t35 * t59 + t36 * t56 - t78)) * t52; Ifges(5,1) * t79 + Ifges(6,1) * t80 + 0.2e1 * t3 * t10 + 0.2e1 * t2 * t11 + 0.2e1 * t39 * t17 - 0.2e1 * t18 * t4 + 0.2e1 * t9 * t22 + 0.2e1 * t8 * t23 + 0.2e1 * t35 * t37 + 0.2e1 * t36 * t38 + Ifges(3,3) + (0.2e1 * Ifges(5,4) * t30 + Ifges(5,2) * t29) * t29 + (0.2e1 * Ifges(6,4) * t16 + Ifges(6,2) * t15) * t15 + (0.2e1 * t27 + 0.2e1 * t28 + t48 + t68 + t81 + t82) * t53 + m(6) * (t18 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t39 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(4) * (t50 * pkin(2) ^ 2 + t35 ^ 2 + t36 ^ 2) + (0.2e1 * pkin(2) * t34 + (Ifges(4,2) * t83 + (Ifges(4,1) * t56 + 0.2e1 * Ifges(4,4) * t59) * t56) * t52 + t85 * t53) * t52; m(6) * (t32 * t15 + t33 * t16) + m(5) * (t29 * t58 + t30 * t55) * pkin(3) - t84; t62 + (t55 * t22 + t58 * t23 + m(5) * (t55 * t9 + t58 * t8)) * pkin(3) + m(6) * (t32 * t2 + t33 * t3) + t32 * t11 + t33 * t10 + t35 * mrSges(4,1) - t36 * mrSges(4,2) + t81; -0.2e1 * t76 + Ifges(4,3) + 0.2e1 * t31 + 0.2e1 * t64 + m(6) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t55 ^ 2 + t58 ^ 2) * pkin(3) ^ 2 + t70; m(6) * (t15 * t57 + t16 * t54) * pkin(4) + t63; (t54 * t10 + t57 * t11 + m(6) * (t2 * t57 + t3 * t54)) * pkin(4) + t62; Ifges(5,3) + t46 + t64 + (m(6) * (t32 * t57 + t33 * t54) - t73) * pkin(4) + t67; -0.2e1 * t69 + 0.2e1 * t46 + m(6) * (t54 ^ 2 + t57 ^ 2) * pkin(4) ^ 2 + t70; t4; t65; t67; Ifges(6,3) + t46 - t69; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
