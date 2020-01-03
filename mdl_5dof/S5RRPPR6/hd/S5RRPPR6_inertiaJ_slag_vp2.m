% Calculate joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:31:38
% DurationCPUTime: 0.58s
% Computational Cost: add. (913->172), mult. (1769->261), div. (0->0), fcn. (1867->8), ass. (0->67)
t70 = cos(qJ(2));
t79 = -qJ(3) - pkin(6);
t52 = t79 * t70;
t64 = sin(pkin(8));
t66 = cos(pkin(8));
t68 = sin(qJ(2));
t74 = t79 * t68;
t32 = -t64 * t52 - t66 * t74;
t89 = t32 ^ 2;
t88 = 0.2e1 * t32;
t57 = -t70 * pkin(2) - pkin(1);
t87 = 0.2e1 * t57;
t65 = cos(pkin(9));
t85 = t65 / 0.2e1;
t54 = t64 * pkin(2) + qJ(4);
t84 = pkin(7) + t54;
t63 = sin(pkin(9));
t83 = Ifges(5,4) * t63;
t82 = Ifges(5,4) * t65;
t45 = t64 * t70 + t66 * t68;
t81 = t45 * t63;
t80 = t45 * t65;
t43 = t64 * t68 - t66 * t70;
t25 = t43 * pkin(3) - t45 * qJ(4) + t57;
t34 = -t66 * t52 + t64 * t74;
t9 = t63 * t25 + t65 * t34;
t24 = mrSges(5,1) * t81 + mrSges(5,2) * t80;
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t44 = -t67 * t63 + t69 * t65;
t46 = t69 * t63 + t67 * t65;
t78 = Ifges(6,5) * t46 + Ifges(6,6) * t44;
t77 = t63 ^ 2 + t65 ^ 2;
t76 = t68 ^ 2 + t70 ^ 2;
t18 = t46 * t45;
t19 = t44 * t45;
t75 = Ifges(6,5) * t19 - Ifges(6,6) * t18 + Ifges(6,3) * t43;
t56 = -t66 * pkin(2) - pkin(3);
t49 = -t65 * mrSges(5,1) + t63 * mrSges(5,2);
t7 = t18 * mrSges(6,1) + t19 * mrSges(6,2);
t8 = t65 * t25 - t63 * t34;
t73 = -t8 * t63 + t9 * t65;
t28 = -t44 * mrSges(6,1) + t46 * mrSges(6,2);
t51 = Ifges(5,1) * t63 + t82;
t50 = Ifges(5,2) * t65 + t83;
t48 = -t65 * pkin(4) + t56;
t39 = t45 * mrSges(4,2);
t38 = t84 * t65;
t37 = t84 * t63;
t30 = Ifges(6,1) * t46 + Ifges(6,4) * t44;
t29 = Ifges(6,4) * t46 + Ifges(6,2) * t44;
t27 = t43 * mrSges(5,1) - mrSges(5,3) * t80;
t26 = -t43 * mrSges(5,2) - mrSges(5,3) * t81;
t23 = -t67 * t37 + t69 * t38;
t22 = -t69 * t37 - t67 * t38;
t14 = pkin(4) * t81 + t32;
t13 = Ifges(5,5) * t43 + (Ifges(5,1) * t65 - t83) * t45;
t12 = Ifges(5,6) * t43 + (-Ifges(5,2) * t63 + t82) * t45;
t11 = t43 * mrSges(6,1) - t19 * mrSges(6,3);
t10 = -t43 * mrSges(6,2) - t18 * mrSges(6,3);
t6 = -pkin(7) * t81 + t9;
t5 = Ifges(6,1) * t19 - Ifges(6,4) * t18 + Ifges(6,5) * t43;
t4 = Ifges(6,4) * t19 - Ifges(6,2) * t18 + Ifges(6,6) * t43;
t3 = t43 * pkin(4) - pkin(7) * t80 + t8;
t2 = t67 * t3 + t69 * t6;
t1 = t69 * t3 - t67 * t6;
t15 = [-0.2e1 * pkin(1) * (-t70 * mrSges(3,1) + t68 * mrSges(3,2)) + t68 * (Ifges(3,1) * t68 + Ifges(3,4) * t70) + t70 * (Ifges(3,4) * t68 + Ifges(3,2) * t70) + t39 * t87 + t19 * t5 + 0.2e1 * t9 * t26 + 0.2e1 * t8 * t27 + t24 * t88 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + 0.2e1 * t14 * t7 - t18 * t4 + Ifges(2,3) + 0.2e1 * t76 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t87 - 0.2e1 * t34 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t43 + t75) * t43 + m(6) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2 + t89) + m(4) * (t34 ^ 2 + t57 ^ 2 + t89) + m(3) * (t76 * pkin(6) ^ 2 + pkin(1) ^ 2) + (mrSges(4,3) * t88 + Ifges(4,1) * t45 - t63 * t12 + t65 * t13 + (Ifges(5,5) * t65 - Ifges(5,6) * t63 - (2 * Ifges(4,4))) * t43) * t45; (m(4) * (-t32 * t66 + t34 * t64) + (-t64 * t43 - t66 * t45) * mrSges(4,3)) * pkin(2) + (Ifges(5,5) * t63 + Ifges(5,6) * t65 + t78) * t43 / 0.2e1 + (t65 * t26 - t63 * t27) * t54 + (t49 - mrSges(4,1)) * t32 + m(6) * (t22 * t1 + t48 * t14 + t23 * t2) + (-t68 * mrSges(3,1) - t70 * mrSges(3,2)) * pkin(6) + (-t1 * t46 + t2 * t44) * mrSges(6,3) + Ifges(3,6) * t70 + t63 * t13 / 0.2e1 + Ifges(3,5) * t68 + t56 * t24 - Ifges(4,6) * t43 + t44 * t4 / 0.2e1 + t46 * t5 / 0.2e1 + t48 * t7 - t34 * mrSges(4,2) + t22 * t11 + t23 * t10 + t14 * t28 - t18 * t29 / 0.2e1 + t19 * t30 / 0.2e1 + m(5) * (t56 * t32 + t54 * t73) + t73 * mrSges(5,3) + t12 * t85 + (Ifges(4,5) - t63 * t50 / 0.2e1 + t51 * t85) * t45; 0.2e1 * t48 * t28 + t44 * t29 + t46 * t30 + 0.2e1 * t56 * t49 + t65 * t50 + t63 * t51 + Ifges(3,3) + Ifges(4,3) + m(6) * (t22 ^ 2 + t23 ^ 2 + t48 ^ 2) + m(5) * (t77 * t54 ^ 2 + t56 ^ 2) + m(4) * (t64 ^ 2 + t66 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t66 * mrSges(4,1) - t64 * mrSges(4,2)) * pkin(2) + 0.2e1 * (-t22 * t46 + t23 * t44) * mrSges(6,3) + 0.2e1 * t77 * t54 * mrSges(5,3); t43 * mrSges(4,1) + t46 * t10 + t44 * t11 + t63 * t26 + t65 * t27 + t39 + m(6) * (t44 * t1 + t46 * t2) + m(5) * (t63 * t9 + t65 * t8) + m(4) * t57; m(6) * (t44 * t22 + t46 * t23); m(4) + m(5) * t77 + m(6) * (t44 ^ 2 + t46 ^ 2); m(5) * t32 + m(6) * t14 + t24 + t7; m(5) * t56 + m(6) * t48 + t28 + t49; 0; m(5) + m(6); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t75; t22 * mrSges(6,1) - t23 * mrSges(6,2) + t78; -t28; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
