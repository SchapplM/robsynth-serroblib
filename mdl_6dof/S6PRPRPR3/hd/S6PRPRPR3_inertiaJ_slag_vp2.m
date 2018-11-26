% Calculate joint inertia matrix for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:56:30
% EndTime: 2018-11-23 14:56:31
% DurationCPUTime: 0.70s
% Computational Cost: add. (519->172), mult. (1170->234), div. (0->0), fcn. (1058->10), ass. (0->78)
t98 = m(5) + m(6);
t52 = sin(qJ(6));
t55 = cos(qJ(6));
t72 = t52 ^ 2 + t55 ^ 2;
t68 = m(7) * t72;
t97 = m(6) + t68;
t53 = sin(qJ(4));
t44 = t53 ^ 2;
t56 = cos(qJ(4));
t46 = t56 ^ 2;
t96 = t44 + t46;
t95 = 0.2e1 * t96;
t94 = m(4) * pkin(2);
t91 = t72 * mrSges(7,3);
t75 = mrSges(6,1) + mrSges(5,3);
t90 = -m(6) * pkin(4) + mrSges(6,2);
t48 = sin(pkin(11));
t49 = sin(pkin(6));
t50 = cos(pkin(11));
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t14 = (t48 * t57 + t50 * t54) * t49;
t51 = cos(pkin(6));
t9 = t56 * t14 + t51 * t53;
t6 = t9 ^ 2;
t12 = (t48 * t54 - t50 * t57) * t49;
t11 = t12 ^ 2;
t89 = -t52 / 0.2e1;
t58 = -pkin(4) - pkin(9);
t5 = t53 * t9;
t88 = t56 * pkin(4);
t7 = t53 * t14 - t51 * t56;
t86 = t7 * t53;
t36 = t48 * pkin(2) + pkin(8);
t85 = pkin(5) + t36;
t84 = Ifges(7,4) * t52;
t83 = Ifges(7,4) * t55;
t78 = t56 * mrSges(7,3);
t23 = -t53 * mrSges(7,2) - t55 * t78;
t82 = t52 * t23;
t81 = t52 * t56;
t80 = t55 * t56;
t79 = t55 * t58;
t77 = t9 * qJ(5);
t76 = -mrSges(5,1) + mrSges(6,2);
t26 = t52 * mrSges(7,1) + t55 * mrSges(7,2);
t74 = t26 + mrSges(6,3);
t73 = t96 * t36 ^ 2;
t69 = -mrSges(5,2) + t74;
t37 = -t50 * pkin(2) - pkin(3);
t67 = t72 * t58;
t38 = t53 * qJ(5);
t66 = t37 - t38;
t1 = -t52 * t12 + t55 * t7;
t2 = t55 * t12 + t52 * t7;
t64 = t55 * t1 + t52 * t2;
t15 = t58 * t56 + t66;
t20 = t85 * t53;
t3 = -t52 * t15 + t55 * t20;
t4 = t55 * t15 + t52 * t20;
t63 = t55 * t3 + t52 * t4;
t62 = -Ifges(7,5) * t52 - Ifges(7,6) * t55;
t61 = (t9 * t56 + t86) * t36;
t59 = qJ(5) ^ 2;
t42 = t51 ^ 2;
t40 = Ifges(7,5) * t55;
t39 = Ifges(7,3) * t53;
t28 = Ifges(7,1) * t55 - t84;
t27 = -Ifges(7,2) * t52 + t83;
t25 = -t56 * mrSges(5,1) + t53 * mrSges(5,2);
t24 = t56 * mrSges(6,2) - t53 * mrSges(6,3);
t22 = t53 * mrSges(7,1) + t52 * t78;
t21 = t85 * t56;
t19 = -mrSges(7,1) * t80 + mrSges(7,2) * t81;
t18 = t66 - t88;
t17 = Ifges(7,5) * t53 + (-Ifges(7,1) * t52 - t83) * t56;
t16 = Ifges(7,6) * t53 + (-Ifges(7,2) * t55 - t84) * t56;
t8 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t6) + m(4) * (t14 ^ 2 + t11 + t42) + m(3) * (t42 + (t54 ^ 2 + t57 ^ 2) * t49 ^ 2) + t98 * (t7 ^ 2 + t11 + t6); -t14 * mrSges(4,2) + t1 * t22 + t2 * t23 + t75 * t86 + (t57 * mrSges(3,1) - t54 * mrSges(3,2)) * t49 + (t75 * t56 - t19) * t9 + (-mrSges(4,1) + t24 + t25) * t12 + m(7) * (t3 * t1 + t4 * t2 + t21 * t9) + m(5) * (t37 * t12 + t61) + m(6) * (t18 * t12 + t61) + (-t12 * t50 + t14 * t48) * t94; 0.2e1 * t18 * t24 - 0.2e1 * t21 * t19 + 0.2e1 * t3 * t22 + 0.2e1 * t4 * t23 + 0.2e1 * t37 * t25 + Ifges(3,3) + Ifges(4,3) + (t39 + (Ifges(6,2) + Ifges(5,1)) * t53) * t53 + m(7) * (t21 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(6) * (t18 ^ 2 + t73) + m(5) * (t37 ^ 2 + t73) + (-t55 * t16 - t52 * t17 + (Ifges(6,3) + Ifges(5,2)) * t56 + ((2 * Ifges(5,4)) + (2 * Ifges(6,6)) + t62) * t53) * t56 + t75 * t36 * t95 + (0.2e1 * t50 * mrSges(4,1) - 0.2e1 * t48 * mrSges(4,2) + (t48 ^ 2 + t50 ^ 2) * t94) * pkin(2); m(4) * t51 + m(7) * (-t56 * t64 + t5) + t98 * (-t56 * t7 + t5); -t53 * t19 + m(7) * (t53 * t21 - t56 * t63) - t23 * t81 - t22 * t80; m(4) + m(7) * (t72 * t46 + t44) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t95; t76 * t7 - t64 * mrSges(7,3) + t69 * t9 + m(7) * (t58 * t64 + t77) + m(6) * (-pkin(4) * t7 + t77); t21 * t26 + t55 * t17 / 0.2e1 + t16 * t89 - qJ(5) * t19 + m(7) * (qJ(5) * t21 + t58 * t63) + t58 * t82 + t22 * t79 - t63 * mrSges(7,3) + (Ifges(7,6) * t89 + t40 / 0.2e1 + Ifges(5,5) - Ifges(6,4) - pkin(4) * mrSges(6,1)) * t53 + (Ifges(5,6) - Ifges(6,5) - t55 * t27 / 0.2e1 + t28 * t89 + qJ(5) * mrSges(6,1)) * t56 + ((m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t56 + (-mrSges(5,1) + t90) * t53) * t36; (-t76 + t91) * t56 + m(6) * (t38 + t88) + m(7) * (-t56 * t67 + t38) + t69 * t53; -0.2e1 * pkin(4) * mrSges(6,2) - t52 * t27 + t55 * t28 + Ifges(6,1) + Ifges(5,3) + m(7) * (t72 * t58 ^ 2 + t59) + m(6) * (pkin(4) ^ 2 + t59) - 0.2e1 * mrSges(7,3) * t67 + 0.2e1 * t74 * qJ(5); m(6) * t7 + m(7) * t64; m(7) * t63 + t82 + t55 * t22 + (m(6) * t36 + mrSges(6,1)) * t53; -t97 * t56; t58 * t68 + t90 - t91; t97; t1 * mrSges(7,1) - t2 * mrSges(7,2); t3 * mrSges(7,1) - t4 * mrSges(7,2) + t56 * t62 + t39; t19; mrSges(7,1) * t79 + t40 + (-mrSges(7,2) * t58 - Ifges(7,6)) * t52; t55 * mrSges(7,1) - t52 * mrSges(7,2); Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
