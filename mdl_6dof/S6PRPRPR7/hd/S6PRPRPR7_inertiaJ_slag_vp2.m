% Calculate joint inertia matrix for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2018-11-23 14:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:58:52
% EndTime: 2018-11-23 14:58:53
% DurationCPUTime: 0.58s
% Computational Cost: add. (438->166), mult. (910->218), div. (0->0), fcn. (706->8), ass. (0->77)
t98 = m(6) + m(5);
t48 = sin(qJ(6));
t51 = cos(qJ(6));
t73 = t48 ^ 2 + t51 ^ 2;
t67 = m(7) * t73;
t97 = m(6) + t67;
t49 = sin(qJ(4));
t41 = t49 ^ 2;
t52 = cos(qJ(4));
t43 = t52 ^ 2;
t72 = -t43 - t41;
t47 = cos(pkin(6));
t46 = sin(pkin(6));
t53 = cos(qJ(2));
t84 = t46 * t53;
t12 = t47 * t52 - t49 * t84;
t5 = t49 * t12;
t10 = t47 * t49 + t52 * t84;
t79 = t52 * t10;
t61 = t5 - t79;
t95 = mrSges(7,3) * t73;
t94 = mrSges(5,1) - mrSges(6,2);
t77 = -mrSges(5,3) - mrSges(6,1);
t35 = qJ(5) * t49;
t93 = m(6) * (pkin(4) * t52 + t35);
t92 = t77 * t72;
t9 = t12 ^ 2;
t91 = -t48 / 0.2e1;
t90 = t51 / 0.2e1;
t55 = -pkin(2) - pkin(8);
t88 = pkin(5) - t55;
t87 = Ifges(7,4) * t48;
t86 = Ifges(7,4) * t51;
t50 = sin(qJ(2));
t85 = t46 * t50;
t81 = t49 * t51;
t18 = -t52 * mrSges(7,2) + mrSges(7,3) * t81;
t83 = t48 * t18;
t82 = t48 * t49;
t54 = -pkin(4) - pkin(9);
t80 = t51 * t54;
t22 = t48 * mrSges(7,1) + t51 * mrSges(7,2);
t76 = t22 + mrSges(6,3);
t75 = t49 * mrSges(5,1) + t52 * mrSges(5,2) + mrSges(4,3);
t74 = t72 * t55 ^ 2;
t71 = t12 * qJ(5);
t69 = -mrSges(5,2) + t76;
t68 = Ifges(7,5) * t82 + Ifges(7,6) * t81 + Ifges(7,3) * t52;
t66 = t73 * t54;
t19 = t49 * pkin(4) - t52 * qJ(5) + qJ(3);
t14 = t49 * pkin(9) + t19;
t21 = t88 * t52;
t1 = -t48 * t14 + t51 * t21;
t2 = t51 * t14 + t48 * t21;
t64 = t51 * t1 + t48 * t2;
t3 = t51 * t10 - t48 * t85;
t4 = t48 * t10 + t51 * t85;
t63 = t51 * t3 + t48 * t4;
t62 = t51 * mrSges(7,1) - t48 * mrSges(7,2);
t17 = t52 * mrSges(7,1) - mrSges(7,3) * t82;
t60 = -t51 * t17 - t83;
t59 = t61 * t55;
t58 = -0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t72;
t57 = qJ(3) ^ 2;
t56 = qJ(5) ^ 2;
t39 = t46 ^ 2;
t37 = Ifges(7,5) * t51;
t32 = t39 * t50 ^ 2;
t28 = qJ(3) * t85;
t26 = Ifges(7,1) * t51 - t87;
t25 = -Ifges(7,2) * t48 + t86;
t23 = -t49 * mrSges(6,2) - t52 * mrSges(6,3);
t20 = t88 * t49;
t15 = t62 * t49;
t7 = Ifges(7,5) * t52 + (Ifges(7,1) * t48 + t86) * t49;
t6 = Ifges(7,6) * t52 + (Ifges(7,2) * t51 + t87) * t49;
t8 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t39 * t53 ^ 2 + t47 ^ 2 + t32) + t98 * (t10 ^ 2 + t32 + t9); t3 * t17 + t4 * t18 - t77 * t79 + (t49 * t77 - t15) * t12 + ((mrSges(3,1) - mrSges(4,2)) * t53 + (-mrSges(3,2) + t23 + t75) * t50) * t46 + m(6) * (t19 * t85 + t59) + m(7) * (t1 * t3 - t20 * t12 + t2 * t4) + m(5) * (t28 + t59) + m(4) * (pkin(2) * t84 + t28); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t1 * t17 + 0.2e1 * t20 * t15 + 0.2e1 * t2 * t18 + 0.2e1 * t19 * t23 + Ifges(4,1) + Ifges(3,3) + ((Ifges(5,1) + Ifges(6,2)) * t52 + t68) * t52 + m(7) * (t1 ^ 2 + t2 ^ 2 + t20 ^ 2) + m(6) * (t19 ^ 2 - t74) + m(5) * (t57 - t74) + m(4) * (pkin(2) ^ 2 + t57) + 0.2e1 * t75 * qJ(3) - 0.2e1 * t55 * t92 + (t48 * t7 + t51 * t6 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t52 + (Ifges(5,2) + Ifges(6,3)) * t49) * t49; -m(4) * t84 + m(7) * (-t52 * t63 + t5) + t98 * t61; -m(4) * pkin(2) - t49 * t15 + mrSges(4,2) + t60 * t52 + m(7) * (-t49 * t20 - t52 * t64) + t55 * t58 - t92; m(4) + m(7) * (t43 * t73 + t41) + t58; -t94 * t10 - t63 * mrSges(7,3) + t69 * t12 + m(6) * (-pkin(4) * t10 + t71) + m(7) * (t54 * t63 + t71); t54 * t83 + t17 * t80 + m(7) * (-qJ(5) * t20 + t54 * t64) - t20 * t22 + t7 * t90 + t6 * t91 - qJ(5) * t15 - t64 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) + Ifges(7,6) * t91 + t37 / 0.2e1 - Ifges(6,4) + Ifges(5,5)) * t52 + (t48 * t26 / 0.2e1 + t25 * t90 - qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6)) * t49 + (t93 + t94 * t52 + (mrSges(6,3) - mrSges(5,2)) * t49) * t55; (t94 + t95) * t52 + t93 + m(7) * (-t52 * t66 + t35) + t69 * t49; -0.2e1 * pkin(4) * mrSges(6,2) - t48 * t25 + t51 * t26 + Ifges(6,1) + Ifges(5,3) + m(7) * (t54 ^ 2 * t73 + t56) + m(6) * (pkin(4) ^ 2 + t56) - 0.2e1 * mrSges(7,3) * t66 + 0.2e1 * t76 * qJ(5); m(6) * t10 + m(7) * t63; m(7) * t64 + (-m(6) * t55 + mrSges(6,1)) * t52 - t60; -t97 * t52; -m(6) * pkin(4) + t54 * t67 + mrSges(6,2) - t95; t97; t3 * mrSges(7,1) - t4 * mrSges(7,2); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t68; -t62 * t52; mrSges(7,1) * t80 + t37 + (-mrSges(7,2) * t54 - Ifges(7,6)) * t48; t62; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
