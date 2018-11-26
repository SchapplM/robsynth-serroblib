% Calculate joint inertia matrix for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2018-11-23 15:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:38
% EndTime: 2018-11-23 15:39:38
% DurationCPUTime: 0.53s
% Computational Cost: add. (1261->197), mult. (2383->289), div. (0->0), fcn. (2547->10), ass. (0->80)
t70 = cos(pkin(10));
t68 = sin(pkin(9));
t56 = t68 * pkin(1) + qJ(3);
t96 = pkin(7) + t56;
t37 = t96 * t70;
t73 = sin(qJ(4));
t67 = sin(pkin(10));
t94 = cos(qJ(4));
t81 = t94 * t67;
t23 = t73 * t37 + t96 * t81;
t102 = t23 ^ 2;
t88 = t73 * t67;
t45 = -t94 * t70 + t88;
t42 = t45 ^ 2;
t65 = t70 ^ 2;
t101 = 0.2e1 * t23;
t71 = cos(pkin(9));
t58 = -t71 * pkin(1) - pkin(2);
t49 = -t70 * pkin(3) + t58;
t100 = 0.2e1 * t49;
t69 = cos(pkin(11));
t98 = t69 / 0.2e1;
t97 = pkin(4) * t45;
t48 = t73 * t70 + t81;
t89 = t48 * t69;
t66 = sin(pkin(11));
t90 = t48 * t66;
t26 = mrSges(6,1) * t90 + mrSges(6,2) * t89;
t72 = sin(qJ(6));
t74 = cos(qJ(6));
t47 = t74 * t66 + t72 * t69;
t21 = t47 * t48;
t44 = -t72 * t66 + t74 * t69;
t22 = t44 * t48;
t9 = t21 * mrSges(7,1) + t22 * mrSges(7,2);
t95 = t26 + t9;
t93 = Ifges(6,4) * t66;
t92 = Ifges(6,4) * t69;
t91 = t23 * t45;
t87 = pkin(8) + qJ(5);
t20 = -t48 * qJ(5) + t49 + t97;
t25 = t94 * t37 - t96 * t88;
t8 = t66 * t20 + t69 * t25;
t86 = Ifges(7,5) * t47 + Ifges(7,6) * t44;
t51 = -t69 * mrSges(6,1) + t66 * mrSges(6,2);
t85 = t51 - mrSges(5,1);
t84 = t66 ^ 2 + t69 ^ 2;
t83 = t67 ^ 2 + t65;
t82 = Ifges(7,5) * t22 - Ifges(7,6) * t21 + Ifges(7,3) * t45;
t80 = -t70 * mrSges(4,1) + t67 * mrSges(4,2);
t7 = t69 * t20 - t66 * t25;
t79 = qJ(5) * t84;
t78 = -t66 * t7 + t69 * t8;
t29 = -t44 * mrSges(7,1) + t47 * mrSges(7,2);
t27 = -t45 * mrSges(6,2) - mrSges(6,3) * t90;
t28 = t45 * mrSges(6,1) - mrSges(6,3) * t89;
t77 = t69 * t27 - t66 * t28;
t59 = -t69 * pkin(5) - pkin(4);
t54 = Ifges(6,1) * t66 + t92;
t53 = Ifges(6,2) * t69 + t93;
t52 = t87 * t69;
t50 = t87 * t66;
t43 = t48 ^ 2;
t38 = t48 * mrSges(5,2);
t33 = -t72 * t50 + t74 * t52;
t32 = -t74 * t50 - t72 * t52;
t31 = Ifges(7,1) * t47 + Ifges(7,4) * t44;
t30 = Ifges(7,4) * t47 + Ifges(7,2) * t44;
t14 = Ifges(6,5) * t45 + (Ifges(6,1) * t69 - t93) * t48;
t13 = Ifges(6,6) * t45 + (-Ifges(6,2) * t66 + t92) * t48;
t12 = pkin(5) * t90 + t23;
t11 = t45 * mrSges(7,1) - t22 * mrSges(7,3);
t10 = -t45 * mrSges(7,2) - t21 * mrSges(7,3);
t6 = Ifges(7,1) * t22 - Ifges(7,4) * t21 + Ifges(7,5) * t45;
t5 = Ifges(7,4) * t22 - Ifges(7,2) * t21 + Ifges(7,6) * t45;
t4 = -pkin(8) * t90 + t8;
t3 = t45 * pkin(5) - pkin(8) * t89 + t7;
t2 = t72 * t3 + t74 * t4;
t1 = t74 * t3 - t72 * t4;
t15 = [Ifges(4,2) * t65 + 0.2e1 * t58 * t80 + t38 * t100 + t26 * t101 + 0.2e1 * t8 * t27 + 0.2e1 * t7 * t28 - t21 * t5 + t22 * t6 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + 0.2e1 * t12 * t9 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t67 + 0.2e1 * Ifges(4,4) * t70) * t67 + (mrSges(5,1) * t100 - 0.2e1 * t25 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t45 + t82) * t45 + (mrSges(5,3) * t101 + Ifges(5,1) * t48 - t66 * t13 + t69 * t14 + (Ifges(6,5) * t69 - Ifges(6,6) * t66 - (2 * Ifges(5,4))) * t45) * t48 + m(5) * (t25 ^ 2 + t49 ^ 2 + t102) + m(4) * (t83 * t56 ^ 2 + t58 ^ 2) + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) + m(6) * (t7 ^ 2 + t8 ^ 2 + t102) + m(3) * (t68 ^ 2 + t71 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t71 * mrSges(3,1) - t68 * mrSges(3,2)) * pkin(1) + 0.2e1 * t83 * t56 * mrSges(4,3); t22 * t10 - t21 * t11 + t77 * t48 + t95 * t45 + m(6) * (t78 * t48 + t91) + m(7) * (-t1 * t21 + t12 * t45 + t2 * t22) + m(5) * (t25 * t48 + t91); m(3) + m(6) * (t84 * t43 + t42) + m(7) * (t21 ^ 2 + t22 ^ 2 + t42) + m(5) * (t43 + t42) + m(4) * t83; t45 * mrSges(5,1) + t47 * t10 + t44 * t11 + t66 * t27 + t69 * t28 + t38 + m(7) * (t44 * t1 + t47 * t2) + m(6) * (t66 * t8 + t69 * t7) + m(5) * t49 + m(4) * t58 + t80; m(7) * (-t44 * t21 + t47 * t22); m(4) + m(5) + m(6) * t84 + m(7) * (t44 ^ 2 + t47 ^ 2); t66 * t14 / 0.2e1 + t13 * t98 + t59 * t9 - Ifges(5,6) * t45 + t47 * t6 / 0.2e1 + t44 * t5 / 0.2e1 - t25 * mrSges(5,2) - pkin(4) * t26 + t12 * t29 - t21 * t30 / 0.2e1 + t22 * t31 / 0.2e1 + t32 * t11 + t33 * t10 + t85 * t23 + t77 * qJ(5) + (-t1 * t47 + t2 * t44) * mrSges(7,3) + t78 * mrSges(6,3) + m(6) * (-pkin(4) * t23 + t78 * qJ(5)) + m(7) * (t32 * t1 + t59 * t12 + t33 * t2) + (-t66 * t53 / 0.2e1 + t54 * t98 + Ifges(5,5)) * t48 + (Ifges(6,5) * t66 + Ifges(6,6) * t69 + t86) * t45 / 0.2e1; -t38 + (t21 * t47 + t22 * t44) * mrSges(7,3) + t84 * t48 * mrSges(6,3) + (t29 + t85) * t45 + m(6) * (t48 * t79 - t97) + m(7) * (-t32 * t21 + t33 * t22 + t59 * t45); m(7) * (t32 * t44 + t33 * t47); -0.2e1 * pkin(4) * t51 + 0.2e1 * t59 * t29 + t44 * t30 + t47 * t31 + t69 * t53 + t66 * t54 + Ifges(5,3) + m(7) * (t32 ^ 2 + t33 ^ 2 + t59 ^ 2) + m(6) * (t84 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t32 * t47 + t33 * t44) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t79; m(6) * t23 + m(7) * t12 + t95; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t45; 0; -m(6) * pkin(4) + m(7) * t59 + t29 + t51; m(6) + m(7); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t82; -t9; -t29; t32 * mrSges(7,1) - t33 * mrSges(7,2) + t86; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
