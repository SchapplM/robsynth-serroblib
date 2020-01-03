% Calculate joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:02
% DurationCPUTime: 0.44s
% Computational Cost: add. (477->135), mult. (868->169), div. (0->0), fcn. (661->6), ass. (0->60)
t53 = cos(qJ(3));
t48 = t53 ^ 2;
t50 = sin(qJ(3));
t94 = t50 ^ 2 + t48;
t93 = t53 * pkin(3) + t50 * qJ(4);
t86 = pkin(7) - pkin(8);
t30 = t86 * t50;
t31 = t86 * t53;
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t6 = t52 * t30 - t49 * t31;
t7 = t49 * t30 + t52 * t31;
t92 = t6 * mrSges(6,1) - t7 * mrSges(6,2);
t51 = sin(qJ(2));
t35 = t51 * pkin(1) + pkin(7);
t79 = -pkin(8) + t35;
t15 = t79 * t50;
t16 = t79 * t53;
t3 = t52 * t15 - t49 * t16;
t4 = t49 * t15 + t52 * t16;
t91 = t3 * mrSges(6,1) - t4 * mrSges(6,2);
t90 = -m(5) * pkin(3) - mrSges(5,1);
t89 = (mrSges(4,3) + mrSges(5,2)) * t94;
t20 = -t50 * t49 - t53 * t52;
t21 = -t53 * t49 + t50 * t52;
t5 = -t20 * mrSges(6,1) + t21 * mrSges(6,2);
t88 = 0.2e1 * t5;
t27 = -t53 * mrSges(5,1) - t50 * mrSges(5,3);
t87 = 0.2e1 * t27;
t85 = m(5) * t50;
t78 = t20 * mrSges(6,3);
t77 = t21 * mrSges(6,3);
t55 = -pkin(3) - pkin(4);
t22 = -t49 * qJ(4) + t52 * t55;
t76 = t22 * mrSges(6,1);
t23 = t52 * qJ(4) + t49 * t55;
t75 = t23 * mrSges(6,2);
t39 = t50 * mrSges(5,2);
t74 = Ifges(6,5) * t21 + Ifges(6,6) * t20;
t73 = t94 * pkin(7) * t35;
t72 = t94 * t35 ^ 2;
t71 = t94 * pkin(7) ^ 2;
t69 = qJ(4) * t53;
t68 = 0.2e1 * mrSges(6,3);
t67 = pkin(2) + t93;
t54 = cos(qJ(2));
t36 = -t54 * pkin(1) - pkin(2);
t65 = t52 * mrSges(6,1) - t49 * mrSges(6,2);
t11 = t36 - t93;
t64 = t49 * t78 - t52 * t77 + t39;
t63 = (t54 * mrSges(3,1) - t51 * mrSges(3,2)) * pkin(1);
t62 = Ifges(6,1) * t21 ^ 2 + Ifges(3,3) + (0.2e1 * Ifges(6,4) * t21 + Ifges(6,2) * t20) * t20 + (Ifges(5,3) + Ifges(4,2)) * t48 + ((Ifges(5,1) + Ifges(4,1)) * t50 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t53) * t50;
t60 = 0.2e1 * t89;
t59 = mrSges(5,2) * t69 - pkin(3) * t39 - t22 * t77 + t23 * t78 - t74 + (Ifges(4,6) - Ifges(5,6)) * t53 + (Ifges(5,4) + Ifges(4,5)) * t50;
t58 = m(5) * t69 + (mrSges(5,3) - mrSges(4,2)) * t53 + (-mrSges(4,1) + t90) * t50;
t45 = t53 * pkin(4);
t28 = -t53 * mrSges(4,1) + t50 * mrSges(4,2);
t12 = t45 + t67;
t9 = -t11 + t45;
t1 = [t62 + t60 * t35 + 0.2e1 * t63 + (t4 * t20 - t3 * t21) * t68 + m(3) * (t51 ^ 2 + t54 ^ 2) * pkin(1) ^ 2 + m(4) * (t36 ^ 2 + t72) + m(6) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t72) + t11 * t87 + 0.2e1 * t36 * t28 + t9 * t88 + Ifges(2,3); t62 + ((-t3 - t6) * t21 + (t4 + t7) * t20) * mrSges(6,3) + m(4) * (-pkin(2) * t36 + t73) + m(6) * (t12 * t9 + t6 * t3 + t7 * t4) + m(5) * (-t11 * t67 + t73) + (t12 + t9) * t5 + (-pkin(2) + t36) * t28 + (-t67 + t11) * t27 + t63 + (pkin(7) + t35) * t89; -0.2e1 * pkin(2) * t28 + t12 * t88 - t67 * t87 + (t7 * t20 - t6 * t21) * t68 + m(6) * (t12 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t67 ^ 2 + t71) + m(4) * (pkin(2) ^ 2 + t71) + t60 * pkin(7) + t62; t59 + t58 * t35 + m(6) * (t22 * t3 + t23 * t4) - t91; t59 + t58 * pkin(7) + m(6) * (t22 * t6 + t23 * t7) - t92; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t76 + 0.2e1 * t75 + 0.2e1 * qJ(4) * mrSges(5,3) + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2); m(6) * (t52 * t3 + t49 * t4) + t35 * t85 + t64; m(6) * (t49 * t7 + t52 * t6) + pkin(7) * t85 + t64; m(6) * (t52 * t22 + t49 * t23) - t65 + t90; m(5) + m(6) * (t49 ^ 2 + t52 ^ 2); t74 + t91; t74 + t92; -Ifges(6,3) - t75 + t76; t65; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
