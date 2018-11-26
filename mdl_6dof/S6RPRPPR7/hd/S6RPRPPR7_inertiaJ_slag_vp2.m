% Calculate joint inertia matrix for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2018-11-23 15:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:55:27
% EndTime: 2018-11-23 15:55:27
% DurationCPUTime: 0.59s
% Computational Cost: add. (783->174), mult. (1334->227), div. (0->0), fcn. (1270->6), ass. (0->70)
t53 = cos(qJ(3));
t92 = t53 ^ 2;
t91 = m(6) + m(5);
t72 = -mrSges(5,3) - mrSges(6,1);
t50 = sin(qJ(6));
t52 = cos(qJ(6));
t70 = t50 ^ 2 + t52 ^ 2;
t29 = m(7) * t70;
t90 = t70 * mrSges(7,3);
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t51 = sin(qJ(3));
t25 = t48 * t53 + t49 * t51;
t27 = -t48 * t51 + t49 * t53;
t40 = t51 * pkin(3) + qJ(2);
t60 = -t27 * qJ(5) + t40;
t8 = t25 * pkin(4) + t60;
t88 = -0.2e1 * t8;
t22 = t25 ^ 2;
t87 = t27 ^ 2;
t82 = t48 * pkin(3);
t36 = qJ(5) + t82;
t86 = t36 ^ 2;
t85 = 0.2e1 * t40;
t83 = m(5) * pkin(3);
t81 = t49 * pkin(3);
t80 = m(6) + t29;
t79 = Ifges(7,4) * t50;
t78 = Ifges(7,4) * t52;
t77 = Ifges(7,6) * t50;
t76 = t25 * t50;
t75 = t25 * t52;
t74 = t36 * t25;
t73 = mrSges(6,2) - mrSges(5,1);
t31 = t50 * mrSges(7,1) + t52 * mrSges(7,2);
t71 = t31 + mrSges(6,3);
t69 = t51 ^ 2 + t92;
t54 = -pkin(1) - pkin(7);
t68 = -qJ(4) + t54;
t67 = Ifges(7,5) * t76 + Ifges(7,6) * t75 + Ifges(7,3) * t27;
t30 = t68 * t51;
t62 = t68 * t53;
t12 = t48 * t30 - t49 * t62;
t14 = t49 * t30 + t48 * t62;
t66 = t12 ^ 2 + t14 ^ 2;
t39 = -pkin(4) - t81;
t65 = m(4) * t69;
t35 = -pkin(8) + t39;
t64 = t70 * t35;
t63 = t69 * mrSges(4,3);
t3 = (pkin(4) + pkin(8)) * t25 + t60;
t6 = t27 * pkin(5) + t12;
t1 = -t50 * t3 + t52 * t6;
t2 = t52 * t3 + t50 * t6;
t61 = t52 * t1 + t50 * t2;
t59 = t52 * mrSges(7,1) - t50 * mrSges(7,2);
t10 = t27 * mrSges(7,1) - mrSges(7,3) * t76;
t11 = -t27 * mrSges(7,2) + mrSges(7,3) * t75;
t58 = -t52 * t10 - t50 * t11;
t55 = qJ(2) ^ 2;
t41 = Ifges(7,5) * t52;
t33 = Ifges(7,1) * t52 - t79;
t32 = -Ifges(7,2) * t50 + t78;
t19 = t27 * mrSges(6,3);
t18 = t27 * mrSges(5,2);
t9 = t59 * t25;
t7 = -t25 * pkin(5) + t14;
t5 = Ifges(7,5) * t27 + (Ifges(7,1) * t50 + t78) * t25;
t4 = Ifges(7,6) * t27 + (Ifges(7,2) * t52 + t79) * t25;
t13 = [Ifges(4,1) * t92 + t18 * t85 + t19 * t88 - 0.2e1 * t7 * t9 + 0.2e1 * t1 * t10 + 0.2e1 * t2 * t11 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t54 * t63 + (mrSges(5,1) * t85 + mrSges(6,2) * t88 + t52 * t4 + t50 * t5 + (Ifges(6,3) + Ifges(5,2)) * t25 + 0.2e1 * t72 * t14) * t25 + m(4) * (t69 * t54 ^ 2 + t55) + m(3) * ((pkin(1) ^ 2) + t55) + m(5) * (t40 ^ 2 + t66) + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(6) * (t8 ^ 2 + t66) + ((Ifges(5,1) + Ifges(6,2)) * t27 + t67 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t25 - 0.2e1 * t72 * t12) * t27 + (-0.2e1 * Ifges(4,4) * t53 + Ifges(4,2) * t51) * t51 + 0.2e1 * (t51 * mrSges(4,1) + t53 * mrSges(4,2) + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) + t58 * t27 - t63 + (t72 * t25 - t9) * t25 + m(7) * (t25 * t7 - t61 * t27) + t54 * t65 + t91 * (-t27 * t12 + t25 * t14) + t72 * t87; m(3) + m(7) * (t70 * t87 + t22) + t65 + t91 * (t22 + t87); t7 * t31 - t36 * t9 + (t54 * mrSges(4,1) + Ifges(4,5)) * t53 + (-t54 * mrSges(4,2) - Ifges(4,6)) * t51 + (-mrSges(5,2) + mrSges(6,3)) * t14 + t73 * t12 + (t35 * t10 - t1 * mrSges(7,3) + t5 / 0.2e1) * t52 + (t35 * t11 - t2 * mrSges(7,3) - t4 / 0.2e1) * t50 + m(7) * (t61 * t35 + t36 * t7) + m(6) * (t39 * t12 + t36 * t14) + (-t12 * t49 + t14 * t48) * t83 + (t39 * mrSges(6,1) - t77 / 0.2e1 + t41 / 0.2e1 + Ifges(5,5) - Ifges(6,4) - mrSges(5,3) * t81) * t27 + (t50 * t33 / 0.2e1 + t52 * t32 / 0.2e1 - t36 * mrSges(6,1) - Ifges(5,6) + Ifges(6,5) - mrSges(5,3) * t82) * t25; t53 * mrSges(4,1) - t51 * mrSges(4,2) + (-mrSges(5,2) + t71) * t25 + (-t73 + t90) * t27 + m(7) * (-t27 * t64 + t74) + m(6) * (-t39 * t27 + t74) + (t25 * t48 + t27 * t49) * t83; 0.2e1 * t39 * mrSges(6,2) - t50 * t32 + t52 * t33 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3) + m(7) * (t70 * t35 ^ 2 + t86) + m(6) * (t39 ^ 2 + t86) + m(5) * (t48 ^ 2 + t49 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t71 * t36 + 0.2e1 * (t49 * mrSges(5,1) - t48 * mrSges(5,2)) * pkin(3) - 0.2e1 * mrSges(7,3) * t64; -t50 * t10 + t52 * t11 + t18 - t19 - t73 * t25 + m(7) * (-t50 * t1 + t52 * t2) + m(6) * t8 + m(5) * t40; 0; 0; m(5) + t80; m(6) * t12 + m(7) * t61 + t27 * mrSges(6,1) - t58; 0.2e1 * (-t29 / 0.2e1 - m(6) / 0.2e1) * t27; m(6) * t39 + t35 * t29 + mrSges(6,2) - t90; 0; t80; t1 * mrSges(7,1) - t2 * mrSges(7,2) + t67; -t59 * t27; t59 * t35 + t41 - t77; -t31; t59; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
