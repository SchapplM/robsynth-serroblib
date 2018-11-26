% Calculate joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 15:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:41:00
% EndTime: 2018-11-23 15:41:00
% DurationCPUTime: 0.52s
% Computational Cost: add. (897->180), mult. (1505->264), div. (0->0), fcn. (1451->8), ass. (0->69)
t54 = sin(pkin(10));
t56 = cos(pkin(10));
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t25 = t54 * t59 - t56 * t61;
t27 = t54 * t61 + t56 * t59;
t58 = sin(qJ(6));
t80 = t27 * t58;
t10 = mrSges(7,2) * t25 + mrSges(7,3) * t80;
t60 = cos(qJ(6));
t79 = t27 * t60;
t11 = -mrSges(7,1) * t25 + mrSges(7,3) * t79;
t91 = t60 * t10 - t58 * t11;
t55 = sin(pkin(9));
t57 = cos(pkin(9));
t62 = -pkin(1) - pkin(2);
t32 = qJ(2) * t57 + t55 * t62;
t30 = -pkin(7) + t32;
t71 = qJ(5) - t30;
t16 = t71 * t61;
t68 = t71 * t59;
t5 = -t16 * t54 - t56 * t68;
t90 = t5 ^ 2;
t17 = t27 * t55;
t89 = t17 ^ 2;
t88 = t25 ^ 2;
t53 = t61 ^ 2;
t87 = m(6) * pkin(4);
t86 = t58 / 0.2e1;
t85 = t17 * t5;
t84 = t25 * t5;
t83 = Ifges(7,4) * t58;
t82 = Ifges(7,4) * t60;
t81 = t25 * t17;
t33 = -mrSges(7,1) * t60 + mrSges(7,2) * t58;
t76 = -mrSges(6,1) + t33;
t75 = Ifges(7,5) * t79 + Ifges(7,3) * t25;
t74 = Ifges(7,5) * t58 + Ifges(7,6) * t60;
t73 = t58 ^ 2 + t60 ^ 2;
t72 = t59 ^ 2 + t53;
t39 = pkin(4) * t54 + pkin(8);
t70 = t73 * t39;
t69 = t72 * mrSges(5,3);
t12 = -t25 * mrSges(6,1) - t27 * mrSges(6,2);
t34 = mrSges(5,1) * t61 - t59 * mrSges(5,2);
t31 = -qJ(2) * t55 + t57 * t62;
t29 = pkin(3) - t31;
t7 = -t56 * t16 + t54 * t68;
t21 = pkin(4) * t61 + t29;
t8 = -pkin(5) * t25 + pkin(8) * t27 + t21;
t1 = -t58 * t7 + t60 * t8;
t2 = t58 * t8 + t60 * t7;
t67 = -t1 * t58 + t2 * t60;
t66 = -mrSges(5,1) * t59 - mrSges(5,2) * t61;
t65 = -mrSges(7,1) * t58 - mrSges(7,2) * t60;
t19 = t25 * t55;
t13 = t19 * t58 - t57 * t60;
t14 = -t19 * t60 - t57 * t58;
t64 = -t13 * t58 + t14 * t60;
t9 = t65 * t27;
t49 = t57 ^ 2;
t48 = t55 ^ 2;
t40 = -pkin(4) * t56 - pkin(5);
t36 = Ifges(7,1) * t58 + t82;
t35 = Ifges(7,2) * t60 + t83;
t24 = t27 ^ 2;
t4 = -Ifges(7,5) * t25 + (-Ifges(7,1) * t60 + t83) * t27;
t3 = -Ifges(7,6) * t25 + (Ifges(7,2) * t58 - t82) * t27;
t6 = [t53 * Ifges(5,2) + (2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t31 * mrSges(4,1) + 0.2e1 * t32 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t1 * t11 + 0.2e1 * t2 * t10 + 0.2e1 * t21 * t12 + 0.2e1 * t29 * t34 + 0.2e1 * t5 * t9 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (Ifges(5,1) * t59 + 0.2e1 * Ifges(5,4) * t61) * t59 - 0.2e1 * t30 * t69 + (0.2e1 * t7 * mrSges(6,3) + Ifges(6,2) * t25 + t75) * t25 + (-0.2e1 * t5 * mrSges(6,3) + Ifges(6,1) * t27 + t58 * t3 - t60 * t4 + (-Ifges(7,6) * t58 - (2 * Ifges(6,4))) * t25) * t27 + m(5) * (t30 ^ 2 * t72 + t29 ^ 2) + m(4) * (t31 ^ 2 + t32 ^ 2) + m(6) * (t21 ^ 2 + t7 ^ 2 + t90) + m(7) * (t1 ^ 2 + t2 ^ 2 + t90) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2); -m(3) * pkin(1) + t14 * t10 + t13 * t11 + t17 * t9 - mrSges(3,1) + (-t17 * t27 - t19 * t25) * mrSges(6,3) + (-mrSges(4,1) - t12 - t34) * t57 + (mrSges(4,2) - t69) * t55 + m(7) * (t1 * t13 + t14 * t2 + t85) + m(6) * (-t19 * t7 - t21 * t57 + t85) + m(5) * (t30 * t55 * t72 - t57 * t29) + m(4) * (t31 * t57 + t32 * t55); m(3) + m(7) * (t13 ^ 2 + t14 ^ 2 + t89) + m(5) * (t48 * t72 + t49) + m(6) * (t19 ^ 2 + t49 + t89) + m(4) * (t48 + t49); t25 * t9 + t91 * t27 + m(7) * (t27 * t67 + t84) + m(6) * (t27 * t7 + t84); m(7) * (t27 * t64 + t81) + m(6) * (-t19 * t27 + t81); m(4) + m(5) * t72 + m(6) * (t24 + t88) + m(7) * (t24 * t73 + t88); Ifges(6,6) * t25 - Ifges(5,5) * t59 - Ifges(5,6) * t61 + t5 * t33 + t4 * t86 + t60 * t3 / 0.2e1 - t7 * mrSges(6,2) - t5 * mrSges(6,1) - t25 * t74 / 0.2e1 + t66 * t30 + t67 * mrSges(7,3) + (-Ifges(6,5) + t35 * t86 - t60 * t36 / 0.2e1) * t27 + (m(6) * (-t5 * t56 + t54 * t7) + (t25 * t54 + t27 * t56) * mrSges(6,3)) * pkin(4) + (m(7) * t5 + t9) * t40 + (m(7) * t67 + t91) * t39; t19 * mrSges(6,2) + t66 * t55 + t76 * t17 + t64 * mrSges(7,3) + m(7) * (t40 * t17 + t39 * t64) + (-t17 * t56 - t19 * t54) * t87; t76 * t25 + (mrSges(7,3) * t73 - mrSges(6,2)) * t27 + m(7) * (t40 * t25 + t27 * t70) + (-t25 * t56 + t27 * t54) * t87 + t34; 0.2e1 * t40 * t33 + t60 * t35 + t58 * t36 + Ifges(5,3) + Ifges(6,3) + m(7) * (t39 ^ 2 * t73 + t40 ^ 2) + m(6) * (t54 ^ 2 + t56 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (t56 * mrSges(6,1) - t54 * mrSges(6,2)) * pkin(4) + 0.2e1 * mrSges(7,3) * t70; t58 * t10 + t60 * t11 + m(7) * (t1 * t60 + t2 * t58) + m(6) * t21 + t12; m(7) * (t13 * t60 + t14 * t58) - m(6) * t57; 0; 0; m(7) * t73 + m(6); mrSges(7,1) * t1 - mrSges(7,2) * t2 + Ifges(7,6) * t80 - t75; mrSges(7,1) * t13 - mrSges(7,2) * t14; t9; t39 * t65 + t74; -t33; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
