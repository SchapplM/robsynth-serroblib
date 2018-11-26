% Calculate joint inertia matrix for
% S6PRPRPR5
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
% Datum: 2018-11-23 14:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:57:52
% EndTime: 2018-11-23 14:57:53
% DurationCPUTime: 0.58s
% Computational Cost: add. (768->178), mult. (1602->241), div. (0->0), fcn. (1681->10), ass. (0->76)
t94 = m(6) + m(5);
t53 = sin(pkin(11));
t55 = cos(pkin(11));
t56 = cos(pkin(6));
t54 = sin(pkin(6));
t59 = sin(qJ(2));
t86 = t54 * t59;
t28 = -t53 * t86 + t56 * t55;
t29 = t56 * t53 + t55 * t86;
t58 = sin(qJ(4));
t91 = cos(qJ(4));
t13 = t58 * t28 + t91 * t29;
t10 = t13 ^ 2;
t49 = t55 ^ 2;
t92 = pkin(4) + pkin(9);
t57 = sin(qJ(6));
t90 = Ifges(7,4) * t57;
t60 = cos(qJ(6));
t89 = Ifges(7,4) * t60;
t34 = t58 * t53 - t91 * t55;
t88 = t34 * t57;
t87 = t34 * t60;
t61 = cos(qJ(2));
t85 = t54 * t61;
t84 = t60 * mrSges(7,1);
t83 = mrSges(6,1) + mrSges(5,3);
t82 = -mrSges(6,2) + mrSges(5,1);
t81 = pkin(8) + qJ(3);
t39 = t57 * mrSges(7,1) + t60 * mrSges(7,2);
t80 = t39 + mrSges(6,3);
t79 = t53 ^ 2 + t49;
t78 = t57 ^ 2 + t60 ^ 2;
t77 = t13 * qJ(5);
t76 = -m(4) - t94;
t35 = t91 * t53 + t58 * t55;
t74 = Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t35;
t38 = t81 * t55;
t72 = t81 * t53;
t21 = t58 * t38 + t91 * t72;
t23 = t91 * t38 - t58 * t72;
t73 = t21 ^ 2 + t23 ^ 2;
t44 = -t55 * pkin(3) - pkin(2);
t36 = m(7) * t78;
t71 = t78 * mrSges(7,3);
t37 = -t55 * mrSges(4,1) + t53 * mrSges(4,2);
t65 = -t35 * qJ(5) + t44;
t5 = t92 * t34 + t65;
t8 = t35 * pkin(5) + t21;
t1 = -t57 * t5 + t60 * t8;
t2 = t60 * t5 + t57 * t8;
t70 = t60 * t1 + t57 * t2;
t11 = -t91 * t28 + t58 * t29;
t3 = t60 * t11 + t57 * t85;
t4 = t57 * t11 - t60 * t85;
t69 = t60 * t3 + t57 * t4;
t68 = -t57 * mrSges(7,2) + t84;
t67 = t21 * t11 + t23 * t13;
t66 = -t28 * t53 + t29 * t55;
t63 = qJ(5) ^ 2;
t48 = t54 ^ 2;
t46 = Ifges(7,5) * t60;
t42 = t48 * t61 ^ 2;
t41 = Ifges(7,1) * t60 - t90;
t40 = -Ifges(7,2) * t57 + t89;
t31 = t35 * mrSges(6,3);
t30 = t35 * mrSges(5,2);
t20 = -t34 * mrSges(6,2) - t31;
t19 = t34 * mrSges(5,1) + t30;
t18 = -t35 * mrSges(7,2) + mrSges(7,3) * t87;
t17 = t35 * mrSges(7,1) - mrSges(7,3) * t88;
t16 = t68 * t34;
t15 = t34 * pkin(4) + t65;
t9 = -t34 * pkin(5) + t23;
t7 = Ifges(7,5) * t35 + (Ifges(7,1) * t57 + t89) * t34;
t6 = Ifges(7,6) * t35 + (Ifges(7,2) * t60 + t90) * t34;
t12 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t10) + m(4) * (t28 ^ 2 + t29 ^ 2 + t42) + m(3) * (t48 * t59 ^ 2 + t56 ^ 2 + t42) + t94 * (t11 ^ 2 + t10 + t42); t3 * t17 + t4 * t18 + t83 * t35 * t11 + t66 * mrSges(4,3) + (-t83 * t34 - t16) * t13 + (-t59 * mrSges(3,2) + (mrSges(3,1) - t19 - t20 - t37) * t61) * t54 + m(7) * (t1 * t3 + t9 * t13 + t2 * t4) + m(6) * (-t15 * t85 + t67) + m(5) * (-t44 * t85 + t67) + m(4) * (pkin(2) * t85 + qJ(3) * t66); Ifges(4,2) * t49 - 0.2e1 * pkin(2) * t37 + 0.2e1 * t1 * t17 + 0.2e1 * t15 * t20 - 0.2e1 * t9 * t16 + 0.2e1 * t2 * t18 + 0.2e1 * t44 * t19 + Ifges(3,3) + (Ifges(4,1) * t53 + 0.2e1 * Ifges(4,4) * t55) * t53 + 0.2e1 * t79 * qJ(3) * mrSges(4,3) + m(4) * (t79 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t44 ^ 2 + t73) + m(6) * (t15 ^ 2 + t73) + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + ((Ifges(5,1) + Ifges(6,2)) * t35 + 0.2e1 * t83 * t21 + t74) * t35 + (t57 * t7 + t60 * t6 + (Ifges(5,2) + Ifges(6,3)) * t34 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t35 - 0.2e1 * t83 * t23) * t34; m(7) * (-t57 * t3 + t60 * t4) + t76 * t85; -m(4) * pkin(2) - t57 * t17 + t60 * t18 + t30 - t31 + t82 * t34 + m(7) * (-t57 * t1 + t60 * t2) + m(6) * t15 + m(5) * t44 + t37; t36 - t76; -t82 * t11 - t69 * mrSges(7,3) + (-mrSges(5,2) + t80) * t13 + m(7) * (-t69 * t92 + t77) + m(6) * (-pkin(4) * t11 + t77); -qJ(5) * t16 + t9 * t39 + (mrSges(6,3) - mrSges(5,2)) * t23 - t82 * t21 + (t7 / 0.2e1 - t92 * t17 - t1 * mrSges(7,3)) * t60 + (-t6 / 0.2e1 - t92 * t18 - t2 * mrSges(7,3)) * t57 + m(7) * (qJ(5) * t9 - t70 * t92) + m(6) * (-pkin(4) * t21 + qJ(5) * t23) + (Ifges(5,5) - Ifges(6,4) - Ifges(7,6) * t57 / 0.2e1 + t46 / 0.2e1 - pkin(4) * mrSges(6,1)) * t35 + (-Ifges(5,6) + Ifges(6,5) - qJ(5) * mrSges(6,1) + t57 * t41 / 0.2e1 + t60 * t40 / 0.2e1) * t34; 0; -0.2e1 * pkin(4) * mrSges(6,2) - t57 * t40 + t60 * t41 + Ifges(6,1) + Ifges(5,3) + m(7) * (t78 * t92 ^ 2 + t63) + m(6) * (pkin(4) ^ 2 + t63) + 0.2e1 * t80 * qJ(5) + 0.2e1 * t92 * t71; m(6) * t11 + m(7) * t69; m(6) * t21 + m(7) * t70 + t35 * mrSges(6,1) + t60 * t17 + t57 * t18; 0; -m(6) * pkin(4) - t36 * t92 + mrSges(6,2) - t71; m(6) + t36; t3 * mrSges(7,1) - t4 * mrSges(7,2); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t74; -t39; -t92 * t84 + t46 + (mrSges(7,2) * t92 - Ifges(7,6)) * t57; t68; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
