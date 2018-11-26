% Calculate joint inertia matrix for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2018-11-23 15:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:06
% EndTime: 2018-11-23 15:38:06
% DurationCPUTime: 0.47s
% Computational Cost: add. (834->159), mult. (1397->228), div. (0->0), fcn. (1362->8), ass. (0->69)
t53 = cos(pkin(10));
t56 = sin(qJ(5));
t51 = sin(pkin(10));
t80 = cos(qJ(5));
t66 = t80 * t51;
t28 = t56 * t53 + t66;
t72 = t56 * t51;
t60 = t53 * t80 - t72;
t12 = mrSges(6,1) * t60 - t28 * mrSges(6,2);
t33 = t53 * mrSges(5,1) - t51 * mrSges(5,2);
t92 = -t12 - t33;
t55 = sin(qJ(6));
t75 = t28 * t55;
t10 = -mrSges(7,2) * t60 + mrSges(7,3) * t75;
t57 = cos(qJ(6));
t74 = t28 * t57;
t11 = mrSges(7,1) * t60 + mrSges(7,3) * t74;
t91 = t57 * t10 - t55 * t11;
t34 = -t57 * mrSges(7,1) + t55 * mrSges(7,2);
t90 = m(7) * pkin(5) + mrSges(6,1) - t34;
t52 = sin(pkin(9));
t54 = cos(pkin(9));
t58 = -pkin(1) - pkin(2);
t32 = t54 * qJ(2) + t52 * t58;
t29 = -qJ(4) + t32;
t81 = pkin(7) - t29;
t16 = t81 * t53;
t5 = -t56 * t16 - t66 * t81;
t89 = t5 ^ 2;
t17 = t28 * t52;
t88 = t17 ^ 2;
t87 = t60 ^ 2;
t47 = t53 ^ 2;
t86 = t55 / 0.2e1;
t85 = -m(5) - m(6);
t84 = t17 * t5;
t83 = t60 * pkin(5);
t82 = t60 * t5;
t79 = Ifges(7,4) * t55;
t78 = Ifges(7,4) * t57;
t76 = t60 * t17;
t70 = Ifges(7,5) * t74 - Ifges(7,3) * t60;
t69 = Ifges(7,5) * t55 + Ifges(7,6) * t57;
t68 = t51 ^ 2 + t47;
t67 = t55 ^ 2 + t57 ^ 2;
t65 = t68 * mrSges(5,3);
t64 = t67 * t28;
t31 = -t52 * qJ(2) + t54 * t58;
t30 = pkin(3) - t31;
t7 = -t16 * t80 + t72 * t81;
t21 = t53 * pkin(4) + t30;
t8 = t28 * pkin(8) + t21 + t83;
t1 = -t55 * t7 + t57 * t8;
t2 = t55 * t8 + t57 * t7;
t63 = -t1 * t55 + t2 * t57;
t62 = -mrSges(7,1) * t55 - mrSges(7,2) * t57;
t19 = t60 * t52;
t13 = -t55 * t19 - t57 * t54;
t14 = t57 * t19 - t55 * t54;
t61 = -t13 * t55 + t14 * t57;
t9 = t62 * t28;
t48 = t54 ^ 2;
t46 = t52 ^ 2;
t36 = Ifges(7,1) * t55 + t78;
t35 = Ifges(7,2) * t57 + t79;
t24 = t28 ^ 2;
t4 = Ifges(7,5) * t60 + (-Ifges(7,1) * t57 + t79) * t28;
t3 = Ifges(7,6) * t60 + (Ifges(7,2) * t55 - t78) * t28;
t6 = [t47 * Ifges(5,2) + (2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t31 * mrSges(4,1) + 0.2e1 * t32 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t1 * t11 + 0.2e1 * t2 * t10 + 0.2e1 * t21 * t12 + 0.2e1 * t30 * t33 + 0.2e1 * t5 * t9 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (Ifges(5,1) * t51 + 0.2e1 * Ifges(5,4) * t53) * t51 - 0.2e1 * t29 * t65 - (0.2e1 * t7 * mrSges(6,3) - Ifges(6,2) * t60 + t70) * t60 + (-0.2e1 * t5 * mrSges(6,3) + Ifges(6,1) * t28 + t55 * t3 - t57 * t4 - (-Ifges(7,6) * t55 - (2 * Ifges(6,4))) * t60) * t28 + m(4) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t29 ^ 2 * t68 + t30 ^ 2) + m(6) * (t21 ^ 2 + t7 ^ 2 + t89) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t89); -m(3) * pkin(1) + t14 * t10 + t13 * t11 + t17 * t9 - mrSges(3,1) + (-t17 * t28 - t19 * t60) * mrSges(6,3) + (-mrSges(4,1) + t92) * t54 + (mrSges(4,2) - t65) * t52 + m(7) * (t13 * t1 + t14 * t2 + t84) + m(6) * (t19 * t7 - t54 * t21 + t84) + m(5) * (t29 * t52 * t68 - t54 * t30) + m(4) * (t54 * t31 + t52 * t32); m(3) + m(4) * (t46 + t48) + m(5) * (t46 * t68 + t48) + m(6) * (t19 ^ 2 + t48 + t88) + m(7) * (t13 ^ 2 + t14 ^ 2 + t88); -t60 * t9 + t91 * t28 + m(7) * (t28 * t63 - t82) + m(6) * (t28 * t7 - t82); m(6) * (t28 * t19 - t76) + m(7) * (t28 * t61 - t76); m(4) + m(5) * t68 + m(6) * (t24 + t87) + m(7) * (t24 * t67 + t87); t55 * t10 + t57 * t11 + m(7) * (t57 * t1 + t55 * t2) + m(6) * t21 + m(5) * t30 - t92; m(7) * (t57 * t13 + t55 * t14) + t85 * t54; 0; m(7) * t67 - t85; t4 * t86 + t57 * t3 / 0.2e1 - pkin(5) * t9 - t7 * mrSges(6,2) + t63 * mrSges(7,3) + (-t57 * t36 / 0.2e1 + t35 * t86 - Ifges(6,5)) * t28 - (-t69 / 0.2e1 + Ifges(6,6)) * t60 - t90 * t5 + (m(7) * t63 + t91) * pkin(8); -t19 * mrSges(6,2) + (m(7) * pkin(8) + mrSges(7,3)) * t61 - t90 * t17; -t60 * t34 + m(7) * (pkin(8) * t64 + t83) + mrSges(7,3) * t64 + t12; 0; Ifges(6,3) + m(7) * (pkin(8) ^ 2 * t67 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t34 + t57 * t35 + t55 * t36 + 0.2e1 * t67 * pkin(8) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,6) * t75 - t70; t13 * mrSges(7,1) - t14 * mrSges(7,2); t9; -t34; pkin(8) * t62 + t69; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
