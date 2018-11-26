% Calculate joint inertia matrix for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2018-11-23 15:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:40:00
% EndTime: 2018-11-23 15:40:00
% DurationCPUTime: 0.48s
% Computational Cost: add. (731->152), mult. (1325->202), div. (0->0), fcn. (1283->8), ass. (0->70)
t88 = m(5) + m(6);
t47 = cos(pkin(10));
t45 = sin(pkin(10));
t50 = sin(qJ(4));
t74 = t50 * t45;
t79 = cos(qJ(4));
t25 = -t79 * t47 + t74;
t87 = 0.2e1 * t25;
t23 = t25 ^ 2;
t64 = t79 * t45;
t27 = t50 * t47 + t64;
t24 = t27 ^ 2;
t41 = t47 ^ 2;
t85 = -2 * mrSges(6,2);
t48 = cos(pkin(9));
t37 = -t48 * pkin(1) - pkin(2);
t30 = -t47 * pkin(3) + t37;
t84 = 0.2e1 * t30;
t83 = pkin(4) + pkin(8);
t82 = pkin(4) * t25;
t49 = sin(qJ(6));
t51 = cos(qJ(6));
t68 = t49 ^ 2 + t51 ^ 2;
t29 = m(7) * t68;
t81 = m(6) + t29;
t46 = sin(pkin(9));
t35 = t46 * pkin(1) + qJ(3);
t80 = pkin(7) + t35;
t78 = Ifges(7,4) * t49;
t77 = Ifges(7,4) * t51;
t76 = t25 * t49;
t75 = t25 * t51;
t73 = t51 * mrSges(7,1);
t72 = mrSges(6,1) + mrSges(5,3);
t71 = -mrSges(6,2) + mrSges(5,1);
t20 = t27 * mrSges(5,2);
t21 = t27 * mrSges(6,3);
t70 = t21 - t20;
t69 = t45 ^ 2 + t41;
t67 = t27 * qJ(5);
t66 = Ifges(7,5) * t76 + Ifges(7,6) * t75 + Ifges(7,3) * t27;
t19 = t80 * t47;
t11 = t79 * t19 - t80 * t74;
t9 = t50 * t19 + t80 * t64;
t65 = t11 ^ 2 + t9 ^ 2;
t63 = t68 * mrSges(7,3);
t62 = t68 * t83;
t61 = -t47 * mrSges(4,1) + t45 * mrSges(4,2);
t55 = t30 - t67;
t3 = t83 * t25 + t55;
t4 = t27 * pkin(5) + t9;
t1 = -t49 * t3 + t51 * t4;
t2 = t51 * t3 + t49 * t4;
t60 = t51 * t1 + t49 * t2;
t58 = mrSges(6,2) - t63;
t57 = -t49 * mrSges(7,2) + t73;
t14 = t27 * mrSges(7,1) - mrSges(7,3) * t76;
t15 = -t27 * mrSges(7,2) + mrSges(7,3) * t75;
t56 = t51 * t14 + t49 * t15;
t53 = qJ(5) ^ 2;
t39 = Ifges(7,5) * t51;
t33 = Ifges(7,1) * t51 - t78;
t32 = -Ifges(7,2) * t49 + t77;
t31 = t49 * mrSges(7,1) + t51 * mrSges(7,2);
t13 = t57 * t25;
t8 = t55 + t82;
t7 = Ifges(7,5) * t27 + (Ifges(7,1) * t49 + t77) * t25;
t6 = Ifges(7,6) * t27 + (Ifges(7,2) * t51 + t78) * t25;
t5 = -t25 * pkin(5) + t11;
t10 = [Ifges(4,2) * t41 + 0.2e1 * t37 * t61 + t20 * t84 - 0.2e1 * t8 * t21 - 0.2e1 * t5 * t13 + 0.2e1 * t1 * t14 + 0.2e1 * t2 * t15 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t45 + 0.2e1 * Ifges(4,4) * t47) * t45 + 0.2e1 * (t48 * mrSges(3,1) - t46 * mrSges(3,2)) * pkin(1) + 0.2e1 * t69 * t35 * mrSges(4,3) + (mrSges(5,1) * t84 + t8 * t85 + t49 * t7 + t51 * t6 + (Ifges(5,2) + Ifges(6,3)) * t25 - 0.2e1 * t72 * t11) * t25 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t8 ^ 2 + t65) + m(5) * (t30 ^ 2 + t65) + m(4) * (t69 * t35 ^ 2 + t37 ^ 2) + m(3) * (t46 ^ 2 + t48 ^ 2) * pkin(1) ^ 2 + ((Ifges(5,1) + Ifges(6,2)) * t27 + t66 + 0.2e1 * t72 * t9 + (-Ifges(5,4) - Ifges(6,6)) * t87) * t27; -t27 * t13 + t56 * t25 + m(7) * (t60 * t25 + t5 * t27) + t88 * (t11 * t27 + t9 * t25); m(3) + m(7) * (t68 * t23 + t24) + m(4) * t69 + t88 * (t23 + t24); -t49 * t14 + t51 * t15 + t71 * t25 + m(7) * (-t49 * t1 + t51 * t2) + m(6) * t8 + m(5) * t30 + m(4) * t37 + t61 - t70; 0; m(4) + m(5) + t81; -qJ(5) * t13 + t5 * t31 - t71 * t9 + (-mrSges(5,2) + mrSges(6,3)) * t11 + (t7 / 0.2e1 - t83 * t14 - t1 * mrSges(7,3)) * t51 + (-t6 / 0.2e1 - t83 * t15 - t2 * mrSges(7,3)) * t49 + m(7) * (qJ(5) * t5 - t60 * t83) + m(6) * (-pkin(4) * t9 + qJ(5) * t11) + (-Ifges(7,6) * t49 / 0.2e1 + t39 / 0.2e1 - Ifges(6,4) + Ifges(5,5) - pkin(4) * mrSges(6,1)) * t27 + (Ifges(6,5) - Ifges(5,6) - qJ(5) * mrSges(6,1) + t49 * t33 / 0.2e1 + t51 * t32 / 0.2e1) * t25; t27 * t31 + (-mrSges(5,1) + t58) * t25 + m(7) * (-t25 * t62 + t67) + m(6) * (t67 - t82) + t70; 0; pkin(4) * t85 - t49 * t32 + t51 * t33 + Ifges(6,1) + Ifges(5,3) + m(7) * (t68 * t83 ^ 2 + t53) + m(6) * (pkin(4) ^ 2 + t53) + 0.2e1 * (t31 + mrSges(6,3)) * qJ(5) + 0.2e1 * t83 * t63; m(6) * t9 + m(7) * t60 + t27 * mrSges(6,1) + t56; (t29 / 0.2e1 + m(6) / 0.2e1) * t87; 0; -m(6) * pkin(4) - m(7) * t62 + t58; t81; t1 * mrSges(7,1) - t2 * mrSges(7,2) + t66; t13; -t31; -t83 * t73 + t39 + (mrSges(7,2) * t83 - Ifges(7,6)) * t49; t57; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
