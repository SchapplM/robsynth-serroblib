% Calculate joint inertia matrix for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:30
% EndTime: 2019-03-09 01:55:31
% DurationCPUTime: 0.51s
% Computational Cost: add. (724->158), mult. (1238->203), div. (0->0), fcn. (1198->6), ass. (0->67)
t89 = m(6) + m(5);
t71 = mrSges(6,1) + mrSges(5,3);
t49 = sin(qJ(6));
t50 = cos(qJ(6));
t67 = t49 ^ 2 + t50 ^ 2;
t29 = m(7) * t67;
t88 = mrSges(7,3) * t67;
t46 = sin(pkin(9));
t51 = cos(qJ(4));
t47 = cos(pkin(9));
t79 = sin(qJ(4));
t62 = t79 * t47;
t25 = t51 * t46 + t62;
t22 = t25 ^ 2;
t73 = t51 * t47;
t27 = -t79 * t46 + t73;
t86 = t27 ^ 2;
t41 = t47 ^ 2;
t85 = -2 * mrSges(6,2);
t34 = t46 * pkin(3) + qJ(2);
t84 = 0.2e1 * t34;
t82 = pkin(4) + pkin(8);
t81 = m(6) + t29;
t48 = -pkin(1) - qJ(3);
t80 = -pkin(7) + t48;
t78 = Ifges(7,4) * t49;
t77 = Ifges(7,4) * t50;
t76 = t25 * t49;
t75 = t25 * t50;
t74 = t50 * mrSges(7,1);
t72 = mrSges(5,1) - mrSges(6,2);
t31 = t49 * mrSges(7,1) + t50 * mrSges(7,2);
t70 = mrSges(6,3) + t31;
t69 = t46 * mrSges(4,1) + t47 * mrSges(4,2);
t68 = t46 ^ 2 + t41;
t66 = qJ(5) * t25;
t65 = Ifges(7,5) * t76 + Ifges(7,6) * t75 + Ifges(7,3) * t27;
t30 = t80 * t46;
t12 = t79 * t30 - t80 * t73;
t14 = t51 * t30 + t80 * t62;
t64 = t12 ^ 2 + t14 ^ 2;
t63 = m(4) * t68;
t61 = t68 * mrSges(4,3);
t60 = t67 * t82;
t58 = -t27 * qJ(5) + t34;
t3 = t82 * t25 + t58;
t6 = t27 * pkin(5) + t12;
t1 = -t49 * t3 + t50 * t6;
t2 = t50 * t3 + t49 * t6;
t59 = t50 * t1 + t49 * t2;
t57 = -t49 * mrSges(7,2) + t74;
t10 = t27 * mrSges(7,1) - mrSges(7,3) * t76;
t11 = -t27 * mrSges(7,2) + mrSges(7,3) * t75;
t56 = -t50 * t10 - t49 * t11;
t54 = qJ(2) ^ 2;
t53 = qJ(5) ^ 2;
t39 = Ifges(7,5) * t50;
t33 = Ifges(7,1) * t50 - t78;
t32 = -Ifges(7,2) * t49 + t77;
t20 = t27 * mrSges(6,3);
t19 = t27 * mrSges(5,2);
t9 = t57 * t25;
t8 = t25 * pkin(4) + t58;
t7 = -t25 * pkin(5) + t14;
t5 = Ifges(7,5) * t27 + (Ifges(7,1) * t49 + t77) * t25;
t4 = Ifges(7,6) * t27 + (Ifges(7,2) * t50 + t78) * t25;
t13 = [Ifges(4,1) * t41 + t19 * t84 - 0.2e1 * t8 * t20 - 0.2e1 * t7 * t9 + 0.2e1 * t1 * t10 + 0.2e1 * t2 * t11 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) + (-0.2e1 * Ifges(4,4) * t47 + Ifges(4,2) * t46) * t46 - 0.2e1 * t48 * t61 + (mrSges(5,1) * t84 + t8 * t85 + t50 * t4 + t49 * t5 + (Ifges(5,2) + Ifges(6,3)) * t25 - 0.2e1 * t71 * t14) * t25 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(6) * (t8 ^ 2 + t64) + m(5) * (t34 ^ 2 + t64) + m(3) * ((pkin(1) ^ 2) + t54) + m(4) * (t68 * t48 ^ 2 + t54) + ((Ifges(5,1) + Ifges(6,2)) * t27 + t65 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t25 + 0.2e1 * t71 * t12) * t27 + 0.2e1 * (t69 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) + t56 * t27 - t61 + (-t25 * t71 - t9) * t25 + m(7) * (t25 * t7 - t27 * t59) + t48 * t63 + t89 * (-t27 * t12 + t25 * t14) - t71 * t86; m(3) + m(7) * (t67 * t86 + t22) + t63 + t89 * (t22 + t86); m(4) * qJ(2) - t49 * t10 + t50 * t11 + t19 - t20 + t72 * t25 + m(7) * (-t49 * t1 + t50 * t2) + m(6) * t8 + m(5) * t34 + t69; 0; m(4) + m(5) + t81; -qJ(5) * t9 + t7 * t31 + (-mrSges(5,2) + mrSges(6,3)) * t14 - t72 * t12 + (t5 / 0.2e1 - t82 * t10 - t1 * mrSges(7,3)) * t50 + (-t4 / 0.2e1 - t82 * t11 - t2 * mrSges(7,3)) * t49 + m(7) * (qJ(5) * t7 - t59 * t82) + m(6) * (-pkin(4) * t12 + qJ(5) * t14) + (-Ifges(7,6) * t49 / 0.2e1 + t39 / 0.2e1 - Ifges(6,4) + Ifges(5,5) - pkin(4) * mrSges(6,1)) * t27 + (Ifges(6,5) - Ifges(5,6) + t49 * t33 / 0.2e1 + t50 * t32 / 0.2e1 - qJ(5) * mrSges(6,1)) * t25; (-mrSges(5,2) + t70) * t25 + (t72 + t88) * t27 + m(6) * (pkin(4) * t27 + t66) + m(7) * (t27 * t60 + t66); 0; pkin(4) * t85 - t49 * t32 + t50 * t33 + Ifges(6,1) + Ifges(5,3) + m(6) * (pkin(4) ^ 2 + t53) + m(7) * (t67 * t82 ^ 2 + t53) + 0.2e1 * mrSges(7,3) * t60 + 0.2e1 * t70 * qJ(5); m(6) * t12 + m(7) * t59 + t27 * mrSges(6,1) - t56; 0.2e1 * (-m(6) / 0.2e1 - t29 / 0.2e1) * t27; 0; -m(6) * pkin(4) - t29 * t82 + mrSges(6,2) - t88; t81; t1 * mrSges(7,1) - t2 * mrSges(7,2) + t65; -t57 * t27; -t31; -t82 * t74 + t39 + (mrSges(7,2) * t82 - Ifges(7,6)) * t49; t57; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
