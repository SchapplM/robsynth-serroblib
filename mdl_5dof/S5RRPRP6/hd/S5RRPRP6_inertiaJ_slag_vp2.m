% Calculate joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:56:48
% DurationCPUTime: 0.57s
% Computational Cost: add. (676->177), mult. (1300->247), div. (0->0), fcn. (1268->6), ass. (0->72)
t91 = Ifges(5,3) + Ifges(6,3);
t90 = 2 * mrSges(6,3);
t60 = sin(qJ(4));
t62 = cos(qJ(4));
t72 = t60 ^ 2 + t62 ^ 2;
t89 = 0.2e1 * t72;
t63 = cos(qJ(2));
t73 = -qJ(3) - pkin(6);
t38 = t73 * t63;
t58 = sin(pkin(8));
t59 = cos(pkin(8));
t61 = sin(qJ(2));
t68 = t73 * t61;
t20 = -t38 * t58 - t59 * t68;
t88 = t20 ^ 2;
t87 = 0.2e1 * t20;
t48 = -pkin(2) * t63 - pkin(1);
t86 = 0.2e1 * t48;
t85 = m(6) * pkin(4);
t84 = t58 * pkin(2);
t83 = t59 * pkin(2);
t82 = Ifges(5,4) * t60;
t81 = Ifges(5,4) * t62;
t80 = Ifges(6,4) * t60;
t79 = Ifges(6,4) * t62;
t33 = t58 * t63 + t59 * t61;
t78 = t33 * t60;
t77 = t33 * t62;
t76 = t60 * mrSges(5,2);
t75 = t60 * mrSges(6,3);
t74 = -Ifges(5,6) - Ifges(6,6);
t32 = t58 * t61 - t59 * t63;
t14 = pkin(3) * t32 - pkin(7) * t33 + t48;
t22 = -t59 * t38 + t58 * t68;
t4 = t60 * t14 + t62 * t22;
t12 = mrSges(6,1) * t78 + mrSges(6,2) * t77;
t71 = t61 ^ 2 + t63 ^ 2;
t70 = qJ(5) * t33;
t46 = pkin(7) + t84;
t69 = qJ(5) + t46;
t47 = -pkin(3) - t83;
t49 = t60 * mrSges(6,2);
t36 = -t62 * mrSges(6,1) + t49;
t3 = t62 * t14 - t22 * t60;
t67 = (Ifges(5,5) + Ifges(6,5)) * t77 + t91 * t32;
t66 = mrSges(5,1) * t60 + mrSges(5,2) * t62;
t53 = Ifges(5,5) * t60;
t52 = Ifges(6,5) * t60;
t51 = Ifges(5,6) * t62;
t50 = Ifges(6,6) * t62;
t42 = Ifges(5,1) * t60 + t81;
t41 = Ifges(6,1) * t60 + t79;
t40 = Ifges(5,2) * t62 + t82;
t39 = Ifges(6,2) * t62 + t80;
t37 = -mrSges(5,1) * t62 + t76;
t35 = -pkin(4) * t62 + t47;
t31 = t69 * t62;
t30 = t69 * t60;
t27 = t33 * mrSges(4,2);
t18 = mrSges(5,1) * t32 - mrSges(5,3) * t77;
t17 = mrSges(6,1) * t32 - mrSges(6,3) * t77;
t16 = -mrSges(5,2) * t32 - mrSges(5,3) * t78;
t15 = -mrSges(6,2) * t32 - t33 * t75;
t13 = t66 * t33;
t9 = pkin(4) * t78 + t20;
t8 = Ifges(5,5) * t32 + (Ifges(5,1) * t62 - t82) * t33;
t7 = Ifges(6,5) * t32 + (Ifges(6,1) * t62 - t80) * t33;
t6 = Ifges(5,6) * t32 + (-Ifges(5,2) * t60 + t81) * t33;
t5 = Ifges(6,6) * t32 + (-Ifges(6,2) * t60 + t79) * t33;
t2 = -t60 * t70 + t4;
t1 = pkin(4) * t32 - t62 * t70 + t3;
t10 = [-0.2e1 * pkin(1) * (-t63 * mrSges(3,1) + t61 * mrSges(3,2)) + t63 * (Ifges(3,4) * t61 + Ifges(3,2) * t63) + t61 * (Ifges(3,1) * t61 + Ifges(3,4) * t63) + t27 * t86 + t13 * t87 + 0.2e1 * t9 * t12 + 0.2e1 * t2 * t15 + 0.2e1 * t4 * t16 + 0.2e1 * t1 * t17 + 0.2e1 * t3 * t18 + Ifges(2,3) + 0.2e1 * t71 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t86 - 0.2e1 * t22 * mrSges(4,3) + Ifges(4,2) * t32 + t67) * t32 + m(6) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t88) + m(4) * (t22 ^ 2 + t48 ^ 2 + t88) + m(3) * (pkin(6) ^ 2 * t71 + pkin(1) ^ 2) + (mrSges(4,3) * t87 + Ifges(4,1) * t33 - 0.2e1 * Ifges(4,4) * t32 + (t7 + t8) * t62 + (t32 * t74 - t5 - t6) * t60) * t33; -t22 * mrSges(4,2) + Ifges(3,5) * t61 + Ifges(3,6) * t63 + t35 * t12 + t47 * t13 + t31 * t15 - t30 * t17 + t9 * t36 + (t52 / 0.2e1 + t50 / 0.2e1 + t53 / 0.2e1 + t51 / 0.2e1 - Ifges(4,6) - mrSges(4,3) * t84) * t32 + (t37 - mrSges(4,1)) * t20 + (-mrSges(3,1) * t61 - mrSges(3,2) * t63) * pkin(6) + (t5 / 0.2e1 + t6 / 0.2e1 + t2 * mrSges(6,3) + t4 * mrSges(5,3) + t46 * t16) * t62 + (t7 / 0.2e1 + t8 / 0.2e1 - t1 * mrSges(6,3) - t3 * mrSges(5,3) - t46 * t18) * t60 + m(6) * (-t1 * t30 + t2 * t31 + t35 * t9) + m(5) * (t47 * t20 + (-t3 * t60 + t4 * t62) * t46) + m(4) * (-t20 * t59 + t22 * t58) * pkin(2) + (-mrSges(4,3) * t83 + Ifges(4,5) + (t41 / 0.2e1 + t42 / 0.2e1) * t62 + (-t39 / 0.2e1 - t40 / 0.2e1) * t60) * t33; 0.2e1 * t35 * t36 + 0.2e1 * t47 * t37 + Ifges(3,3) + Ifges(4,3) + (t31 * t90 + t39 + t40) * t62 + (t30 * t90 + t41 + t42) * t60 + m(6) * (t30 ^ 2 + t31 ^ 2 + t35 ^ 2) + m(5) * (t46 ^ 2 * t72 + t47 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t59 * mrSges(4,1) - t58 * mrSges(4,2)) * pkin(2) + t46 * mrSges(5,3) * t89; t32 * mrSges(4,1) + t27 + (t17 + t18) * t62 + (t15 + t16) * t60 + m(5) * (t3 * t62 + t4 * t60) + m(6) * (t1 * t62 + t2 * t60) + m(4) * t48; m(6) * (-t30 * t62 + t31 * t60); m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t89; t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + t74 * t78 + (m(6) * t1 + t17) * pkin(4) + t67; -t30 * mrSges(6,1) - t31 * mrSges(6,2) + t50 + t51 + t52 + t53 - t66 * t46 + (-m(6) * t30 - t75) * pkin(4); -t76 - t49 + (mrSges(5,1) + mrSges(6,1) + t85) * t62; (0.2e1 * mrSges(6,1) + t85) * pkin(4) + t91; m(6) * t9 + t12; m(6) * t35 + t36; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
