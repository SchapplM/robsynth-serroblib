% Calculate joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:35
% DurationCPUTime: 0.50s
% Computational Cost: add. (615->182), mult. (1208->264), div. (0->0), fcn. (1063->6), ass. (0->78)
t64 = (-pkin(1) - pkin(6));
t91 = -2 * t64;
t60 = sin(qJ(3));
t62 = cos(qJ(4));
t63 = cos(qJ(3));
t81 = t62 * t63;
t90 = Ifges(5,5) * t81 + Ifges(5,3) * t60;
t89 = 2 * qJ(2);
t88 = t62 / 0.2e1;
t87 = -pkin(8) - pkin(7);
t59 = sin(qJ(4));
t86 = Ifges(5,4) * t59;
t85 = Ifges(5,4) * t62;
t84 = Ifges(5,6) * t60;
t83 = t59 * t63;
t82 = t60 * t64;
t39 = -t62 * mrSges(5,1) + t59 * mrSges(5,2);
t80 = mrSges(4,1) - t39;
t38 = t60 * pkin(3) - t63 * pkin(7) + qJ(2);
t18 = t59 * t38 + t62 * t82;
t79 = t59 ^ 2 + t62 ^ 2;
t54 = t60 ^ 2;
t56 = t63 ^ 2;
t78 = t56 + t54;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t35 = t58 * t62 + t61 * t59;
t25 = t35 * t63;
t34 = -t58 * t59 + t61 * t62;
t27 = t34 * t63;
t77 = Ifges(6,5) * t27 - Ifges(6,6) * t25 + Ifges(6,3) * t60;
t76 = t79 * mrSges(5,3);
t75 = t78 * mrSges(4,3);
t24 = t35 * t60;
t26 = t34 * t60;
t74 = -t24 * mrSges(6,1) - t26 * mrSges(6,2);
t73 = t59 * mrSges(5,1) + t62 * mrSges(5,2);
t32 = t62 * t38;
t17 = -t59 * t82 + t32;
t72 = -t17 * t59 + t18 * t62;
t36 = -t60 * mrSges(5,2) - mrSges(5,3) * t83;
t37 = t60 * mrSges(5,1) - mrSges(5,3) * t81;
t71 = t62 * t36 - t59 * t37;
t42 = t87 * t59;
t43 = t87 * t62;
t13 = t61 * t42 + t58 * t43;
t14 = t58 * t42 - t61 * t43;
t29 = Ifges(6,6) * t34;
t30 = Ifges(6,5) * t35;
t70 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t29 + t30;
t7 = -pkin(8) * t81 + t32 + (-t59 * t64 + pkin(4)) * t60;
t8 = -pkin(8) * t83 + t18;
t2 = -t58 * t8 + t61 * t7;
t3 = t58 * t7 + t61 * t8;
t69 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t77;
t68 = (t61 * mrSges(6,1) - t58 * mrSges(6,2)) * pkin(4);
t65 = qJ(2) ^ 2;
t57 = t64 ^ 2;
t52 = Ifges(5,5) * t59;
t51 = Ifges(5,6) * t62;
t48 = t56 * t64;
t47 = t56 * t57;
t46 = -t62 * pkin(4) - pkin(3);
t41 = Ifges(5,1) * t59 + t85;
t40 = Ifges(5,2) * t62 + t86;
t33 = (pkin(4) * t59 - t64) * t63;
t28 = t73 * t63;
t23 = Ifges(5,5) * t60 + (Ifges(5,1) * t62 - t86) * t63;
t22 = t84 + (-Ifges(5,2) * t59 + t85) * t63;
t16 = t60 * mrSges(6,1) - t27 * mrSges(6,3);
t15 = -t60 * mrSges(6,2) - t25 * mrSges(6,3);
t11 = Ifges(6,1) * t35 + Ifges(6,4) * t34;
t10 = Ifges(6,4) * t35 + Ifges(6,2) * t34;
t9 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t6 = t25 * mrSges(6,1) + t27 * mrSges(6,2);
t5 = Ifges(6,1) * t27 - Ifges(6,4) * t25 + Ifges(6,5) * t60;
t4 = Ifges(6,4) * t27 - Ifges(6,2) * t25 + Ifges(6,6) * t60;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t89) + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + 0.2e1 * t17 * t37 + 0.2e1 * t18 * t36 - t25 * t4 + t27 * t5 + 0.2e1 * t33 * t6 + Ifges(3,1) + Ifges(2,3) + t75 * t91 + (mrSges(4,1) * t89 + Ifges(4,2) * t60 + t77 + t90) * t60 + ((mrSges(4,2) * t89) + Ifges(4,1) * t63 - 0.2e1 * Ifges(4,4) * t60 + t62 * t23 + t28 * t91 + (-t22 - t84) * t59) * t63 + m(6) * (t2 ^ 2 + t3 ^ 2 + t33 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t47) + (m(3) * (pkin(1) ^ 2 + t65)) + m(4) * (t54 * t57 + t47 + t65); -(m(3) * pkin(1)) + t26 * t15 - t24 * t16 + mrSges(3,2) + (-t28 - t6) * t63 + t71 * t60 - t75 + m(6) * (-t24 * t2 + t26 * t3 - t63 * t33) + m(5) * (t60 * t72 + t48) + m(4) * (t54 * t64 + t48); m(3) + m(4) * t78 + m(5) * (t54 * t79 + t56) + m(6) * (t24 ^ 2 + t26 ^ 2 + t56); t22 * t88 + t59 * t23 / 0.2e1 + t34 * t4 / 0.2e1 + t35 * t5 / 0.2e1 + t46 * t6 - t25 * t10 / 0.2e1 + t27 * t11 / 0.2e1 - pkin(3) * t28 + t33 * t9 + t14 * t15 + t13 * t16 + m(6) * (t13 * t2 + t14 * t3 + t46 * t33) + (t52 / 0.2e1 + t51 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1 - Ifges(4,6) - (t64 * mrSges(4,2))) * t60 + (-t2 * t35 + t3 * t34) * mrSges(6,3) + t72 * mrSges(5,3) + (m(5) * t72 + t71) * pkin(7) + (Ifges(4,5) - t59 * t40 / 0.2e1 + t41 * t88 + (m(5) * pkin(3) + t80) * t64) * t63; (t24 * t35 + t26 * t34) * mrSges(6,3) + (-t9 + t80) * t63 + (-mrSges(4,2) + t76) * t60 + m(5) * (pkin(7) * t60 * t79 + pkin(3) * t63) + m(6) * (-t13 * t24 + t14 * t26 - t46 * t63); -0.2e1 * pkin(3) * t39 + t34 * t10 + t35 * t11 + t62 * t40 + t59 * t41 + 0.2e1 * t46 * t9 + Ifges(4,3) + m(6) * (t13 ^ 2 + t14 ^ 2 + t46 ^ 2) + m(5) * (pkin(7) ^ 2 * t79 + pkin(3) ^ 2) + 0.2e1 * (-t13 * t35 + t14 * t34) * mrSges(6,3) + 0.2e1 * pkin(7) * t76; -Ifges(5,6) * t83 + t17 * mrSges(5,1) - t18 * mrSges(5,2) + (m(6) * (t2 * t61 + t3 * t58) + t58 * t15 + t61 * t16) * pkin(4) + t69 + t90; -t73 * t60 + m(6) * (-t24 * t61 + t26 * t58) * pkin(4) + t74; t51 + t52 - t73 * pkin(7) + (m(6) * (t13 * t61 + t14 * t58) + (t58 * t34 - t61 * t35) * mrSges(6,3)) * pkin(4) + t70; Ifges(5,3) + Ifges(6,3) + m(6) * (t58 ^ 2 + t61 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t68; t69; t74; t70; Ifges(6,3) + t68; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
