% Calculate joint inertia matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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

function Mq = S6RPPPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:05
% EndTime: 2018-11-23 15:39:05
% DurationCPUTime: 0.46s
% Computational Cost: add. (519->163), mult. (849->223), div. (0->0), fcn. (651->6), ass. (0->71)
t47 = sin(qJ(6));
t49 = cos(qJ(6));
t22 = -mrSges(7,1) * t49 + mrSges(7,2) * t47;
t83 = -mrSges(6,1) + t22;
t48 = sin(qJ(5));
t79 = m(7) * t48;
t82 = m(7) * pkin(8) + mrSges(7,3);
t50 = cos(qJ(5));
t81 = -t50 * mrSges(6,2) - pkin(5) * t79 + t48 * t83;
t80 = t49 / 0.2e1;
t77 = t50 * pkin(5);
t43 = sin(pkin(9));
t37 = t43 ^ 2;
t44 = cos(pkin(9));
t38 = t44 ^ 2;
t76 = m(4) + m(5) * (t38 + t37);
t75 = Ifges(7,4) * t47;
t74 = Ifges(7,4) * t49;
t73 = Ifges(7,6) * t50;
t45 = qJ(2) + pkin(3);
t46 = -pkin(1) - qJ(3);
t18 = t43 * t45 + t44 * t46;
t16 = pkin(7) + t18;
t42 = t50 ^ 2;
t72 = t16 * t42;
t40 = t48 ^ 2;
t71 = t40 * t16;
t70 = t43 * t44;
t69 = t47 * t48;
t68 = t47 * t50;
t65 = t48 * t49;
t64 = t49 * t50;
t23 = -t50 * mrSges(6,1) + t48 * mrSges(6,2);
t62 = mrSges(5,1) - t23;
t61 = t47 ^ 2 + t49 ^ 2;
t60 = t40 + t42;
t59 = t61 * t48;
t58 = t60 * mrSges(6,3);
t17 = -t43 * t46 + t44 * t45;
t15 = -pkin(4) - t17;
t3 = -pkin(8) * t48 + t15 - t77;
t1 = -t16 * t68 + t3 * t49;
t2 = t16 * t64 + t3 * t47;
t57 = -t1 * t47 + t2 * t49;
t11 = t43 * t64 - t44 * t47;
t9 = -t43 * t68 - t44 * t49;
t56 = t11 * t49 - t47 * t9;
t10 = t43 * t49 - t44 * t68;
t12 = t43 * t47 + t44 * t64;
t55 = -t10 * t47 + t12 * t49;
t20 = mrSges(7,2) * t50 - mrSges(7,3) * t69;
t21 = -mrSges(7,1) * t50 - mrSges(7,3) * t65;
t54 = t49 * t20 - t47 * t21;
t14 = -mrSges(7,1) * t69 - mrSges(7,2) * t65;
t53 = -t48 * t14 - mrSges(5,2) + t58;
t51 = (qJ(2) ^ 2);
t36 = Ifges(7,5) * t47;
t35 = Ifges(7,6) * t49;
t34 = t40 * t38;
t33 = t40 * t37;
t29 = Ifges(7,5) * t65;
t26 = t40 * t70;
t25 = Ifges(7,1) * t47 + t74;
t24 = Ifges(7,2) * t49 + t75;
t13 = t16 ^ 2;
t8 = t40 * t13;
t7 = -Ifges(7,5) * t50 + (Ifges(7,1) * t49 - t75) * t48;
t6 = -t73 + (-Ifges(7,2) * t47 + t74) * t48;
t5 = t44 * t71;
t4 = t43 * t71;
t19 = [0.2e1 * t17 * mrSges(5,1) - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t18 * mrSges(5,2) - 0.2e1 * t46 * mrSges(4,3) + 0.2e1 * t1 * t21 + 0.2e1 * t15 * t23 + 0.2e1 * t2 * t20 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + Ifges(5,3) + (-t29 + (Ifges(7,3) + Ifges(6,2)) * t50) * t50 + (Ifges(6,1) * t48 + 0.2e1 * Ifges(6,4) * t50 - 0.2e1 * t16 * t14 + t49 * t7 + (-t6 + t73) * t47) * t48 + m(7) * (t1 ^ 2 + t2 ^ 2 + t8) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(6) * (t13 * t42 + t15 ^ 2 + t8) + (m(3) * (pkin(1) ^ 2 + t51)) + m(4) * (t46 ^ 2 + t51) + (2 * (mrSges(3,3) + mrSges(4,1)) * qJ(2)) + 0.2e1 * t16 * t58; -(m(3) * pkin(1)) + t10 * t21 + t12 * t20 + mrSges(3,2) - mrSges(4,3) - t62 * t43 + t53 * t44 + m(7) * (t1 * t10 + t12 * t2 + t5) + m(6) * (t15 * t43 + t44 * t72 + t5) + m(5) * (-t17 * t43 + t18 * t44) + m(4) * t46; m(3) + m(6) * (t38 * t42 + t34 + t37) + m(7) * (t10 ^ 2 + t12 ^ 2 + t34) + t76; m(4) * qJ(2) + t11 * t20 + t9 * t21 + mrSges(4,1) + t62 * t44 + t53 * t43 + m(7) * (t1 * t9 + t11 * t2 + t4) + m(6) * (-t15 * t44 + t43 * t72 + t4) + m(5) * (t17 * t44 + t18 * t43); m(6) * (t26 + (t42 - 0.1e1) * t70) + m(7) * (t10 * t9 + t11 * t12 + t26); m(6) * (t37 * t42 + t33 + t38) + m(7) * (t11 ^ 2 + t9 ^ 2 + t33) + t76; t50 * t14 + (m(7) * (-t16 * t50 + t57) + t54) * t48; (-t44 * t50 + t55) * t79; (-t43 * t50 + t56) * t79; m(5) + m(7) * (t40 * t61 + t42) + m(6) * t60; pkin(5) * t14 + t47 * t7 / 0.2e1 + t6 * t80 + (-t16 * mrSges(6,2) - t36 / 0.2e1 - t35 / 0.2e1 + Ifges(6,6)) * t50 + t57 * mrSges(7,3) + (m(7) * t57 + t54) * pkin(8) + (t25 * t80 - t47 * t24 / 0.2e1 + Ifges(6,5) + (-m(7) * pkin(5) + t83) * t16) * t48; t44 * t81 + t55 * t82; t43 * t81 + t56 * t82; m(7) * (pkin(8) * t59 + t77) - t50 * t22 + mrSges(7,3) * t59 - t23; Ifges(6,3) - 0.2e1 * pkin(5) * t22 + m(7) * (pkin(8) ^ 2 * t61 + pkin(5) ^ 2) + t47 * t25 + t49 * t24 + 0.2e1 * t61 * pkin(8) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,6) * t69 - Ifges(7,3) * t50 + t29; mrSges(7,1) * t10 - mrSges(7,2) * t12; mrSges(7,1) * t9 - mrSges(7,2) * t11; t14; t36 + t35 + (-mrSges(7,1) * t47 - mrSges(7,2) * t49) * pkin(8); Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
