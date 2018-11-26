% Calculate joint inertia matrix for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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

function Mq = S6RPPPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:31
% EndTime: 2018-11-23 15:38:31
% DurationCPUTime: 0.46s
% Computational Cost: add. (480->148), mult. (799->207), div. (0->0), fcn. (582->6), ass. (0->66)
t39 = sin(qJ(6));
t41 = cos(qJ(6));
t21 = -t41 * mrSges(7,1) + t39 * mrSges(7,2);
t76 = mrSges(6,1) - t21;
t42 = cos(qJ(5));
t46 = -mrSges(7,1) * t39 - mrSges(7,2) * t41;
t10 = t46 * t42;
t40 = sin(qJ(5));
t34 = t40 ^ 2;
t36 = t42 ^ 2;
t54 = t36 + t34;
t50 = t54 * mrSges(6,3);
t75 = t42 * t10 + mrSges(5,2) - t50;
t74 = -t40 * mrSges(6,2) + t76 * t42;
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t43 = -pkin(1) - pkin(2);
t20 = t38 * qJ(2) + t37 * t43;
t13 = -qJ(4) + t20;
t73 = t13 ^ 2;
t72 = t39 / 0.2e1;
t71 = m(7) * t42;
t70 = pkin(5) * t42;
t69 = t40 * pkin(5);
t68 = m(6) * t54 + m(5);
t67 = Ifges(7,4) * t39;
t66 = Ifges(7,4) * t41;
t65 = Ifges(7,5) * t41;
t64 = t36 * t38;
t63 = t37 * t13;
t62 = t39 * t40;
t60 = t40 * t41;
t58 = t42 * mrSges(7,3);
t22 = -t40 * mrSges(6,1) - t42 * mrSges(6,2);
t56 = t22 - mrSges(5,3);
t55 = t39 ^ 2 + t41 ^ 2;
t53 = t55 * t40;
t52 = t55 * t42;
t51 = t54 * t38;
t19 = -t37 * qJ(2) + t38 * t43;
t15 = pkin(3) - t19;
t12 = pkin(7) + t15;
t3 = t42 * pkin(8) + t13 - t69;
t1 = -t12 * t62 + t41 * t3;
t2 = t12 * t60 + t39 * t3;
t48 = -t1 * t39 + t2 * t41;
t8 = t41 * t37 + t38 * t62;
t9 = t39 * t37 - t38 * t60;
t47 = -t39 * t8 + t41 * t9;
t17 = t40 * mrSges(7,2) + t39 * t58;
t18 = -t40 * mrSges(7,1) + t41 * t58;
t45 = t41 * t17 - t39 * t18;
t32 = t38 ^ 2;
t31 = t37 ^ 2;
t30 = Ifges(7,5) * t39;
t29 = Ifges(7,6) * t41;
t27 = t36 * t32;
t26 = Ifges(7,6) * t39 * t42;
t24 = Ifges(7,1) * t39 + t66;
t23 = Ifges(7,2) * t41 + t67;
t11 = t12 ^ 2;
t7 = t36 * t12;
t6 = t36 * t11;
t5 = -Ifges(7,5) * t40 + (-Ifges(7,1) * t41 + t67) * t42;
t4 = -Ifges(7,6) * t40 + (Ifges(7,2) * t39 - t66) * t42;
t14 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t19 * mrSges(4,1) + 0.2e1 * t20 * mrSges(4,2) - 0.2e1 * t15 * mrSges(5,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t1 * t18 + 0.2e1 * t2 * t17 + Ifges(5,1) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (-t26 + (Ifges(7,3) + Ifges(6,2)) * t40) * t40 + (Ifges(6,1) * t42 - 0.2e1 * t12 * t10 + t39 * t4 - t41 * t5 + (-(2 * Ifges(6,4)) + t65) * t40) * t42 + m(6) * (t34 * t11 + t6 + t73) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t6) + m(5) * (t15 ^ 2 + t73) + m(4) * (t19 ^ 2 + t20 ^ 2) + 0.2e1 * t12 * t50 + 0.2e1 * t56 * t13; -m(3) * pkin(1) + t9 * t17 + t8 * t18 - mrSges(3,1) + (mrSges(4,2) + t56) * t37 + (-mrSges(4,1) + t75) * t38 + m(7) * (t8 * t1 - t12 * t64 + t9 * t2) + m(6) * (-t12 * t51 + t63) + m(5) * (-t38 * t15 + t63) + m(4) * (t38 * t19 + t37 * t20); m(3) + m(6) * (t34 * t32 + t27 + t31) + m(7) * (t8 ^ 2 + t9 ^ 2 + t27) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t31 + t32); t40 * t10 + (m(7) * (-t12 * t40 + t48) + t45) * t42; (t38 * t40 + t47) * t71; m(4) + m(7) * (t55 * t36 + t34) + t68; t45 * t40 + m(7) * (t48 * t40 + t7) + m(6) * (t34 * t12 + t7) + m(5) * t15 - t75; -m(5) * t38 - m(6) * t51 + m(7) * (t47 * t40 - t64); (-0.1e1 + t55) * t40 * t71; m(7) * (t55 * t34 + t36) + t68; t5 * t72 + t41 * t4 / 0.2e1 - pkin(5) * t10 + (-t12 * mrSges(6,2) - t30 / 0.2e1 - t29 / 0.2e1 + Ifges(6,6)) * t40 + t48 * mrSges(7,3) + (m(7) * t48 + t45) * pkin(8) + (t23 * t72 - t41 * t24 / 0.2e1 - Ifges(6,5) + (m(7) * pkin(5) + t76) * t12) * t42; (m(7) * pkin(8) + mrSges(7,3)) * t47 + (-m(7) * t70 - t74) * t38; t40 * t21 + m(7) * (pkin(8) * t52 - t69) + mrSges(7,3) * t52 + t22; m(7) * (pkin(8) * t53 + t70) + mrSges(7,3) * t53 + t74; Ifges(6,3) - 0.2e1 * pkin(5) * t21 + t41 * t23 + t39 * t24 + m(7) * (t55 * pkin(8) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t55 * pkin(8) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2) - Ifges(7,3) * t40 - t42 * t65 + t26; t8 * mrSges(7,1) - t9 * mrSges(7,2); t10; t46 * t40; t46 * pkin(8) + t29 + t30; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
