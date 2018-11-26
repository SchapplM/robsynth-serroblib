% Calculate joint inertia matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:08
% EndTime: 2018-11-23 15:37:08
% DurationCPUTime: 0.30s
% Computational Cost: add. (321->123), mult. (572->164), div. (0->0), fcn. (375->6), ass. (0->52)
t34 = sin(qJ(6));
t36 = cos(qJ(6));
t12 = -t36 * mrSges(7,1) + t34 * mrSges(7,2);
t63 = m(7) * pkin(5) + mrSges(6,1) - t12;
t33 = cos(pkin(9));
t24 = -t33 * pkin(1) - pkin(2);
t17 = qJ(4) - t24;
t62 = t17 ^ 2;
t32 = sin(pkin(9));
t22 = t32 * pkin(1) + qJ(3);
t61 = t22 ^ 2;
t60 = 0.2e1 * t17;
t59 = t36 / 0.2e1;
t35 = sin(qJ(5));
t29 = t35 ^ 2;
t37 = cos(qJ(5));
t31 = t37 ^ 2;
t47 = t31 + t29;
t58 = m(6) * t47 + m(5);
t57 = t35 * pkin(5);
t56 = Ifges(7,4) * t34;
t55 = Ifges(7,4) * t36;
t54 = Ifges(7,6) * t35;
t16 = -pkin(7) + t22;
t53 = t16 * t35;
t52 = t34 * t37;
t51 = t36 * t37;
t50 = t37 * mrSges(7,3);
t49 = Ifges(7,5) * t51 + Ifges(7,3) * t35;
t48 = t34 ^ 2 + t36 ^ 2;
t45 = t48 * t37;
t44 = t47 * mrSges(6,3);
t5 = -t37 * pkin(8) + t17 + t57;
t1 = -t34 * t53 + t36 * t5;
t2 = t34 * t5 + t36 * t53;
t43 = -t1 * t34 + t2 * t36;
t8 = -t35 * mrSges(7,2) - t34 * t50;
t9 = t35 * mrSges(7,1) - t36 * t50;
t42 = -t34 * t9 + t36 * t8;
t41 = -t35 * mrSges(6,1) - t37 * mrSges(6,2);
t40 = -mrSges(7,1) * t34 - mrSges(7,2) * t36;
t27 = Ifges(7,5) * t34;
t26 = Ifges(7,6) * t36;
t15 = t16 ^ 2;
t14 = Ifges(7,1) * t34 + t55;
t13 = Ifges(7,2) * t36 + t56;
t11 = t31 * t16;
t10 = t31 * t15;
t6 = -mrSges(7,1) * t52 - mrSges(7,2) * t51;
t4 = Ifges(7,5) * t35 + (Ifges(7,1) * t36 - t56) * t37;
t3 = t54 + (-Ifges(7,2) * t34 + t55) * t37;
t7 = [0.2e1 * t24 * mrSges(4,2) + mrSges(5,3) * t60 + 0.2e1 * t1 * t9 + 0.2e1 * t2 * t8 + Ifges(4,1) + Ifges(5,1) + Ifges(2,3) + Ifges(3,3) + (mrSges(6,1) * t60 + Ifges(6,2) * t35 + t49) * t35 + (mrSges(6,2) * t60 + Ifges(6,1) * t37 - 0.2e1 * Ifges(6,4) * t35 + 0.2e1 * t16 * t6 + t36 * t4 + (-t3 - t54) * t34) * t37 + m(7) * (t1 ^ 2 + t2 ^ 2 + t10) + m(5) * (t61 + t62) + m(4) * (t24 ^ 2 + t61) + m(6) * (t29 * t15 + t10 + t62) + m(3) * (t32 ^ 2 + t33 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t22 + 0.2e1 * (t33 * mrSges(3,1) - t32 * mrSges(3,2)) * pkin(1) - 0.2e1 * t16 * t44; -t35 * t6 + (m(7) * (t43 - t53) + t42) * t37; m(7) * (t48 * t31 + t29) + m(4) + m(3) + t58; -t34 * t8 - t36 * t9 + mrSges(4,2) - mrSges(5,3) + m(7) * (-t36 * t1 - t34 * t2) + m(4) * t24 + (-m(6) / 0.2e1 - m(5) / 0.2e1) * t60 + t41; 0; m(7) * t48 + m(4) + m(5) + m(6); t37 * t6 + mrSges(5,2) + t42 * t35 - t44 + m(7) * (t35 * t43 + t11) + m(6) * (t29 * t16 + t11) + m(5) * t22; m(7) * (-0.1e1 + t48) * t37 * t35; 0; m(7) * (t48 * t29 + t31) + t58; t34 * t4 / 0.2e1 + t3 * t59 + pkin(5) * t6 + (-t16 * mrSges(6,2) + t27 / 0.2e1 + t26 / 0.2e1 - Ifges(6,6)) * t35 + t43 * mrSges(7,3) + (m(7) * t43 + t42) * pkin(8) + (t14 * t59 - t34 * t13 / 0.2e1 + Ifges(6,5) + t63 * t16) * t37; t35 * t12 + m(7) * (pkin(8) * t45 - t57) + mrSges(7,3) * t45 + t41; 0; t63 * t37 + (-mrSges(6,2) + (m(7) * pkin(8) + mrSges(7,3)) * t48) * t35; Ifges(6,3) + t34 * t14 + t36 * t13 - 0.2e1 * pkin(5) * t12 + m(7) * (t48 * pkin(8) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t48 * pkin(8) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2) - Ifges(7,6) * t52 + t49; t6; t12; t40 * t35; pkin(8) * t40 + t26 + t27; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
