% Calculate joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:48
% EndTime: 2022-01-23 09:31:49
% DurationCPUTime: 0.59s
% Computational Cost: add. (557->138), mult. (1135->185), div. (0->0), fcn. (1056->6), ass. (0->60)
t50 = sin(qJ(4));
t51 = sin(qJ(3));
t52 = cos(qJ(4));
t53 = cos(qJ(3));
t36 = t50 * t53 + t52 * t51;
t48 = sin(pkin(8));
t25 = t36 * t48;
t68 = t48 * t51;
t34 = pkin(3) * t68 + t48 * qJ(2);
t10 = t25 * pkin(4) + t34;
t77 = 0.2e1 * t10;
t76 = 0.2e1 * t34;
t75 = Ifges(5,6) + Ifges(6,6);
t35 = -t50 * t51 + t52 * t53;
t26 = t35 * t48;
t74 = (-Ifges(5,5) - Ifges(6,5)) * t26;
t73 = 0.2e1 * qJ(2);
t65 = -mrSges(5,2) - mrSges(6,2);
t72 = (t52 * mrSges(5,1) + t65 * t50) * pkin(3);
t71 = 2 * mrSges(6,1);
t70 = m(6) * pkin(4);
t69 = t52 * pkin(3);
t49 = cos(pkin(8));
t37 = -pkin(2) * t49 - t48 * pkin(6) - pkin(1);
t19 = t53 * t49 * qJ(2) + t51 * t37;
t12 = -pkin(7) * t68 + t19;
t31 = t53 * t37;
t62 = qJ(2) * t51;
t66 = t53 * t48;
t9 = -pkin(7) * t66 + t31 + (-pkin(3) - t62) * t49;
t6 = t52 * t12 + t50 * t9;
t64 = Ifges(5,3) + Ifges(6,3);
t14 = t49 * mrSges(6,2) - t25 * mrSges(6,3);
t15 = t49 * mrSges(5,2) - t25 * mrSges(5,3);
t63 = t14 + t15;
t61 = Ifges(4,3) + t64;
t16 = -t49 * mrSges(6,1) - t26 * mrSges(6,3);
t5 = -t50 * t12 + t52 * t9;
t2 = -t49 * pkin(4) - t26 * qJ(5) + t5;
t60 = m(6) * t2 + t16;
t58 = t75 * t25 + t74;
t57 = t65 * t36 + (mrSges(5,1) + mrSges(6,1)) * t35;
t3 = -t25 * qJ(5) + t6;
t56 = t5 * mrSges(5,1) + t2 * mrSges(6,1) - t6 * mrSges(5,2) - t3 * mrSges(6,2) - t58;
t55 = pkin(3) ^ 2;
t54 = qJ(2) ^ 2;
t47 = t49 ^ 2;
t46 = t48 ^ 2;
t45 = t50 ^ 2 * t55;
t44 = t48 * mrSges(3,2);
t42 = t46 * t54;
t41 = pkin(4) + t69;
t39 = Ifges(4,5) * t66;
t33 = -t49 * mrSges(4,1) - mrSges(4,3) * t66;
t32 = t49 * mrSges(4,2) - mrSges(4,3) * t68;
t27 = t50 * pkin(3) * t36;
t20 = t26 * mrSges(6,2);
t18 = -t49 * t62 + t31;
t17 = -t49 * mrSges(5,1) - t26 * mrSges(5,3);
t1 = [-0.2e1 * pkin(1) * t44 + t20 * t77 + 0.2e1 * t3 * t14 + 0.2e1 * t6 * t15 + 0.2e1 * t2 * t16 + 0.2e1 * t5 * t17 + 0.2e1 * t18 * t33 + 0.2e1 * t19 * t32 + Ifges(2,3) + (mrSges(5,2) * t76 + (Ifges(6,1) + Ifges(5,1)) * t26) * t26 + m(3) * (pkin(1) ^ 2 + t47 * t54 + t42) + m(4) * (t18 ^ 2 + t19 ^ 2 + t42) + m(5) * (t34 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(6) * (t10 ^ 2 + t2 ^ 2 + t3 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t39 + (Ifges(3,2) + t61) * t49 + t74 + t58) * t49 + ((Ifges(3,1) + (mrSges(4,2) * t73 + Ifges(4,1) * t53) * t53 + (mrSges(4,1) * t73 - 0.2e1 * Ifges(4,4) * t53 + Ifges(4,2) * t51) * t51) * t48 + (-t53 * Ifges(4,5) + 0.2e1 * Ifges(4,6) * t51 + (2 * Ifges(3,4))) * t49) * t48 + (t46 + t47) * mrSges(3,3) * t73 + (mrSges(5,1) * t76 + mrSges(6,1) * t77 + t75 * t49 - 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t26 + (Ifges(5,2) + Ifges(6,2)) * t25) * t25; -m(3) * pkin(1) - t49 * mrSges(3,1) + t51 * t32 + t53 * t33 + t44 + t63 * t36 + (t16 + t17) * t35 + m(6) * (t35 * t2 + t36 * t3) + m(5) * (t35 * t5 + t36 * t6) + m(4) * (t53 * t18 + t51 * t19); m(3) + m(4) * (t51 ^ 2 + t53 ^ 2) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t35 ^ 2 + t36 ^ 2); -Ifges(4,6) * t68 + t18 * mrSges(4,1) - t19 * mrSges(4,2) + t39 + t60 * t41 + ((m(5) * t5 + t17) * t52 + (m(5) * t6 + m(6) * t3 + t63) * t50) * pkin(3) - t61 * t49 + t56; t53 * mrSges(4,1) - t51 * mrSges(4,2) + m(5) * (t35 * t69 + t27) + m(6) * (t41 * t35 + t27) + t57; t41 * t71 + m(6) * (t41 ^ 2 + t45) + m(5) * (t52 ^ 2 * t55 + t45) + 0.2e1 * t72 + t61; t60 * pkin(4) - t64 * t49 + t56; t35 * t70 + t57; t41 * t70 + (pkin(4) + t41) * mrSges(6,1) + t72 + t64; (t71 + t70) * pkin(4) + t64; m(6) * t10 + t25 * mrSges(6,1) + t20; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
