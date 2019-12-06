% Calculate time derivative of joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:19
% EndTime: 2019-12-05 18:30:20
% DurationCPUTime: 0.31s
% Computational Cost: add. (662->90), mult. (1626->138), div. (0->0), fcn. (1120->8), ass. (0->70)
t83 = Ifges(6,1) - Ifges(6,2);
t53 = cos(qJ(2));
t40 = t53 * pkin(1) + pkin(2);
t47 = cos(pkin(9));
t46 = sin(pkin(9));
t50 = sin(qJ(2));
t68 = t46 * t50;
t58 = -pkin(1) * t68 + t47 * t40;
t22 = pkin(3) + t58;
t67 = t47 * t50;
t28 = pkin(1) * t67 + t46 * t40;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t11 = t49 * t22 + t52 * t28;
t39 = t47 * pkin(2) + pkin(3);
t79 = pkin(2) * t46;
t29 = t49 * t39 + t52 * t79;
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t66 = t51 * mrSges(6,1);
t34 = t48 * mrSges(6,2) - t66;
t64 = pkin(1) * qJD(2);
t25 = (-t46 * t53 - t67) * t64;
t26 = (t47 * t53 - t68) * t64;
t6 = t11 * qJD(4) - t52 * t25 + t49 * t26;
t1 = t6 * t34;
t10 = t52 * t22 - t49 * t28;
t5 = t10 * qJD(4) + t49 * t25 + t52 * t26;
t44 = t48 ^ 2;
t76 = mrSges(6,3) * t44;
t2 = t5 * t76;
t45 = t51 ^ 2;
t75 = mrSges(6,3) * t45;
t3 = t5 * t75;
t4 = t6 * mrSges(5,1);
t62 = qJD(5) * t51;
t63 = qJD(5) * t48;
t33 = -mrSges(6,1) * t63 - mrSges(6,2) * t62;
t8 = -pkin(4) - t10;
t7 = t8 * t33;
t82 = t1 + t2 + t3 - t4 - t7;
t81 = 2 * m(5);
t80 = 2 * m(6);
t78 = pkin(4) * t33;
t77 = t5 * mrSges(5,2);
t72 = Ifges(6,6) * t48;
t27 = t52 * t39 - t49 * t79;
t20 = t27 * qJD(4);
t71 = t20 * mrSges(5,2);
t70 = t25 * mrSges(4,1);
t69 = t26 * mrSges(4,2);
t65 = t44 + t45;
t61 = t65 * t5;
t60 = t65 * t20;
t24 = pkin(8) + t29;
t59 = t65 * t24;
t57 = -mrSges(6,1) * t48 - mrSges(6,2) * t51;
t56 = (-0.2e1 * Ifges(6,4) * t48 + t83 * t51) * t63 + (0.2e1 * Ifges(6,4) * t51 + t83 * t48) * t62;
t55 = (-mrSges(3,1) * t50 - mrSges(3,2) * t53) * t64;
t23 = -pkin(4) - t27;
t12 = t23 * t33;
t21 = t29 * qJD(4);
t13 = t21 * t34;
t14 = t20 * t76;
t15 = t20 * t75;
t16 = t21 * mrSges(5,1);
t54 = -t12 + t13 + t14 + t15 - t16 + t56;
t43 = Ifges(6,5) * t62;
t9 = pkin(8) + t11;
t17 = [0.2e1 * t70 - 0.2e1 * t69 - 0.2e1 * t77 + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 - 0.2e1 * t7 + 0.2e1 * t55 + (t8 * t6 + t9 * t61) * t80 + (-t10 * t6 + t11 * t5) * t81 + 0.2e1 * m(4) * (t58 * t25 + t28 * t26) + t56; t54 + m(6) * (t21 * t8 + t23 * t6 + t5 * t59 + t9 * t60) + m(5) * (-t21 * t10 + t20 * t11 - t27 * t6 + t29 * t5) + (-t20 - t5) * mrSges(5,2) + m(4) * (t25 * t47 + t26 * t46) * pkin(2) + t55 + t70 - t69 + t82; -0.2e1 * t71 - 0.2e1 * t12 + 0.2e1 * t13 + 0.2e1 * t14 + 0.2e1 * t15 - 0.2e1 * t16 + (t20 * t59 + t23 * t21) * t80 + (t29 * t20 - t27 * t21) * t81 + t56; 0; 0; 0; t78 + m(6) * (-pkin(4) * t6 + pkin(8) * t61) - t77 + t56 + t82; t78 + m(6) * (-pkin(4) * t21 + pkin(8) * t60) - t71 + t54; 0; t56 + 0.2e1 * t78; t43 + t57 * t5 + (-t9 * t66 + (mrSges(6,2) * t9 - Ifges(6,6)) * t48) * qJD(5); t43 + t57 * t20 + (t34 * t24 - t72) * qJD(5); t33; t43 + (t34 * pkin(8) - t72) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
