% Calculate time derivative of joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:37
% DurationCPUTime: 0.52s
% Computational Cost: add. (676->123), mult. (1496->180), div. (0->0), fcn. (1083->8), ass. (0->69)
t34 = cos(pkin(8)) * pkin(1) + pkin(2);
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t78 = pkin(1) * sin(pkin(8));
t54 = t45 * t34 - t43 * t78;
t17 = t54 * qJD(3);
t14 = qJD(4) + t17;
t67 = t43 * t34 + t45 * t78;
t19 = qJ(4) + t67;
t87 = qJ(4) * t14 + qJD(4) * t19;
t39 = sin(pkin(9));
t35 = t39 ^ 2;
t41 = cos(pkin(9));
t36 = t41 ^ 2;
t86 = (t35 + t36) * mrSges(5,3);
t85 = qJD(4) + t14;
t84 = 2 * m(5);
t83 = 2 * m(6);
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t61 = qJD(5) * t39;
t20 = (mrSges(6,1) * t44 - mrSges(6,2) * t42) * t61;
t82 = 0.2e1 * t20;
t51 = mrSges(6,1) * t42 + mrSges(6,2) * t44;
t24 = t51 * t39;
t81 = 0.2e1 * t24;
t76 = mrSges(6,3) * t39;
t25 = t41 * mrSges(6,2) - t42 * t76;
t80 = 0.2e1 * t25;
t26 = -t41 * mrSges(6,1) - t44 * t76;
t79 = 0.2e1 * t26;
t77 = t87 * t35;
t75 = Ifges(6,4) * t42;
t74 = Ifges(6,4) * t44;
t73 = t17 * mrSges(4,2);
t72 = t19 * t14;
t71 = t41 * t42;
t70 = t41 * t44;
t69 = t42 * (-Ifges(6,2) * t44 - t75) * t61;
t64 = qJ(4) * t41;
t62 = qJD(4) * t41;
t60 = qJD(5) * t42;
t59 = qJD(5) * t44;
t58 = qJ(4) * qJD(4);
t57 = 0.2e1 * mrSges(6,3);
t18 = t67 * qJD(3);
t56 = (-t41 * mrSges(5,1) + t39 * mrSges(5,2) - mrSges(4,1)) * t18;
t21 = (-Ifges(6,5) * t42 - Ifges(6,6) * t44) * t61;
t55 = -t41 * t21 + t35 * (-Ifges(6,1) * t42 - t74) * t59;
t53 = t25 * t59 - t26 * t60;
t28 = -t41 * pkin(4) - t39 * pkin(7) - pkin(3);
t8 = t28 - t54;
t3 = -t19 * t71 + t44 * t8;
t4 = t19 * t70 + t42 * t8;
t52 = t3 * t42 - t4 * t44;
t12 = t44 * t28 - t42 * t64;
t13 = t42 * t28 + t44 * t64;
t50 = t12 * t42 - t13 * t44;
t49 = -(-Ifges(6,6) * t41 + (-Ifges(6,2) * t42 + t74) * t39) * t44 - (-Ifges(6,5) * t41 + (Ifges(6,1) * t44 - t75) * t39) * t42;
t48 = 0.2e1 * t86;
t47 = -t25 * t60 - t26 * t59;
t46 = -t41 * t20 + (-t42 ^ 2 - t44 ^ 2) * t35 * qJD(5) * mrSges(6,3);
t33 = t35 * t58;
t7 = -qJD(5) * t13 - t42 * t62;
t6 = qJD(5) * t12 + t44 * t62;
t5 = t35 * t72;
t2 = -qJD(5) * t4 - t14 * t71 + t44 * t18;
t1 = qJD(5) * t3 + t14 * t70 + t42 * t18;
t9 = [-0.2e1 * t73 + t1 * t80 + t2 * t79 + 0.2e1 * t56 + t14 * t48 + (t4 * t1 + t3 * t2 + t5) * t83 + (t36 * t72 + t5 + (-pkin(3) - t54) * t18) * t84 + 0.2e1 * m(4) * (t67 * t17 - t54 * t18) + (t14 * t81 + t19 * t82 - t69 + (t52 * t57 + t49) * qJD(5)) * t39 + t55; (m(6) * (t1 * t44 - t14 * t41 - t2 * t42 - t3 * t59 - t4 * t60) + t47) * t39 + t46; 0; -t73 + (t7 + t2) * t26 + (t6 + t1) * t25 + t56 + m(5) * (-pkin(3) * t18 + t87 * t36 + t77) + m(6) * (t13 * t1 + t12 * t2 + t7 * t3 + t6 * t4 + t77) + (-t69 + t85 * t24 + (qJ(4) + t19) * t20 + (((-t13 - t4) * t44 + (t12 + t3) * t42) * mrSges(6,3) + t49) * qJD(5)) * t39 + t55 + t85 * t86; (m(6) * (-t12 * t59 - t13 * t60 - t42 * t7 + t44 * t6 - t62) + t47) * t39 + t46; t6 * t80 + t7 * t79 + qJD(4) * t48 + (t12 * t7 + t13 * t6 + t33) * t83 + (t36 * t58 + t33) * t84 + (qJ(4) * t82 + qJD(4) * t81 - t69 + (t50 * t57 + t49) * qJD(5)) * t39 + t55; m(6) * (-t52 * qJD(5) + t42 * t1 + t44 * t2) + m(5) * t18 + t53; 0; m(6) * (-qJD(5) * t50 + t42 * t6 + t44 * t7) + t53; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t21; -t20; t7 * mrSges(6,1) - t6 * mrSges(6,2) + t21; -t51 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
