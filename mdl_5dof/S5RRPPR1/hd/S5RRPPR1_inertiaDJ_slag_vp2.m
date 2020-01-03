% Calculate time derivative of joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:44
% EndTime: 2020-01-03 11:55:46
% DurationCPUTime: 0.41s
% Computational Cost: add. (677->106), mult. (1536->162), div. (0->0), fcn. (1204->8), ass. (0->75)
t41 = sin(pkin(9));
t43 = cos(pkin(9));
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t29 = t47 * t41 + t45 * t43;
t25 = t29 * qJD(5);
t28 = -t45 * t41 + t47 * t43;
t88 = t28 * t25;
t24 = t28 * qJD(5);
t87 = t29 * t24;
t39 = t41 ^ 2;
t40 = t43 ^ 2;
t67 = t39 + t40;
t86 = mrSges(5,3) * t67 * qJD(4);
t84 = 2 * m(5);
t83 = 2 * m(6);
t82 = t43 * pkin(4);
t13 = t25 * mrSges(6,1) + t24 * mrSges(6,2);
t48 = cos(qJ(2));
t37 = t48 * pkin(1) + pkin(2);
t42 = sin(pkin(8));
t44 = cos(pkin(8));
t46 = sin(qJ(2));
t52 = -t42 * t46 * pkin(1) + t44 * t37;
t49 = -pkin(3) - t52;
t17 = t49 - t82;
t81 = t17 * t13;
t66 = pkin(1) * qJD(2);
t57 = t48 * t66;
t58 = t46 * t66;
t21 = -t42 * t58 + t44 * t57;
t18 = qJD(4) + t21;
t80 = t18 * mrSges(5,3);
t70 = t44 * t46;
t20 = (t42 * t48 + t70) * t66;
t79 = t20 * mrSges(4,1);
t78 = t20 * (-t28 * mrSges(6,1) + t29 * mrSges(6,2));
t77 = t20 * (-t43 * mrSges(5,1) + t41 * mrSges(5,2));
t76 = t21 * mrSges(4,2);
t75 = t24 * mrSges(6,3);
t74 = t25 * mrSges(6,3);
t73 = t28 * mrSges(6,3);
t72 = t29 * mrSges(6,3);
t56 = -t44 * pkin(2) - pkin(3);
t30 = t56 - t82;
t71 = t30 * t13;
t69 = Ifges(6,5) * t24 - Ifges(6,6) * t25;
t68 = pkin(1) * t70 + t42 * t37;
t19 = qJ(4) + t68;
t15 = (-pkin(7) - t19) * t41;
t38 = t43 * pkin(7);
t16 = t43 * t19 + t38;
t3 = t47 * t15 - t45 * t16;
t1 = qJD(5) * t3 + t18 * t28;
t64 = t1 * t73;
t4 = t45 * t15 + t47 * t16;
t2 = -qJD(5) * t4 - t18 * t29;
t63 = t2 * t72;
t62 = t3 * t75;
t61 = t4 * t74;
t60 = t39 * t80;
t59 = t40 * t80;
t55 = 0.2e1 * Ifges(6,1) * t87 - 0.2e1 * Ifges(6,2) * t88 + 0.2e1 * (t28 * t24 - t25 * t29) * Ifges(6,4);
t54 = t67 * t19;
t36 = t42 * pkin(2) + qJ(4);
t53 = t67 * t36;
t51 = mrSges(3,1) * t58;
t50 = mrSges(3,2) * t57;
t26 = (-pkin(7) - t36) * t41;
t27 = t43 * t36 + t38;
t11 = t47 * t26 - t45 * t27;
t12 = t45 * t26 + t47 * t27;
t6 = -qJD(4) * t29 - qJD(5) * t12;
t5 = qJD(4) * t28 + qJD(5) * t11;
t7 = [t55 - 0.2e1 * t50 - 0.2e1 * t51 + (t18 * t54 + t20 * t49) * t84 + 0.2e1 * m(4) * (-t20 * t52 + t21 * t68) + (t4 * t1 + t17 * t20 + t3 * t2) * t83 + 0.2e1 * t59 - 0.2e1 * t61 + 0.2e1 * t64 - 0.2e1 * t63 + 0.2e1 * t60 + 0.2e1 * t77 + 0.2e1 * t81 + 0.2e1 * t78 - 0.2e1 * t79 - 0.2e1 * t76 - 0.2e1 * t62; t55 - t50 - t51 + m(5) * (qJD(4) * t54 + t18 * t53 + t20 * t56) + m(4) * (-t44 * t20 + t21 * t42) * pkin(2) + m(6) * (t12 * t1 + t11 * t2 + t30 * t20 + t6 * t3 + t5 * t4) + t59 - t12 * t74 - t61 + t5 * t73 + t64 - t6 * t72 - t63 + t60 + t71 + t77 + t81 + t78 - t79 - t76 - t11 * t75 - t62 + t86; 0.2e1 * t71 + (t11 * t6 + t12 * t5) * t83 + t53 * t84 * qJD(4) + t55 + 0.2e1 * t86 + 0.2e1 * (-t11 * t24 - t12 * t25 + t5 * t28 - t6 * t29) * mrSges(6,3); m(6) * (t1 * t29 + t2 * t28 + t4 * t24 - t3 * t25); m(6) * (-t11 * t25 + t12 * t24 + t6 * t28 + t5 * t29); (t87 - t88) * t83; 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t20 + t13; t13; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t69; t6 * mrSges(6,1) - t5 * mrSges(6,2) + t69; -t13; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
