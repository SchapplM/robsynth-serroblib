% Calculate time derivative of joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:49
% EndTime: 2019-12-05 17:11:53
% DurationCPUTime: 1.06s
% Computational Cost: add. (1796->187), mult. (4539->306), div. (0->0), fcn. (4065->8), ass. (0->79)
t67 = sin(qJ(3));
t96 = -pkin(7) - pkin(6);
t59 = t96 * t67;
t71 = cos(qJ(3));
t60 = t96 * t71;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t39 = t70 * t59 + t60 * t66;
t68 = sin(qJ(2));
t76 = t66 * t67 - t70 * t71;
t43 = t76 * t68;
t40 = t66 * t59 - t70 * t60;
t90 = t67 ^ 2 + t71 ^ 2;
t98 = qJD(3) + qJD(4);
t97 = 2 * m(6);
t61 = pkin(3) * t70 + pkin(4);
t69 = cos(qJ(5));
t85 = qJD(5) * t69;
t65 = sin(qJ(5));
t86 = qJD(5) * t65;
t92 = t65 * t66;
t30 = t61 * t85 + (-t66 * t86 + (t69 * t70 - t92) * qJD(4)) * pkin(3);
t94 = t30 * mrSges(6,2);
t91 = t66 * t69;
t72 = cos(qJ(2));
t89 = qJD(2) * t72;
t88 = qJD(3) * t67;
t87 = qJD(3) * t71;
t84 = 0.2e1 * t67;
t83 = pkin(3) * t88;
t82 = t68 * t89;
t52 = t66 * t71 + t67 * t70;
t37 = t98 * t52;
t19 = -t37 * t68 - t76 * t89;
t20 = t43 * t98 - t52 * t89;
t42 = t52 * t68;
t24 = -t42 * t69 + t43 * t65;
t6 = t24 * qJD(5) + t19 * t69 + t20 * t65;
t25 = -t42 * t65 - t43 * t69;
t7 = -t25 * qJD(5) - t19 * t65 + t20 * t69;
t81 = t7 * mrSges(6,1) - t6 * mrSges(6,2);
t62 = -pkin(3) * t71 - pkin(2);
t80 = qJD(3) * t96;
t31 = -t61 * t86 + (-t66 * t85 + (-t65 * t70 - t91) * qJD(4)) * pkin(3);
t29 = t31 * mrSges(6,1);
t79 = t29 - t94;
t32 = -t52 * t65 - t69 * t76;
t36 = t98 * t76;
t10 = t32 * qJD(5) - t36 * t69 - t37 * t65;
t33 = t52 * t69 - t65 * t76;
t11 = -t33 * qJD(5) + t36 * t65 - t37 * t69;
t57 = t67 * t80;
t58 = t71 * t80;
t22 = t39 * qJD(4) + t70 * t57 + t66 * t58;
t12 = -pkin(8) * t37 + t22;
t23 = -qJD(4) * t40 - t57 * t66 + t70 * t58;
t13 = pkin(8) * t36 + t23;
t26 = -pkin(8) * t52 + t39;
t27 = -pkin(8) * t76 + t40;
t14 = t26 * t69 - t27 * t65;
t2 = t14 * qJD(5) + t12 * t69 + t13 * t65;
t15 = t26 * t65 + t27 * t69;
t3 = -t15 * qJD(5) - t12 * t65 + t13 * t69;
t78 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t11;
t77 = -t71 * mrSges(4,1) + t67 * mrSges(4,2);
t75 = t20 * mrSges(5,1) - t19 * mrSges(5,2) + t81;
t74 = (-mrSges(5,1) * t66 - mrSges(5,2) * t70) * qJD(4) * pkin(3);
t73 = t23 * mrSges(5,1) - t22 * mrSges(5,2) - Ifges(5,5) * t36 - Ifges(5,6) * t37 + t78;
t56 = (mrSges(4,1) * t67 + mrSges(4,2) * t71) * qJD(3);
t49 = (-t65 * mrSges(6,1) - t69 * mrSges(6,2)) * qJD(5) * pkin(4);
t45 = pkin(3) * t91 + t61 * t65;
t44 = -pkin(3) * t92 + t61 * t69;
t41 = pkin(4) * t76 + t62;
t38 = mrSges(5,1) * t76 + mrSges(5,2) * t52;
t28 = pkin(4) * t37 + t83;
t17 = mrSges(5,1) * t37 - mrSges(5,2) * t36;
t16 = -mrSges(6,1) * t32 + mrSges(6,2) * t33;
t4 = -mrSges(6,1) * t11 + mrSges(6,2) * t10;
t1 = [0.2e1 * m(6) * (t24 * t7 + t25 * t6 - t82) + 0.2e1 * m(5) * (-t43 * t19 - t42 * t20 - t82) + 0.2e1 * m(4) * (-0.1e1 + t90) * t82; (-t17 - t4 - t56) * t72 + m(6) * (t14 * t7 + t15 * t6 + t2 * t25 + t3 * t24 - t28 * t72) + m(5) * (t40 * t19 + t39 * t20 - t22 * t43 - t23 * t42 - t72 * t83) + (-t24 * t10 + t25 * t11 + t6 * t32 - t7 * t33) * mrSges(6,3) + (-t19 * t76 - t20 * t52 - t42 * t36 + t43 * t37) * mrSges(5,3) + ((-m(4) * pkin(2) + m(5) * t62 + m(6) * t41 - mrSges(3,1) + t16 + t38 + t77) * t68 + (-mrSges(3,2) + (m(4) * pkin(6) + mrSges(4,3)) * t90) * t72) * qJD(2); -0.2e1 * t36 * t52 * Ifges(5,1) + 0.2e1 * t33 * t10 * Ifges(6,1) + 0.2e1 * t76 * Ifges(5,2) * t37 + 0.2e1 * t11 * Ifges(6,2) * t32 - 0.2e1 * pkin(2) * t56 + 0.2e1 * t28 * t16 + 0.2e1 * t62 * t17 + 0.2e1 * t41 * t4 + (-Ifges(4,4) * t67 + pkin(3) * t38) * qJD(3) * t84 + 0.2e1 * m(5) * (t22 * t40 + t23 * t39 + t62 * t83) + (t14 * t3 + t15 * t2 + t28 * t41) * t97 + (0.2e1 * Ifges(4,4) * t71 + (Ifges(4,1) - Ifges(4,2)) * t84) * t87 + 0.2e1 * (t32 * t10 + t11 * t33) * Ifges(6,4) + 0.2e1 * (t36 * t76 - t37 * t52) * Ifges(5,4) + 0.2e1 * (-t14 * t10 + t15 * t11 + t2 * t32 - t3 * t33) * mrSges(6,3) + 0.2e1 * (-t22 * t76 - t23 * t52 + t39 * t36 - t40 * t37) * mrSges(5,3); m(6) * (t24 * t31 + t25 * t30 + t44 * t7 + t45 * t6) + (t68 * t88 - t71 * t89) * mrSges(4,2) + (-t67 * t89 - t68 * t87) * mrSges(4,1) + m(5) * (t19 * t66 + t20 * t70 + (t42 * t66 - t43 * t70) * qJD(4)) * pkin(3) + t75; m(6) * (t14 * t31 + t15 * t30 + t2 * t45 + t3 * t44) + (Ifges(4,5) * t71 - Ifges(4,6) * t67 + t77 * pkin(6)) * qJD(3) + (-t44 * t10 + t45 * t11 + t30 * t32 - t31 * t33) * mrSges(6,3) + (m(5) * (t22 * t66 + t23 * t70 + (-t39 * t66 + t40 * t70) * qJD(4)) + (t70 * t36 - t66 * t37 + (t52 * t66 - t70 * t76) * qJD(4)) * mrSges(5,3)) * pkin(3) + t73; (t30 * t45 + t31 * t44) * t97 - 0.2e1 * t94 + 0.2e1 * t29 + 0.2e1 * t74; m(6) * (t6 * t65 + t69 * t7 + (-t24 * t65 + t25 * t69) * qJD(5)) * pkin(4) + t75; (m(6) * (t2 * t65 + t3 * t69 + (-t14 * t65 + t15 * t69) * qJD(5)) + (-t69 * t10 + t65 * t11 + (t32 * t69 + t33 * t65) * qJD(5)) * mrSges(6,3)) * pkin(4) + t73; t74 + (m(6) * (t30 * t65 + t31 * t69 - t44 * t86 + t45 * t85) - mrSges(6,2) * t85 - mrSges(6,1) * t86) * pkin(4) + t79; 0.2e1 * t49; t81; t78; t79; t49; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
