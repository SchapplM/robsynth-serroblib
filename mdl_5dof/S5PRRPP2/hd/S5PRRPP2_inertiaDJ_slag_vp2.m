% Calculate time derivative of joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:44
% EndTime: 2019-12-05 16:08:47
% DurationCPUTime: 0.62s
% Computational Cost: add. (473->111), mult. (1297->166), div. (0->0), fcn. (1020->6), ass. (0->62)
t69 = (mrSges(6,2) + mrSges(5,3));
t82 = 2 * t69;
t81 = m(5) + m(6);
t39 = sin(pkin(8));
t42 = cos(qJ(3));
t40 = sin(qJ(3));
t66 = cos(pkin(8));
t53 = t66 * t40;
t26 = t39 * t42 + t53;
t23 = t26 * qJD(3);
t52 = t66 * t42;
t70 = t39 * t40;
t48 = t52 - t70;
t24 = t48 * qJD(3);
t10 = t23 * mrSges(6,1) - t24 * mrSges(6,3);
t11 = t23 * mrSges(5,1) + t24 * mrSges(5,2);
t80 = -t10 - t11;
t67 = t40 ^ 2 + t42 ^ 2;
t41 = sin(qJ(2));
t63 = qJD(3) * t41;
t43 = cos(qJ(2));
t65 = qJD(2) * t43;
t79 = t40 * t63 - t42 * t65;
t73 = t39 * pkin(3);
t33 = qJ(5) + t73;
t78 = m(6) * t33 + mrSges(6,3);
t54 = t66 * pkin(3);
t35 = -t54 - pkin(4);
t75 = m(5) * pkin(3);
t77 = m(6) * t35 - t66 * t75 - mrSges(5,1) - mrSges(6,1);
t76 = t39 * t75 - mrSges(5,2) + t78;
t72 = mrSges(4,1) * t40;
t71 = t24 * mrSges(6,2);
t68 = -qJ(4) - pkin(6);
t64 = qJD(3) * t40;
t62 = qJD(3) * t42;
t61 = 0.2e1 * t40;
t60 = pkin(3) * t64;
t59 = mrSges(4,1) * t62;
t56 = t41 * t65;
t36 = -t42 * pkin(3) - pkin(2);
t32 = t68 * t42;
t15 = -t39 * t32 - t68 * t53;
t16 = -t66 * t32 + t68 * t70;
t51 = qJD(3) * t68;
t22 = t42 * qJD(4) + t40 * t51;
t45 = -t40 * qJD(4) + t42 * t51;
t8 = t39 * t22 - t66 * t45;
t9 = t66 * t22 + t39 * t45;
t55 = t15 * t8 + t16 * t9;
t50 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t18 = t26 * t41;
t19 = t48 * t41;
t6 = t79 * t39 - t52 * t63 - t53 * t65;
t7 = -t26 * t63 + t48 * t65;
t49 = -t15 * t6 + t16 * t7 + t8 * t18 + t9 * t19;
t28 = (mrSges(4,2) * t42 + t72) * qJD(3);
t14 = -mrSges(5,1) * t48 + t26 * mrSges(5,2);
t13 = -mrSges(6,1) * t48 - t26 * mrSges(6,3);
t12 = -pkin(4) * t48 - t26 * qJ(5) + t36;
t5 = t23 * pkin(4) - t24 * qJ(5) - t26 * qJD(5) + t60;
t1 = [0.2e1 * m(4) * (-0.1e1 + t67) * t56 + 0.2e1 * t81 * (-t18 * t6 + t19 * t7 - t56); (-t28 + t80) * t43 + m(5) * (-t43 * t60 + t49) + m(6) * (-t5 * t43 + t49) + ((-m(4) * pkin(2) + m(5) * t36 + m(6) * t12 - t42 * mrSges(4,1) + t40 * mrSges(4,2) - mrSges(3,1) + t13 + t14) * t41 + (-mrSges(3,2) + (m(4) * pkin(6) + mrSges(4,3)) * t67) * t43) * qJD(2) + t69 * (t18 * t24 - t19 * t23 - t6 * t26 + t48 * t7); -0.2e1 * pkin(2) * t28 + 0.2e1 * t12 * t10 + 0.2e1 * t36 * t11 + 0.2e1 * t5 * t13 + (-Ifges(4,4) * t40 + pkin(3) * t14) * qJD(3) * t61 + 0.2e1 * m(5) * (t36 * t60 + t55) + 0.2e1 * m(6) * (t12 * t5 + t55) + (-t48 * t50 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t26 + t15 * t82) * t24 + (t26 * t50 - 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t48 - 0.2e1 * t69 * t16) * t23 + (0.2e1 * Ifges(4,4) * t42 + (Ifges(4,1) - Ifges(4,2)) * t61) * t62 + (t8 * t26 + t48 * t9) * t82; m(6) * t19 * qJD(5) + t79 * mrSges(4,2) - t41 * t59 - t77 * t6 - t65 * t72 + t76 * t7; Ifges(4,5) * t62 - Ifges(4,6) * t64 + t35 * t71 + (m(6) * t16 + mrSges(6,2) * t48) * qJD(5) + t76 * t9 + t77 * t8 + (mrSges(4,2) * t64 - t59) * pkin(6) + (-mrSges(5,3) * t54 + Ifges(6,4) + Ifges(5,5)) * t24 + (-t33 * mrSges(6,2) - mrSges(5,3) * t73 - Ifges(5,6) + Ifges(6,6)) * t23; 0.2e1 * t78 * qJD(5); t81 * t41 * qJD(2); m(5) * t60 + m(6) * t5 - t80; 0; 0; -m(6) * t6; m(6) * t8 + t71; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
