% Calculate time derivative of joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:35
% DurationCPUTime: 0.66s
% Computational Cost: add. (461->156), mult. (1002->247), div. (0->0), fcn. (630->4), ass. (0->72)
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t20 = -t36 * mrSges(6,1) + t34 * mrSges(6,2);
t76 = -mrSges(5,1) + t20;
t61 = t34 ^ 2 + t36 ^ 2;
t75 = 2 * m(6);
t74 = 2 * qJD(3);
t73 = -t34 / 0.2e1;
t35 = sin(qJ(4));
t72 = t35 * pkin(4);
t37 = cos(qJ(4));
t71 = t37 * pkin(7);
t70 = Ifges(6,4) * t34;
t69 = Ifges(6,4) * t36;
t68 = Ifges(6,5) * t34;
t67 = Ifges(6,6) * t34;
t66 = Ifges(6,6) * t35;
t65 = Ifges(6,6) * t36;
t32 = qJ(2) - pkin(6);
t64 = t32 * t35;
t63 = t37 * mrSges(6,3);
t33 = (pkin(1) + qJ(3));
t60 = qJD(4) * t35;
t50 = t34 * t60;
t59 = qJD(4) * t37;
t62 = Ifges(6,6) * t50 + Ifges(6,3) * t59;
t58 = qJD(5) * t34;
t57 = qJD(5) * t35;
t56 = qJD(5) * t36;
t55 = qJD(5) * t37;
t29 = t35 ^ 2;
t54 = t29 * qJD(2);
t31 = t37 ^ 2;
t27 = t31 * qJD(2);
t53 = t33 * qJD(3);
t52 = t32 * t59;
t51 = t36 * t55;
t49 = t34 * t59;
t48 = t36 * t59;
t47 = -Ifges(6,5) * t36 + (2 * Ifges(5,4));
t19 = t33 - t71 + t72;
t38 = t35 * qJD(2) + qJD(5) * t19 + t52;
t44 = -t32 * t57 + qJD(3) + (pkin(4) * t37 + pkin(7) * t35) * qJD(4);
t1 = t44 * t34 + t38 * t36;
t6 = t36 * t19 - t34 * t64;
t46 = -qJD(5) * t6 + t1;
t45 = m(6) * pkin(4) - t76;
t43 = mrSges(6,1) * t34 + mrSges(6,2) * t36;
t42 = Ifges(6,1) * t36 - t70;
t22 = Ifges(6,1) * t34 + t69;
t41 = -Ifges(6,2) * t34 + t69;
t21 = Ifges(6,2) * t36 + t70;
t40 = t34 * t55 + t36 * t60;
t39 = t50 - t51;
t26 = Ifges(6,5) * t56;
t24 = t32 * t27;
t18 = t35 * mrSges(6,1) - t36 * t63;
t17 = -t35 * mrSges(6,2) - t34 * t63;
t16 = t42 * qJD(5);
t15 = t41 * qJD(5);
t14 = t43 * qJD(5);
t12 = t43 * t37;
t11 = Ifges(6,5) * t35 + t42 * t37;
t10 = t41 * t37 + t66;
t9 = -mrSges(6,2) * t59 + t39 * mrSges(6,3);
t8 = mrSges(6,1) * t59 + t40 * mrSges(6,3);
t7 = t34 * t19 + t36 * t64;
t5 = -t39 * mrSges(6,1) - t40 * mrSges(6,2);
t4 = -t22 * t55 + (Ifges(6,5) * t37 - t42 * t35) * qJD(4);
t3 = -t21 * t55 + (Ifges(6,6) * t37 - t41 * t35) * qJD(4);
t2 = -t38 * t34 + t44 * t36;
t13 = [(mrSges(4,3) * t74) + 0.2e1 * t1 * t17 + 0.2e1 * t2 * t18 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3) + (-t29 - t31) * mrSges(5,3)) * qJD(2) + (t7 * t1 + t6 * t2 + t24) * t75 + 0.2e1 * m(5) * (t32 * t54 + t24 + t53) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t53) + (mrSges(5,1) * t74 + (-(2 * t33 * mrSges(5,2)) + t34 * t10 - t36 * t11 + 0.2e1 * t32 * t12 + t47 * t35) * qJD(4) + t62) * t35 + ((mrSges(5,2) * t74) - 0.2e1 * qJD(2) * t12 - t34 * t3 - 0.2e1 * t32 * t5 + t36 * t4 + (t35 * (-t65 - t68) - t36 * t10 - t34 * t11) * qJD(5) + (0.2e1 * t33 * mrSges(5,1) + (-t47 - t67) * t37 + (-0.2e1 * m(6) * t32 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(6,3)) * t35) * qJD(4)) * t37; m(6) * (-t34 * t1 - t36 * t2 + (t34 * t6 - t36 * t7) * qJD(5)) - t17 * t56 - t34 * t9 + t18 * t58 - t36 * t8 + (-t37 * mrSges(5,1) + t35 * mrSges(5,2)) * qJD(4) + (-m(5) - m(4)) * qJD(3); 0; m(4) * qJD(2) + (-t5 + (t36 * t17 - t34 * t18) * qJD(4)) * t37 + m(6) * (t7 * t48 - t6 * t49 + t27) + m(5) * (t27 + t54) + (m(6) * (t1 * t36 - t2 * t34 - t6 * t56 - t7 * t58 - 0.2e1 * t52) - t17 * t58 + t36 * t9 - t18 * t56 - t34 * t8 + qJD(4) * t12) * t35; 0; (-0.1e1 + t61) * t35 * t59 * t75; -pkin(4) * t5 + (t26 / 0.2e1 - qJD(2) * mrSges(5,2) + (-t45 * t32 - Ifges(5,5)) * qJD(4)) * t35 + (t4 / 0.2e1 + t21 * t60 / 0.2e1 - t2 * mrSges(6,3) + (-t66 / 0.2e1 - t7 * mrSges(6,3) - t10 / 0.2e1) * qJD(5) + (m(6) * (-qJD(5) * t7 - t2) - qJD(5) * t17 - t8) * pkin(7)) * t34 + (t3 / 0.2e1 - t22 * t60 / 0.2e1 + qJD(5) * t11 / 0.2e1 + t46 * mrSges(6,3) + (m(6) * t46 - qJD(5) * t18 + t9) * pkin(7)) * t36 + (t36 * t16 / 0.2e1 + t15 * t73 - t32 * t14 + (-t36 * t21 / 0.2e1 + t22 * t73) * qJD(5) + (-t32 * mrSges(5,2) + t68 / 0.2e1 + t65 / 0.2e1 - Ifges(5,6)) * qJD(4) + t45 * qJD(2)) * t37; 0; -t37 * t14 + (-t37 * mrSges(5,2) + m(6) * (t61 * t71 - t72) + t61 * t63 + t76 * t35) * qJD(4); -0.2e1 * pkin(4) * t14 + t36 * t15 + t34 * t16 + (-t21 * t34 + t22 * t36) * qJD(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) - t40 * Ifges(6,5) - Ifges(6,6) * t51 + t62; t14; (t34 * t57 - t48) * mrSges(6,2) + (-t35 * t56 - t49) * mrSges(6,1); t26 + (t20 * pkin(7) - t67) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
