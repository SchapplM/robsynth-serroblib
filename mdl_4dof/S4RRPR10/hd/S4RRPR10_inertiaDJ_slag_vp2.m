% Calculate time derivative of joint inertia matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:12
% DurationCPUTime: 0.56s
% Computational Cost: add. (393->132), mult. (908->202), div. (0->0), fcn. (577->4), ass. (0->70)
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t37 = cos(qJ(2));
t55 = qJD(4) * t37;
t35 = sin(qJ(2));
t57 = qJD(2) * t35;
t40 = t34 * t55 + t36 * t57;
t48 = pkin(2) * t57 - t35 * qJD(3);
t56 = qJD(2) * t37;
t79 = -0.2e1 * qJ(3) * t56 + 0.2e1 * t48;
t71 = mrSges(5,3) * t37;
t19 = t35 * mrSges(5,1) + t34 * t71;
t20 = -t35 * mrSges(5,2) - t36 * t71;
t38 = -pkin(2) - pkin(6);
t58 = t35 * qJ(3);
t15 = t38 * t37 - pkin(1) - t58;
t72 = pkin(3) + pkin(5);
t27 = t72 * t35;
t6 = -t34 * t15 + t36 * t27;
t7 = t36 * t15 + t34 * t27;
t46 = t34 * t6 - t36 * t7;
t78 = -m(5) * t46 - t34 * t19 + t36 * t20;
t77 = -0.2e1 * pkin(1);
t41 = -t37 * pkin(2) - t58;
t23 = -pkin(1) + t41;
t76 = -0.2e1 * t23;
t75 = -t34 / 0.2e1;
t74 = -t36 / 0.2e1;
t73 = t36 / 0.2e1;
t70 = Ifges(5,4) * t34;
t69 = Ifges(5,4) * t36;
t68 = Ifges(5,5) * t35;
t44 = Ifges(5,1) * t34 + t69;
t12 = -t44 * t37 + t68;
t67 = t34 * t12;
t26 = Ifges(5,1) * t36 - t70;
t65 = t34 * t26;
t64 = t34 * t38;
t63 = t35 * mrSges(4,3);
t43 = Ifges(5,2) * t36 + t70;
t11 = Ifges(5,6) * t35 - t43 * t37;
t62 = t36 * t11;
t25 = -Ifges(5,2) * t34 + t69;
t60 = t36 * t25;
t59 = t36 * t38;
t28 = t72 * t37;
t53 = t36 * t55;
t52 = t34 * t57;
t50 = m(4) * pkin(5) + mrSges(4,1);
t49 = Ifges(5,5) * t52 + t40 * Ifges(5,6) + Ifges(5,3) * t56;
t10 = (pkin(6) * t35 - qJ(3) * t37) * qJD(2) + t48;
t22 = qJD(2) * t28;
t1 = t6 * qJD(4) + t36 * t10 + t34 * t22;
t2 = -t7 * qJD(4) - t34 * t10 + t36 * t22;
t47 = -t34 * t1 - t36 * t2;
t45 = mrSges(5,1) * t36 - mrSges(5,2) * t34;
t24 = t34 * mrSges(5,1) + t36 * mrSges(5,2);
t42 = -Ifges(5,5) * t34 - Ifges(5,6) * t36;
t39 = t52 - t53;
t21 = t72 * t57;
t18 = t44 * qJD(4);
t17 = t43 * qJD(4);
t16 = t45 * qJD(4);
t14 = t45 * t37;
t9 = mrSges(5,1) * t56 - t39 * mrSges(5,3);
t8 = -mrSges(5,2) * t56 + t40 * mrSges(5,3);
t5 = -t40 * mrSges(5,1) + t39 * mrSges(5,2);
t4 = -t26 * t55 + (Ifges(5,5) * t37 + t44 * t35) * qJD(2);
t3 = -t25 * t55 + (Ifges(5,6) * t37 + t43 * t35) * qJD(2);
t13 = [0.2e1 * m(5) * (t7 * t1 + t6 * t2 - t28 * t21) + t35 * t49 + 0.2e1 * t2 * t19 + 0.2e1 * t1 * t20 - 0.2e1 * t21 * t14 + 0.2e1 * t28 * t5 + 0.2e1 * t7 * t8 + 0.2e1 * t6 * t9 + (m(4) * t23 - t63) * t79 + (mrSges(4,2) * t79 - t36 * t3 - t34 * t4 + (t11 * t34 + (-t12 - t68) * t36) * qJD(4)) * t37 + ((mrSges(3,1) * t77 + mrSges(4,2) * t76 + t62 + t67 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6)) * t35) * t35 + (mrSges(3,2) * t77 + mrSges(4,3) * t76 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t42) * t37 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3)) * t35) * t37) * qJD(2); t8 * t64 + m(5) * (-qJ(3) * t21 + qJD(3) * t28 + t1 * t64 + t2 * t59) + t9 * t59 + t4 * t73 + t3 * t75 - t21 * t24 + t28 * t16 + qJD(3) * t14 + qJ(3) * t5 + t47 * mrSges(5,3) + (t50 * qJD(3) - t17 * t74 - t18 * t75) * t37 + (-t67 / 0.2e1 - t62 / 0.2e1 + t35 * t42 / 0.2e1 + (t34 * t25 / 0.2e1 + t26 * t74) * t37 + t46 * mrSges(5,3) + t78 * t38) * qJD(4) + ((Ifges(4,5) - Ifges(3,6) + t65 / 0.2e1 - qJ(3) * mrSges(4,1) + t60 / 0.2e1) * t35 + (m(4) * t41 + t35 * mrSges(3,2) - t63) * pkin(5) + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t73 + Ifges(5,6) * t75 - Ifges(4,4) + Ifges(3,5) + (-mrSges(3,1) + mrSges(4,2)) * pkin(5)) * t37) * qJD(2); 0.2e1 * qJ(3) * t16 + t34 * t17 - t36 * t18 + (-t60 - t65) * qJD(4) + 0.2e1 * (mrSges(4,3) + t24 + (m(4) + m(5)) * qJ(3)) * qJD(3); -m(5) * t47 + t78 * qJD(4) + t34 * t8 + t36 * t9 + t50 * t56; 0; 0; t2 * mrSges(5,1) - t1 * mrSges(5,2) - Ifges(5,5) * t53 + t49; ((-mrSges(5,2) * t38 - Ifges(5,6)) * t36 + (-mrSges(5,1) * t38 - Ifges(5,5)) * t34) * qJD(4); -t24 * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1), t13(2), t13(4), t13(7); t13(2), t13(3), t13(5), t13(8); t13(4), t13(5), t13(6), t13(9); t13(7), t13(8), t13(9), t13(10);];
Mq = res;
