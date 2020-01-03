% Calculate time derivative of joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:55
% EndTime: 2019-12-31 19:28:57
% DurationCPUTime: 0.86s
% Computational Cost: add. (1056->146), mult. (2300->225), div. (0->0), fcn. (1972->6), ass. (0->61)
t66 = (mrSges(5,2) + mrSges(4,3));
t74 = 2 * t66;
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t73 = -t51 * mrSges(6,1) - t53 * mrSges(6,2);
t72 = 2 * m(6);
t54 = cos(qJ(2));
t48 = -pkin(2) * t54 - pkin(1);
t71 = 0.2e1 * t48;
t50 = cos(pkin(8));
t47 = -pkin(2) * t50 - pkin(3);
t44 = -pkin(4) + t47;
t49 = sin(pkin(8));
t45 = pkin(2) * t49 + qJ(4);
t24 = t44 * t53 - t45 * t51;
t20 = qJD(4) * t53 + qJD(5) * t24;
t70 = t20 * mrSges(6,2);
t25 = t44 * t51 + t45 * t53;
t21 = -qJD(4) * t51 - qJD(5) * t25;
t69 = t21 * mrSges(6,1);
t65 = -qJ(3) - pkin(6);
t52 = sin(qJ(2));
t61 = qJD(2) * t65;
t33 = qJD(3) * t54 + t52 * t61;
t56 = -t52 * qJD(3) + t54 * t61;
t14 = t50 * t33 + t49 * t56;
t42 = t65 * t52;
t43 = t65 * t54;
t23 = t49 * t42 - t50 * t43;
t64 = 0.2e1 * t54;
t63 = pkin(2) * qJD(2) * t52;
t13 = t33 * t49 - t50 * t56;
t22 = -t50 * t42 - t43 * t49;
t62 = t13 * t22 + t23 * t14;
t60 = -2 * Ifges(4,4) + 2 * Ifges(5,5);
t37 = t49 * t52 - t50 * t54;
t38 = t49 * t54 + t50 * t52;
t18 = t37 * t53 - t38 * t51;
t34 = t38 * qJD(2);
t35 = t37 * qJD(2);
t5 = qJD(5) * t18 + t34 * t51 - t35 * t53;
t19 = t37 * t51 + t38 * t53;
t6 = -qJD(5) * t19 + t34 * t53 + t35 * t51;
t59 = -t6 * mrSges(6,1) + t5 * mrSges(6,2);
t15 = -pkin(7) * t38 + t22;
t16 = pkin(7) * t37 + t23;
t3 = t15 * t53 - t16 * t51;
t4 = t15 * t51 + t16 * t53;
t58 = t38 * qJ(4) - t48;
t11 = t34 * pkin(3) + t35 * qJ(4) - t38 * qJD(4) + t63;
t10 = pkin(7) * t34 + t14;
t9 = pkin(7) * t35 + t13;
t1 = qJD(5) * t3 + t10 * t53 + t51 * t9;
t2 = -qJD(5) * t4 - t10 * t51 + t53 * t9;
t55 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t5 + Ifges(6,6) * t6;
t31 = t35 * mrSges(4,2);
t30 = t34 * mrSges(5,1);
t17 = t37 * pkin(3) - t58;
t12 = (-pkin(3) - pkin(4)) * t37 + t58;
t7 = pkin(4) * t34 + t11;
t8 = [-t31 * t71 + 0.2e1 * t11 * (t37 * mrSges(5,1) - t38 * mrSges(5,3)) + 0.2e1 * t17 * t30 + 0.2e1 * t18 * Ifges(6,2) * t6 + 0.2e1 * t5 * t19 * Ifges(6,1) - 0.2e1 * t7 * (-t18 * mrSges(6,1) + t19 * mrSges(6,2)) + 0.2e1 * t12 * t59 + 0.2e1 * m(4) * t62 + (t1 * t4 - t12 * t7 + t2 * t3) * t72 + 0.2e1 * m(5) * (t11 * t17 + t62) + 0.2e1 * (t18 * t5 + t6 * t19) * Ifges(6,4) + 0.2e1 * (t1 * t18 - t2 * t19 - t3 * t5 + t4 * t6) * mrSges(6,3) - (-0.2e1 * t17 * mrSges(5,3) + t37 * t60 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t38 + t22 * t74) * t35 + (mrSges(4,1) * t71 + t38 * t60 + 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t37 - 0.2e1 * t66 * t23) * t34 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t54) * t64 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t71 + 0.2e1 * pkin(2) * (mrSges(4,1) * t37 + mrSges(4,2) * t38) - 0.2e1 * Ifges(3,4) * t52 + (Ifges(3,1) - Ifges(3,2)) * t64) * t52) * qJD(2) + (t13 * t38 - t14 * t37) * t74; -(Ifges(5,4) + Ifges(4,5)) * t35 + (Ifges(5,6) - Ifges(4,6)) * t34 + (-mrSges(4,2) + mrSges(5,3)) * t14 + (-mrSges(4,1) - mrSges(5,1)) * t13 + m(5) * (qJD(4) * t23 + t13 * t47 + t14 * t45) + m(6) * (t1 * t25 + t2 * t24 + t20 * t4 + t21 * t3) + (-qJD(4) * t37 - t34 * t45 - t35 * t47) * mrSges(5,2) + (Ifges(3,5) * t54 - Ifges(3,6) * t52 + (-mrSges(3,1) * t54 + mrSges(3,2) * t52) * pkin(6)) * qJD(2) + (m(4) * (-t13 * t50 + t14 * t49) + (-t34 * t49 + t35 * t50) * mrSges(4,3)) * pkin(2) + (t18 * t20 - t19 * t21 - t24 * t5 + t25 * t6) * mrSges(6,3) - t55; (t20 * t25 + t21 * t24) * t72 + 0.2e1 * t70 - 0.2e1 * t69 + 0.2e1 * (m(5) * t45 + mrSges(5,3)) * qJD(4); m(4) * t63 + m(5) * t11 + m(6) * t7 + t34 * mrSges(4,1) + t35 * mrSges(5,3) + t30 - t31 - t59; 0; 0; -t35 * mrSges(5,2) + m(6) * (t1 * t51 + t2 * t53 + (-t3 * t51 + t4 * t53) * qJD(5)) + m(5) * t13 + (-t53 * t5 + t51 * t6 + (t18 * t53 + t19 * t51) * qJD(5)) * mrSges(6,3); m(6) * (t20 * t51 + t21 * t53) + (m(6) * (-t24 * t51 + t25 * t53) - t73) * qJD(5); 0; 0; t55; t69 - t70; 0; t73 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
