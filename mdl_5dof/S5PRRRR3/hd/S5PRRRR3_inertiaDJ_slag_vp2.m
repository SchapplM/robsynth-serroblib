% Calculate time derivative of joint inertia matrix for
% S5PRRRR3
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:05
% EndTime: 2019-12-05 17:06:06
% DurationCPUTime: 0.30s
% Computational Cost: add. (310->72), mult. (863->119), div. (0->0), fcn. (489->6), ass. (0->54)
t32 = sin(qJ(5));
t30 = t32 ^ 2;
t35 = cos(qJ(5));
t31 = t35 ^ 2;
t65 = t30 + t31;
t64 = Ifges(6,1) - Ifges(6,2);
t36 = cos(qJ(4));
t63 = t65 * t36;
t18 = -t35 * mrSges(6,1) + t32 * mrSges(6,2);
t33 = sin(qJ(4));
t49 = qJD(4) * t33;
t15 = pkin(3) * t18 * t49;
t48 = qJD(4) * t36;
t45 = pkin(3) * t48;
t43 = mrSges(6,3) * t45;
t20 = t30 * t43;
t21 = t31 * t43;
t46 = qJD(5) * t35;
t47 = qJD(5) * t32;
t17 = -mrSges(6,1) * t47 - mrSges(6,2) * t46;
t25 = -t36 * pkin(3) - pkin(4);
t8 = t25 * t17;
t62 = t15 + t20 + t21 - t8;
t61 = pkin(4) * t17;
t37 = cos(qJ(3));
t26 = t37 * pkin(2) + pkin(3);
t34 = sin(qJ(3));
t52 = t33 * t34;
t5 = t26 * t48 + (-t34 * t49 + (t36 * t37 - t52) * qJD(3)) * pkin(2);
t60 = t30 * t5;
t59 = t31 * t5;
t58 = t5 * mrSges(5,2);
t55 = Ifges(6,6) * t32;
t51 = t34 * t36;
t12 = pkin(2) * t51 + t33 * t26;
t50 = pkin(3) * qJD(4);
t44 = t65 * t5;
t42 = -t33 * mrSges(5,1) - t36 * mrSges(5,2);
t41 = -mrSges(6,1) * t32 - mrSges(6,2) * t35;
t11 = -pkin(2) * t52 + t36 * t26;
t40 = (-0.2e1 * Ifges(6,4) * t32 + t64 * t35) * t47 + (0.2e1 * Ifges(6,4) * t35 + t64 * t32) * t46;
t39 = (-mrSges(4,1) * t34 - mrSges(4,2) * t37) * qJD(3) * pkin(2);
t6 = t26 * t49 + (t34 * t48 + (t33 * t37 + t51) * qJD(3)) * pkin(2);
t1 = t6 * t18;
t2 = mrSges(6,3) * t60;
t3 = mrSges(6,3) * t59;
t4 = t6 * mrSges(5,1);
t9 = -pkin(4) - t11;
t7 = t9 * t17;
t38 = t1 + t2 + t3 - t4 + t40 - t7 - t58;
t29 = Ifges(6,5) * t46;
t24 = t33 * pkin(3) + pkin(8);
t10 = pkin(8) + t12;
t13 = [0; 0; -0.2e1 * t58 + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 - 0.2e1 * t7 + 0.2e1 * t39 + 0.2e1 * m(6) * (t10 * t44 + t9 * t6) + 0.2e1 * m(5) * (-t11 * t6 + t12 * t5) + t40; 0; t38 + m(6) * (t25 * t6 + (t59 + t60) * t24) + (m(5) * (t33 * t5 - t36 * t6) + (m(6) * (t63 * t10 + t33 * t9) + m(5) * (-t11 * t33 + t12 * t36) + t42) * qJD(4)) * pkin(3) + t39 + t62; 0.2e1 * t15 + 0.2e1 * t20 + 0.2e1 * t21 - 0.2e1 * t8 + 0.2e1 * (m(6) * (t63 * t24 + t25 * t33) + t42) * t50 + t40; 0; t61 + m(6) * (-pkin(4) * t6 + pkin(8) * t44) + t38; t61 + (m(6) * (-pkin(4) * t33 + t63 * pkin(8)) + t42) * t50 + t40 + t62; t40 + 0.2e1 * t61; t17; t29 + t41 * t5 + (t18 * t10 - t55) * qJD(5); t29 + t41 * t45 + (t18 * t24 - t55) * qJD(5); t29 + (t18 * pkin(8) - t55) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
