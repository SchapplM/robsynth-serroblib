% Calculate time derivative of joint inertia matrix for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:38
% EndTime: 2019-12-05 15:12:39
% DurationCPUTime: 0.31s
% Computational Cost: add. (327->59), mult. (879->105), div. (0->0), fcn. (791->8), ass. (0->38)
t29 = sin(qJ(5));
t32 = cos(qJ(5));
t44 = t29 ^ 2 + t32 ^ 2;
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t17 = -t31 * t27 + t34 * t28;
t12 = t17 * qJD(3);
t18 = t34 * t27 + t31 * t28;
t13 = t18 * qJD(3);
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t37 = t33 * t17 - t30 * t18;
t5 = t37 * qJD(4) + t33 * t12 - t30 * t13;
t63 = t44 * t5;
t20 = -t32 * mrSges(6,1) + t29 * mrSges(6,2);
t62 = -mrSges(5,1) + t20;
t61 = Ifges(6,1) - Ifges(6,2);
t59 = t44 * t33;
t58 = (mrSges(6,3) * t44 - mrSges(5,2)) * t33 + t62 * t30;
t57 = m(5) / 0.2e1;
t9 = t30 * t17 + t33 * t18;
t6 = qJD(4) * t9 + t30 * t12 + t33 * t13;
t56 = t37 * t6;
t53 = t30 * t37;
t50 = Ifges(6,6) * t29;
t43 = pkin(3) * qJD(4);
t42 = qJD(5) * t29;
t41 = qJD(5) * t32;
t38 = mrSges(6,1) * t29 + mrSges(6,2) * t32;
t19 = t38 * qJD(5);
t36 = -t5 * mrSges(5,2) + mrSges(6,3) * t63 - t37 * t19 + t62 * t6;
t35 = (-0.2e1 * Ifges(6,4) * t29 + t61 * t32) * t42 + (0.2e1 * Ifges(6,4) * t32 + t61 * t29) * t41;
t24 = Ifges(6,5) * t41;
t23 = -t33 * pkin(3) - pkin(4);
t22 = t30 * pkin(3) + pkin(7);
t1 = [0.2e1 * m(6) * (t63 * t9 - t56) + 0.2e1 * m(5) * (t9 * t5 - t56) + 0.2e1 * m(4) * (t18 * t12 - t17 * t13); 0; 0; m(6) * (t22 * t63 + t23 * t6) - t12 * mrSges(4,2) - t13 * mrSges(4,1) + 0.2e1 * ((t30 * t5 - t33 * t6) * t57 + (m(6) * (t59 * t9 - t53) / 0.2e1 + (t33 * t9 - t53) * t57) * qJD(4)) * pkin(3) + t36; 0; 0.2e1 * t23 * t19 + 0.2e1 * (m(6) * (t59 * t22 + t23 * t30) + t58) * t43 + t35; m(6) * (-pkin(4) * t6 + pkin(7) * t63) + t36; 0; (t23 - pkin(4)) * t19 + (m(6) * (-pkin(4) * t30 + t59 * pkin(7)) + t58) * t43 + t35; -0.2e1 * pkin(4) * t19 + t35; (-t32 * t5 + t9 * t42) * mrSges(6,2) + (-t29 * t5 - t9 * t41) * mrSges(6,1); -t19; t24 - t38 * t33 * t43 + (t20 * t22 - t50) * qJD(5); t24 + (t20 * pkin(7) - t50) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
