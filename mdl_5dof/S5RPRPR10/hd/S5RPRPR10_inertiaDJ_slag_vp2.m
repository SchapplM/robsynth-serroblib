% Calculate time derivative of joint inertia matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:25:56
% DurationCPUTime: 0.38s
% Computational Cost: add. (486->86), mult. (946->138), div. (0->0), fcn. (622->6), ass. (0->45)
t34 = sin(pkin(8));
t35 = cos(pkin(8));
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t15 = t34 * t39 + t35 * t37;
t10 = t15 * qJD(3);
t14 = t34 * t37 - t35 * t39;
t11 = t14 * qJD(3);
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t44 = mrSges(6,1) * t36 + mrSges(6,2) * t38;
t16 = t44 * qJD(5);
t22 = -t38 * mrSges(6,1) + t36 * mrSges(6,2);
t55 = mrSges(5,1) - t22;
t54 = t36 ^ 2 + t38 ^ 2;
t49 = t54 * mrSges(6,3);
t64 = mrSges(5,2) - t49;
t65 = t64 * t11 - (t37 * mrSges(4,1) + t39 * mrSges(4,2)) * qJD(3) - t55 * t10 + t14 * t16;
t40 = -pkin(1) - pkin(2);
t20 = -t37 * qJ(2) + t39 * t40;
t8 = t39 * qJD(2) + t20 * qJD(3);
t21 = t39 * qJ(2) + t37 * t40;
t9 = -t37 * qJD(2) - t21 * qJD(3);
t1 = t34 * t8 - t35 * t9;
t63 = t55 * t1;
t51 = qJD(5) * t38;
t52 = qJD(5) * t36;
t41 = -(t36 * (Ifges(6,4) * t36 + Ifges(6,2) * t38) - t38 * (Ifges(6,1) * t36 + Ifges(6,4) * t38)) * qJD(5) + t38 * (Ifges(6,4) * t51 - Ifges(6,2) * t52) + t36 * (Ifges(6,1) * t51 - Ifges(6,4) * t52);
t61 = m(5) * pkin(3);
t60 = t14 * t1;
t59 = t8 * mrSges(4,2);
t58 = t9 * mrSges(4,1);
t57 = t14 * t10;
t19 = -pkin(3) + t20;
t6 = t34 * t19 + t35 * t21;
t4 = -pkin(7) + t6;
t50 = t54 * t4;
t48 = t54 * t15;
t26 = t34 * pkin(3) + pkin(7);
t47 = t54 * t26;
t5 = t35 * t19 - t34 * t21;
t27 = -t35 * pkin(3) - pkin(4);
t3 = pkin(4) - t5;
t2 = t34 * t9 + t35 * t8;
t7 = [-0.2e1 * t58 + 0.2e1 * t59 - 0.2e1 * t3 * t16 + 0.2e1 * t2 * mrSges(5,2) + 0.2e1 * m(6) * (t3 * t1 + t2 * t50) + 0.2e1 * m(5) * (-t5 * t1 + t6 * t2) + 0.2e1 * m(4) * (t20 * t9 + t21 * t8) + t41 + 0.2e1 * t63 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) - 0.2e1 * t49 * t2; m(6) * (t10 * t3 - t11 * t50 + t2 * t48 + t60) + m(5) * (-t10 * t5 - t11 * t6 + t15 * t2 + t60) + m(4) * (t37 * t8 + t39 * t9 + (-t20 * t37 + t21 * t39) * qJD(3)) - t65; 0.2e1 * m(5) * (-t15 * t11 + t57) + 0.2e1 * m(6) * (-t11 * t48 + t57); t58 - t59 + (-t27 + t3) * t16 - t63 - t64 * t2 + m(6) * (t27 * t1 + t2 * t47) + (-t1 * t35 + t2 * t34) * t61 - t41; m(6) * (t27 * t10 - t11 * t47) + (-t10 * t35 - t11 * t34) * t61 + t65; 0.2e1 * t27 * t16 + t41; 0; 0; 0; 0; -t44 * t2 + ((-mrSges(6,1) * t4 - Ifges(6,5)) * t38 + (mrSges(6,2) * t4 + Ifges(6,6)) * t36) * qJD(5); (t11 * t38 + t15 * t52) * mrSges(6,2) + (t11 * t36 - t15 * t51) * mrSges(6,1); (Ifges(6,5) * t38 - Ifges(6,6) * t36 + t22 * t26) * qJD(5); -t16; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
