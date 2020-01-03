% Calculate time derivative of joint inertia matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:56
% EndTime: 2020-01-03 11:46:59
% DurationCPUTime: 0.82s
% Computational Cost: add. (805->112), mult. (1754->174), div. (0->0), fcn. (1413->6), ass. (0->52)
t45 = cos(qJ(3));
t53 = -cos(pkin(8)) * pkin(1) - pkin(2);
t38 = -t45 * pkin(3) + t53;
t62 = 0.2e1 * t38;
t70 = 2 * mrSges(5,3);
t69 = 2 * mrSges(6,3);
t44 = cos(qJ(4));
t68 = t44 * pkin(3);
t42 = sin(qJ(4));
t43 = sin(qJ(3));
t37 = t42 * t45 + t44 * t43;
t63 = qJD(3) + qJD(4);
t22 = t63 * t37;
t36 = -t42 * t43 + t44 * t45;
t67 = t36 * t22;
t21 = t63 * t36;
t66 = t37 * t21;
t65 = -mrSges(5,1) - mrSges(6,1);
t39 = sin(pkin(8)) * pkin(1) + pkin(6);
t60 = pkin(7) + t39;
t34 = t60 * t43;
t35 = t60 * t45;
t12 = -t42 * t34 + t44 * t35;
t64 = pkin(3) * qJD(4);
t61 = t22 * pkin(4);
t56 = qJD(4) * t44;
t59 = (t21 * t42 + t37 * t56) * pkin(3);
t58 = -t22 * mrSges(6,1) - t21 * mrSges(6,2);
t57 = qJD(4) * t42;
t55 = 0.2e1 * t45;
t54 = t36 * t57;
t52 = (-mrSges(5,2) - mrSges(6,2)) * t44;
t11 = -t44 * t34 - t42 * t35;
t51 = qJD(3) * t60;
t50 = 2 * Ifges(5,4) + 2 * Ifges(6,4);
t32 = t43 * t51;
t33 = t45 * t51;
t7 = -t12 * qJD(4) + t42 * t32 - t44 * t33;
t3 = -t21 * qJ(5) - t37 * qJD(5) + t7;
t49 = m(6) * t3 - t21 * mrSges(6,3);
t16 = t22 * mrSges(5,1);
t48 = -t21 * mrSges(5,2) - t16 + t58;
t6 = -t44 * t32 - t42 * t33 - t34 * t56 - t35 * t57;
t47 = -t42 * t22 + (t36 * t44 + t37 * t42) * qJD(4);
t2 = -t22 * qJ(5) + t36 * qJD(5) + t6;
t46 = t7 * mrSges(5,1) + t3 * mrSges(6,1) - t6 * mrSges(5,2) - t2 * mrSges(6,2) - (Ifges(5,6) + Ifges(6,6)) * t22 + (Ifges(5,5) + Ifges(6,5)) * t21;
t40 = pkin(4) + t68;
t23 = -t36 * pkin(4) + t38;
t10 = qJD(3) * t43 * pkin(3) + t61;
t9 = t36 * qJ(5) + t12;
t8 = -t37 * qJ(5) + t11;
t1 = [0.2e1 * t10 * (-t36 * mrSges(6,1) + t37 * mrSges(6,2)) + t16 * t62 - 0.2e1 * t23 * t58 + 0.2e1 * m(5) * (t11 * t7 + t12 * t6) + 0.2e1 * m(6) * (t23 * t10 + t9 * t2 + t8 * t3) + (mrSges(5,2) * t62 - 0.2e1 * t11 * mrSges(5,3) - 0.2e1 * t8 * mrSges(6,3) + t36 * t50) * t21 - (t12 * t70 + t37 * t50 + t9 * t69) * t22 + ((t53 * mrSges(4,2) + Ifges(4,4) * t45) * t55 + (0.2e1 * t53 * mrSges(4,1) + 0.2e1 * pkin(3) * (-t36 * mrSges(5,1) + t37 * mrSges(5,2)) + m(5) * pkin(3) * t62 - 0.2e1 * Ifges(4,4) * t43 + (Ifges(4,1) - Ifges(4,2)) * t55) * t43) * qJD(3) + (t2 * t36 - t3 * t37) * t69 + (t6 * t36 - t7 * t37) * t70 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t66 - 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t67; m(5) * (-t11 * t22 + t12 * t21 + t7 * t36 + t6 * t37) + m(6) * (t2 * t37 + t9 * t21 - t8 * t22 + t3 * t36); 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t66 - t67); t49 * t40 + (Ifges(4,5) * t45 - Ifges(4,6) * t43 + (-mrSges(4,1) * t45 + mrSges(4,2) * t43) * t39) * qJD(3) + (m(6) * (t2 * t42 + t9 * t56 - t8 * t57) + m(5) * (-t11 * t57 + t12 * t56 + t42 * t6 + t44 * t7) + t47 * mrSges(6,3) + (-t44 * t21 + t47) * mrSges(5,3)) * pkin(3) + t46; (-t43 * mrSges(4,1) - t45 * mrSges(4,2)) * qJD(3) + m(5) * ((-t22 * t44 - t54) * pkin(3) + t59) + m(6) * (-pkin(3) * t54 - t40 * t22 + t59) + t48; 0.2e1 * (t52 + ((-t40 + t68) * m(6) + t65) * t42) * t64; t49 * pkin(4) + t46; -m(6) * t61 + t48; (t52 + (-m(6) * pkin(4) + t65) * t42) * t64; 0; m(6) * t10 - t58; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
