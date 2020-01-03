% Calculate time derivative of joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:13
% DurationCPUTime: 0.67s
% Computational Cost: add. (972->127), mult. (2163->195), div. (0->0), fcn. (1945->6), ass. (0->57)
t59 = (mrSges(5,2) + mrSges(4,3));
t69 = 2 * t59;
t43 = cos(pkin(8));
t57 = -pkin(2) * t43 - pkin(1);
t67 = 0.2e1 * t57;
t42 = sin(pkin(8));
t58 = pkin(6) + qJ(2);
t36 = t58 * t42;
t37 = t58 * t43;
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t21 = -t45 * t36 + t47 * t37;
t44 = sin(qJ(5));
t46 = cos(qJ(5));
t66 = -t44 * mrSges(6,1) - t46 * mrSges(6,2);
t65 = 2 * m(6);
t48 = -pkin(3) - pkin(4);
t34 = -t44 * qJ(4) + t46 * t48;
t22 = t46 * qJD(4) + qJD(5) * t34;
t64 = t22 * mrSges(6,2);
t35 = t46 * qJ(4) + t44 * t48;
t23 = -t44 * qJD(4) - qJD(5) * t35;
t63 = t23 * mrSges(6,1);
t30 = t47 * t36;
t60 = t47 * t43;
t13 = -qJD(3) * t30 + qJD(2) * t60 + (-qJD(2) * t42 - qJD(3) * t37) * t45;
t33 = t42 * t47 + t43 * t45;
t14 = t33 * qJD(2) + qJD(3) * t21;
t20 = t37 * t45 + t30;
t56 = t21 * t13 + t14 * t20;
t54 = 2 * Ifges(5,5) - 2 * Ifges(4,4);
t32 = t42 * t45 - t60;
t18 = t32 * t46 - t33 * t44;
t26 = t32 * qJD(3);
t27 = t33 * qJD(3);
t5 = t18 * qJD(5) - t26 * t46 + t27 * t44;
t19 = t32 * t44 + t33 * t46;
t6 = -t19 * qJD(5) + t26 * t44 + t27 * t46;
t53 = -t6 * mrSges(6,1) + t5 * mrSges(6,2);
t15 = -pkin(7) * t33 + t20;
t16 = pkin(7) * t32 + t21;
t3 = t15 * t46 - t16 * t44;
t4 = t15 * t44 + t16 * t46;
t52 = -qJ(4) * t26 + qJD(4) * t33;
t51 = qJ(4) * t33 - t57;
t10 = pkin(7) * t26 + t14;
t9 = pkin(7) * t27 + t13;
t1 = qJD(5) * t3 + t10 * t44 + t46 * t9;
t2 = -qJD(5) * t4 + t10 * t46 - t44 * t9;
t49 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t5 + Ifges(6,6) * t6;
t25 = t27 * mrSges(5,1);
t24 = t26 * mrSges(4,2);
t17 = pkin(3) * t32 - t51;
t12 = t48 * t32 + t51;
t11 = pkin(3) * t27 - t52;
t8 = t48 * t27 + t52;
t7 = [0.2e1 * t17 * t25 - t24 * t67 + 0.2e1 * t11 * (t32 * mrSges(5,1) - t33 * mrSges(5,3)) + 0.2e1 * t18 * Ifges(6,2) * t6 + 0.2e1 * t8 * (-t18 * mrSges(6,1) + t19 * mrSges(6,2)) + 0.2e1 * t5 * t19 * Ifges(6,1) + 0.2e1 * t12 * t53 + 0.2e1 * (t18 * t5 + t6 * t19) * Ifges(6,4) + 0.2e1 * (t1 * t18 - t19 * t2 - t3 * t5 + t4 * t6) * mrSges(6,3) + (mrSges(4,1) * t67 + t33 * t54 + 0.2e1 * (Ifges(5,3) + Ifges(4,2)) * t32 - 0.2e1 * t59 * t21) * t27 - (-0.2e1 * mrSges(5,3) * t17 + t32 * t54 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t33 + t20 * t69) * t26 + 0.2e1 * m(5) * (t11 * t17 + t56) + 0.2e1 * m(4) * t56 + (t1 * t4 + t12 * t8 + t2 * t3) * t65 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t42 ^ 2 + t43 ^ 2) * qJD(2) + (-t13 * t32 + t14 * t33) * t69; m(5) * t11 - m(6) * t8 + t27 * mrSges(4,1) + t26 * mrSges(5,3) - t24 + t25 - t53; 0; (Ifges(5,6) - Ifges(4,6)) * t27 - (Ifges(5,4) + Ifges(4,5)) * t26 + (-mrSges(4,1) - mrSges(5,1)) * t14 + (-mrSges(4,2) + mrSges(5,3)) * t13 + m(6) * (t1 * t35 + t2 * t34 + t22 * t4 + t23 * t3) + m(5) * (-pkin(3) * t14 + qJ(4) * t13 + qJD(4) * t21) + (pkin(3) * t26 - qJ(4) * t27 - qJD(4) * t32) * mrSges(5,2) + (t18 * t22 - t19 * t23 - t34 * t5 + t35 * t6) * mrSges(6,3) - t49; 0; (t22 * t35 + t23 * t34) * t65 + 0.2e1 * t64 - 0.2e1 * t63 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); -t26 * mrSges(5,2) + m(6) * (t1 * t44 + t2 * t46 + (-t3 * t44 + t4 * t46) * qJD(5)) + m(5) * t14 + (t44 * t6 - t46 * t5 + (t18 * t46 + t19 * t44) * qJD(5)) * mrSges(6,3); 0; m(6) * (t22 * t44 + t23 * t46) + (m(6) * (-t34 * t44 + t35 * t46) - t66) * qJD(5); 0; t49; 0; t63 - t64; t66 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
