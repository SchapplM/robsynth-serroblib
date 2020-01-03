% Calculate time derivative of joint inertia matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:20
% EndTime: 2019-12-31 17:42:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (446->81), mult. (1195->133), div. (0->0), fcn. (986->8), ass. (0->53)
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t57 = t39 ^ 2 + t42 ^ 2;
t51 = t57 * mrSges(6,3);
t71 = -mrSges(5,2) + t51;
t27 = -t42 * mrSges(6,1) + t39 * mrSges(6,2);
t70 = -mrSges(5,1) + t27;
t69 = Ifges(6,1) - Ifges(6,2);
t68 = qJD(2) + qJD(3);
t48 = mrSges(6,1) * t39 + mrSges(6,2) * t42;
t26 = t48 * qJD(5);
t67 = 0.2e1 * t26;
t66 = m(5) * pkin(3);
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t43 = cos(qJ(3));
t44 = cos(qJ(2));
t24 = -t40 * t41 + t43 * t44;
t12 = t68 * t24;
t25 = t40 * t44 + t43 * t41;
t13 = t68 * t25;
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t5 = t37 * t12 + t38 * t13;
t9 = -t38 * t24 + t37 * t25;
t65 = t9 * t5;
t56 = pkin(2) * qJD(3);
t58 = t38 * t40;
t16 = (t37 * t43 + t58) * t56;
t63 = t16 * t9;
t60 = Ifges(6,6) * t39;
t59 = t37 * t40;
t33 = t43 * pkin(2) + pkin(3);
t19 = pkin(2) * t58 + t37 * t33;
t55 = qJD(5) * t39;
t54 = qJD(5) * t42;
t6 = t38 * t12 - t37 * t13;
t53 = t57 * t6;
t52 = t70 * t16;
t17 = (t38 * t43 - t59) * t56;
t50 = t57 * t17;
t31 = t37 * pkin(3) + pkin(7);
t49 = t57 * t31;
t18 = -pkin(2) * t59 + t38 * t33;
t47 = (-0.2e1 * Ifges(6,4) * t39 + t42 * t69) * t55 + (0.2e1 * Ifges(6,4) * t42 + t39 * t69) * t54;
t46 = (-mrSges(4,1) * t40 - mrSges(4,2) * t43) * t56;
t45 = -t13 * mrSges(4,1) - t12 * mrSges(4,2) + t9 * t26 + t5 * t70 + t71 * t6;
t34 = Ifges(6,5) * t54;
t32 = -t38 * pkin(3) - pkin(4);
t15 = pkin(7) + t19;
t14 = -pkin(4) - t18;
t10 = t37 * t24 + t38 * t25;
t1 = [0.2e1 * m(6) * (t10 * t53 + t65) + 0.2e1 * m(4) * (t25 * t12 - t24 * t13) + 0.2e1 * m(5) * (t10 * t6 + t65); (-t41 * mrSges(3,1) - t44 * mrSges(3,2)) * qJD(2) + m(6) * (t10 * t50 + t14 * t5 + t15 * t53 + t63) + m(5) * (t17 * t10 - t18 * t5 + t19 * t6 + t63) + m(4) * (t12 * t40 - t13 * t43 + (-t24 * t40 + t25 * t43) * qJD(3)) * pkin(2) + t45; t14 * t67 - 0.2e1 * t17 * mrSges(5,2) + 0.2e1 * m(6) * (t14 * t16 + t15 * t50) + 0.2e1 * m(5) * (-t18 * t16 + t19 * t17) + t47 + 0.2e1 * t51 * t17 + 0.2e1 * t46 + 0.2e1 * t52; m(6) * (t32 * t5 + t6 * t49) + (t37 * t6 - t38 * t5) * t66 + t45; (t14 + t32) * t26 + t52 + t46 + t71 * t17 + m(6) * (t32 * t16 + t17 * t49) + (-t16 * t38 + t17 * t37) * t66 + t47; t32 * t67 + t47; 0; 0; 0; 0; (t10 * t55 - t42 * t6) * mrSges(6,2) + (-t10 * t54 - t39 * t6) * mrSges(6,1); t34 - t48 * t17 + (t27 * t15 - t60) * qJD(5); t34 + (t27 * t31 - t60) * qJD(5); -t26; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
