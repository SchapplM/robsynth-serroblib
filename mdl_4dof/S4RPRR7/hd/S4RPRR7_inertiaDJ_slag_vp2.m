% Calculate time derivative of joint inertia matrix for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:53:44
% DurationCPUTime: 0.51s
% Computational Cost: add. (683->125), mult. (1605->196), div. (0->0), fcn. (1419->6), ass. (0->66)
t75 = -2 * mrSges(4,3);
t39 = sin(pkin(7));
t61 = pkin(5) + qJ(2);
t31 = t61 * t39;
t40 = cos(pkin(7));
t32 = t61 * t40;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t47 = -t44 * t31 - t42 * t32;
t73 = -0.2e1 * t47;
t27 = t44 * t39 + t42 * t40;
t72 = -t27 / 0.2e1;
t43 = cos(qJ(4));
t41 = sin(qJ(4));
t69 = Ifges(5,4) * t41;
t33 = Ifges(5,2) * t43 + t69;
t71 = -t33 / 0.2e1;
t70 = mrSges(5,3) * t27;
t68 = Ifges(5,4) * t43;
t24 = t27 * qJD(3);
t67 = Ifges(5,5) * t24;
t66 = Ifges(5,6) * t24;
t26 = t42 * t39 - t44 * t40;
t65 = Ifges(5,6) * t26;
t64 = Ifges(5,6) * t41;
t19 = -t42 * t31 + t44 * t32;
t13 = t27 * qJD(2) + t19 * qJD(3);
t63 = t47 * t13;
t23 = t26 * qJD(3);
t62 = t43 * t23;
t60 = -Ifges(5,5) * t62 + Ifges(5,3) * t24;
t53 = -t40 * pkin(2) - pkin(1);
t15 = t26 * pkin(3) - t27 * pkin(6) + t53;
t6 = t43 * t15 - t41 * t19;
t58 = qJD(4) * t6;
t7 = t41 * t15 + t43 * t19;
t57 = qJD(4) * t7;
t56 = qJD(4) * t41;
t55 = qJD(4) * t43;
t54 = t27 * t56;
t52 = t24 * mrSges(4,1) - t23 * mrSges(4,2);
t51 = -(2 * Ifges(4,4)) - t64;
t50 = mrSges(5,1) * t41 + mrSges(5,2) * t43;
t49 = Ifges(5,1) * t43 - t69;
t48 = -Ifges(5,2) * t41 + t68;
t46 = -t41 * t23 + t27 * t55;
t45 = t54 + t62;
t36 = Ifges(5,5) * t55;
t34 = Ifges(5,1) * t41 + t68;
t30 = t49 * qJD(4);
t29 = t48 * qJD(4);
t28 = t50 * qJD(4);
t17 = t26 * mrSges(5,1) - t43 * t70;
t16 = -t26 * mrSges(5,2) - t41 * t70;
t14 = t24 * pkin(3) + t23 * pkin(6);
t12 = -t26 * qJD(2) + t47 * qJD(3);
t11 = Ifges(5,5) * t26 + t49 * t27;
t10 = t48 * t27 + t65;
t9 = -t24 * mrSges(5,2) - t46 * mrSges(5,3);
t8 = t24 * mrSges(5,1) + t45 * mrSges(5,3);
t5 = t46 * mrSges(5,1) - t45 * mrSges(5,2);
t4 = -t45 * Ifges(5,1) - t46 * Ifges(5,4) + t67;
t3 = -t45 * Ifges(5,4) - t46 * Ifges(5,2) + t66;
t2 = -t41 * t12 + t43 * t14 - t57;
t1 = t43 * t12 + t41 * t14 + t58;
t18 = [t19 * t24 * t75 + 0.2e1 * t53 * t52 + 0.2e1 * t1 * t16 + 0.2e1 * t2 * t17 + t5 * t73 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + 0.2e1 * m(5) * (t7 * t1 + t6 * t2 - t63) + 0.2e1 * m(4) * (t19 * t12 - t63) - (mrSges(4,3) * t73 - t41 * t10 + t43 * t11) * t23 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t39 ^ 2 + t40 ^ 2) + (t12 * t75 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t24 - t51 * t23 + t60) * t26 + (t43 * t4 - t41 * t3 - 0.2e1 * Ifges(4,1) * t23 + (Ifges(5,5) * t43 + t51) * t24 + (t26 * (-Ifges(5,5) * t41 - Ifges(5,6) * t43) - t43 * t10 - t41 * t11) * qJD(4) + 0.2e1 * (mrSges(4,3) + t50) * t13) * t27; m(5) * (t41 * t1 + t43 * t2 + (-t41 * t6 + t43 * t7) * qJD(4)) + t16 * t55 + t41 * t9 - t17 * t56 + t43 * t8 + t52; 0; -t47 * t28 + t26 * t36 / 0.2e1 - Ifges(4,5) * t23 - Ifges(4,6) * t24 - t12 * mrSges(4,2) - t13 * mrSges(4,1) + (-m(5) * t13 - t5) * pkin(3) + (t29 * t72 - t23 * t71 - t2 * mrSges(5,3) + t4 / 0.2e1 + t67 / 0.2e1 + t13 * mrSges(5,2) + (-t10 / 0.2e1 - t65 / 0.2e1 + t34 * t72 - t7 * mrSges(5,3)) * qJD(4) + (m(5) * (-t2 - t57) - t8 - qJD(4) * t16) * pkin(6)) * t41 + (t27 * t30 / 0.2e1 - t23 * t34 / 0.2e1 + t1 * mrSges(5,3) + t66 / 0.2e1 - t13 * mrSges(5,1) + t3 / 0.2e1 + (t11 / 0.2e1 + t27 * t71 - t6 * mrSges(5,3)) * qJD(4) + (m(5) * (t1 - t58) + t9 - qJD(4) * t17) * pkin(6)) * t43; 0; -0.2e1 * pkin(3) * t28 + t43 * t29 + t41 * t30 + (-t41 * t33 + t43 * t34) * qJD(4); t2 * mrSges(5,1) - t1 * mrSges(5,2) - Ifges(5,5) * t54 - t46 * Ifges(5,6) + t60; -t28; t36 + (-t64 + (-t43 * mrSges(5,1) + t41 * mrSges(5,2)) * pkin(6)) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t18(1), t18(2), t18(4), t18(7); t18(2), t18(3), t18(5), t18(8); t18(4), t18(5), t18(6), t18(9); t18(7), t18(8), t18(9), t18(10);];
Mq = res;
