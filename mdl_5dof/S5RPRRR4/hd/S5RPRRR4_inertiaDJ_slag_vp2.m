% Calculate time derivative of joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:20
% EndTime: 2019-12-05 18:14:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (628->77), mult. (1415->123), div. (0->0), fcn. (967->8), ass. (0->58)
t39 = sin(qJ(5));
t36 = t39 ^ 2;
t42 = cos(qJ(5));
t37 = t42 ^ 2;
t73 = t36 + t37;
t72 = Ifges(6,1) - Ifges(6,2);
t30 = cos(pkin(9)) * pkin(1) + pkin(2);
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t69 = pkin(1) * sin(pkin(9));
t49 = t44 * t30 - t41 * t69;
t18 = pkin(3) + t49;
t19 = t41 * t30 + t44 * t69;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t11 = t40 * t18 + t43 * t19;
t71 = t73 * t43;
t54 = qJD(5) * t42;
t55 = qJD(5) * t39;
t24 = -mrSges(6,1) * t55 - mrSges(6,2) * t54;
t32 = -t43 * pkin(3) - pkin(4);
t17 = t32 * t24;
t57 = t42 * mrSges(6,1);
t25 = t39 * mrSges(6,2) - t57;
t56 = pkin(3) * qJD(4);
t22 = t40 * t25 * t56;
t52 = t43 * t56;
t50 = mrSges(6,3) * t52;
t28 = t36 * t50;
t29 = t37 * t50;
t70 = t22 + t28 + t29 - t17;
t68 = pkin(4) * t24;
t10 = t43 * t18 - t40 * t19;
t15 = t49 * qJD(3);
t16 = t19 * qJD(3);
t5 = qJD(4) * t10 + t43 * t15 - t40 * t16;
t67 = t36 * t5;
t66 = t37 * t5;
t65 = t5 * mrSges(5,2);
t62 = Ifges(6,6) * t39;
t61 = t15 * mrSges(4,2);
t60 = t16 * mrSges(4,1);
t51 = t73 * t5;
t48 = -t40 * mrSges(5,1) - t43 * mrSges(5,2);
t47 = -mrSges(6,1) * t39 - mrSges(6,2) * t42;
t46 = (-0.2e1 * Ifges(6,4) * t39 + t72 * t42) * t55 + (0.2e1 * Ifges(6,4) * t42 + t72 * t39) * t54;
t6 = qJD(4) * t11 + t40 * t15 + t43 * t16;
t1 = t6 * t25;
t2 = mrSges(6,3) * t67;
t3 = mrSges(6,3) * t66;
t4 = t6 * mrSges(5,1);
t8 = -pkin(4) - t10;
t7 = t8 * t24;
t45 = t1 + t2 + t3 - t4 + t46 - t7 - t65;
t35 = Ifges(6,5) * t54;
t31 = t40 * pkin(3) + pkin(8);
t9 = pkin(8) + t11;
t12 = [-0.2e1 * t60 - 0.2e1 * t61 - 0.2e1 * t65 + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 - 0.2e1 * t7 + 0.2e1 * m(6) * (t9 * t51 + t8 * t6) + 0.2e1 * m(5) * (-t10 * t6 + t11 * t5) + 0.2e1 * m(4) * (t19 * t15 - t49 * t16) + t46; 0; 0; t45 + m(6) * (t32 * t6 + (t66 + t67) * t31) + (m(5) * (t40 * t5 - t43 * t6) + (m(6) * (t40 * t8 + t71 * t9) + m(5) * (-t10 * t40 + t11 * t43) + t48) * qJD(4)) * pkin(3) - t61 - t60 + t70; 0; -0.2e1 * t17 + 0.2e1 * t22 + 0.2e1 * t28 + 0.2e1 * t29 + 0.2e1 * (m(6) * (t71 * t31 + t32 * t40) + t48) * t56 + t46; t68 + m(6) * (-pkin(4) * t6 + pkin(8) * t51) + t45; 0; t68 + (m(6) * (-pkin(4) * t40 + t71 * pkin(8)) + t48) * t56 + t46 + t70; t46 + 0.2e1 * t68; t35 + t47 * t5 + (-t9 * t57 + (mrSges(6,2) * t9 - Ifges(6,6)) * t39) * qJD(5); t24; t35 + t47 * t52 + (t25 * t31 - t62) * qJD(5); t35 + (t25 * pkin(8) - t62) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;
