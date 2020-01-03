% Calculate time derivative of joint inertia matrix for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:11
% DurationCPUTime: 0.73s
% Computational Cost: add. (1051->152), mult. (2238->237), div. (0->0), fcn. (1912->6), ass. (0->66)
t83 = m(4) * pkin(6) + mrSges(4,2);
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t82 = -t61 * pkin(2) - t59 * qJ(3);
t81 = 2 * m(6);
t80 = -2 * pkin(1);
t47 = -pkin(1) + t82;
t79 = 0.2e1 * t47;
t62 = -pkin(2) - pkin(3);
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t45 = -qJ(3) * t56 + t57 * t62;
t44 = -pkin(4) + t45;
t46 = qJ(3) * t57 + t56 * t62;
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t19 = t44 * t60 - t46 * t58;
t37 = -t56 * t58 + t57 * t60;
t10 = qJD(3) * t37 + qJD(5) * t19;
t77 = t10 * mrSges(6,2);
t20 = t44 * t58 + t46 * t60;
t40 = t56 * t60 + t57 * t58;
t11 = -qJD(3) * t40 - qJD(5) * t20;
t76 = t11 * mrSges(6,1);
t75 = Ifges(3,4) - Ifges(4,5);
t74 = pkin(6) - qJ(4);
t72 = qJD(2) * t59;
t29 = -qJD(4) * t61 - t72 * t74;
t49 = t74 * t61;
t31 = qJD(2) * t49 - qJD(4) * t59;
t13 = t57 * t29 + t56 * t31;
t48 = t74 * t59;
t23 = t56 * t48 + t57 * t49;
t71 = qJD(2) * t61;
t73 = qJ(3) * t71 + t59 * qJD(3);
t69 = 0.2e1 * t61;
t39 = -t56 * t61 + t57 * t59;
t64 = t56 * t59 + t57 * t61;
t17 = -t39 * t58 - t60 * t64;
t34 = t39 * qJD(2);
t35 = t64 * qJD(2);
t6 = qJD(5) * t17 + t34 * t58 + t35 * t60;
t18 = t39 * t60 - t58 * t64;
t7 = -qJD(5) * t18 + t34 * t60 - t35 * t58;
t68 = -t7 * mrSges(6,1) + t6 * mrSges(6,2);
t67 = -t34 * mrSges(5,1) + t35 * mrSges(5,2);
t12 = -t29 * t56 + t57 * t31;
t22 = t57 * t48 - t49 * t56;
t36 = t61 * pkin(3) - t47;
t66 = -t61 * mrSges(4,1) - t59 * mrSges(4,3);
t32 = t37 * qJD(5);
t33 = t40 * qJD(5);
t65 = -t33 * mrSges(6,1) - t32 * mrSges(6,2);
t14 = -pkin(7) * t39 + t22;
t15 = -pkin(7) * t64 + t23;
t3 = t14 * t60 - t15 * t58;
t4 = t14 * t58 + t15 * t60;
t24 = t62 * t72 + t73;
t8 = -pkin(7) * t35 + t12;
t9 = pkin(7) * t34 + t13;
t1 = qJD(5) * t3 + t58 * t8 + t60 * t9;
t2 = -qJD(5) * t4 - t58 * t9 + t60 * t8;
t63 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t6 + Ifges(6,6) * t7;
t21 = pkin(4) * t64 + t36;
t16 = -pkin(4) * t34 + t24;
t5 = [0.2e1 * t36 * t67 - 0.2e1 * t64 * Ifges(5,2) * t34 + 0.2e1 * t39 * t35 * Ifges(5,1) + 0.2e1 * t24 * (mrSges(5,1) * t64 + t39 * mrSges(5,2)) + 0.2e1 * t21 * t68 + 0.2e1 * t17 * Ifges(6,2) * t7 + 0.2e1 * t16 * (-t17 * mrSges(6,1) + t18 * mrSges(6,2)) + 0.2e1 * t6 * t18 * Ifges(6,1) + 0.2e1 * m(5) * (t12 * t22 + t13 * t23 + t24 * t36) + (t1 * t4 + t16 * t21 + t2 * t3) * t81 + (m(4) * t79 + 0.2e1 * t66) * (pkin(2) * t72 - t73) + 0.2e1 * (t17 * t6 + t18 * t7) * Ifges(6,4) + 0.2e1 * (t34 * t39 - t35 * t64) * Ifges(5,4) + 0.2e1 * (t1 * t17 - t18 * t2 - t3 * t6 + t4 * t7) * mrSges(6,3) + 0.2e1 * (-t12 * t39 - t13 * t64 - t22 * t35 + t23 * t34) * mrSges(5,3) + (((mrSges(3,2) * t80) - 0.2e1 * t47 * mrSges(4,3) + t69 * t75) * t61 + ((mrSges(3,1) * t80) + mrSges(4,1) * t79 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t69 - 0.2e1 * t75 * t59) * t59) * qJD(2); -t12 * mrSges(5,1) + t13 * mrSges(5,2) - Ifges(5,5) * t35 - Ifges(5,6) * t34 + m(5) * (t12 * t45 + t13 * t46 + (-t22 * t56 + t23 * t57) * qJD(3)) + m(6) * (t1 * t20 + t10 * t4 + t11 * t3 + t19 * t2) + (t10 * t17 - t11 * t18 - t19 * t6 + t20 * t7) * mrSges(6,3) + (t46 * t34 - t45 * t35 + (t39 * t56 - t57 * t64) * qJD(3)) * mrSges(5,3) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t61 + (-mrSges(4,2) * qJ(3) - Ifges(3,6) + Ifges(4,6)) * t59 + (m(4) * t82 - t61 * mrSges(3,1) + t59 * mrSges(3,2) + t66) * pkin(6)) * qJD(2) - t63 + t83 * qJD(3) * t61; (t10 * t20 + t11 * t19) * t81 + 0.2e1 * t77 - 0.2e1 * t76 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3) + t57 * mrSges(5,2) + t56 * mrSges(5,1) + m(5) * (-t45 * t56 + t46 * t57)) * qJD(3); t83 * t71 + (t34 * t56 - t35 * t57) * mrSges(5,3) + m(6) * (t1 * t40 + t2 * t37 - t3 * t33 + t32 * t4) + m(5) * (t12 * t57 + t13 * t56) + (t17 * t32 + t18 * t33 - t37 * t6 + t40 * t7) * mrSges(6,3); m(6) * (t10 * t40 + t11 * t37 - t19 * t33 + t20 * t32) - t65; (t32 * t40 - t33 * t37) * t81; m(5) * t24 + m(6) * t16 + t67 + t68; 0; 0; 0; t63; t76 - t77; t65; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
