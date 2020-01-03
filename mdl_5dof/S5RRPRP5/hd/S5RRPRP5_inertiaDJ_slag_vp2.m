% Calculate time derivative of joint inertia matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:02
% EndTime: 2019-12-31 19:54:04
% DurationCPUTime: 0.92s
% Computational Cost: add. (1515->149), mult. (3374->226), div. (0->0), fcn. (3023->6), ass. (0->70)
t76 = (mrSges(6,2) + mrSges(5,3));
t90 = 2 * t76;
t89 = -mrSges(5,1) - mrSges(6,1);
t60 = cos(pkin(8));
t55 = pkin(2) * t60 + pkin(3);
t61 = sin(qJ(4));
t82 = cos(qJ(4));
t59 = sin(pkin(8));
t84 = pkin(2) * t59;
t41 = t61 * t55 + t82 * t84;
t88 = 2 * m(5);
t87 = 2 * m(6);
t62 = sin(qJ(2));
t63 = cos(qJ(2));
t77 = t60 * t63;
t45 = -t59 * t62 + t77;
t56 = -pkin(2) * t63 - pkin(1);
t31 = -t45 * pkin(3) + t56;
t86 = 0.2e1 * t31;
t85 = 0.2e1 * t56;
t35 = t41 * qJD(4);
t75 = -qJ(3) - pkin(6);
t51 = t75 * t62;
t52 = t75 * t63;
t27 = t60 * t51 + t52 * t59;
t46 = t59 * t63 + t60 * t62;
t23 = -pkin(7) * t46 + t27;
t47 = t59 * t51;
t28 = -t60 * t52 + t47;
t24 = pkin(7) * t45 + t28;
t81 = t24 * t61;
t8 = -t82 * t23 + t81;
t83 = t35 * t8;
t26 = t61 * t45 + t82 * t46;
t80 = t26 * t35;
t70 = qJD(4) * t82;
t73 = t61 * t84;
t34 = -qJD(4) * t73 + t55 * t70;
t30 = qJD(5) + t34;
t79 = t30 * mrSges(6,3);
t78 = t34 * mrSges(5,2);
t69 = qJD(2) * t75;
t22 = t60 * (qJD(3) * t63 + t62 * t69) + t59 * (-t62 * qJD(3) + t63 * t69);
t74 = 0.2e1 * t63;
t57 = qJD(2) * t62 * pkin(2);
t43 = t46 * qJD(2);
t19 = -pkin(7) * t43 + t22;
t21 = -t46 * qJD(3) + (t75 * t77 - t47) * qJD(2);
t44 = t45 * qJD(2);
t64 = -t44 * pkin(7) + t21;
t4 = -qJD(4) * t81 + t82 * t19 + t23 * t70 + t61 * t64;
t9 = t61 * t23 + t82 * t24;
t5 = t9 * qJD(4) + t61 * t19 - t82 * t64;
t72 = t9 * t4 + t5 * t8;
t29 = pkin(3) * t43 + t57;
t71 = t43 * mrSges(4,1) + t44 * mrSges(4,2);
t68 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t67 = t82 * t45 - t61 * t46;
t40 = t82 * t55 - t73;
t16 = t67 * qJD(4) - t61 * t43 + t82 * t44;
t17 = t26 * qJD(4) + t82 * t43 + t61 * t44;
t65 = t89 * t5 + (-mrSges(5,2) + mrSges(6,3)) * t4 + (-Ifges(5,6) + Ifges(6,6)) * t17 + (Ifges(6,4) + Ifges(5,5)) * t16;
t58 = qJD(5) * mrSges(6,3);
t38 = -pkin(4) - t40;
t37 = qJ(5) + t41;
t11 = t17 * mrSges(6,1);
t10 = t16 * mrSges(5,2);
t7 = -pkin(4) * t67 - t26 * qJ(5) + t31;
t6 = pkin(4) * t17 - qJ(5) * t16 - qJD(5) * t26 + t29;
t1 = [0.2e1 * t7 * t11 - 0.2e1 * t45 * Ifges(4,2) * t43 + 0.2e1 * t44 * t46 * Ifges(4,1) + t71 * t85 + 0.2e1 * t6 * (-mrSges(6,1) * t67 - mrSges(6,3) * t26) + 0.2e1 * t29 * (-mrSges(5,1) * t67 + mrSges(5,2) * t26) + t10 * t86 + (t29 * t31 + t72) * t88 + (t6 * t7 + t72) * t87 + 0.2e1 * m(4) * (t21 * t27 + t22 * t28) + 0.2e1 * (-t43 * t46 + t44 * t45) * Ifges(4,4) + 0.2e1 * (-t21 * t46 + t22 * t45 - t27 * t44 - t28 * t43) * mrSges(4,3) + (mrSges(5,1) * t86 + t26 * t68 - 0.2e1 * t76 * t9 - 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t67) * t17 + (-0.2e1 * t7 * mrSges(6,3) - t67 * t68 + t8 * t90 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t26) * t16 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t63) * t74 + (-0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (-mrSges(4,1) * t45 + mrSges(4,2) * t46) + m(4) * pkin(2) * t85 - 0.2e1 * Ifges(3,4) * t62 + (Ifges(3,1) - Ifges(3,2)) * t74) * t62) * qJD(2) + (t5 * t26 + t4 * t67) * t90; (Ifges(3,5) * t63 - Ifges(3,6) * t62 + (-mrSges(3,1) * t63 + mrSges(3,2) * t62) * pkin(6)) * qJD(2) + (m(4) * (t21 * t60 + t22 * t59) + (-t43 * t59 - t44 * t60) * mrSges(4,3)) * pkin(2) + m(6) * (t30 * t9 + t37 * t4 + t38 * t5 + t83) + m(5) * (t34 * t9 + t4 * t41 - t40 * t5 + t83) + (-t16 * t40 - t17 * t41 + t34 * t67 + t80) * mrSges(5,3) + (t16 * t38 - t17 * t37 + t30 * t67 + t80) * mrSges(6,2) + t65 - Ifges(4,6) * t43 + Ifges(4,5) * t44 + t21 * mrSges(4,1) - t22 * mrSges(4,2); -0.2e1 * t78 + 0.2e1 * t79 + 0.2e1 * t89 * t35 + (t30 * t37 + t35 * t38) * t87 + (t34 * t41 - t35 * t40) * t88; m(4) * t57 + m(5) * t29 + m(6) * t6 + t17 * mrSges(5,1) - t16 * mrSges(6,3) + t10 + t11 + t71; 0; 0; m(6) * (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t9) + (-pkin(4) * t16 - qJ(5) * t17 + qJD(5) * t67) * mrSges(6,2) + t65; m(6) * (qJ(5) * t30 + qJD(5) * t37) + t79 + t58 - t78 + (-m(6) * pkin(4) + t89) * t35; 0; 0.2e1 * m(6) * qJ(5) * qJD(5) + 0.2e1 * t58; m(6) * t5 + t16 * mrSges(6,2); m(6) * t35; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
