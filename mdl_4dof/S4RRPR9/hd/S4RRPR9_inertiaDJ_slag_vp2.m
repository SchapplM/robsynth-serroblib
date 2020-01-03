% Calculate time derivative of joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:22
% EndTime: 2019-12-31 17:09:24
% DurationCPUTime: 0.91s
% Computational Cost: add. (749->176), mult. (1913->293), div. (0->0), fcn. (1532->6), ass. (0->85)
t69 = sin(pkin(7));
t70 = cos(pkin(7));
t102 = mrSges(4,1) * t69 + mrSges(4,2) * t70;
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t75 = t69 * t71 - t70 * t73;
t46 = t75 * qJD(4);
t101 = 2 * m(4);
t100 = 2 * m(5);
t99 = -2 * pkin(1);
t98 = 2 * pkin(5);
t97 = -t75 / 0.2e1;
t54 = t69 * t73 + t70 * t71;
t96 = t54 / 0.2e1;
t95 = t70 / 0.2e1;
t92 = Ifges(4,4) * t69;
t91 = Ifges(4,4) * t70;
t72 = sin(qJ(2));
t90 = t69 * t72;
t74 = cos(qJ(2));
t89 = t69 * t74;
t88 = t70 * t72;
t87 = t70 * t74;
t86 = pkin(6) + qJ(3);
t47 = t54 * qJD(4);
t85 = -Ifges(5,5) * t46 - Ifges(5,6) * t47;
t82 = qJD(2) * t74;
t42 = t102 * t82;
t45 = -t72 * qJD(3) + (pkin(2) * t72 - qJ(3) * t74) * qJD(2);
t83 = qJD(2) * t72;
t80 = pkin(5) * t83;
t27 = t70 * t45 + t69 * t80;
t58 = -pkin(2) * t74 - t72 * qJ(3) - pkin(1);
t37 = pkin(5) * t87 + t69 * t58;
t18 = -t72 * t47 - t75 * t82;
t19 = t72 * t46 - t54 * t82;
t81 = Ifges(5,5) * t18 + Ifges(5,6) * t19 + Ifges(5,3) * t83;
t79 = pkin(3) * t69 + pkin(5);
t5 = -t19 * mrSges(5,1) + t18 * mrSges(5,2);
t78 = t70 * Ifges(4,1) - t92;
t77 = -t69 * Ifges(4,2) + t91;
t76 = -Ifges(4,5) * t70 + Ifges(4,6) * t69;
t52 = t70 * t58;
t24 = -pkin(6) * t88 + t52 + (-pkin(5) * t69 - pkin(3)) * t74;
t29 = -pkin(6) * t90 + t37;
t6 = t24 * t73 - t29 * t71;
t7 = t24 * t71 + t29 * t73;
t59 = t86 * t69;
t60 = t86 * t70;
t30 = -t59 * t73 - t60 * t71;
t31 = -t59 * t71 + t60 * t73;
t65 = -pkin(3) * t70 - pkin(2);
t57 = t79 * t72;
t56 = -mrSges(4,1) * t74 - mrSges(4,3) * t88;
t55 = mrSges(4,2) * t74 - mrSges(4,3) * t90;
t50 = t79 * t82;
t49 = (mrSges(4,1) * t72 - mrSges(4,3) * t87) * qJD(2);
t48 = (-mrSges(4,2) * t72 - mrSges(4,3) * t89) * qJD(2);
t41 = t75 * t72;
t40 = t54 * t72;
t38 = t69 * t45;
t36 = -pkin(5) * t89 + t52;
t35 = (Ifges(4,5) * t72 + t78 * t74) * qJD(2);
t34 = (Ifges(4,6) * t72 + t77 * t74) * qJD(2);
t33 = -mrSges(5,1) * t74 + t41 * mrSges(5,3);
t32 = mrSges(5,2) * t74 - t40 * mrSges(5,3);
t28 = -t70 * t80 + t38;
t26 = Ifges(5,1) * t54 - Ifges(5,4) * t75;
t25 = Ifges(5,4) * t54 - Ifges(5,2) * t75;
t23 = -Ifges(5,1) * t46 - Ifges(5,4) * t47;
t22 = -Ifges(5,4) * t46 - Ifges(5,2) * t47;
t21 = mrSges(5,1) * t47 - mrSges(5,2) * t46;
t20 = t38 + (-pkin(5) * t88 - pkin(6) * t89) * qJD(2);
t14 = (pkin(3) * t72 - pkin(6) * t87) * qJD(2) + t27;
t13 = -Ifges(5,1) * t41 - Ifges(5,4) * t40 - Ifges(5,5) * t74;
t12 = -Ifges(5,4) * t41 - Ifges(5,2) * t40 - Ifges(5,6) * t74;
t11 = -t54 * qJD(3) - t31 * qJD(4);
t10 = -t75 * qJD(3) + t30 * qJD(4);
t9 = -mrSges(5,2) * t83 + mrSges(5,3) * t19;
t8 = mrSges(5,1) * t83 - mrSges(5,3) * t18;
t4 = Ifges(5,1) * t18 + Ifges(5,4) * t19 + Ifges(5,5) * t83;
t3 = Ifges(5,4) * t18 + Ifges(5,2) * t19 + Ifges(5,6) * t83;
t2 = -t7 * qJD(4) + t14 * t73 - t20 * t71;
t1 = t6 * qJD(4) + t14 * t71 + t20 * t73;
t15 = [-t74 * t81 + 0.2e1 * t27 * t56 + 0.2e1 * t57 * t5 + 0.2e1 * t37 * t48 + 0.2e1 * t36 * t49 + 0.2e1 * t50 * (t40 * mrSges(5,1) - t41 * mrSges(5,2)) + 0.2e1 * t28 * t55 + 0.2e1 * t1 * t32 + 0.2e1 * t2 * t33 - t40 * t3 - t41 * t4 + t19 * t12 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + t18 * t13 + (t36 * t27 + t37 * t28) * t101 + (t1 * t7 + t2 * t6 + t50 * t57) * t100 + (-t69 * t34 + t70 * t35 + t42 * t98) * t72 + (((mrSges(3,1) * t99) - Ifges(5,5) * t41 - Ifges(5,6) * t40 + (-(2 * Ifges(3,4)) - t76) * t72) * t72 + ((mrSges(3,2) * t99) + 0.2e1 * (Ifges(3,4) + t76) * t74 + ((pkin(5) ^ 2 * t101) + t102 * t98 - t69 * t77 + t70 * t78 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t72) * t74) * qJD(2); -t74 * t85 / 0.2e1 + t57 * t21 + t65 * t5 - t46 * t13 / 0.2e1 - t47 * t12 / 0.2e1 + t3 * t97 + t50 * (mrSges(5,1) * t75 + t54 * mrSges(5,2)) + t4 * t96 + t10 * t32 + t11 * t33 - t40 * t22 / 0.2e1 - t41 * t23 / 0.2e1 - pkin(2) * t42 + t19 * t25 / 0.2e1 + t18 * t26 / 0.2e1 + t30 * t8 + t31 * t9 + m(5) * (t1 * t31 + t10 * t7 + t11 * t6 + t2 * t30 + t50 * t65) + (m(4) * (qJ(3) * t28 + qJD(3) * t37) + t34 / 0.2e1 + qJD(3) * t55 + t28 * mrSges(4,3) + qJ(3) * t48) * t70 + (m(4) * (-qJ(3) * t27 - qJD(3) * t36) + t35 / 0.2e1 - qJD(3) * t56 - t27 * mrSges(4,3) - qJ(3) * t49) * t69 + (-t1 * t75 - t2 * t54 + t46 * t6 - t47 * t7) * mrSges(5,3) + (((pkin(5) * mrSges(3,2)) + Ifges(4,5) * t69 / 0.2e1 + Ifges(4,6) * t95 + Ifges(5,5) * t96 + Ifges(5,6) * t97 - Ifges(3,6)) * t72 + ((Ifges(4,1) * t69 + t91) * t95 - t69 * (Ifges(4,2) * t70 + t92) / 0.2e1 + Ifges(3,5) + (-m(4) * pkin(2) - t70 * mrSges(4,1) + t69 * mrSges(4,2) - mrSges(3,1)) * pkin(5)) * t74) * qJD(2); (t10 * t31 + t11 * t30) * t100 - t46 * t26 + t54 * t23 + 0.2e1 * t65 * t21 - t47 * t25 - t75 * t22 + 0.2e1 * (-t10 * t75 - t11 * t54 + t30 * t46 - t31 * t47) * mrSges(5,3) + (qJ(3) * t101 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t69 ^ 2 + t70 ^ 2); m(4) * pkin(5) * t82 + m(5) * t50 + t42 + t5; t21; 0; mrSges(5,1) * t2 - mrSges(5,2) * t1 + t81; mrSges(5,1) * t11 - mrSges(5,2) * t10 + t85; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1), t15(2), t15(4), t15(7); t15(2), t15(3), t15(5), t15(8); t15(4), t15(5), t15(6), t15(9); t15(7), t15(8), t15(9), t15(10);];
Mq = res;
