% Calculate time derivative of joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:27
% EndTime: 2019-12-05 16:02:31
% DurationCPUTime: 0.85s
% Computational Cost: add. (612->187), mult. (1627->304), div. (0->0), fcn. (1267->8), ass. (0->92)
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t32 = t45 * mrSges(5,1) + t48 * mrSges(5,2);
t102 = mrSges(4,3) + t32;
t44 = sin(qJ(5));
t47 = cos(qJ(5));
t78 = t44 ^ 2 + t47 ^ 2;
t43 = cos(pkin(5));
t42 = sin(pkin(5));
t49 = cos(qJ(2));
t83 = t42 * t49;
t18 = t43 * t45 + t48 * t83;
t67 = t45 * t83;
t19 = t43 * t48 - t67;
t46 = sin(qJ(2));
t84 = t42 * t46;
t66 = qJD(2) * t84;
t76 = qJD(4) * t18;
t8 = t45 * t66 - t76;
t73 = qJD(4) * t48;
t9 = -qJD(4) * t67 + t43 * t73 - t48 * t66;
t91 = t48 * t9;
t101 = (t18 * t45 + t19 * t48) * qJD(4) + t45 * t8 - t91;
t31 = -mrSges(6,1) * t47 + mrSges(6,2) * t44;
t100 = -m(6) * pkin(4) - mrSges(5,1) + t31;
t99 = 0.2e1 * t45;
t50 = (-pkin(2) - pkin(7));
t98 = 2 * t50;
t97 = m(5) / 0.2e1;
t96 = m(6) / 0.2e1;
t95 = -t44 / 0.2e1;
t94 = pkin(8) * t48;
t93 = t18 * t9;
t92 = t45 * pkin(4);
t90 = Ifges(6,4) * t44;
t89 = Ifges(6,4) * t47;
t88 = Ifges(6,5) * t44;
t87 = Ifges(6,6) * t44;
t86 = Ifges(6,6) * t45;
t85 = Ifges(6,6) * t47;
t82 = t45 * t50;
t81 = t48 * mrSges(6,3);
t77 = qJD(2) * t49;
t65 = t42 * t77;
t80 = qJ(3) * t65 + qJD(3) * t84;
t75 = qJD(4) * t45;
t62 = t44 * t75;
t79 = Ifges(6,6) * t62 + Ifges(6,3) * t73;
t74 = qJD(4) * t47;
t28 = qJ(3) + t92 - t94;
t14 = t47 * t28 - t44 * t82;
t72 = qJD(5) * t14;
t15 = t44 * t28 + t47 * t82;
t71 = qJD(5) * t15;
t70 = qJD(5) * t44;
t69 = qJD(5) * t47;
t68 = qJD(5) * t48;
t64 = t50 * t73;
t63 = t47 * t68;
t61 = -Ifges(6,5) * t47 + (2 * Ifges(5,4));
t21 = qJD(3) + (pkin(4) * t48 + pkin(8) * t45) * qJD(4);
t3 = t44 * t21 + t47 * t64 + t72;
t60 = t3 - t72;
t11 = t19 * t47 + t44 * t84;
t1 = -qJD(5) * t11 - t44 * t8 + t47 * t65;
t10 = -t19 * t44 + t47 * t84;
t2 = qJD(5) * t10 + t44 * t65 + t47 * t8;
t59 = -t1 * t44 + t2 * t47;
t57 = mrSges(6,1) * t44 + mrSges(6,2) * t47;
t56 = Ifges(6,1) * t47 - t90;
t34 = Ifges(6,1) * t44 + t89;
t55 = -Ifges(6,2) * t44 + t89;
t33 = Ifges(6,2) * t47 + t90;
t53 = t44 * t68 + t45 * t74;
t52 = t62 - t63;
t39 = Ifges(6,5) * t69;
t27 = mrSges(6,1) * t45 - t47 * t81;
t26 = -mrSges(6,2) * t45 - t44 * t81;
t25 = t56 * qJD(5);
t24 = t55 * qJD(5);
t23 = (mrSges(5,1) * t48 - mrSges(5,2) * t45) * qJD(4);
t22 = t57 * qJD(5);
t20 = t57 * t48;
t17 = Ifges(6,5) * t45 + t56 * t48;
t16 = t55 * t48 + t86;
t13 = -mrSges(6,2) * t73 + t52 * mrSges(6,3);
t12 = mrSges(6,1) * t73 + t53 * mrSges(6,3);
t7 = -t52 * mrSges(6,1) - t53 * mrSges(6,2);
t6 = -t34 * t68 + (Ifges(6,5) * t48 - t56 * t45) * qJD(4);
t5 = -t33 * t68 + (Ifges(6,6) * t48 - t55 * t45) * qJD(4);
t4 = t47 * t21 - t44 * t64 - t71;
t29 = [0.2e1 * m(5) * (t42 ^ 2 * t46 * t77 + t19 * t8 + t93) + 0.2e1 * m(6) * (t1 * t10 + t11 * t2 + t93); t1 * t27 + t10 * t12 + t11 * t13 + t18 * t7 + t2 * t26 + t9 * t20 + (t46 * t23 + ((-mrSges(3,1) + mrSges(4,2)) * t46 + (-mrSges(3,2) + t102) * t49) * qJD(2)) * t42 + m(4) * (-pkin(2) * t66 + t80) + m(5) * t80 + m(6) * (t14 * t1 + t4 * t10 + t3 * t11 + t15 * t2) + (t101 * t97 + (t18 * t75 - t91) * t96) * t98 - t101 * mrSges(5,3); 0.2e1 * qJ(3) * t23 + 0.2e1 * t4 * t27 + 0.2e1 * t14 * t12 + 0.2e1 * t3 * t26 + 0.2e1 * t15 * t13 + 0.2e1 * m(6) * (t14 * t4 + t15 * t3) + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t102) * qJD(3) + ((t44 * t16 - t47 * t17 + t20 * t98 + t61 * t45) * qJD(4) + t79) * t45 + (-t44 * t5 + t47 * t6 - 0.2e1 * t50 * t7 + (t45 * (-t85 - t88) - t44 * t17 - t47 * t16) * qJD(5) + ((-t61 - t87) * t48 + (-0.2e1 * m(6) * (t50 ^ 2) - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(6,3)) * t45) * qJD(4)) * t48; m(4) * t66 + 0.2e1 * ((qJD(4) * t19 - t9) * t97 + (-qJD(4) * t10 * t44 + t11 * t74 - t9) * t96) * t48 + ((t8 + t76) * t97 + (-t10 * t69 - t11 * t70 + t59 + t76) * t96) * t99; (-t7 + (m(6) * (-t14 * t44 + t15 * t47) + t47 * t26 - t44 * t27) * qJD(4)) * t48 + (m(6) * (-t14 * t69 - t15 * t70 + t3 * t47 - t4 * t44 - 0.2e1 * t64) - t26 * t70 + t47 * t13 - t27 * t69 - t44 * t12 + qJD(4) * t20) * t45; m(6) * (-0.1e1 + t78) * t73 * t99; -t8 * mrSges(5,2) + t18 * t22 + (m(6) * pkin(8) + mrSges(6,3)) * ((-t10 * t47 - t11 * t44) * qJD(5) + t59) + t100 * t9; -pkin(4) * t7 + (t39 / 0.2e1 + (t100 * t50 - Ifges(5,5)) * qJD(4)) * t45 + (t6 / 0.2e1 - t4 * mrSges(6,3) + t33 * t75 / 0.2e1 + (-t86 / 0.2e1 - t16 / 0.2e1 - t15 * mrSges(6,3)) * qJD(5) + (m(6) * (-t4 - t71) - qJD(5) * t26 - t12) * pkin(8)) * t44 + (t5 / 0.2e1 + qJD(5) * t17 / 0.2e1 - t34 * t75 / 0.2e1 + t60 * mrSges(6,3) + (m(6) * t60 - qJD(5) * t27 + t13) * pkin(8)) * t47 + (t47 * t25 / 0.2e1 + t24 * t95 - t50 * t22 + (-t47 * t33 / 0.2e1 + t34 * t95) * qJD(5) + (t88 / 0.2e1 + t85 / 0.2e1 - Ifges(5,6) - t50 * mrSges(5,2)) * qJD(4)) * t48; -t48 * t22 + (t45 * t31 + m(6) * (t78 * t94 - t92) + t78 * t81 - t32) * qJD(4); -0.2e1 * pkin(4) * t22 + t24 * t47 + t25 * t44 + (-t33 * t44 + t34 * t47) * qJD(5); mrSges(6,1) * t1 - mrSges(6,2) * t2; mrSges(6,1) * t4 - mrSges(6,2) * t3 - t53 * Ifges(6,5) - Ifges(6,6) * t63 + t79; (t45 * t70 - t47 * t73) * mrSges(6,2) + (-t44 * t73 - t45 * t69) * mrSges(6,1); t39 + (t31 * pkin(8) - t87) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t29(1), t29(2), t29(4), t29(7), t29(11); t29(2), t29(3), t29(5), t29(8), t29(12); t29(4), t29(5), t29(6), t29(9), t29(13); t29(7), t29(8), t29(9), t29(10), t29(14); t29(11), t29(12), t29(13), t29(14), t29(15);];
Mq = res;
