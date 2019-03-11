% Calculate time derivative of joint inertia matrix for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:47
% EndTime: 2019-03-09 02:58:50
% DurationCPUTime: 0.93s
% Computational Cost: add. (872->214), mult. (1651->298), div. (0->0), fcn. (1037->4), ass. (0->100)
t107 = -2 * mrSges(6,3);
t45 = cos(qJ(3));
t106 = 0.2e1 * t45;
t105 = 2 * qJ(2);
t44 = cos(qJ(6));
t94 = -t44 / 0.2e1;
t46 = -pkin(3) - pkin(4);
t38 = -pkin(8) + t46;
t43 = sin(qJ(3));
t104 = t38 * t43;
t42 = sin(qJ(6));
t91 = mrSges(7,3) * t43;
t23 = -mrSges(7,2) * t45 - t42 * t91;
t24 = mrSges(7,1) * t45 - t44 * t91;
t103 = (t42 * t23 + t44 * t24) * qJD(6);
t76 = t45 * qJ(4) - qJ(2);
t10 = pkin(5) * t45 + t104 + t76;
t47 = (-pkin(1) - pkin(7));
t75 = qJ(5) + t47;
t26 = t75 * t45;
t3 = t10 * t44 + t26 * t42;
t4 = t10 * t42 - t26 * t44;
t62 = qJD(3) * t75;
t15 = -t45 * qJD(5) + t43 * t62;
t41 = qJ(4) + pkin(5);
t65 = t45 * qJD(4) - qJD(2);
t7 = (t38 * t45 - t41 * t43) * qJD(3) + t65;
t1 = t3 * qJD(6) + t15 * t44 + t42 * t7;
t2 = -t4 * qJD(6) - t15 * t42 + t44 * t7;
t60 = t1 * t44 - t2 * t42;
t102 = (t3 * t44 + t4 * t42) * qJD(6) - t60;
t101 = 2 * mrSges(6,1);
t16 = t43 * qJD(5) + t45 * t62;
t100 = 0.2e1 * t16;
t25 = t75 * t43;
t99 = 0.2e1 * t25;
t98 = 0.2e1 * pkin(3) * t43 - 0.2e1 * t76;
t97 = 0.2e1 * t41;
t96 = 0.2e1 * t43;
t95 = -t42 / 0.2e1;
t93 = m(5) + m(6);
t92 = t3 * t42;
t90 = Ifges(7,4) * t42;
t89 = Ifges(7,4) * t44;
t88 = Ifges(7,5) * t42;
t87 = Ifges(7,5) * t44;
t86 = Ifges(7,6) * t44;
t85 = t16 * t25;
t29 = -Ifges(7,2) * t44 - t90;
t84 = t42 * t29;
t30 = -Ifges(7,1) * t42 - t89;
t83 = t44 * t30;
t82 = -mrSges(4,1) - mrSges(5,1);
t81 = (mrSges(5,3) - mrSges(4,2));
t28 = mrSges(7,1) * t44 - mrSges(7,2) * t42;
t80 = t28 + mrSges(6,1);
t36 = qJD(4) * t43;
t71 = qJD(3) * t45;
t79 = qJ(4) * t71 + t36;
t78 = t42 ^ 2 + t44 ^ 2;
t77 = qJ(4) * t43;
t74 = qJD(3) * t25;
t73 = qJD(3) * t43;
t72 = qJD(3) * t44;
t70 = qJD(4) * t25;
t69 = qJD(6) * t42;
t68 = qJD(6) * t43;
t67 = qJD(6) * t44;
t66 = qJD(6) * t45;
t64 = t42 * t71;
t63 = t44 * t71;
t61 = m(5) * t47 - mrSges(5,2) + mrSges(6,3);
t58 = mrSges(7,1) * t42 + mrSges(7,2) * t44;
t57 = Ifges(7,1) * t44 - t90;
t56 = -Ifges(7,2) * t42 + t89;
t55 = -t86 - t88;
t50 = -t42 * t68 + t63;
t11 = -mrSges(7,1) * t73 - t50 * mrSges(7,3);
t49 = t43 * t67 + t64;
t12 = mrSges(7,2) * t73 - t49 * mrSges(7,3);
t54 = t42 * t11 - t44 * t12;
t53 = t44 * t23 - t42 * t24;
t51 = -Ifges(7,6) * t42 - (2 * Ifges(4,4)) - (2 * Ifges(6,4)) + (2 * Ifges(5,5));
t48 = -m(7) * t102 - t54;
t34 = Ifges(7,6) * t69;
t33 = mrSges(6,2) * t71;
t31 = Ifges(7,5) * t63;
t22 = t57 * qJD(6);
t21 = t56 * qJD(6);
t20 = t58 * qJD(6);
t19 = t58 * t43;
t18 = t46 * t43 + t76;
t17 = (pkin(3) * t45 + t77) * qJD(3) - t65;
t14 = t45 * Ifges(7,5) + t57 * t43;
t13 = t45 * Ifges(7,6) + t56 * t43;
t9 = (t46 * t45 - t77) * qJD(3) + t65;
t8 = t49 * mrSges(7,1) + t50 * mrSges(7,2);
t6 = t30 * t68 + (-Ifges(7,5) * t43 + t57 * t45) * qJD(3);
t5 = t29 * t68 + (-Ifges(7,6) * t43 + t56 * t45) * qJD(3);
t27 = [m(5) * t17 * t98 + 0.2e1 * t1 * t23 + 0.2e1 * t3 * t11 + 0.2e1 * t4 * t12 + t19 * t100 + 0.2e1 * t18 * t33 + 0.2e1 * t2 * t24 + t8 * t99 + 0.2e1 * m(7) * (t1 * t4 + t2 * t3 + t85) + 0.2e1 * m(6) * (-t15 * t26 + t18 * t9 + t85) + (-0.2e1 * t17 * mrSges(5,3) + t9 * t101 + t15 * t107 + t31) * t45 + (0.2e1 * t17 * mrSges(5,1) + 0.2e1 * t9 * mrSges(6,2) + mrSges(6,3) * t100 - t42 * t5 + t44 * t6 + (-t44 * t13 - t42 * t14 + t45 * t55) * qJD(6)) * t43 + (mrSges(4,1) * t96 + mrSges(4,2) * t106 + (2 * mrSges(3,3)) + ((m(3) + m(4)) * t105)) * qJD(2) + ((mrSges(4,1) * t105 + mrSges(5,1) * t98 + mrSges(6,3) * t99 - t42 * t13 + t44 * t14 + t51 * t45) * t45 + (-0.2e1 * t18 * mrSges(6,1) - 0.2e1 * qJ(2) * mrSges(4,2) + mrSges(5,3) * t98 + t26 * t107 + (-t51 - t87) * t43 + (-Ifges(4,1) - Ifges(5,1) + Ifges(6,1) + Ifges(4,2) - Ifges(6,2) + Ifges(5,3) - Ifges(7,3)) * t106) * t43) * qJD(3); (t8 + t53 * qJD(3) + m(7) * (-qJD(3) * t92 + t4 * t72 + t16) + m(6) * (-qJD(3) * t26 + t16)) * t43 + (qJD(3) * t19 + t103 + m(7) * (t3 * t67 + t4 * t69 - t60 + t74) + m(6) * (-t15 + t74) + t54) * t45; m(7) * (0.1e1 - t78) * t71 * t96; -t14 * t67 / 0.2e1 + t5 * t94 + t45 * (-Ifges(7,5) * t67 + t34) / 0.2e1 + t41 * t8 - t25 * t20 + t15 * mrSges(6,2) + qJD(4) * t19 + (qJD(6) * t13 / 0.2e1 - t6 / 0.2e1) * t42 + t80 * t16 + m(6) * (qJ(4) * t16 + t15 * t46 + t70) + m(7) * (t16 * t41 + t70) + t102 * mrSges(7,3) + (-t23 * t69 - t24 * t67 + t48) * t38 + (t22 * t94 - t21 * t95 + (t29 * t94 + t30 * t95) * qJD(6) + t61 * qJD(4)) * t43 + ((t88 / 0.2e1 + t86 / 0.2e1 - Ifges(6,6) - Ifges(5,4) - Ifges(4,5) + t46 * mrSges(6,3) + pkin(3) * mrSges(5,2) + (-m(5) * pkin(3) + t82) * t47) * t43 + (-Ifges(6,5) + Ifges(5,6) - Ifges(4,6) + t83 / 0.2e1 - t84 / 0.2e1 + (t81 * t47) + t61 * qJ(4)) * t45) * qJD(3); -t43 * t20 + ((t80 + t81) * t45 + (-t78 * mrSges(7,3) + mrSges(6,2) + t82) * t43) * qJD(3) + m(6) * (t46 * t73 + t79) + m(7) * (t36 + (t78 * t104 + t41 * t45) * qJD(3)) + m(5) * (-pkin(3) * t73 + t79); -t20 * t97 + t21 * t44 + t22 * t42 + (-t83 + t84) * qJD(6) + (m(7) * t97 + 0.2e1 * qJ(4) * t93 + 0.2e1 * mrSges(5,3) + t101 + 0.2e1 * t28) * qJD(4); m(6) * t15 + t61 * t73 - t103 + t48; (m(7) * t78 + t93) * t73; 0; 0; -mrSges(6,1) * t73 + t44 * t11 + t42 * t12 + t33 + t53 * qJD(6) + m(7) * (t1 * t42 + t2 * t44 + (t4 * t44 - t92) * qJD(6)) + m(6) * t9; 0; 0; 0; 0; -Ifges(7,6) * t64 + mrSges(7,1) * t2 - mrSges(7,2) * t1 + t31 + (-Ifges(7,3) * qJD(3) + t55 * qJD(6)) * t43; (-t42 * t66 - t43 * t72) * mrSges(7,2) + (-t42 * t73 + t44 * t66) * mrSges(7,1); t34 + (-t28 * t38 - t87) * qJD(6); -t28 * qJD(6); -t20; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t27(1) t27(2) t27(4) t27(7) t27(11) t27(16); t27(2) t27(3) t27(5) t27(8) t27(12) t27(17); t27(4) t27(5) t27(6) t27(9) t27(13) t27(18); t27(7) t27(8) t27(9) t27(10) t27(14) t27(19); t27(11) t27(12) t27(13) t27(14) t27(15) t27(20); t27(16) t27(17) t27(18) t27(19) t27(20) t27(21);];
Mq  = res;
