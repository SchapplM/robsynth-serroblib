% Calculate time derivative of joint inertia matrix for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:44
% EndTime: 2019-03-09 01:50:46
% DurationCPUTime: 0.97s
% Computational Cost: add. (720->197), mult. (1402->289), div. (0->0), fcn. (864->4), ass. (0->89)
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t79 = t45 ^ 2 + t47 ^ 2;
t44 = sin(qJ(6));
t46 = cos(qJ(6));
t68 = qJD(6) * t46;
t71 = qJD(4) * t47;
t50 = t44 * t71 + t45 * t68;
t83 = mrSges(6,3) - mrSges(5,2);
t109 = m(6) * qJ(5) + t83;
t108 = 2 * qJ(5);
t107 = m(6) + m(7);
t105 = 2 * m(7);
t104 = 2 * mrSges(6,3);
t102 = -t44 / 0.2e1;
t101 = t44 / 0.2e1;
t100 = t46 / 0.2e1;
t99 = pkin(4) + pkin(8);
t42 = qJ(2) - pkin(7);
t98 = pkin(5) - t42;
t97 = Ifges(7,4) * t44;
t96 = Ifges(7,4) * t46;
t95 = Ifges(7,5) * t44;
t94 = Ifges(7,6) * t46;
t93 = Ifges(7,6) * t47;
t53 = Ifges(7,1) * t44 + t96;
t12 = Ifges(7,5) * t47 + t45 * t53;
t92 = t44 * t12;
t29 = Ifges(7,1) * t46 - t97;
t91 = t44 * t29;
t90 = t44 * t45;
t89 = t44 * t99;
t88 = t45 * t46;
t52 = Ifges(7,2) * t46 + t97;
t11 = t45 * t52 + t93;
t87 = t46 * t11;
t28 = -Ifges(7,2) * t44 + t96;
t86 = t46 * t28;
t85 = t46 * t99;
t84 = -mrSges(5,1) + mrSges(6,2);
t43 = pkin(1) + qJ(3);
t82 = t79 * qJD(2) * t42;
t80 = t44 ^ 2 + t46 ^ 2;
t78 = qJ(5) * t45;
t76 = qJD(2) * t47;
t75 = qJD(3) * t43;
t74 = qJD(4) * t44;
t73 = qJD(4) * t45;
t72 = qJD(4) * t46;
t70 = qJD(6) * t44;
t69 = qJD(6) * t45;
t67 = qJD(6) * t47;
t66 = 0.2e1 * t47;
t62 = t46 * t71;
t65 = t50 * Ifges(7,5) + Ifges(7,6) * t62;
t61 = m(7) * t80;
t60 = qJD(4) * t98;
t59 = -qJD(5) * t47 + qJD(3);
t58 = m(6) * pkin(4) - t84;
t13 = -t45 * t60 - t76;
t23 = t45 * pkin(4) - qJ(5) * t47 + t43;
t16 = pkin(8) * t45 + t23;
t25 = t98 * t47;
t6 = -t16 * t44 + t25 * t46;
t8 = (t99 * t47 + t78) * qJD(4) + t59;
t1 = qJD(6) * t6 + t13 * t44 + t46 * t8;
t7 = t16 * t46 + t25 * t44;
t2 = -qJD(6) * t7 + t13 * t46 - t44 * t8;
t56 = t1 * t44 + t2 * t46;
t55 = t6 * t44 - t7 * t46;
t54 = mrSges(7,1) * t46 - mrSges(7,2) * t44;
t27 = mrSges(7,1) * t44 + mrSges(7,2) * t46;
t51 = -t44 * t69 + t62;
t10 = -mrSges(7,1) * t73 - mrSges(7,3) * t50;
t21 = mrSges(7,1) * t47 - mrSges(7,3) * t90;
t22 = -mrSges(7,2) * t47 + mrSges(7,3) * t88;
t9 = mrSges(7,2) * t73 + mrSges(7,3) * t51;
t49 = -t46 * t10 + t21 * t70 - t22 * t68 - t44 * t9;
t24 = t98 * t45;
t20 = t53 * qJD(6);
t19 = t52 * qJD(6);
t18 = t54 * qJD(6);
t17 = t54 * t45;
t15 = (pkin(4) * t47 + t78) * qJD(4) + t59;
t14 = qJD(2) * t45 - t47 * t60;
t5 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t4 = t29 * t69 + (-Ifges(7,5) * t45 + t47 * t53) * qJD(4);
t3 = t28 * t69 + (-Ifges(7,6) * t45 + t47 * t52) * qJD(4);
t26 = [t47 * t65 + t3 * t88 + t4 * t90 + 0.2e1 * t15 * (-t45 * mrSges(6,2) - t47 * mrSges(6,3)) - 0.2e1 * t14 * t17 + 0.2e1 * t2 * t21 + 0.2e1 * t1 * t22 - 0.2e1 * t24 * t5 + 0.2e1 * t7 * t9 + 0.2e1 * t6 * t10 + (t12 * t46 + (-t11 - t93) * t44) * t69 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3) - (mrSges(6,1) + mrSges(5,3)) * t79) * qJD(2) + 0.2e1 * m(5) * (t75 + t82) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t75) + 0.2e1 * m(6) * (t15 * t23 + t82) + (t1 * t7 - t14 * t24 + t2 * t6) * t105 + ((0.2e1 * t43 * mrSges(5,1) - 0.2e1 * t23 * mrSges(6,2) + t87 + t92 + (-Ifges(6,6) - Ifges(5,4)) * t66) * t47 + (-0.2e1 * t43 * mrSges(5,2) + t23 * t104 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) - t94 - t95) * t45 + (-Ifges(5,1) + Ifges(5,2) - Ifges(6,2) + Ifges(6,3) - Ifges(7,3)) * t66) * t45) * qJD(4) + 0.2e1 * (t45 * mrSges(5,1) + t47 * mrSges(5,2) + mrSges(4,3)) * qJD(3); t44 * t10 - t46 * t9 + (t46 * t21 + t44 * t22) * qJD(6) + (-m(5) - m(4)) * qJD(3) + m(7) * (-t1 * t46 + t2 * t44 + (t44 * t7 + t46 * t6) * qJD(6)) - m(6) * t15 + (-t83 * t45 + t84 * t47) * qJD(4); 0; (m(7) * (t6 * t72 + t7 * t74 + t14) + t22 * t74 + t21 * t72 + t5) * t45 + (m(7) * (-qJD(4) * t24 + t6 * t70 - t7 * t68 - t56) - qJD(4) * t17 + t49) * t47 + (m(4) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t79) * qJD(2); 0; (0.1e1 - t80) * t45 * t71 * t105; m(7) * (qJ(5) * t14 - qJD(5) * t24 - t1 * t89 - t2 * t85) - t10 * t85 - t9 * t89 + t4 * t100 + t3 * t102 + t14 * t27 - qJD(5) * t17 - t24 * t18 + qJ(5) * t5 - t56 * mrSges(7,3) + (-t87 / 0.2e1 - t92 / 0.2e1 + t55 * mrSges(7,3) - (-m(7) * t55 - t44 * t21 + t46 * t22) * t99) * qJD(6) + ((-t95 / 0.2e1 - t94 / 0.2e1) * qJD(6) + t58 * qJD(2) + (-Ifges(5,6) + Ifges(6,5) + t86 / 0.2e1 + t91 / 0.2e1 - qJ(5) * mrSges(6,1) + t109 * t42) * qJD(4)) * t47 + (-t19 * t100 - t20 * t101 + (t29 * t100 + t28 * t102) * qJD(6) + (-Ifges(7,5) * t46 / 0.2e1 + Ifges(7,6) * t101 + Ifges(6,4) - Ifges(5,5) + pkin(4) * mrSges(6,1) - t58 * t42) * qJD(4) + (m(6) * t42 - mrSges(6,1)) * qJD(5) + t109 * qJD(2)) * t45; 0; t45 * t18 + ((t27 + t83) * t47 + (-t80 * mrSges(7,3) - t61 * t99 - t58) * t45) * qJD(4) + t107 * (qJ(5) * t71 + qJD(5) * t45); t18 * t108 + t19 * t44 - t20 * t46 + (-t86 - t91) * qJD(6) + (t107 * t108 + t104 + 0.2e1 * t27) * qJD(5); -mrSges(6,1) * t73 + m(7) * (-qJD(6) * t55 + t56) + (t42 * t73 - t76) * m(6) - t49; 0; (m(6) + t61) * t73; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + (-Ifges(7,6) * t70 - Ifges(7,3) * qJD(4)) * t45 + t65; t18; (-t44 * t73 + t46 * t67) * mrSges(7,2) + (t44 * t67 + t45 * t72) * mrSges(7,1); ((mrSges(7,2) * t99 - Ifges(7,6)) * t46 + (mrSges(7,1) * t99 - Ifges(7,5)) * t44) * qJD(6); -t27 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t26(1) t26(2) t26(4) t26(7) t26(11) t26(16); t26(2) t26(3) t26(5) t26(8) t26(12) t26(17); t26(4) t26(5) t26(6) t26(9) t26(13) t26(18); t26(7) t26(8) t26(9) t26(10) t26(14) t26(19); t26(11) t26(12) t26(13) t26(14) t26(15) t26(20); t26(16) t26(17) t26(18) t26(19) t26(20) t26(21);];
Mq  = res;
