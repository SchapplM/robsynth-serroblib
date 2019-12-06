% Calculate time derivative of joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:47
% EndTime: 2019-12-05 18:44:50
% DurationCPUTime: 1.21s
% Computational Cost: add. (2584->180), mult. (5841->282), div. (0->0), fcn. (5214->6), ass. (0->81)
t100 = -pkin(7) - pkin(6);
t74 = sin(qJ(2));
t65 = t100 * t74;
t77 = cos(qJ(2));
t66 = t100 * t77;
t73 = sin(qJ(3));
t76 = cos(qJ(3));
t47 = t76 * t65 + t66 * t73;
t108 = 2 * mrSges(5,3);
t107 = 2 * mrSges(6,3);
t72 = sin(qJ(4));
t106 = pkin(3) * t72;
t75 = cos(qJ(4));
t99 = pkin(3) * t75;
t95 = -mrSges(5,1) - mrSges(6,1);
t59 = t73 * t77 + t76 * t74;
t34 = -t59 * pkin(8) + t47;
t48 = t73 * t65 - t76 * t66;
t58 = -t73 * t74 + t76 * t77;
t35 = t58 * pkin(8) + t48;
t23 = t72 * t34 + t75 * t35;
t105 = qJD(2) + qJD(3);
t104 = 2 * m(5);
t103 = 2 * m(6);
t70 = -pkin(2) * t77 - pkin(1);
t102 = 0.2e1 * t70;
t101 = m(6) * pkin(4);
t97 = t72 * t73;
t96 = t73 * t75;
t94 = -mrSges(5,2) - mrSges(6,2);
t69 = pkin(2) * t76 + pkin(3);
t90 = qJD(4) * t75;
t91 = qJD(4) * t72;
t39 = t69 * t90 + (-t73 * t91 + (t75 * t76 - t97) * qJD(3)) * pkin(2);
t53 = pkin(2) * t96 + t69 * t72;
t93 = pkin(3) * t53 * t90 + t39 * t106;
t92 = pkin(3) * qJD(4);
t89 = 0.2e1 * t77;
t46 = t105 * t59;
t36 = qJD(2) * t74 * pkin(2) + t46 * pkin(3);
t88 = qJD(2) * t100;
t87 = t94 * t75;
t22 = t75 * t34 - t35 * t72;
t86 = 0.2e1 * t94;
t85 = 2 * Ifges(5,4) + 2 * Ifges(6,4);
t52 = -pkin(2) * t97 + t75 * t69;
t41 = t58 * t75 - t59 * t72;
t45 = t105 * t58;
t15 = qJD(4) * t41 + t45 * t75 - t46 * t72;
t42 = t58 * t72 + t59 * t75;
t63 = t74 * t88;
t64 = t77 * t88;
t25 = t47 * qJD(3) + t76 * t63 + t73 * t64;
t20 = -t46 * pkin(8) + t25;
t26 = -t48 * qJD(3) - t63 * t73 + t76 * t64;
t21 = -t45 * pkin(8) + t26;
t6 = -t23 * qJD(4) - t20 * t72 + t75 * t21;
t3 = -t15 * qJ(5) - t42 * qJD(5) + t6;
t84 = m(6) * t3 - t15 * mrSges(6,3);
t50 = -t58 * pkin(3) + t70;
t40 = -t69 * t91 + (-t73 * t90 + (-t72 * t76 - t96) * qJD(3)) * pkin(2);
t37 = t40 * mrSges(6,1);
t38 = t40 * mrSges(5,1);
t83 = t39 * t94 + t37 + t38;
t5 = t75 * t20 + t72 * t21 + t34 * t90 - t35 * t91;
t16 = -qJD(4) * t42 - t45 * t72 - t46 * t75;
t82 = t16 * t53 + t39 * t41 - t40 * t42;
t81 = (-mrSges(4,1) * t73 - mrSges(4,2) * t76) * qJD(3) * pkin(2);
t80 = t72 * t16 + (t41 * t75 + t42 * t72) * qJD(4);
t2 = t16 * qJ(5) + t41 * qJD(5) + t5;
t79 = t6 * mrSges(5,1) + t3 * mrSges(6,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(5,6) + Ifges(6,6)) * t16 + (Ifges(5,5) + Ifges(6,5)) * t15;
t78 = t26 * mrSges(4,1) - t25 * mrSges(4,2) + Ifges(4,5) * t45 - Ifges(4,6) * t46 + t79;
t68 = pkin(4) + t99;
t51 = pkin(4) + t52;
t28 = t53 * t39;
t27 = -t41 * pkin(4) + t50;
t10 = t15 * mrSges(6,2);
t9 = qJ(5) * t41 + t23;
t8 = -t42 * qJ(5) + t22;
t7 = -t16 * pkin(4) + t36;
t1 = [0.2e1 * t59 * t45 * Ifges(4,1) - 0.2e1 * t46 * Ifges(4,2) * t58 + (mrSges(4,1) * t46 + mrSges(4,2) * t45) * t102 + 0.2e1 * t36 * (-mrSges(5,1) * t41 + mrSges(5,2) * t42) + 0.2e1 * t7 * (-mrSges(6,1) * t41 + mrSges(6,2) * t42) + 0.2e1 * t27 * t10 + (t22 * t6 + t23 * t5 + t36 * t50) * t104 + (t2 * t9 + t27 * t7 + t3 * t8) * t103 + 0.2e1 * m(4) * (t48 * t25 + t47 * t26) + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t77) * t89 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t102 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t58 + mrSges(4,2) * t59) - 0.2e1 * Ifges(3,4) * t74 + (Ifges(3,1) - Ifges(3,2)) * t89) * t74) * qJD(2) + (t2 * t41 - t3 * t42) * t107 + (t5 * t41 - t6 * t42) * t108 + 0.2e1 * (t45 * t58 - t46 * t59) * Ifges(4,4) + 0.2e1 * (t25 * t58 - t26 * t59 - t45 * t47 - t46 * t48) * mrSges(4,3) + (-0.2e1 * mrSges(5,1) * t50 - 0.2e1 * mrSges(6,1) * t27 + t9 * t107 + t23 * t108 + t42 * t85 + 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t41) * t16 + (0.2e1 * mrSges(5,2) * t50 - 0.2e1 * mrSges(5,3) * t22 - 0.2e1 * mrSges(6,3) * t8 + t41 * t85 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t42) * t15; t78 + (Ifges(3,5) * t77 - Ifges(3,6) * t74 + (-mrSges(3,1) * t77 + mrSges(3,2) * t74) * pkin(6)) * qJD(2) + m(5) * (t22 * t40 + t23 * t39 + t5 * t53 + t52 * t6) + (m(4) * (t25 * t73 + t26 * t76 + (-t47 * t73 + t48 * t76) * qJD(3)) + (-t76 * t45 - t73 * t46 + (t58 * t76 + t59 * t73) * qJD(3)) * mrSges(4,3)) * pkin(2) + m(6) * (t2 * t53 + t3 * t51 + t39 * t9 + t40 * t8) + (-t15 * t52 + t82) * mrSges(5,3) + (-t15 * t51 + t82) * mrSges(6,3); 0.2e1 * t37 + 0.2e1 * t38 + t39 * t86 + 0.2e1 * t81 + (t40 * t51 + t28) * t103 + (t40 * t52 + t28) * t104; t78 + t84 * t68 + (m(6) * (t2 * t72 - t8 * t91 + t9 * t90) + m(5) * (-t22 * t91 + t23 * t90 + t5 * t72 + t6 * t75) + t80 * mrSges(6,3) + (-t75 * t15 + t80) * mrSges(5,3)) * pkin(3); t81 + (t72 * t95 + t87) * t92 + m(6) * (-pkin(3) * t51 * t91 + t40 * t68 + t93) + m(5) * ((t40 * t75 - t52 * t91) * pkin(3) + t93) + t83; (t86 * t99 + 0.2e1 * ((-t68 + t99) * m(6) + t95) * t106) * qJD(4); pkin(4) * t84 + t79; t101 * t40 + t83; (t87 + (t95 - t101) * t72) * t92; 0; m(6) * t7 - t16 * mrSges(6,1) + t10; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
