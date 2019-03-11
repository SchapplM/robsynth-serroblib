% Calculate time derivative of joint inertia matrix for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:14
% EndTime: 2019-03-09 02:56:18
% DurationCPUTime: 1.82s
% Computational Cost: add. (1748->217), mult. (3492->321), div. (0->0), fcn. (3038->6), ass. (0->102)
t87 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t153 = m(6) + m(5);
t62 = cos(qJ(6));
t103 = qJD(6) * t62;
t106 = sin(pkin(9));
t107 = cos(pkin(9));
t129 = sin(qJ(3));
t130 = cos(qJ(3));
t40 = -t106 * t129 + t107 * t130;
t36 = t40 * qJD(3);
t39 = t106 * t130 + t107 * t129;
t61 = sin(qJ(6));
t70 = t39 * t103 + t61 * t36;
t104 = qJD(6) * t61;
t69 = t39 * t104 - t62 * t36;
t63 = -pkin(1) - pkin(7);
t142 = -qJ(4) + t63;
t45 = t142 * t130;
t32 = qJD(3) * t45 - qJD(4) * t129;
t44 = t142 * t129;
t64 = -qJD(3) * t44 - qJD(4) * t130;
t17 = t106 * t32 - t107 * t64;
t18 = t106 * t64 + t107 * t32;
t25 = t106 * t44 - t107 * t45;
t26 = t106 * t45 + t107 * t44;
t37 = t39 * qJD(3);
t152 = -t17 * t40 + t18 * t39 + t25 * t37 + t26 * t36;
t55 = t129 * pkin(3) + qJ(2);
t81 = -qJ(5) * t40 + t55;
t21 = pkin(4) * t39 + t81;
t151 = 0.2e1 * t21;
t150 = -0.2e1 * t55;
t132 = pkin(4) + pkin(8);
t14 = t132 * t39 + t81;
t19 = t40 * pkin(5) + t25;
t3 = -t14 * t61 + t19 * t62;
t91 = qJD(3) * t130;
t49 = pkin(3) * t91 + qJD(2);
t68 = qJ(5) * t37 - qJD(5) * t40 + t49;
t8 = t132 * t36 + t68;
t9 = -t37 * pkin(5) + t17;
t1 = qJD(6) * t3 + t61 * t9 + t62 * t8;
t4 = t14 * t62 + t19 * t61;
t2 = -qJD(6) * t4 - t61 * t8 + t62 * t9;
t137 = -t1 * t61 - t2 * t62 + (t3 * t61 - t4 * t62) * qJD(6);
t149 = m(7) * t137;
t108 = t61 ^ 2 + t62 ^ 2;
t148 = m(7) * t108;
t147 = mrSges(6,2) - mrSges(5,1);
t46 = mrSges(7,1) * t61 + mrSges(7,2) * t62;
t146 = -mrSges(6,3) - t46;
t145 = -Ifges(4,1) + Ifges(4,2);
t90 = qJD(3) * t129;
t144 = -mrSges(4,1) * t90 - mrSges(4,2) * t91;
t143 = -t106 * t36 + t107 * t37;
t54 = -pkin(3) * t107 - pkin(4);
t52 = pkin(3) * t106 + qJ(5);
t71 = qJD(5) * t39 + t36 * t52;
t140 = t37 * t54 + t71;
t12 = mrSges(7,2) * t37 - mrSges(7,3) * t69;
t13 = -mrSges(7,1) * t37 - mrSges(7,3) * t70;
t75 = -t61 * t12 - t62 * t13;
t126 = Ifges(7,4) * t62;
t47 = -Ifges(7,2) * t61 + t126;
t127 = Ifges(7,4) * t61;
t48 = Ifges(7,1) * t62 - t127;
t139 = t62 * t47 + t61 * t48;
t116 = t61 * mrSges(7,3);
t23 = mrSges(7,1) * t40 - t116 * t39;
t112 = t62 * mrSges(7,3);
t24 = -mrSges(7,2) * t40 + t112 * t39;
t138 = qJD(6) * (t61 * t23 - t62 * t24) + t75;
t135 = 0.2e1 * qJD(2);
t134 = m(5) * pkin(3);
t27 = t39 * t36;
t119 = t39 * t61;
t118 = t39 * t62;
t117 = t40 * t37;
t95 = t118 / 0.2e1;
t83 = t3 * t62 + t4 * t61;
t78 = Ifges(7,1) * t61 + t126;
t77 = Ifges(7,2) * t62 + t127;
t76 = Ifges(7,5) * t61 + Ifges(7,6) * t62;
t74 = t17 * t25 + t18 * t26;
t73 = t62 * t23 + t61 * t24;
t67 = Ifges(7,5) * t70 - Ifges(7,6) * t69 - Ifges(7,3) * t37;
t51 = -pkin(8) + t54;
t43 = t78 * qJD(6);
t42 = t77 * qJD(6);
t41 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t34 = t37 * mrSges(5,2);
t33 = t37 * mrSges(6,3);
t22 = (-mrSges(7,1) * t62 + mrSges(7,2) * t61) * t39;
t20 = -t39 * pkin(5) + t26;
t16 = t40 * Ifges(7,5) + t39 * t78;
t15 = t40 * Ifges(7,6) + t39 * t77;
t11 = pkin(4) * t36 + t68;
t10 = -t36 * pkin(5) + t18;
t7 = mrSges(7,1) * t69 + mrSges(7,2) * t70;
t6 = Ifges(7,1) * t70 - Ifges(7,4) * t69 - t37 * Ifges(7,5);
t5 = Ifges(7,4) * t70 - Ifges(7,2) * t69 - t37 * Ifges(7,6);
t28 = [((m(4) + m(3)) * t135 + 0.2e1 * (mrSges(4,1) * t130 - mrSges(4,2) * t129) * qJD(3)) * qJ(2) + 0.2e1 * t49 * (mrSges(5,1) * t39 + mrSges(5,2) * t40) + 0.2e1 * t11 * (-mrSges(6,2) * t39 - mrSges(6,3) * t40) + (-0.2e1 * Ifges(4,4) * t130 + t129 * t145) * t91 + (0.2e1 * Ifges(4,4) * t129 + t130 * t145) * t90 + 0.2e1 * t20 * t7 + 0.2e1 * t10 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t1 * t24 - t69 * t15 + t70 * t16 + 0.2e1 * t4 * t12 + 0.2e1 * t3 * t13 + 0.2e1 * m(7) * (t1 * t4 + t10 * t20 + t2 * t3) + 0.2e1 * m(5) * (t49 * t55 + t74) + 0.2e1 * m(6) * (t11 * t21 + t74) + t40 * t67 - t152 * t87 - (mrSges(5,1) * t150 + mrSges(6,2) * t151 + 0.2e1 * (Ifges(5,4) + Ifges(6,6)) * t40 + 0.2e1 * (-Ifges(5,2) - Ifges(6,3)) * t39) * t36 + (mrSges(4,1) * t129 + mrSges(4,2) * t130 + mrSges(3,3)) * t135 + t34 * t150 + t33 * t151 + ((-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t40 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) - t76) * t39) * t37 + t5 * t118 + t6 * t119; t39 * t7 + t73 * t37 - (t39 * t87 - t22) * t36 + (t37 * t87 + t138) * t40 + m(7) * (t10 * t39 + t137 * t40 + t20 * t36 + t83 * t37) + t153 * t152; 0.2e1 * m(7) * (-t108 * t117 + t27) + 0.2e1 * t153 * (t27 - t117); (-Ifges(7,5) * t62 / 0.2e1 + Ifges(7,6) * t61 / 0.2e1 - Ifges(5,5) + Ifges(6,4)) * t37 + (m(6) * t52 + t106 * t134 - mrSges(5,2) + mrSges(6,3)) * t18 + (m(7) * t52 + t46) * t10 + (t48 * t95 - t40 * t76 / 0.2e1) * qJD(6) + (m(6) * t26 + m(7) * t20 + t22) * qJD(5) + (-t103 * t4 + t104 * t3) * mrSges(7,3) - t61 * t5 / 0.2e1 + t62 * t6 / 0.2e1 + t52 * t7 + (t103 * t24 - t104 * t23 - t149 - t75) * t51 - t20 * t41 + (m(6) * t54 - t107 * t134 + t147) * t17 + t144 * t63 - t140 * mrSges(6,1) - t43 * t119 / 0.2e1 - t1 * t116 - t2 * t112 - t15 * t103 / 0.2e1 - Ifges(4,5) * t90 - Ifges(4,6) * t91 - t42 * t95 - (t39 * t47 + t16) * t104 / 0.2e1 + t143 * mrSges(5,3) * pkin(3) - (Ifges(5,6) - Ifges(6,5) - t139 / 0.2e1) * t36; m(6) * t140 + m(7) * t71 - t39 * t41 - t143 * t134 - (mrSges(5,2) + t146) * t36 + t144 + (-mrSges(7,3) * t108 + t148 * t51 + t147) * t37; -0.2e1 * t41 * t52 + t42 * t61 - t43 * t62 - t139 * qJD(6) + 0.2e1 * ((m(6) + m(7)) * t52 - t146) * qJD(5); t62 * t12 - t61 * t13 + t33 - t34 - t147 * t36 - t73 * qJD(6) + m(7) * (-qJD(6) * t83 + t1 * t62 - t2 * t61) + m(6) * t11 + m(5) * t49; 0; 0; 0; m(6) * t17 - t37 * mrSges(6,1) - t138 - t149; 0.2e1 * (m(6) / 0.2e1 + t148 / 0.2e1) * t37; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t67; (t103 * t40 - t37 * t61) * mrSges(7,2) + (t104 * t40 + t37 * t62) * mrSges(7,1); (-t46 * t51 - t76) * qJD(6); t41; -t46 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t28(1) t28(2) t28(4) t28(7) t28(11) t28(16); t28(2) t28(3) t28(5) t28(8) t28(12) t28(17); t28(4) t28(5) t28(6) t28(9) t28(13) t28(18); t28(7) t28(8) t28(9) t28(10) t28(14) t28(19); t28(11) t28(12) t28(13) t28(14) t28(15) t28(20); t28(16) t28(17) t28(18) t28(19) t28(20) t28(21);];
Mq  = res;
