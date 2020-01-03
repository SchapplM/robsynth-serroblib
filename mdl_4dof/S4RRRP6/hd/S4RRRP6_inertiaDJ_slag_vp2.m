% Calculate time derivative of joint inertia matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:13
% DurationCPUTime: 1.17s
% Computational Cost: add. (564->205), mult. (1496->306), div. (0->0), fcn. (964->4), ass. (0->98)
t99 = Ifges(4,5) + Ifges(5,5);
t122 = -Ifges(4,3) - Ifges(5,3);
t66 = sin(qJ(3));
t68 = cos(qJ(3));
t69 = cos(qJ(2));
t90 = qJD(2) * t69;
t82 = t68 * t90;
t67 = sin(qJ(2));
t88 = qJD(3) * t67;
t71 = -t66 * t88 + t82;
t85 = t66 * t90;
t87 = qJD(3) * t68;
t70 = t67 * t87 + t85;
t98 = Ifges(4,6) + Ifges(5,6);
t121 = t99 * t66 + t98 * t68;
t103 = Ifges(5,4) * t68;
t73 = -Ifges(5,2) * t66 + t103;
t105 = Ifges(4,4) * t68;
t74 = -Ifges(4,2) * t66 + t105;
t118 = (t73 + t74) * qJD(3);
t104 = Ifges(5,4) * t66;
t75 = Ifges(5,1) * t68 - t104;
t106 = Ifges(4,4) * t66;
t76 = Ifges(4,1) * t68 - t106;
t117 = (t75 + t76) * qJD(3);
t44 = -pkin(2) * t69 - t67 * pkin(6) - pkin(1);
t107 = pkin(5) * t69;
t58 = t68 * t107;
t19 = t66 * t44 + t58;
t116 = 2 * m(4);
t115 = 2 * m(5);
t114 = -0.2e1 * pkin(1);
t113 = 0.2e1 * pkin(5);
t112 = -2 * mrSges(5,3);
t111 = m(5) * pkin(3);
t108 = pkin(5) * t66;
t102 = t66 * t67;
t101 = t67 * t68;
t97 = -qJ(4) - pkin(6);
t22 = -Ifges(5,5) * t69 + t75 * t67;
t23 = -Ifges(4,5) * t69 + t76 * t67;
t96 = t22 + t23;
t42 = (pkin(2) * t67 - pkin(6) * t69) * qJD(2);
t95 = t66 * t42 + t44 * t87;
t91 = qJD(2) * t67;
t94 = t91 * t108 + t68 * t42;
t89 = qJD(3) * t66;
t32 = mrSges(5,1) * t89 + mrSges(5,2) * t87;
t93 = qJ(4) * t67;
t92 = qJ(4) * t68;
t86 = qJD(4) * t68;
t48 = Ifges(5,2) * t68 + t104;
t49 = Ifges(4,2) * t68 + t106;
t81 = -t48 / 0.2e1 - t49 / 0.2e1;
t50 = Ifges(5,1) * t66 + t103;
t51 = Ifges(4,1) * t66 + t105;
t80 = -t50 / 0.2e1 - t51 / 0.2e1;
t79 = qJD(3) * t97;
t78 = t122 * t91 - t99 * t82;
t77 = mrSges(4,1) * t66 + mrSges(4,2) * t68;
t20 = -t69 * Ifges(5,6) + t73 * t67;
t21 = -t69 * Ifges(4,6) + t74 * t67;
t72 = t98 * t69 - t20 - t21;
t9 = t70 * mrSges(5,1) + t71 * mrSges(5,2);
t65 = Ifges(4,5) * t87;
t64 = Ifges(5,5) * t87;
t59 = -pkin(3) * t68 - pkin(2);
t47 = t97 * t68;
t46 = -mrSges(5,1) * t68 + mrSges(5,2) * t66;
t45 = t97 * t66;
t43 = (pkin(3) * t66 + pkin(5)) * t67;
t41 = -mrSges(4,1) * t69 - mrSges(4,3) * t101;
t40 = -mrSges(5,1) * t69 - mrSges(5,3) * t101;
t39 = mrSges(4,2) * t69 - mrSges(4,3) * t102;
t38 = mrSges(5,2) * t69 - mrSges(5,3) * t102;
t33 = t77 * qJD(3);
t31 = t68 * t44;
t28 = (mrSges(5,1) * t66 + mrSges(5,2) * t68) * t67;
t25 = -qJD(4) * t66 + t68 * t79;
t24 = t66 * t79 + t86;
t18 = -t66 * t107 + t31;
t17 = t70 * pkin(3) + pkin(5) * t90;
t16 = -mrSges(4,2) * t91 - t70 * mrSges(4,3);
t15 = -mrSges(5,2) * t91 - t70 * mrSges(5,3);
t14 = mrSges(4,1) * t91 - t71 * mrSges(4,3);
t13 = mrSges(5,1) * t91 - t71 * mrSges(5,3);
t12 = -t66 * t93 + t19;
t11 = -t67 * t92 + t31 + (-pkin(3) - t108) * t69;
t10 = t70 * mrSges(4,1) + t71 * mrSges(4,2);
t8 = -t51 * t88 + (Ifges(4,5) * t67 + t76 * t69) * qJD(2);
t7 = -t50 * t88 + (Ifges(5,5) * t67 + t75 * t69) * qJD(2);
t6 = -t49 * t88 + (Ifges(4,6) * t67 + t74 * t69) * qJD(2);
t5 = -t48 * t88 + (Ifges(5,6) * t67 + t73 * t69) * qJD(2);
t4 = -qJD(3) * t19 + t94;
t3 = (-t68 * t91 - t69 * t89) * pkin(5) + t95;
t2 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t101 + (-qJD(4) * t67 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t69) * t66 + t95;
t1 = -t67 * t86 + (pkin(3) * t67 - t69 * t92) * qJD(2) + (-t58 + (-t44 + t93) * t66) * qJD(3) + t94;
t26 = [0.2e1 * t1 * t40 + 0.2e1 * t11 * t13 + 0.2e1 * t12 * t15 + 0.2e1 * t18 * t14 + 0.2e1 * t19 * t16 + 0.2e1 * t17 * t28 + 0.2e1 * t2 * t38 + 0.2e1 * t3 * t39 + 0.2e1 * t4 * t41 + 0.2e1 * t43 * t9 + (t18 * t4 + t19 * t3) * t116 + (t1 * t11 + t12 * t2 + t17 * t43) * t115 + ((mrSges(3,2) * t114 + 0.2e1 * Ifges(3,4) * t69 + t72 * t66 + t96 * t68) * qJD(2) + t78) * t69 + (t10 * t113 + (t7 + t8) * t68 + (-t5 - t6) * t66 + (t72 * t68 + (t99 * t69 - t96) * t66) * qJD(3) + (mrSges(3,1) * t114 + (-t98 * t66 + t99 * t68 - 0.2e1 * Ifges(3,4)) * t67 + (pkin(5) ^ 2 * t116 + t77 * t113 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + t122) * t69) * qJD(2)) * t67; t45 * t13 + t17 * t46 - t47 * t15 + t59 * t9 + t24 * t38 + t25 * t40 + t43 * t32 - pkin(2) * t10 + m(5) * (t1 * t45 + t11 * t25 + t12 * t24 + t17 * t59 - t2 * t47) + (-t65 / 0.2e1 - t64 / 0.2e1 + (Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1)) * pkin(5)) * qJD(2)) * t69 + (t7 / 0.2e1 + t8 / 0.2e1 - t4 * mrSges(4,3) - t1 * mrSges(5,3) + (pkin(5) * mrSges(4,2) + t81) * t90 + (-t12 * mrSges(5,3) - t19 * mrSges(4,3) + pkin(3) * t28 - t20 / 0.2e1 - t21 / 0.2e1 + (Ifges(4,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t69 + t43 * t111) * qJD(3) + (-m(4) * t4 - t14 + (-m(4) * t19 - t39) * qJD(3)) * pkin(6)) * t66 + (t2 * mrSges(5,3) + t3 * mrSges(4,3) + t5 / 0.2e1 + t6 / 0.2e1 + (-pkin(5) * mrSges(4,1) - t80) * t90 + (t22 / 0.2e1 + t23 / 0.2e1 - t18 * mrSges(4,3) - t11 * mrSges(5,3)) * qJD(3) + (t16 + m(4) * (-t18 * qJD(3) + t3) - qJD(3) * t41) * pkin(6)) * t68 + (pkin(5) * t33 + (t80 * t66 + t81 * t68) * qJD(3) - t118 * t66 / 0.2e1 + t117 * t68 / 0.2e1 + (-Ifges(3,6) + pkin(5) * mrSges(3,2) + t121 / 0.2e1) * qJD(2)) * t67; 0.2e1 * t59 * t32 + (-t24 * t47 + t25 * t45) * t115 - 0.2e1 * pkin(2) * t33 + (t25 * t112 + (-t47 * t112 - t48 - t49 + 0.2e1 * (m(5) * t59 + t46) * pkin(3)) * qJD(3) + t117) * t66 + (0.2e1 * t24 * mrSges(5,3) + (t45 * t112 + t50 + t51) * qJD(3) + t118) * t68; t4 * mrSges(4,1) + t1 * mrSges(5,1) - t3 * mrSges(4,2) - t2 * mrSges(5,2) - t98 * t85 + (m(5) * t1 + t13) * pkin(3) - t121 * t88 - t78; -mrSges(5,2) * t24 + t64 + t65 + (mrSges(5,1) + t111) * t25 + ((-mrSges(4,1) * pkin(6) - (mrSges(5,3) * pkin(3))) * t68 + (mrSges(4,2) * pkin(6) - t98) * t66) * qJD(3); 0; m(5) * t17 + t9; t89 * t111 + t32; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t26(1), t26(2), t26(4), t26(7); t26(2), t26(3), t26(5), t26(8); t26(4), t26(5), t26(6), t26(9); t26(7), t26(8), t26(9), t26(10);];
Mq = res;
