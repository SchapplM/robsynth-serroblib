% Calculate time derivative of joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:30
% DurationCPUTime: 0.76s
% Computational Cost: add. (1084->147), mult. (2111->204), div. (0->0), fcn. (1562->6), ass. (0->85)
t61 = sin(qJ(5));
t62 = sin(qJ(4));
t64 = cos(qJ(5));
t65 = cos(qJ(4));
t37 = -t61 * t65 - t64 * t62;
t94 = t64 * t65;
t38 = -t61 * t62 + t94;
t117 = t37 * t64 + t38 * t61;
t113 = qJD(4) + qJD(5);
t87 = qJD(5) * t61;
t89 = qJD(4) * t62;
t21 = t113 * t94 - t61 * t89 - t62 * t87;
t101 = t21 * t37;
t22 = t113 * t37;
t99 = t22 * t38;
t116 = t99 - t101;
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t91 = pkin(1) * qJD(2);
t92 = t62 ^ 2 + t65 ^ 2;
t115 = (-mrSges(3,2) * t66 + (-t92 * mrSges(5,3) - mrSges(3,1) + mrSges(4,2)) * t63) * t91;
t44 = mrSges(5,1) * t62 + mrSges(5,2) * t65;
t114 = t44 + mrSges(4,3);
t100 = t21 * t61;
t90 = pkin(4) * qJD(5);
t112 = (-pkin(4) * t100 + t117 * t90) * mrSges(6,3);
t110 = 2 * m(6);
t9 = mrSges(6,1) * t21 + mrSges(6,2) * t22;
t109 = 0.2e1 * t9;
t23 = -mrSges(6,1) * t37 + mrSges(6,2) * t38;
t108 = 0.2e1 * t23;
t75 = mrSges(5,1) * t65 - mrSges(5,2) * t62;
t39 = t75 * qJD(4);
t107 = 0.2e1 * t39;
t78 = -pkin(1) * t66 - pkin(2);
t51 = -pkin(7) + t78;
t106 = -pkin(8) + t51;
t67 = -pkin(2) - pkin(7);
t105 = -pkin(8) + t67;
t103 = Ifges(5,4) * t62;
t102 = Ifges(5,4) * t65;
t98 = t22 * t64;
t47 = t66 * t91 + qJD(3);
t52 = pkin(1) * t63 + qJ(3);
t95 = t52 * t47;
t93 = Ifges(6,5) * t22 - Ifges(6,6) * t21;
t88 = qJD(4) * t65;
t86 = qJD(5) * t64;
t84 = 0.2e1 * mrSges(6,3);
t82 = mrSges(6,3) * t98;
t81 = -0.2e1 * t116 * mrSges(6,3);
t79 = t63 * t91;
t32 = t106 * t65;
t42 = t105 * t65;
t76 = t22 * mrSges(6,1) - t21 * mrSges(6,2);
t31 = t106 * t62;
t16 = t31 * t64 + t32 * t61;
t15 = -t31 * t61 + t32 * t64;
t41 = t105 * t62;
t25 = t41 * t64 + t42 * t61;
t24 = -t41 * t61 + t42 * t64;
t56 = pkin(8) * t89;
t26 = -t51 * t89 + t65 * t79 + t56;
t27 = qJD(4) * t32 + t62 * t79;
t2 = qJD(5) * t15 + t26 * t61 + t27 * t64;
t3 = -qJD(5) * t16 + t26 * t64 - t27 * t61;
t74 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t93;
t34 = -t67 * t89 + t56;
t35 = qJD(4) * t42;
t7 = qJD(5) * t24 + t34 * t61 + t35 * t64;
t8 = -qJD(5) * t25 + t34 * t64 - t35 * t61;
t73 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + t93;
t72 = t92 * t79;
t71 = qJ(3) * t47 + qJD(3) * t52;
t70 = t15 * t22 + t16 * t21 - t2 * t37 + t3 * t38;
t69 = t21 * t25 + t22 * t24 - t37 * t7 + t38 * t8;
t68 = 0.2e1 * Ifges(6,1) * t99 + (-Ifges(5,1) * t62 - t102) * t88 - 0.2e1 * Ifges(6,2) * t101 + (-t65 * (-Ifges(5,2) * t62 + t102) - t62 * (Ifges(5,1) * t65 - t103)) * qJD(4) - (-Ifges(5,2) * t65 - t103) * t89 + 0.2e1 * (-t21 * t38 + t22 * t37) * Ifges(6,4);
t58 = t62 * pkin(4);
t57 = pkin(4) * t88;
t53 = qJ(3) + t58;
t48 = qJD(3) + t57;
t43 = t52 + t58;
t36 = t47 + t57;
t33 = (-mrSges(6,1) * t61 - mrSges(6,2) * t64) * t90;
t1 = [t36 * t108 + t52 * t107 + t43 * t109 + 0.2e1 * t114 * t47 + 0.2e1 * t115 + 0.2e1 * m(5) * (t51 * t72 + t95) + 0.2e1 * m(4) * (t78 * t79 + t95) + (t15 * t3 + t16 * t2 + t36 * t43) * t110 - t70 * t84 + t68; m(4) * (-pkin(2) * t79 + t71) + m(6) * (t15 * t8 + t16 * t7 + t2 * t25 + t24 * t3 + t36 * t53 + t43 * t48) + (t53 + t43) * t9 + (t52 + qJ(3)) * t39 + (t48 + t36) * t23 + t68 + m(5) * (t67 * t72 + t71) + t115 + ((-t3 - t8) * t38 + (t2 + t7) * t37 + (-t15 - t24) * t22 + (-t16 - t25) * t21) * mrSges(6,3) + t114 * (t47 + qJD(3)); t48 * t108 + t53 * t109 + qJ(3) * t107 + (t24 * t8 + t25 * t7 + t48 * t53) * t110 + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t114) * qJD(3) - t69 * t84 + t68; m(6) * t70 + (m(5) * t92 + m(4)) * t79 + t81; m(6) * t69 + t81; t116 * t110; t75 * t79 + (-Ifges(5,5) * t62 - Ifges(5,6) * t65 - t44 * t51) * qJD(4) + (-t82 + m(6) * (-t15 * t87 + t16 * t86 + t2 * t61 + t3 * t64)) * pkin(4) + t74 + t112; ((-mrSges(5,2) * t67 - Ifges(5,6)) * t65 + (-mrSges(5,1) * t67 - Ifges(5,5)) * t62) * qJD(4) + (-t82 + m(6) * (-t24 * t87 + t25 * t86 + t61 * t7 + t64 * t8)) * pkin(4) + t73 + t112; -t44 * qJD(4) + m(6) * (-t117 * qJD(5) + t100 + t98) * pkin(4) + t76; 0.2e1 * t33; t74; t73; t76; t33; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
