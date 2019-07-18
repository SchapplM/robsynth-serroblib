% Calculate time derivative of joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:46
% EndTime: 2019-07-18 17:20:50
% DurationCPUTime: 1.06s
% Computational Cost: add. (1346->179), mult. (3077->276), div. (0->0), fcn. (2494->6), ass. (0->96)
t78 = cos(qJ(2));
t134 = 0.2e1 * t78;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t56 = -mrSges(6,1) * t76 + mrSges(6,2) * t73;
t141 = t56 - mrSges(5,1);
t108 = t73 ^ 2 + t76 ^ 2;
t140 = t108 * mrSges(6,3) - mrSges(5,2);
t112 = pkin(3) + qJ(3);
t75 = sin(qJ(2));
t57 = t112 * t75;
t58 = t112 * t78;
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t87 = -t77 * t57 - t58 * t74;
t139 = t87 * qJD(4);
t103 = qJD(5) * t76;
t135 = qJD(2) + qJD(4);
t50 = t74 * t75 - t77 * t78;
t34 = t135 * t50;
t51 = t74 * t78 + t75 * t77;
t85 = t51 * t103 - t34 * t73;
t93 = qJD(2) * t112;
t45 = t78 * qJD(3) - t75 * t93;
t82 = -t75 * qJD(3) - t78 * t93;
t13 = t77 * t45 + t74 * t82 + t139;
t42 = -t57 * t74 + t58 * t77;
t79 = pkin(2) + pkin(1);
t61 = t79 * t78;
t43 = -pkin(4) * t51 - t61;
t15 = -t42 * t73 + t43 * t76;
t107 = qJD(2) * t75;
t68 = pkin(1) * t107;
t55 = pkin(2) * t107 + t68;
t23 = pkin(4) * t34 + t55;
t2 = t15 * qJD(5) + t13 * t76 + t23 * t73;
t16 = t42 * t76 + t43 * t73;
t3 = -t16 * qJD(5) - t13 * t73 + t23 * t76;
t138 = t2 * t76 - t3 * t73;
t137 = m(4) * qJ(3) + mrSges(4,3);
t136 = -t15 * t73 + t16 * t76;
t133 = -2 * mrSges(5,3);
t14 = t42 * qJD(4) + t45 * t74 - t77 * t82;
t132 = 0.2e1 * t14;
t131 = -0.2e1 * t87;
t130 = 0.2e1 * t55;
t125 = Ifges(6,4) * t73;
t124 = Ifges(6,4) * t76;
t123 = Ifges(6,6) * t73;
t122 = t14 * t87;
t118 = t51 * t73;
t117 = t51 * t76;
t116 = t74 * t79;
t115 = t76 * t34;
t113 = Ifges(3,4) + Ifges(4,4);
t35 = t135 * t51;
t111 = -Ifges(6,5) * t115 + Ifges(6,3) * t35;
t106 = qJD(2) * t78;
t110 = mrSges(4,1) * t107 + mrSges(4,2) * t106;
t105 = qJD(4) * t77;
t104 = qJD(5) * t73;
t102 = t51 * t104;
t100 = -m(4) * pkin(1) - mrSges(4,1);
t99 = m(6) * t108;
t98 = t141 * t74;
t96 = -t104 / 0.2e1;
t95 = t35 * mrSges(5,1) - t34 * mrSges(5,2);
t94 = -(2 * Ifges(5,4)) - t123;
t92 = mrSges(6,1) * t73 + mrSges(6,2) * t76;
t91 = Ifges(6,1) * t76 - t125;
t90 = -Ifges(6,2) * t73 + t124;
t89 = Ifges(6,5) * t73 + Ifges(6,6) * t76;
t28 = -mrSges(6,2) * t50 - mrSges(6,3) * t118;
t29 = mrSges(6,1) * t50 - mrSges(6,3) * t117;
t88 = t76 * t28 - t73 * t29;
t86 = -t139 * t74 - t14 * t77;
t84 = t102 + t115;
t53 = t90 * qJD(5);
t54 = t91 * qJD(5);
t59 = Ifges(6,2) * t76 + t125;
t60 = Ifges(6,1) * t73 + t124;
t83 = t60 * t103 - t59 * t104 + t76 * t53 + t73 * t54;
t10 = -mrSges(6,2) * t35 - t85 * mrSges(6,3);
t9 = mrSges(6,1) * t35 + t84 * mrSges(6,3);
t81 = -t29 * t103 - t28 * t104 + m(6) * (-t15 * t103 - t16 * t104 + t138) + t76 * t10 - t73 * t9;
t19 = t50 * Ifges(6,6) + t90 * t51;
t20 = t50 * Ifges(6,5) + t91 * t51;
t52 = t92 * qJD(5);
t6 = -t84 * Ifges(6,4) - t85 * Ifges(6,2) + t35 * Ifges(6,6);
t66 = Ifges(6,5) * t103;
t7 = -t84 * Ifges(6,1) - t85 * Ifges(6,4) + t35 * Ifges(6,5);
t80 = t19 * t96 + t20 * t103 / 0.2e1 - t87 * t52 - Ifges(5,5) * t34 - t53 * t118 / 0.2e1 + t54 * t117 / 0.2e1 + t73 * t7 / 0.2e1 + t50 * (-Ifges(6,6) * t104 + t66) / 0.2e1 + t76 * t6 / 0.2e1 - t13 * mrSges(5,2) + (-t115 / 0.2e1 + t51 * t96) * t60 + (t89 / 0.2e1 - Ifges(5,6)) * t35 + t141 * t14 - t85 * t59 / 0.2e1 + ((-t15 * t76 - t16 * t73) * qJD(5) + t138) * mrSges(6,3);
t62 = pkin(4) + t116;
t26 = t92 * t51;
t8 = t85 * mrSges(6,1) - t84 * mrSges(6,2);
t1 = [-0.2e1 * t61 * t95 + t8 * t131 + t26 * t132 + 0.2e1 * t2 * t28 + 0.2e1 * t3 * t29 + 0.2e1 * t15 * t9 + 0.2e1 * t16 * t10 + t42 * t35 * t133 + 0.2e1 * m(5) * (t13 * t42 - t55 * t61 - t122) + 0.2e1 * m(6) * (t15 * t3 + t16 * t2 - t122) - (mrSges(5,3) * t131 - t73 * t19 + t76 * t20) * t34 + (-pkin(1) * t110 + t106 * t113) * t134 + 0.2e1 * t137 * qJD(3) * (t75 ^ 2 + t78 ^ 2) + (mrSges(5,1) * t130 + t13 * t133 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t35 - t94 * t34 + t111) * t50 + (mrSges(5,2) * t130 + mrSges(5,3) * t132 - 0.2e1 * Ifges(5,1) * t34 - t73 * t6 + t76 * t7 + (Ifges(6,5) * t76 + t94) * t35 + (-t76 * t19 - t73 * t20 - t50 * t89) * qJD(5)) * t51 + (0.2e1 * (mrSges(4,2) * pkin(1) - t113) * t75 + (pkin(1) * t100 + Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,2)) * t134) * t107; t80 + ((mrSges(4,2) * qJ(3) - Ifges(3,6) - Ifges(4,6)) * t75 + (-mrSges(4,1) * qJ(3) - t137 * pkin(1) + Ifges(3,5) + Ifges(4,5)) * t78) * qJD(2) + (-t78 * mrSges(4,2) + t100 * t75) * qJD(3) + t81 * t62 + (-t77 * t8 + (t74 * t26 + t88 * t77) * qJD(4) + m(5) * (t42 * t105 + t13 * t74 + t86) + m(6) * (t136 * t105 + t86) + (t77 * t34 - t74 * t35 + (-t50 * t77 + t51 * t74) * qJD(4)) * mrSges(5,3)) * t79; t83 + (0.2e1 * t98 * qJD(4) + (-0.2e1 * t52 + 0.2e1 * (-m(6) * t116 + t62 * t99 + t140) * qJD(4)) * t77) * t79; m(4) * t68 + t73 * t10 + t76 * t9 + t88 * qJD(5) + m(6) * (t136 * qJD(5) + t2 * t73 + t3 * t76) + m(5) * t55 + t95 + t110; 0; 0; t81 * pkin(4) + t80; (-t77 * t52 + (t98 + (pkin(4) * t99 + t140) * t77) * qJD(4)) * t79 + t83; 0; t83; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t102 - t85 * Ifges(6,6) + t111; t66 - t92 * t79 * t105 + (t56 * t62 - t123) * qJD(5); -t52; t66 + (t56 * pkin(4) - t123) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
