% Calculate time derivative of joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:25
% EndTime: 2019-12-31 20:49:27
% DurationCPUTime: 0.78s
% Computational Cost: add. (1032->155), mult. (2400->198), div. (0->0), fcn. (1790->6), ass. (0->87)
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t132 = t72 ^ 2 + t74 ^ 2;
t129 = (mrSges(5,3) + mrSges(6,2));
t131 = 2 * t129;
t97 = qJD(3) * t74;
t98 = qJD(3) * t72;
t128 = -mrSges(4,1) * t97 + mrSges(4,2) * t98;
t127 = m(6) * qJD(5);
t75 = cos(qJ(2));
t126 = t132 * t75;
t71 = sin(pkin(8));
t113 = pkin(3) * t71;
t58 = qJ(5) + t113;
t125 = m(6) * t58 + mrSges(6,3);
t73 = sin(qJ(2));
t124 = (t132 * mrSges(4,3) - mrSges(3,2)) * t75 + (-mrSges(4,1) * t74 + mrSges(4,2) * t72 - mrSges(3,1)) * t73;
t116 = m(5) * pkin(3);
t99 = cos(pkin(8));
t91 = t99 * pkin(3);
t60 = -t91 - pkin(4);
t123 = m(6) * t60 - t99 * t116 - mrSges(5,1) - mrSges(6,1);
t122 = t71 * t116 - mrSges(5,2) + t125;
t121 = 0.2e1 * m(5);
t120 = 0.2e1 * m(6);
t90 = t99 * t72;
t50 = t71 * t74 + t90;
t43 = t50 * qJD(3);
t106 = t71 * t72;
t84 = t99 * t74 - t106;
t44 = t84 * qJD(3);
t22 = t43 * mrSges(6,1) - t44 * mrSges(6,3);
t119 = 0.2e1 * t22;
t23 = t43 * mrSges(5,1) + t44 * mrSges(5,2);
t118 = 0.2e1 * t23;
t26 = -mrSges(6,1) * t84 - mrSges(6,3) * t50;
t117 = 0.2e1 * t26;
t115 = pkin(1) * t75;
t27 = -mrSges(5,1) * t84 + mrSges(5,2) * t50;
t114 = pkin(3) * t27;
t112 = mrSges(4,1) * t72;
t111 = Ifges(4,4) * t72;
t36 = t44 * mrSges(6,2);
t65 = pkin(3) * t98;
t101 = pkin(1) * qJD(2);
t66 = t73 * t101;
t53 = t66 + t65;
t109 = t53 * t27;
t102 = -qJ(4) - pkin(7);
t61 = pkin(1) * t73 + pkin(7);
t100 = -qJ(4) - t61;
t96 = t75 * t101;
t63 = -t74 * pkin(3) - pkin(2);
t68 = t74 * qJ(4);
t47 = t61 * t74 + t68;
t20 = -t100 * t90 + t47 * t71;
t21 = t100 * t106 + t99 * t47;
t67 = t74 * qJD(4);
t87 = t74 * t96;
t88 = qJD(3) * t100;
t28 = t72 * t88 + t67 + t87;
t78 = (-qJD(4) - t96) * t72 + t74 * t88;
t5 = t28 * t71 - t99 * t78;
t6 = t99 * t28 + t71 * t78;
t93 = t20 * t5 + t21 * t6;
t89 = qJD(3) * t102;
t42 = t72 * t89 + t67;
t82 = -t72 * qJD(4) + t74 * t89;
t17 = t42 * t71 - t99 * t82;
t18 = t99 * t42 + t71 * t82;
t56 = pkin(7) * t74 + t68;
t29 = -t102 * t90 + t56 * t71;
t30 = t102 * t106 + t99 * t56;
t92 = t29 * t17 + t30 * t18;
t85 = t17 * t20 + t18 * t21 + t29 * t5 + t30 * t6;
t83 = t22 + t23;
t24 = -pkin(4) * t84 - t50 * qJ(5) + t63;
t8 = t43 * pkin(4) - t44 * qJ(5) - t50 * qJD(5) + t65;
t81 = (Ifges(4,1) * t74 - t111) * t98 + (0.2e1 * Ifges(4,4) * t74 + (Ifges(4,1) - Ifges(4,2)) * t72) * t97 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t50 * t44 + 0.2e1 * (-Ifges(5,2) - Ifges(6,3)) * t84 * t43 + 0.2e1 * (Ifges(6,5) - Ifges(5,4)) * (t50 * t43 - t44 * t84);
t76 = qJD(5) * t84 * mrSges(6,2) + Ifges(4,5) * t97 - Ifges(4,6) * t98 + t60 * t36 + (-mrSges(5,3) * t91 + Ifges(6,4) + Ifges(5,5)) * t44 + (-t58 * mrSges(6,2) - mrSges(5,3) * t113 - Ifges(5,6) + Ifges(6,6)) * t43;
t62 = -pkin(2) - t115;
t57 = Ifges(4,2) * t74 + t111;
t54 = t63 - t115;
t52 = (mrSges(4,2) * t74 + t112) * qJD(3);
t19 = t24 - t115;
t7 = t66 + t8;
t1 = [(t53 * t54 + t93) * t121 + (t19 * t7 + t93) * t120 + t81 - t57 * t98 + t54 * t118 + 0.2e1 * t62 * t52 + 0.2e1 * t109 + t19 * t119 + t7 * t117 + 0.2e1 * (m(4) * (t126 * t61 + t62 * t73) + t124) * t101 + (t20 * t44 - t21 * t43 + t5 * t50 + t6 * t84) * t131; (t8 + t7) * t26 + (t54 + t63) * t23 + (t19 + t24) * t22 + t81 + t109 + (t62 - pkin(2)) * t52 + (-t57 + t114) * t98 + m(5) * (t53 * t63 + t54 * t65 + t85) + m(6) * (t19 * t8 + t24 * t7 + t85) + (m(4) * (-pkin(2) * t73 + t126 * pkin(7)) + t124) * t101 + t129 * ((t17 + t5) * t50 - (-t18 - t6) * t84 + (t20 + t29) * t44 + (-t21 - t30) * t43); t81 + (t63 * t65 + t92) * t121 + (t24 * t8 + t92) * t120 + (-t57 + 0.2e1 * t114) * t98 + t63 * t118 - 0.2e1 * pkin(2) * t52 + t24 * t119 + t8 * t117 + (t17 * t50 + t18 * t84 + t29 * t44 - t30 * t43) * t131; -mrSges(4,2) * t87 - t96 * t112 + t122 * t6 + t123 * t5 + t21 * t127 + t128 * t61 + t76; t128 * pkin(7) + t122 * t18 + t123 * t17 + t30 * t127 + t76; 0.2e1 * t125 * qJD(5); m(5) * t53 + m(6) * t7 + t83; m(5) * t65 + m(6) * t8 + t83; 0; 0; m(6) * t5 + t36; m(6) * t17 + t36; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
