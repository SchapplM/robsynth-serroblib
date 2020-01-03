% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:25
% DurationCPUTime: 2.04s
% Computational Cost: add. (861->205), mult. (2340->269), div. (0->0), fcn. (1441->4), ass. (0->90)
t66 = sin(pkin(7));
t64 = t66 ^ 2;
t67 = cos(pkin(7));
t65 = t67 ^ 2;
t131 = (t64 + t65) * (mrSges(4,2) + mrSges(3,3));
t97 = -mrSges(5,1) - mrSges(6,1);
t128 = Ifges(5,1) + Ifges(6,1);
t126 = Ifges(6,4) + Ifges(5,5);
t127 = -Ifges(5,4) + Ifges(6,5);
t69 = sin(qJ(4));
t70 = cos(qJ(4));
t41 = t66 * t69 + t67 * t70;
t34 = t41 * qJD(4);
t29 = qJD(1) * t34;
t118 = -t29 / 0.2e1;
t125 = Ifges(6,6) - Ifges(5,6);
t36 = t41 * qJD(1);
t107 = Ifges(6,5) * t36;
t32 = Ifges(5,4) * t36;
t91 = qJD(1) * t67;
t82 = t69 * t91;
t87 = t66 * qJD(1);
t37 = t70 * t87 - t82;
t124 = t126 * qJD(4) + t128 * t37 + t107 - t32;
t104 = t36 * mrSges(5,3);
t27 = -t36 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t95 = -qJD(4) * mrSges(5,2) - t104 + t27;
t103 = t37 * mrSges(5,3);
t94 = mrSges(6,2) * t37 + t97 * qJD(4) + t103;
t122 = -t66 * t70 + t67 * t69;
t96 = -pkin(6) + qJ(2);
t45 = t96 * t66;
t46 = t96 * t67;
t22 = t45 * t69 + t46 * t70;
t88 = qJD(4) * t70;
t83 = t66 * t88;
t30 = qJD(1) * t83 - qJD(4) * t82;
t51 = qJ(2) * t87 + qJD(3);
t39 = -pkin(6) * t87 + t51;
t44 = qJD(1) * t46;
t20 = t39 * t69 + t44 * t70;
t86 = qJD(1) * qJD(2);
t79 = t70 * t86;
t80 = t69 * t86;
t6 = qJD(4) * t20 - t66 * t79 + t67 * t80;
t121 = -t122 * t6 - t22 * t30;
t102 = t44 * t69;
t19 = t39 * t70 - t102;
t120 = -t19 + qJD(5);
t116 = -t36 / 0.2e1;
t115 = t36 / 0.2e1;
t113 = t37 / 0.2e1;
t72 = t70 * t45 - t46 * t69;
t111 = t72 * t6;
t109 = t6 * t70;
t108 = Ifges(5,4) * t37;
t105 = t29 * mrSges(6,2);
t99 = -qJD(4) / 0.2e1;
t98 = qJD(4) / 0.2e1;
t71 = qJD(1) ^ 2;
t92 = qJ(2) * t71;
t90 = qJD(3) * t66;
t89 = qJD(4) * t69;
t85 = -t67 * pkin(2) - t66 * qJ(3) - pkin(1);
t84 = t39 * t88 + t66 * t80 + t67 * t79;
t78 = qJD(3) * t87;
t77 = qJ(2) * t86;
t38 = t67 * pkin(3) - t85;
t33 = -qJD(1) * pkin(1) - pkin(2) * t91 - qJ(3) * t87 + qJD(2);
t74 = t65 * t77;
t23 = pkin(3) * t91 - t33;
t43 = (-mrSges(4,1) * t67 - mrSges(4,3) * t66) * qJD(1);
t58 = t65 * t92;
t52 = t64 * t77;
t35 = -t67 * t89 + t83;
t31 = Ifges(6,5) * t37;
t18 = mrSges(5,1) * t36 + mrSges(5,2) * t37;
t17 = mrSges(6,1) * t36 - mrSges(6,3) * t37;
t16 = pkin(4) * t37 + qJ(5) * t36;
t15 = qJD(4) * qJ(5) + t20;
t14 = -qJD(4) * pkin(4) + t120;
t11 = -Ifges(5,2) * t36 + Ifges(5,6) * qJD(4) + t108;
t10 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t36 + t31;
t9 = pkin(4) * t41 + qJ(5) * t122 + t38;
t5 = -t44 * t89 + t84;
t4 = pkin(4) * t35 + qJ(5) * t34 + qJD(5) * t122 + t90;
t3 = pkin(4) * t36 - qJ(5) * t37 + t23;
t2 = (qJD(5) - t102) * qJD(4) + t84;
t1 = pkin(4) * t30 + qJ(5) * t29 - qJD(5) * t37 + t78;
t7 = [m(5) * (t22 * t5 - t111) + m(6) * (t1 * t9 + t2 * t22 + t3 * t4 - t111) + m(4) * (0.2e1 * t74 + t52 + (t51 * qJD(2) + (-qJD(1) * t85 - t33) * qJD(3)) * t66) + t38 * (t30 * mrSges(5,1) - t29 * mrSges(5,2)) + t9 * (t30 * mrSges(6,1) + t29 * mrSges(6,3)) + t4 * t17 + 0.2e1 * m(3) * (t52 + t74) + t72 * t105 + (m(5) * (qJD(1) * t38 + t23) - 0.2e1 * t43 + t18) * t90 + (-m(5) * t19 + m(6) * t14 + t94) * (qJD(2) * t122 + qJD(4) * t22) + (m(5) * t20 + m(6) * t15 + t95) * (qJD(2) * t41 + qJD(4) * t72) - (t127 * t30 - t128 * t29) * t122 / 0.2e1 + t121 * mrSges(6,2) + (t29 * t72 + t121) * mrSges(5,3) - (mrSges(5,2) * t78 - t1 * mrSges(6,3) + (Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t30 + t128 * t118) * t122 + (mrSges(5,1) * t78 + t1 * mrSges(6,1) - t2 * mrSges(6,2) - t5 * mrSges(5,3) + (Ifges(5,2) + Ifges(6,3)) * t30 + 0.2e1 * t127 * t118) * t41 + (t3 * mrSges(6,1) + t10 / 0.2e1 + t23 * mrSges(5,1) - t11 / 0.2e1 + Ifges(6,3) * t115 - Ifges(5,2) * t116 + t125 * t98 + t127 * t113 - t15 * mrSges(6,2) - t20 * mrSges(5,3)) * t35 + (-t23 * mrSges(5,2) + t3 * mrSges(6,3) - Ifges(5,4) * t116 - Ifges(6,5) * t115 - t126 * t98 - t124 / 0.2e1 - t128 * t113 - t14 * mrSges(6,2) + t19 * mrSges(5,3)) * t34 + 0.2e1 * t86 * t131; t94 * t37 - t95 * t36 + t97 * t30 - (-mrSges(5,2) + mrSges(6,3)) * t29 + (-m(4) - m(5)) * t78 - m(3) * (t64 * t92 + t58) - m(4) * (t51 * t87 + t58) - m(5) * (t19 * t37 + t20 * t36) - t71 * t131 + (t14 * t37 - t15 * t36 - t1) * m(6); (t94 * t69 + t95 * t70) * qJD(4) + m(5) * (t5 * t69 - t109 + (-t19 * t69 + t20 * t70) * qJD(4)) + m(6) * (t2 * t69 - t109 + (t14 * t69 + t15 * t70) * qJD(4)) + (-m(5) * t23 - m(6) * t3 - t17 - t18 + t43 + (qJD(2) + t33) * m(4)) * t87 + (mrSges(6,2) + mrSges(5,3)) * (t29 * t70 - t30 * t69); -t23 * (mrSges(5,1) * t37 - mrSges(5,2) * t36) + t11 * t113 + (Ifges(6,3) * t37 - t107) * t116 - t3 * (t37 * mrSges(6,1) + t36 * mrSges(6,3)) + qJD(5) * t27 - t16 * t17 - t5 * mrSges(5,2) + t2 * mrSges(6,3) + t97 * t6 + t125 * t30 - t126 * t29 + (-t94 + t103) * t20 + (-t95 - t104) * t19 + (pkin(4) * t29 - qJ(5) * t30 + t14 * t36 + t15 * t37) * mrSges(6,2) + (t125 * t37 - t126 * t36) * t99 + (-t6 * pkin(4) + t2 * qJ(5) + t120 * t15 - t14 * t20 - t3 * t16) * m(6) + (-Ifges(5,2) * t37 + t124 - t32) * t115 - (-t128 * t36 + t10 - t108 + t31) * t37 / 0.2e1; -t105 - qJD(4) * t27 + t37 * t17 + 0.2e1 * (t6 / 0.2e1 + t15 * t99 + t3 * t113) * m(6);];
tauc = t7(:);
