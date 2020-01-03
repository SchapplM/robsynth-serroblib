% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:53:46
% DurationCPUTime: 1.44s
% Computational Cost: add. (1559->204), mult. (4303->295), div. (0->0), fcn. (2955->6), ass. (0->102)
t73 = cos(qJ(4));
t116 = t73 / 0.2e1;
t71 = sin(qJ(4));
t117 = -t71 / 0.2e1;
t69 = sin(pkin(7));
t70 = cos(pkin(7));
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t60 = t69 * t74 + t70 * t72;
t55 = t60 * qJD(1);
t45 = qJD(3) * t71 + t55 * t73;
t120 = t45 / 0.2e1;
t111 = Ifges(5,4) * t45;
t44 = qJD(3) * t73 - t55 * t71;
t59 = t69 * t72 - t74 * t70;
t54 = t59 * qJD(1);
t53 = qJD(4) + t54;
t13 = Ifges(5,2) * t44 + Ifges(5,6) * t53 + t111;
t43 = Ifges(5,4) * t44;
t14 = Ifges(5,1) * t45 + Ifges(5,5) * t53 + t43;
t103 = pkin(5) + qJ(2);
t64 = t103 * t69;
t61 = qJD(1) * t64;
t65 = t103 * t70;
t62 = qJD(1) * t65;
t39 = -t61 * t74 - t62 * t72;
t34 = -qJD(3) * pkin(3) - t39;
t85 = Ifges(5,5) * t73 - Ifges(5,6) * t71;
t109 = Ifges(5,4) * t73;
t87 = -Ifges(5,2) * t71 + t109;
t110 = Ifges(5,4) * t71;
t89 = Ifges(5,1) * t73 - t110;
t90 = mrSges(5,1) * t71 + mrSges(5,2) * t73;
t98 = -pkin(2) * t70 - pkin(1);
t63 = t98 * qJD(1) + qJD(2);
t28 = pkin(3) * t54 - pkin(6) * t55 + t63;
t40 = -t61 * t72 + t62 * t74;
t35 = qJD(3) * pkin(6) + t40;
t8 = t28 * t73 - t35 * t71;
t9 = t28 * t71 + t35 * t73;
t93 = t9 * t71 + t8 * t73;
t130 = -t93 * mrSges(5,3) + t53 * t85 / 0.2e1 + t89 * t120 + t44 * t87 / 0.2e1 + t34 * t90 + t14 * t116 + t13 * t117;
t52 = Ifges(4,4) * t54;
t126 = Ifges(4,1) * t55 / 0.2e1 - t52 / 0.2e1;
t129 = t63 * mrSges(4,2) + Ifges(4,5) * qJD(3) + t126 + t130;
t79 = t59 * qJD(2);
t19 = -qJD(1) * t79 + t39 * qJD(3);
t56 = t59 * qJD(3);
t50 = qJD(1) * t56;
t57 = t60 * qJD(3);
t51 = qJD(1) * t57;
t31 = pkin(3) * t51 + pkin(6) * t50;
t1 = t8 * qJD(4) + t19 * t73 + t31 * t71;
t2 = -t9 * qJD(4) - t19 * t71 + t31 * t73;
t128 = t1 * t73 - t2 * t71;
t24 = qJD(4) * t44 - t50 * t73;
t25 = -qJD(4) * t45 + t50 * t71;
t127 = t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t24 + Ifges(5,6) * t25;
t125 = (m(3) * qJ(2) + mrSges(3,3)) * (t69 ^ 2 + t70 ^ 2);
t124 = t24 / 0.2e1;
t123 = t25 / 0.2e1;
t122 = -t44 / 0.2e1;
t121 = -t45 / 0.2e1;
t119 = t51 / 0.2e1;
t118 = -t53 / 0.2e1;
t107 = Ifges(4,2) * t54;
t80 = t60 * qJD(2);
t20 = qJD(1) * t80 + t40 * qJD(3);
t82 = -t74 * t64 - t65 * t72;
t104 = t20 * t82;
t102 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t44 + mrSges(5,2) * t45 + t55 * mrSges(4,3);
t97 = t51 * mrSges(4,1) - t50 * mrSges(4,2);
t94 = -t1 * t71 - t2 * t73;
t92 = t8 * t71 - t9 * t73;
t91 = t73 * mrSges(5,1) - t71 * mrSges(5,2);
t88 = Ifges(5,1) * t71 + t109;
t86 = Ifges(5,2) * t73 + t110;
t84 = Ifges(5,5) * t71 + Ifges(5,6) * t73;
t26 = -mrSges(5,2) * t53 + mrSges(5,3) * t44;
t27 = mrSges(5,1) * t53 - mrSges(5,3) * t45;
t83 = t73 * t26 - t71 * t27;
t38 = pkin(3) * t59 - pkin(6) * t60 + t98;
t42 = -t64 * t72 + t65 * t74;
t17 = t38 * t73 - t42 * t71;
t18 = t38 * t71 + t42 * t73;
t78 = t9 * mrSges(5,2) - Ifges(5,3) * t53 - Ifges(5,6) * t44 - Ifges(5,5) * t45 - t107 / 0.2e1 + Ifges(4,6) * qJD(3) + Ifges(4,4) * t55 - t63 * mrSges(4,1) - t8 * mrSges(5,1);
t49 = Ifges(5,3) * t51;
t46 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t54;
t37 = pkin(3) * t57 + pkin(6) * t56;
t36 = pkin(3) * t55 + pkin(6) * t54;
t30 = t42 * qJD(3) + t80;
t29 = t82 * qJD(3) - t79;
t16 = -mrSges(5,2) * t51 + mrSges(5,3) * t25;
t15 = mrSges(5,1) * t51 - mrSges(5,3) * t24;
t11 = t36 * t71 + t39 * t73;
t10 = t36 * t73 - t39 * t71;
t7 = -mrSges(5,1) * t25 + mrSges(5,2) * t24;
t6 = Ifges(5,1) * t24 + Ifges(5,4) * t25 + Ifges(5,5) * t51;
t5 = Ifges(5,4) * t24 + Ifges(5,2) * t25 + Ifges(5,6) * t51;
t4 = -t18 * qJD(4) - t29 * t71 + t37 * t73;
t3 = t17 * qJD(4) + t29 * t73 + t37 * t71;
t12 = [t98 * t97 + t29 * t46 - t82 * t7 + t3 * t26 + t4 * t27 + t17 * t15 + t18 * t16 + (t107 / 0.2e1 - t78) * t57 - (t126 + t129) * t56 + t102 * t30 + m(5) * (t1 * t18 + t17 * t2 + t3 * t9 + t30 * t34 + t4 * t8 - t104) + m(4) * (t19 * t42 + t29 * t40 - t30 * t39 - t104) + 0.2e1 * t125 * qJD(2) * qJD(1) + (t39 * t56 - t40 * t57 - t42 * t51 + t50 * t82) * mrSges(4,3) + (Ifges(4,4) * t50 + t49 / 0.2e1 - t19 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t51 + t127) * t59 + (t85 * t119 + t89 * t124 + t87 * t123 - Ifges(4,4) * t51 - Ifges(4,1) * t50 + t6 * t116 + t5 * t117 + (mrSges(4,3) + t90) * t20 + t94 * mrSges(5,3) + (t34 * t91 + t84 * t118 + t88 * t121 + t86 * t122 + t14 * t117 - t73 * t13 / 0.2e1 + t92 * mrSges(5,3)) * qJD(4)) * t60; t73 * t15 + t71 * t16 - t102 * t55 + t83 * qJD(4) - (-t46 - t83) * t54 - m(4) * (-t39 * t55 - t40 * t54) + t97 - t125 * qJD(1) ^ 2 + (-t34 * t55 - t53 * t92 - t94) * m(5); t84 * t119 + t88 * t124 + t86 * t123 + t5 * t116 + t71 * t6 / 0.2e1 - Ifges(4,6) * t51 - t39 * t46 - Ifges(4,5) * t50 - t11 * t26 - t10 * t27 - t19 * mrSges(4,2) - pkin(3) * t7 - t102 * t40 + (-mrSges(4,1) - t91) * t20 + t128 * mrSges(5,3) + (t40 * mrSges(4,3) + t78) * t55 - (t52 / 0.2e1 + t39 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t55 - t129) * t54 + t130 * qJD(4) + (-pkin(3) * t20 - t10 * t8 - t11 * t9 - t34 * t40) * m(5) + (-t71 * t15 + t73 * t16 + m(5) * t128 + (-m(5) * t93 - t71 * t26 - t73 * t27) * qJD(4)) * pkin(6); t49 - t34 * (mrSges(5,1) * t45 + mrSges(5,2) * t44) + (Ifges(5,1) * t44 - t111) * t121 + t13 * t120 + (Ifges(5,5) * t44 - Ifges(5,6) * t45) * t118 - t8 * t26 + t9 * t27 + (t44 * t8 + t45 * t9) * mrSges(5,3) + (-Ifges(5,2) * t45 + t14 + t43) * t122 + t127;];
tauc = t12(:);
