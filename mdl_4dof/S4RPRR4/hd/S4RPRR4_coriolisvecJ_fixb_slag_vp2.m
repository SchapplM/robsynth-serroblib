% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:18
% DurationCPUTime: 1.20s
% Computational Cost: add. (903->190), mult. (2304->277), div. (0->0), fcn. (1221->6), ass. (0->99)
t117 = qJD(3) / 0.2e1;
t53 = sin(qJ(3));
t88 = qJD(1) * t53;
t116 = t88 / 0.2e1;
t73 = Ifges(4,5) * t117;
t48 = sin(pkin(7)) * pkin(1) + pkin(5);
t42 = t48 * qJD(1);
t55 = cos(qJ(3));
t31 = qJD(2) * t55 - t53 * t42;
t87 = qJD(1) * t55;
t49 = Ifges(4,4) * t87;
t54 = cos(qJ(4));
t105 = t54 / 0.2e1;
t52 = sin(qJ(4));
t106 = -t52 / 0.2e1;
t86 = qJD(3) * t52;
t39 = t54 * t88 + t86;
t108 = t39 / 0.2e1;
t102 = Ifges(5,4) * t39;
t84 = qJD(3) * t54;
t38 = -t52 * t88 + t84;
t47 = qJD(4) - t87;
t11 = Ifges(5,2) * t38 + Ifges(5,6) * t47 + t102;
t37 = Ifges(5,4) * t38;
t12 = Ifges(5,1) * t39 + Ifges(5,5) * t47 + t37;
t26 = -qJD(3) * pkin(3) - t31;
t100 = Ifges(5,4) * t54;
t63 = -Ifges(5,2) * t52 + t100;
t101 = Ifges(5,4) * t52;
t65 = Ifges(5,1) * t54 - t101;
t66 = mrSges(5,1) * t52 + mrSges(5,2) * t54;
t32 = t53 * qJD(2) + t42 * t55;
t27 = qJD(3) * pkin(6) + t32;
t76 = -cos(pkin(7)) * pkin(1) - pkin(2);
t36 = -pkin(3) * t55 - t53 * pkin(6) + t76;
t30 = t36 * qJD(1);
t8 = -t27 * t52 + t30 * t54;
t9 = t27 * t54 + t30 * t52;
t68 = t9 * t52 + t8 * t54;
t96 = Ifges(5,6) * t52;
t98 = Ifges(5,5) * t54;
t56 = -t68 * mrSges(5,3) + t47 * (-t96 + t98) / 0.2e1 + t38 * t63 / 0.2e1 + t65 * t108 + t26 * t66 + t12 * t105 + t11 * t106;
t115 = Ifges(4,1) * t116 + t49 / 0.2e1 + t73 - t31 * mrSges(4,3) + t56;
t83 = qJD(3) * t55;
t28 = qJD(3) * t31;
t70 = pkin(3) * t53 - pkin(6) * t55;
t41 = t70 * qJD(3);
t35 = qJD(1) * t41;
t1 = t8 * qJD(4) + t28 * t54 + t35 * t52;
t2 = -t9 * qJD(4) - t28 * t52 + t35 * t54;
t69 = t1 * t54 - t2 * t52;
t80 = qJD(3) * qJD(4);
t82 = qJD(4) * t52;
t22 = t54 * t80 + (-t53 * t82 + t54 * t83) * qJD(1);
t81 = qJD(4) * t54;
t23 = -t52 * t80 + (-t52 * t83 - t53 * t81) * qJD(1);
t113 = -t2 * mrSges(5,1) + t1 * mrSges(5,2) - Ifges(5,5) * t22 - Ifges(5,6) * t23;
t112 = t22 / 0.2e1;
t111 = t23 / 0.2e1;
t110 = -t38 / 0.2e1;
t109 = -t39 / 0.2e1;
t107 = -t47 / 0.2e1;
t44 = t76 * qJD(1);
t93 = t44 * mrSges(4,2);
t92 = t48 * t55;
t78 = mrSges(4,3) * t88;
t91 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t38 - mrSges(5,2) * t39 - t78;
t89 = Ifges(4,6) * qJD(3);
t29 = qJD(3) * t32;
t85 = qJD(3) * t53;
t79 = t48 * t85;
t77 = mrSges(4,3) * t87;
t75 = m(4) * t48 + mrSges(4,3);
t74 = qJD(1) * t85;
t72 = -t89 / 0.2e1;
t67 = mrSges(5,1) * t54 - mrSges(5,2) * t52;
t64 = Ifges(5,1) * t52 + t100;
t62 = Ifges(5,2) * t54 + t101;
t61 = Ifges(5,5) * t52 + Ifges(5,6) * t54;
t15 = mrSges(5,1) * t74 - mrSges(5,3) * t22;
t16 = -mrSges(5,2) * t74 + mrSges(5,3) * t23;
t60 = -t52 * t15 + t54 * t16;
t24 = -mrSges(5,2) * t47 + mrSges(5,3) * t38;
t25 = mrSges(5,1) * t47 - mrSges(5,3) * t39;
t59 = -t52 * t24 - t54 * t25;
t18 = t52 * t36 + t54 * t92;
t17 = t54 * t36 - t52 * t92;
t58 = t9 * mrSges(5,2) - Ifges(5,3) * t47 - Ifges(5,6) * t38 - Ifges(5,5) * t39 + t89 / 0.2e1 + (Ifges(4,4) * t53 + Ifges(4,2) * t55) * qJD(1) / 0.2e1 - t44 * mrSges(4,1) - t8 * mrSges(5,1);
t46 = Ifges(5,3) * t74;
t45 = -qJD(3) * mrSges(4,2) + t77;
t40 = t70 * qJD(1);
t14 = t31 * t54 + t40 * t52;
t13 = -t31 * t52 + t40 * t54;
t7 = -mrSges(5,1) * t23 + mrSges(5,2) * t22;
t6 = Ifges(5,1) * t22 + Ifges(5,4) * t23 + Ifges(5,5) * t74;
t5 = Ifges(5,4) * t22 + Ifges(5,2) * t23 + Ifges(5,6) * t74;
t4 = -t18 * qJD(4) + t54 * t41 + t52 * t79;
t3 = t17 * qJD(4) + t52 * t41 - t54 * t79;
t10 = [m(5) * (t1 * t18 + t2 * t17 + t9 * t3 + t8 * t4) + t17 * t15 + t18 * t16 + t3 * t24 + t4 * t25 + (-t46 / 0.2e1 + t75 * t28 + (0.3e1 / 0.2e1 * t49 + t73 + 0.2e1 * t93 + t115) * qJD(3) + t113) * t55 + (t65 * t112 + t63 * t111 + t5 * t106 + t6 * t105 + (-t1 * t52 - t2 * t54) * mrSges(5,3) + (mrSges(4,3) + t66) * t29 + (t26 * t67 + t61 * t107 + t62 * t110 + t64 * t109 - t54 * t11 / 0.2e1 + t12 * t106 + (t8 * t52 - t9 * t54) * mrSges(5,3)) * qJD(4) + (t72 - t75 * t32 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t98 / 0.2e1 - t96 / 0.2e1) * t53 + t76 * mrSges(4,1) + (-0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,1)) * t55) * qJD(1) - t58) * qJD(3)) * t53 + ((-m(4) * t31 + m(5) * t26 - t91) * t83 + (t7 + (m(4) + m(5)) * t29 - qJD(3) * t45) * t53) * t48; (-t7 + (t54 * t24 - t52 * t25 + t45 - t77) * qJD(3) + m(5) * (-t8 * t86 + t9 * t84 - t29)) * t55 + (t59 * qJD(4) + (-t78 - t91) * qJD(3) + m(5) * (qJD(3) * t26 - t8 * t81 - t9 * t82 + t69) + t60) * t53; t64 * t112 + t62 * t111 + t5 * t105 + t52 * t6 / 0.2e1 - t31 * t45 - t28 * mrSges(4,2) - t14 * t24 - t13 * t25 - pkin(3) * t7 + t91 * t32 + (-mrSges(4,1) - t67) * t29 + t69 * mrSges(5,3) + t56 * qJD(4) + ((t32 * mrSges(4,3) + Ifges(4,4) * t116 + t61 * t117 + t58 + t72) * t53 + (t73 - t93 - t49 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t88 - t115) * t55) * qJD(1) + (-pkin(3) * t29 - t13 * t8 - t14 * t9 - t26 * t32) * m(5) + (t60 + m(5) * t69 + (-m(5) * t68 + t59) * qJD(4)) * pkin(6); t46 - t26 * (mrSges(5,1) * t39 + mrSges(5,2) * t38) + (Ifges(5,1) * t38 - t102) * t109 + t11 * t108 + (Ifges(5,5) * t38 - Ifges(5,6) * t39) * t107 - t8 * t24 + t9 * t25 + (t38 * t8 + t39 * t9) * mrSges(5,3) + (-Ifges(5,2) * t39 + t12 + t37) * t110 - t113;];
tauc = t10(:);
