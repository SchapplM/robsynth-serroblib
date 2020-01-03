% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:32
% DurationCPUTime: 1.42s
% Computational Cost: add. (1559->193), mult. (4411->291), div. (0->0), fcn. (3206->6), ass. (0->89)
t81 = sin(pkin(7));
t84 = sin(qJ(3));
t82 = cos(pkin(7));
t86 = cos(qJ(3));
t98 = t82 * t86;
t69 = -t81 * t84 + t98;
t61 = t69 * qJD(1);
t70 = t81 * t86 + t82 * t84;
t62 = t70 * qJD(1);
t83 = sin(qJ(4));
t85 = cos(qJ(4));
t43 = t61 * t83 + t62 * t85;
t102 = Ifges(5,4) * t43;
t91 = t85 * t61 - t62 * t83;
t37 = Ifges(5,4) * t91;
t80 = qJD(3) + qJD(4);
t13 = Ifges(5,1) * t43 + Ifges(5,5) * t80 + t37;
t63 = t69 * qJD(3);
t56 = qJD(1) * t63;
t64 = t70 * qJD(3);
t57 = qJD(1) * t64;
t17 = t91 * qJD(4) + t56 * t85 - t57 * t83;
t18 = -qJD(4) * t43 - t56 * t83 - t57 * t85;
t97 = pkin(5) + qJ(2);
t74 = t97 * t81;
t71 = qJD(1) * t74;
t75 = t97 * t82;
t72 = qJD(1) * t75;
t77 = qJD(2) * t98;
t94 = qJD(1) * qJD(2);
t95 = qJD(3) * t86;
t26 = -t71 * t95 + qJD(1) * t77 + (-qJD(3) * t72 - t81 * t94) * t84;
t22 = -pkin(6) * t57 + t26;
t47 = -t71 * t84 + t72 * t86;
t88 = t70 * qJD(2);
t27 = -qJD(1) * t88 - qJD(3) * t47;
t23 = -pkin(6) * t56 + t27;
t30 = pkin(6) * t61 + t47;
t101 = t30 * t83;
t46 = -t86 * t71 - t72 * t84;
t29 = -pkin(6) * t62 + t46;
t28 = qJD(3) * pkin(3) + t29;
t6 = t28 * t85 - t101;
t2 = qJD(4) * t6 + t22 * t85 + t23 * t83;
t100 = t30 * t85;
t7 = t28 * t83 + t100;
t3 = -qJD(4) * t7 - t22 * t83 + t23 * t85;
t93 = -pkin(2) * t82 - pkin(1);
t73 = qJD(1) * t93 + qJD(2);
t48 = -pkin(3) * t61 + t73;
t119 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t17 + Ifges(5,6) * t18 - (Ifges(5,5) * t91 - Ifges(5,6) * t43) * t80 / 0.2e1 + (t43 * t7 + t6 * t91) * mrSges(5,3) - (-Ifges(5,2) * t43 + t13 + t37) * t91 / 0.2e1 - t48 * (mrSges(5,1) * t43 + mrSges(5,2) * t91) - (Ifges(5,1) * t91 - t102) * t43 / 0.2e1;
t12 = Ifges(5,2) * t91 + Ifges(5,6) * t80 + t102;
t117 = t12 / 0.2e1;
t50 = -t84 * t74 + t86 * t75;
t112 = (m(3) * qJ(2) + mrSges(3,3)) * (-t81 ^ 2 - t82 ^ 2);
t110 = t91 / 0.2e1;
t108 = t43 / 0.2e1;
t106 = t63 / 0.2e1;
t105 = -t64 / 0.2e1;
t103 = Ifges(4,4) * t62;
t92 = -t18 * mrSges(5,1) + t17 * mrSges(5,2);
t49 = -t86 * t74 - t75 * t84;
t35 = -pkin(6) * t70 + t49;
t36 = pkin(6) * t69 + t50;
t10 = t35 * t85 - t36 * t83;
t11 = t35 * t83 + t36 * t85;
t44 = t69 * t85 - t70 * t83;
t45 = t69 * t83 + t70 * t85;
t33 = -t74 * t95 + t77 + (-qJD(2) * t81 - qJD(3) * t75) * t84;
t34 = -t50 * qJD(3) - t88;
t58 = Ifges(4,4) * t61;
t54 = t56 * mrSges(4,2);
t53 = -pkin(3) * t69 + t93;
t52 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t62;
t51 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t61;
t39 = Ifges(4,1) * t62 + Ifges(4,5) * qJD(3) + t58;
t38 = Ifges(4,2) * t61 + Ifges(4,6) * qJD(3) + t103;
t32 = mrSges(5,1) * t80 - mrSges(5,3) * t43;
t31 = -mrSges(5,2) * t80 + mrSges(5,3) * t91;
t25 = -pkin(6) * t63 + t34;
t24 = -pkin(6) * t64 + t33;
t21 = -qJD(4) * t45 - t63 * t83 - t64 * t85;
t20 = qJD(4) * t44 + t63 * t85 - t64 * t83;
t19 = -mrSges(5,1) * t91 + mrSges(5,2) * t43;
t9 = t29 * t85 - t101;
t8 = -t29 * t83 - t100;
t5 = -qJD(4) * t11 - t24 * t83 + t25 * t85;
t4 = qJD(4) * t10 + t24 * t85 + t25 * t83;
t1 = [(-t10 * t17 + t11 * t18 + t2 * t44 - t20 * t6 + t21 * t7 - t3 * t45) * mrSges(5,3) + qJD(3) * (Ifges(4,5) * t63 - Ifges(4,6) * t64) / 0.2e1 + t73 * (mrSges(4,1) * t64 + mrSges(4,2) * t63) + (t26 * t69 - t27 * t70 - t46 * t63 - t47 * t64 - t49 * t56 - t50 * t57) * mrSges(4,3) + m(5) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + (t48 * t64 + t53 * t57) * pkin(3)) + (t64 * t19 + t57 * (-mrSges(5,1) * t44 + mrSges(5,2) * t45)) * pkin(3) + (t105 * t61 - t69 * t57) * Ifges(4,2) + (t105 * t62 + t106 * t61 + t69 * t56 - t57 * t70) * Ifges(4,4) + t93 * (t57 * mrSges(4,1) + t54) - 0.2e1 * t112 * t94 + t80 * (Ifges(5,5) * t20 + Ifges(5,6) * t21) / 0.2e1 + t48 * (-t21 * mrSges(5,1) + t20 * mrSges(5,2)) + t33 * t51 + t34 * t52 + t4 * t31 + t5 * t32 + t20 * t13 / 0.2e1 + (t106 * t62 + t56 * t70) * Ifges(4,1) + (t108 * t20 + t17 * t45) * Ifges(5,1) + (t110 * t21 + t18 * t44) * Ifges(5,2) + (t108 * t21 + t110 * t20 + t17 * t44 + t18 * t45) * Ifges(5,4) + t53 * t92 + m(4) * (t26 * t50 + t27 * t49 + t33 * t47 + t34 * t46) + t38 * t105 + t39 * t106 + t21 * t117; -t91 * t31 + t43 * t32 - t61 * t51 + t62 * t52 + t54 - (-m(5) * pkin(3) - mrSges(4,1)) * t57 - m(4) * (-t46 * t62 + t47 * t61) - m(5) * (-t43 * t6 + t7 * t91) + t112 * qJD(1) ^ 2 + t92; -m(5) * (t6 * t8 + t7 * t9) - (-Ifges(4,2) * t62 + t39 + t58) * t61 / 0.2e1 + t43 * t117 + (t46 * t61 + t47 * t62) * mrSges(4,3) - t73 * (mrSges(4,1) * t62 + mrSges(4,2) * t61) - qJD(3) * (Ifges(4,5) * t61 - Ifges(4,6) * t62) / 0.2e1 + t62 * t38 / 0.2e1 + Ifges(4,5) * t56 - Ifges(4,6) * t57 - t46 * t51 + t47 * t52 - t9 * t31 - t8 * t32 - t26 * mrSges(4,2) + t27 * mrSges(4,1) - t62 * (Ifges(4,1) * t61 - t103) / 0.2e1 + (-t62 * t19 + (t31 * t85 - t32 * t83) * qJD(4) + (-t17 * t85 + t18 * t83) * mrSges(5,3) + (t2 * t83 + t3 * t85 - t48 * t62 + (-t6 * t83 + t7 * t85) * qJD(4)) * m(5)) * pkin(3) + t119; t12 * t108 - t6 * t31 + t7 * t32 + t119;];
tauc = t1(:);
