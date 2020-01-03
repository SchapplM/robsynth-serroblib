% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:07
% DurationCPUTime: 0.71s
% Computational Cost: add. (843->145), mult. (1854->227), div. (0->0), fcn. (1010->4), ass. (0->83)
t103 = ((m(3) + m(4)) * qJ(2));
t61 = qJD(3) + qJD(4);
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t64 = cos(qJ(4));
t65 = cos(qJ(3));
t42 = -t62 * t65 - t64 * t63;
t25 = t61 * t42;
t102 = t25 / 0.2e1;
t38 = t42 * qJD(1);
t101 = -t38 / 0.2e1;
t80 = qJD(1) * t65;
t81 = qJD(1) * t63;
t39 = -t62 * t81 + t64 * t80;
t100 = -t39 / 0.2e1;
t99 = t39 / 0.2e1;
t98 = t63 / 0.2e1;
t97 = -t65 / 0.2e1;
t66 = -pkin(1) - pkin(5);
t95 = pkin(6) - t66;
t94 = mrSges(5,3) * t38;
t93 = Ifges(4,4) * t63;
t92 = Ifges(4,4) * t65;
t91 = Ifges(5,4) * t39;
t52 = t66 * qJD(1) + qJD(2);
t33 = -pkin(6) * t81 + t52 * t63;
t90 = t33 * t62;
t89 = t33 * t64;
t88 = t39 * mrSges(5,3);
t84 = t64 * t65;
t71 = t62 * t63 - t84;
t21 = t61 * t71 * qJD(1);
t87 = t42 * t21;
t20 = t25 * qJD(1);
t86 = t71 * t20;
t83 = Ifges(4,5) * qJD(3);
t82 = Ifges(4,6) * qJD(3);
t79 = qJD(3) * t63;
t78 = qJD(3) * t65;
t77 = qJD(4) * t62;
t76 = qJD(4) * t64;
t48 = t95 * t65;
t23 = -mrSges(5,1) * t38 + mrSges(5,2) * t39;
t57 = pkin(3) * t63 + qJ(2);
t51 = t57 * qJD(1);
t75 = -m(5) * t51 - t23;
t74 = pkin(6) * qJD(1) - t52;
t34 = -pkin(6) * t80 + t65 * t52;
t53 = pkin(3) * t78 + qJD(2);
t73 = mrSges(4,1) * t63 + mrSges(4,2) * t65;
t30 = qJD(3) * pkin(3) + t34;
t11 = t30 * t64 - t90;
t12 = t30 * t62 + t89;
t47 = t95 * t63;
t27 = -t47 * t64 - t48 * t62;
t26 = t47 * t62 - t48 * t64;
t49 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t81;
t50 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t80;
t72 = t65 * t49 - t63 * t50;
t70 = qJ(2) * (mrSges(4,1) * t65 - mrSges(4,2) * t63);
t31 = t74 * t79;
t32 = t74 * t78;
t2 = qJD(4) * t11 + t31 * t62 - t32 * t64;
t24 = t61 * t84 - t62 * t79 - t63 * t77;
t3 = -qJD(4) * t12 + t31 * t64 + t32 * t62;
t69 = t11 * t25 + t12 * t24 - t2 * t42 - t3 * t71;
t16 = Ifges(5,2) * t38 + Ifges(5,6) * t61 + t91;
t35 = Ifges(5,4) * t38;
t17 = Ifges(5,1) * t39 + Ifges(5,5) * t61 + t35;
t68 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + t11 * t94 + t16 * t99 - t51 * (mrSges(5,1) * t39 + mrSges(5,2) * t38) + (Ifges(5,1) * t38 - t91) * t100 - t61 * (Ifges(5,5) * t38 - Ifges(5,6) * t39) / 0.2e1 + Ifges(5,6) * t21 + Ifges(5,5) * t20 + (-Ifges(5,2) * t39 + t17 + t35) * t101;
t46 = t53 * qJD(1);
t44 = qJD(1) * t73;
t41 = qJD(3) * t48;
t40 = t95 * t79;
t37 = t83 + (Ifges(4,1) * t65 - t93) * qJD(1);
t36 = t82 + (-Ifges(4,2) * t63 + t92) * qJD(1);
t29 = mrSges(5,1) * t61 - t88;
t28 = -mrSges(5,2) * t61 + t94;
t14 = t34 * t64 - t90;
t13 = -t34 * t62 - t89;
t5 = -t27 * qJD(4) + t40 * t64 + t41 * t62;
t4 = t26 * qJD(4) + t40 * t62 - t41 * t64;
t1 = [t57 * (-t21 * mrSges(5,1) + t20 * mrSges(5,2)) + t61 * (Ifges(5,5) * t25 - Ifges(5,6) * t24) / 0.2e1 + qJD(2) * t44 + t46 * (-t42 * mrSges(5,1) - mrSges(5,2) * t71) + t51 * (t24 * mrSges(5,1) + t25 * mrSges(5,2)) + t53 * t23 + t5 * t29 - t24 * t16 / 0.2e1 + t17 * t102 + t4 * t28 + m(5) * (t11 * t5 + t12 * t4 + t2 * t27 + t26 * t3 + t46 * t57 + t51 * t53) + (t24 * t101 + t87) * Ifges(5,2) + (t25 * t99 - t86) * Ifges(5,1) + ((t66 * t49 - t36 / 0.2e1 - t82 / 0.2e1) * t65 + (-t66 * t50 - t37 / 0.2e1 - t83 / 0.2e1) * t63) * qJD(3) + (-t20 * t26 + t21 * t27 - t69) * mrSges(5,3) + (t24 * t100 + t38 * t102 + t20 * t42 - t21 * t71) * Ifges(5,4) + (((2 * mrSges(3,3)) + t73 + (2 * t103)) * qJD(2) + (0.2e1 * t70 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1)) * t63 * t65 + (0.3e1 / 0.2e1 * t63 ^ 2 - 0.3e1 / 0.2e1 * t65 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1); m(5) * t69 + t24 * t28 + t25 * t29 + t72 * qJD(3) + (t86 - t87) * mrSges(5,3) + (-t44 + t75 + (-mrSges(3,3) - t103) * qJD(1)) * qJD(1); (-t73 * qJD(3) - t72) * t52 - m(5) * (t11 * t13 + t12 * t14) + t68 - t13 * t29 - t14 * t28 + (t37 * t98 + t65 * t36 / 0.2e1 + ((-Ifges(4,1) * t63 - t92) * t97 + (-Ifges(4,2) * t65 - t93) * t98 - t70) * qJD(1) + t75 * t65 * pkin(3) + (-Ifges(4,5) * t63 / 0.2e1 + Ifges(4,6) * t97) * qJD(3)) * qJD(1) + (t28 * t76 - t29 * t77 + m(5) * (-t11 * t77 + t12 * t76 + t2 * t62 + t3 * t64) + (-t20 * t64 + t21 * t62) * mrSges(5,3)) * pkin(3) + t12 * t88; t68 + (t29 + t88) * t12 - t11 * t28;];
tauc = t1(:);
