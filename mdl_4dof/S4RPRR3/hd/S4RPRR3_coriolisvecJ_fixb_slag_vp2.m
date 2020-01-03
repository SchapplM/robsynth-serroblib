% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR3
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:07
% DurationCPUTime: 0.77s
% Computational Cost: add. (889->153), mult. (2236->232), div. (0->0), fcn. (1305->6), ass. (0->84)
t60 = qJD(3) + qJD(4);
t66 = cos(qJ(3));
t56 = sin(pkin(7)) * pkin(1) + pkin(5);
t52 = t56 * qJD(1);
t70 = pkin(6) * qJD(1) + t52;
t64 = sin(qJ(3));
t76 = t64 * qJD(2);
t32 = t70 * t66 + t76;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t49 = -t63 * t64 + t65 * t66;
t44 = t49 * qJD(1);
t96 = t44 / 0.2e1;
t50 = t63 * t66 + t65 * t64;
t45 = t50 * qJD(1);
t95 = t45 / 0.2e1;
t75 = -cos(pkin(7)) * pkin(1) - pkin(2);
t51 = -pkin(3) * t66 + t75;
t46 = qJD(1) * t51;
t94 = m(5) * t46;
t93 = pkin(6) + t56;
t92 = mrSges(5,3) * t44;
t91 = Ifges(4,4) * t64;
t90 = Ifges(5,4) * t45;
t89 = t32 * t63;
t88 = t32 * t65;
t87 = t45 * mrSges(5,3);
t26 = t60 * t49;
t22 = t26 * qJD(1);
t86 = t49 * t22;
t27 = t60 * t50;
t23 = t27 * qJD(1);
t85 = t50 * t23;
t84 = t52 * t64;
t54 = t75 * qJD(1);
t83 = t54 * mrSges(4,2);
t82 = Ifges(4,5) * qJD(3);
t81 = Ifges(4,6) * qJD(3);
t80 = qJD(1) * t64;
t79 = qJD(1) * t66;
t78 = qJD(4) * t63;
t77 = qJD(4) * t65;
t59 = t66 * qJD(2);
t58 = Ifges(4,4) * t79;
t74 = m(4) * t56 + mrSges(4,3);
t73 = t82 / 0.2e1;
t72 = -t81 / 0.2e1;
t71 = qJD(3) * t93;
t69 = t81 / 0.2e1 + (Ifges(4,2) * t66 + t91) * qJD(1) / 0.2e1 - t54 * mrSges(4,1);
t68 = t70 * t64;
t31 = t59 - t68;
t30 = qJD(3) * pkin(3) + t31;
t7 = t30 * t65 - t89;
t8 = t30 * t63 + t88;
t47 = t93 * t64;
t48 = t93 * t66;
t20 = -t47 * t65 - t48 * t63;
t21 = -t47 * t63 + t48 * t65;
t38 = t52 * t66 + t76;
t16 = Ifges(5,2) * t44 + Ifges(5,6) * t60 + t90;
t41 = Ifges(5,4) * t44;
t17 = Ifges(5,1) * t45 + Ifges(5,5) * t60 + t41;
t57 = qJD(3) * t59;
t28 = -qJD(3) * t68 + t57;
t29 = t32 * qJD(3);
t2 = qJD(4) * t7 + t28 * t65 - t29 * t63;
t3 = -qJD(4) * t8 - t28 * t63 - t29 * t65;
t67 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + t16 * t95 - t46 * (mrSges(5,1) * t45 + mrSges(5,2) * t44) + t7 * t92 - t45 * (Ifges(5,1) * t44 - t90) / 0.2e1 - t60 * (Ifges(5,5) * t44 - Ifges(5,6) * t45) / 0.2e1 - Ifges(5,6) * t23 + Ifges(5,5) * t22 - (-Ifges(5,2) * t45 + t17 + t41) * t44 / 0.2e1;
t55 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t79;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t80;
t43 = Ifges(4,1) * t80 + t58 + t82;
t40 = t66 * t71;
t39 = t64 * t71;
t37 = t59 - t84;
t36 = t38 * qJD(3);
t35 = -qJD(3) * t84 + t57;
t34 = mrSges(5,1) * t60 - t87;
t33 = -mrSges(5,2) * t60 + t92;
t25 = -mrSges(5,1) * t44 + mrSges(5,2) * t45;
t12 = t31 * t65 - t89;
t11 = -t31 * t63 - t88;
t5 = -t21 * qJD(4) + t39 * t63 - t40 * t65;
t4 = t20 * qJD(4) - t39 * t65 - t40 * t63;
t1 = [m(5) * (t2 * t21 + t20 * t3 + t4 * t8 + t5 * t7) + t60 * (Ifges(5,5) * t26 - Ifges(5,6) * t27) / 0.2e1 + t51 * (t23 * mrSges(5,1) + t22 * mrSges(5,2)) + t46 * (t27 * mrSges(5,1) + t26 * mrSges(5,2)) - t27 * t16 / 0.2e1 + t4 * t33 + t5 * t34 + t26 * t17 / 0.2e1 + (-t23 * t49 - t27 * t96) * Ifges(5,2) + (t50 * t22 + t26 * t95) * Ifges(5,1) + (t2 * t49 - t20 * t22 - t21 * t23 - t26 * t7 - t27 * t8 - t3 * t50) * mrSges(5,3) + (t74 * t35 + (t43 / 0.2e1 - t56 * t53 + 0.3e1 / 0.2e1 * t58 + t73 - t74 * t37 + 0.2e1 * t83) * qJD(3)) * t66 + (t26 * t96 - t27 * t95 - t85 + t86) * Ifges(5,4) + (t74 * t36 + (-t56 * t55 + t72 - t74 * t38 + (t75 * mrSges(4,1) - 0.3e1 / 0.2e1 * t91 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t66) * qJD(1) + (0.2e1 * t94 + t25 + qJD(1) * (-mrSges(5,1) * t49 + mrSges(5,2) * t50)) * pkin(3) - t69) * qJD(3)) * t64; t26 * t33 - t27 * t34 + (-t85 - t86) * mrSges(5,3) + m(4) * (t35 * t64 - t36 * t66) + m(5) * (t2 * t50 + t26 * t8 - t27 * t7 + t3 * t49) + (m(4) * (-t37 * t64 + t38 * t66) + t66 * t55 - t64 * t53 + (-t64 ^ 2 - t66 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3); -m(5) * (t11 * t7 + t12 * t8) + (m(5) * (t2 * t63 + t3 * t65 - t7 * t78 + t8 * t77) + t33 * t77 - t34 * t78 + (-t22 * t65 - t23 * t63) * mrSges(5,3)) * pkin(3) + ((-t83 + t73 - t43 / 0.2e1 - t58 / 0.2e1 + t37 * mrSges(4,3)) * t66 + (t72 + t38 * mrSges(4,3) + (t91 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t66) * qJD(1) + (-t25 - t94) * pkin(3) + t69) * t64) * qJD(1) + t38 * t53 - t37 * t55 - t35 * mrSges(4,2) - t36 * mrSges(4,1) - t12 * t33 - t11 * t34 + t8 * t87 + t67; (t34 + t87) * t8 - t33 * t7 + t67;];
tauc = t1(:);
