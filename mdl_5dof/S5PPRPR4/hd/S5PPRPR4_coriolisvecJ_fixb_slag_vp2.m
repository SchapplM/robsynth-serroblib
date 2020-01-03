% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:18
% DurationCPUTime: 0.79s
% Computational Cost: add. (577->126), mult. (1510->195), div. (0->0), fcn. (997->6), ass. (0->65)
t53 = cos(qJ(3));
t64 = qJD(2) * t53;
t60 = qJD(4) - t64;
t48 = sin(pkin(8));
t69 = pkin(6) + qJ(4);
t39 = t69 * t48;
t49 = cos(pkin(8));
t40 = t69 * t49;
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t13 = -t39 * t52 - t40 * t50;
t34 = t48 * t50 - t52 * t49;
t82 = t13 * qJD(5) - t60 * t34;
t14 = -t39 * t50 + t40 * t52;
t35 = t48 * t52 + t49 * t50;
t56 = t35 * t53;
t81 = qJD(2) * t56 - t35 * qJD(4) - t14 * qJD(5);
t51 = sin(qJ(3));
t25 = t34 * t51;
t67 = t48 ^ 2 + t49 ^ 2;
t80 = t67 * mrSges(5,3);
t79 = t35 * qJD(5);
t78 = m(5) / 0.2e1;
t31 = t34 * qJD(5);
t76 = -t31 / 0.2e1;
t75 = -t79 / 0.2e1;
t30 = t35 * qJD(3);
t74 = Ifges(6,4) * t30;
t29 = t34 * qJD(3);
t22 = qJD(5) * t29;
t73 = t34 * t22;
t23 = qJD(3) * t79;
t72 = t35 * t23;
t59 = -mrSges(5,1) * t49 + mrSges(5,2) * t48;
t68 = mrSges(6,1) * t29 + mrSges(6,2) * t30 + qJD(3) * t59;
t66 = pkin(6) * qJD(3);
t65 = qJD(1) * t49;
t63 = t51 * qJD(2);
t45 = -pkin(4) * t49 - pkin(3);
t38 = (qJD(4) + t64) * qJD(3);
t62 = t67 * t38;
t61 = qJD(3) * t63;
t9 = t23 * mrSges(6,1) - t22 * mrSges(6,2);
t42 = qJD(3) * qJ(4) + t63;
t27 = -qJD(1) * t48 + t49 * t42;
t15 = -t65 + (-t42 - t66) * t48;
t16 = t49 * t66 + t27;
t3 = t15 * t52 - t16 * t50;
t4 = t15 * t50 + t16 * t52;
t58 = -(-t42 * t48 - t65) * t48 + t27 * t49;
t57 = t58 * t53;
t54 = qJD(3) ^ 2;
t41 = -qJD(3) * pkin(3) + t60;
t33 = t45 * qJD(3) + t60;
t28 = Ifges(6,4) * t29;
t24 = t35 * t51;
t18 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t30;
t17 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t29;
t11 = Ifges(6,1) * t30 + Ifges(6,5) * qJD(5) - t28;
t10 = -Ifges(6,2) * t29 + Ifges(6,6) * qJD(5) + t74;
t8 = -qJD(3) * t56 + qJD(5) * t25;
t7 = -t53 * t29 - t51 * t79;
t2 = -t4 * qJD(5) - t35 * t38;
t1 = t3 * qJD(5) - t34 * t38;
t5 = [t31 * t17 + t79 * t18 + m(6) * (-t1 * t35 + t2 * t34 + t3 * t79 + t31 * t4) + (t72 + t73) * mrSges(6,3); t7 * t17 + t8 * t18 - t53 * t9 + t68 * t51 * qJD(3) + (-t22 * t24 + t23 * t25) * mrSges(6,3) + m(6) * (-t1 * t25 - t2 * t24 + t3 * t8 + t4 * t7) + m(5) * t51 * t62 + 0.2e1 * (t57 * t78 + ((t41 - t64) * t78 + m(6) * (t33 - t64) / 0.2e1) * t51) * qJD(3) + (-t51 * mrSges(4,1) + (-mrSges(4,2) + t80) * t53) * t54; t45 * t9 + t11 * t76 + qJD(5) * (-Ifges(6,5) * t31 - Ifges(6,6) * t79) / 0.2e1 + t10 * t75 + t33 * (mrSges(6,1) * t79 - t31 * mrSges(6,2)) + t81 * t18 + t82 * t17 + (t34 * t23 - t29 * t75) * Ifges(6,2) + (-t22 * t35 + t30 * t76) * Ifges(6,1) + ((mrSges(6,1) * t34 + mrSges(6,2) * t35 + t59) * qJD(3) - t68) * t63 + (-t1 * t34 + t13 * t22 - t14 * t23 - t2 * t35 + t3 * t31 - t4 * t79) * mrSges(6,3) + (qJD(3) * t60 * t67 + t62) * mrSges(5,3) + (-t29 * t76 + t30 * t75 - t72 + t73) * Ifges(6,4) + (t1 * t14 + t13 * t2 + t81 * t3 - t33 * t63 + t82 * t4 + t45 * t61) * m(6) + (-pkin(3) * t61 + qJ(4) * t62 + t58 * qJD(4) - (t41 * t51 + t57) * qJD(2)) * m(5); -m(6) * (-t29 * t4 - t3 * t30) + t29 * t17 + t30 * t18 - t54 * t80 + (m(6) * t63 + (-t58 + t63) * m(5)) * qJD(3) + t9; -Ifges(6,5) * t22 - Ifges(6,6) * t23 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t33 * (mrSges(6,1) * t30 - mrSges(6,2) * t29) - t30 * (-Ifges(6,1) * t29 - t74) / 0.2e1 + t30 * t10 / 0.2e1 - qJD(5) * (-Ifges(6,5) * t29 - Ifges(6,6) * t30) / 0.2e1 - t3 * t17 + t4 * t18 + (-t29 * t3 + t30 * t4) * mrSges(6,3) + (-Ifges(6,2) * t30 + t11 - t28) * t29 / 0.2e1;];
tauc = t5(:);
