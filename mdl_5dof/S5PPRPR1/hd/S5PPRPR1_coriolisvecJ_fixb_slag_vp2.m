% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:00:56
% DurationCPUTime: 0.88s
% Computational Cost: add. (766->137), mult. (2108->201), div. (0->0), fcn. (1550->8), ass. (0->74)
t60 = sin(pkin(8));
t64 = sin(qJ(3));
t79 = t60 * t64;
t53 = qJD(1) * t79;
t62 = cos(pkin(8));
t66 = cos(qJ(3));
t78 = t62 * t66;
t73 = qJD(1) * t78;
t38 = -t53 + t73;
t68 = qJD(4) - t38;
t59 = sin(pkin(9));
t77 = pkin(6) + qJ(4);
t50 = t77 * t59;
t61 = cos(pkin(9));
t51 = t77 * t61;
t63 = sin(qJ(5));
t65 = cos(qJ(5));
t23 = -t50 * t63 + t51 * t65;
t47 = t59 * t65 + t61 * t63;
t88 = -t23 * qJD(5) - t68 * t47;
t22 = -t50 * t65 - t51 * t63;
t69 = t59 * t63 - t61 * t65;
t87 = t22 * qJD(5) - t68 * t69;
t75 = t59 ^ 2 + t61 ^ 2;
t86 = t75 * mrSges(5,3);
t41 = t69 * qJD(5);
t84 = -t41 / 0.2e1;
t42 = t47 * qJD(5);
t83 = -t42 / 0.2e1;
t39 = t47 * qJD(3);
t82 = Ifges(6,4) * t39;
t48 = t60 * t66 + t62 * t64;
t44 = t48 * qJD(3);
t35 = qJD(1) * t44;
t46 = -t78 + t79;
t19 = t35 * t46;
t32 = qJD(3) * t41;
t81 = t69 * t32;
t33 = qJD(3) * t42;
t80 = t47 * t33;
t37 = t69 * qJD(3);
t71 = -t61 * mrSges(5,1) + t59 * mrSges(5,2);
t76 = mrSges(6,1) * t37 + mrSges(6,2) * t39 + t71 * qJD(3);
t40 = t48 * qJD(1);
t31 = qJD(3) * qJ(4) + t40;
t21 = t59 * qJD(2) + t61 * t31;
t74 = pkin(6) * qJD(3);
t54 = -pkin(4) * t61 - pkin(3);
t52 = qJD(3) * t73;
t26 = t52 + (qJD(4) - t53) * qJD(3);
t72 = t75 * t26;
t11 = t33 * mrSges(6,1) - t32 * mrSges(6,2);
t56 = t61 * qJD(2);
t17 = t56 + (-t31 - t74) * t59;
t18 = t61 * t74 + t21;
t5 = t17 * t65 - t18 * t63;
t6 = t17 * t63 + t18 * t65;
t70 = -(-t31 * t59 + t56) * t59 + t21 * t61;
t43 = t46 * qJD(3);
t36 = Ifges(6,4) * t37;
t34 = -qJD(3) * t53 + t52;
t30 = -qJD(3) * pkin(3) + t68;
t28 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t39;
t27 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t37;
t24 = t54 * qJD(3) + t68;
t15 = Ifges(6,1) * t39 + Ifges(6,5) * qJD(5) - t36;
t14 = -Ifges(6,2) * t37 + Ifges(6,6) * qJD(5) + t82;
t13 = t69 * t48;
t12 = t47 * t48;
t4 = t48 * t41 + t47 * t43;
t3 = -t48 * t42 + t69 * t43;
t2 = -t6 * qJD(5) - t47 * t26;
t1 = t5 * qJD(5) - t69 * t26;
t7 = [t46 * t11 + t3 * t27 + t4 * t28 + t76 * t44 + (-t12 * t32 + t13 * t33) * mrSges(6,3) + m(4) * (t34 * t48 - t38 * t44 - t40 * t43 + t19) + m(5) * (t30 * t44 - t70 * t43 + t48 * t72 + t19) + m(6) * (-t1 * t13 - t12 * t2 + t24 * t44 + t3 * t6 + t4 * t5 + t19) + (-t44 * mrSges(4,1) - (-mrSges(4,2) + t86) * t43) * qJD(3); m(6) * (t1 * t47 - t2 * t69 - t41 * t6 - t42 * t5) - t41 * t27 - t42 * t28 + (-t80 - t81) * mrSges(6,3); t54 * t11 + t15 * t84 + t24 * (t42 * mrSges(6,1) - t41 * mrSges(6,2)) + qJD(5) * (-Ifges(6,5) * t41 - Ifges(6,6) * t42) / 0.2e1 + t14 * t83 + t88 * t28 + t87 * t27 + (t33 * t69 - t37 * t83) * Ifges(6,2) + (-t32 * t47 + t39 * t84) * Ifges(6,1) + (t38 * qJD(3) - t34) * mrSges(4,2) + (qJD(3) * mrSges(4,1) - t76) * t40 + (mrSges(6,1) * t69 + t47 * mrSges(6,2) - mrSges(4,1) + t71) * t35 + (-t1 * t69 - t2 * t47 + t22 * t32 - t23 * t33 + t41 * t5 - t42 * t6) * mrSges(6,3) + (t68 * qJD(3) * t75 + t72) * mrSges(5,3) + (-t37 * t84 + t39 * t83 - t80 + t81) * Ifges(6,4) + (t1 * t23 + t2 * t22 - t24 * t40 + t35 * t54 + t88 * t5 + t87 * t6) * m(6) + (-pkin(3) * t35 + qJ(4) * t72 - t30 * t40 + t68 * t70) * m(5); -qJD(3) ^ 2 * t86 + t37 * t27 + t39 * t28 + t11 + (t6 * t37 + t5 * t39 + t35) * m(6) + (-t70 * qJD(3) + t35) * m(5); -Ifges(6,5) * t32 - Ifges(6,6) * t33 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t24 * (mrSges(6,1) * t39 - mrSges(6,2) * t37) - t39 * (-Ifges(6,1) * t37 - t82) / 0.2e1 + t39 * t14 / 0.2e1 - qJD(5) * (-Ifges(6,5) * t37 - Ifges(6,6) * t39) / 0.2e1 - t5 * t27 + t6 * t28 + (-t37 * t5 + t39 * t6) * mrSges(6,3) + (-Ifges(6,2) * t39 + t15 - t36) * t37 / 0.2e1;];
tauc = t7(:);
