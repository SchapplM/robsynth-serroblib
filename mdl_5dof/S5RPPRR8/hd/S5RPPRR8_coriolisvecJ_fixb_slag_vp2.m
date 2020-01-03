% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:06
% DurationCPUTime: 0.88s
% Computational Cost: add. (1434->146), mult. (2489->221), div. (0->0), fcn. (1234->6), ass. (0->78)
t104 = qJD(1) - qJD(4);
t48 = sin(pkin(8));
t49 = cos(pkin(8));
t103 = (m(3) * qJ(2) + t48 * mrSges(4,1) + t49 * mrSges(4,2) + mrSges(3,3)) * qJD(1);
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t31 = t48 * t53 + t49 * t51;
t78 = t104 * t31;
t30 = t48 * t51 - t49 * t53;
t102 = t104 * t30;
t54 = -pkin(1) - pkin(2);
t75 = qJ(2) * t48;
t66 = t49 * t54 - t75;
t32 = -pkin(3) + t66;
t39 = qJ(2) * t49 + t48 * t54;
t76 = t51 * t32 + t53 * t39;
t44 = t54 * qJD(1) + qJD(2);
t40 = t49 * t44;
t15 = t40 + (-pkin(3) - t75) * qJD(1);
t70 = qJD(1) * qJ(2);
t22 = t44 * t48 + t49 * t70;
t11 = t15 * t53 - t22 * t51;
t57 = t30 * qJD(2);
t3 = -qJD(1) * t57 + t11 * qJD(4);
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t12 = t15 * t51 + t22 * t53;
t8 = -pkin(7) * t104 + t12;
t5 = qJD(3) * t52 - t50 * t8;
t1 = qJD(5) * t5 + t3 * t52;
t6 = qJD(3) * t50 + t52 * t8;
t2 = -qJD(5) * t6 - t3 * t50;
t101 = t1 * t52 - t2 * t50;
t84 = t50 * mrSges(6,3);
t37 = qJD(5) * mrSges(6,1) + t104 * t84;
t81 = t52 * mrSges(6,3);
t38 = -qJD(5) * mrSges(6,2) - t104 * t81;
t61 = -t5 * t50 + t52 * t6;
t99 = m(6) * t61 - t50 * t37 + t52 * t38;
t58 = t31 * qJD(2);
t4 = qJD(1) * t58 + t12 * qJD(4);
t95 = t4 * t30;
t94 = t5 * mrSges(6,3);
t93 = t6 * mrSges(6,3);
t91 = Ifges(6,1) * t50;
t90 = Ifges(6,4) * t50;
t89 = Ifges(6,4) * t52;
t88 = t104 * (Ifges(6,2) * t52 + t90);
t87 = t104 * (t89 + t91);
t36 = (Ifges(6,1) * t52 - t90) * qJD(5);
t83 = t50 * t36;
t35 = (-Ifges(6,2) * t50 + t89) * qJD(5);
t80 = t52 * t35;
t74 = Ifges(6,5) * qJD(5);
t73 = Ifges(6,6) * qJD(5);
t72 = qJD(5) * t50;
t71 = qJD(5) * t52;
t16 = t73 - t88;
t65 = -t16 / 0.2e1 - t93;
t45 = t104 * t89;
t17 = -t104 * t91 - t45 + t74;
t64 = t17 / 0.2e1 - t94;
t33 = (mrSges(6,1) * t50 + mrSges(6,2) * t52) * qJD(5);
t7 = pkin(4) * t104 - t11;
t63 = -t3 * mrSges(5,2) + t7 * t33;
t62 = -t5 * t52 - t6 * t50;
t41 = -mrSges(6,1) * t52 + mrSges(6,2) * t50;
t60 = -(-t48 * t70 + t40) * t48 + t22 * t49;
t59 = t32 * t53 - t39 * t51;
t56 = t62 * qJD(5) + t101;
t34 = qJD(5) * (Ifges(6,5) * t52 - Ifges(6,6) * t50);
t23 = t41 * t104;
t20 = t104 * t36;
t19 = t104 * t35;
t18 = t104 * t33;
t13 = pkin(4) - t59;
t10 = t76 * qJD(4) + t58;
t9 = [-t1 * t81 + m(5) * (-t11 * t10 + t3 * t76 - t4 * t59) - t63 + m(6) * (t10 * t7 + t13 * t4) + t50 * t20 / 0.2e1 + t52 * t19 / 0.2e1 - qJD(5) * t34 / 0.2e1 - t4 * t41 - t13 * t18 - t10 * t23 + t2 * t84 + t72 * t93 + t71 * t94 + (t83 + t80) * t104 / 0.2e1 + (t16 - t88) * t72 / 0.2e1 - (t17 - t87) * t71 / 0.2e1 + (m(6) * t56 - t71 * t37 - t72 * t38) * (-pkin(7) + t76) + (t10 * t104 + t4) * mrSges(5,1) + (m(5) * t12 + mrSges(5,2) * t104 + t99) * (t59 * qJD(4) - t57) + (m(4) * ((t49 * t39 - t48 * t66) * qJD(1) + t60) + 0.2e1 * t103) * qJD(2); -t30 * t18 + t78 * t23 + (t102 * t52 - t31 * t72) * t38 + (-t102 * t50 - t31 * t71) * t37 - (t78 * mrSges(5,1) - mrSges(5,2) * t102) * t104 + (-m(4) * t60 - t103) * qJD(1) + (t102 * t61 + t56 * t31 - t78 * t7 + t95) * m(6) + (t102 * t12 + t78 * t11 + t3 * t31 + t95) * m(5); m(6) * (t1 * t50 + t2 * t52) + (-(-t50 ^ 2 - t52 ^ 2) * t104 * mrSges(6,3) + t99) * qJD(5); pkin(4) * t18 + t12 * t23 + (t41 - mrSges(5,1)) * t4 + (-t11 * t38 + t1 * mrSges(6,3) - t19 / 0.2e1) * t52 + (t11 * t37 - t2 * mrSges(6,3) - t20 / 0.2e1) * t50 - m(6) * (t61 * t11 + t7 * t12) + m(6) * (-t4 * pkin(4) + t101 * pkin(7)) - (t80 / 0.2e1 + t83 / 0.2e1 + t11 * mrSges(5,2) + t12 * mrSges(5,1)) * t104 + (t34 / 0.2e1 + (-t87 / 0.2e1 + t64) * t52 + (t88 / 0.2e1 + t65) * t50 + (m(6) * t62 - t52 * t37 - t50 * t38) * pkin(7)) * qJD(5) + t63; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t6 * t37 - t5 * t38 - ((t74 / 0.2e1 - t7 * mrSges(6,2) + t45 / 0.2e1 - t64) * t52 + (-t73 / 0.2e1 - t7 * mrSges(6,1) - (t90 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t52) * t104 - t65) * t50) * t104;];
tauc = t9(:);
