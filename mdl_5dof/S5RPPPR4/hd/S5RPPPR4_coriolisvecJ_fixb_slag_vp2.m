% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:08
% DurationCPUTime: 0.58s
% Computational Cost: add. (658->116), mult. (1409->169), div. (0->0), fcn. (842->6), ass. (0->60)
t41 = sin(pkin(7)) * pkin(1) + qJ(3);
t82 = qJD(1) * t41;
t83 = m(4) * t82;
t49 = sin(pkin(8));
t51 = cos(pkin(8));
t53 = sin(qJ(5));
t54 = cos(qJ(5));
t81 = -t53 * t49 + t54 * t51;
t29 = t81 * qJD(5);
t66 = t49 ^ 2 + t51 ^ 2;
t80 = qJD(1) * t66;
t78 = -t29 / 0.2e1;
t32 = -t54 * t49 - t53 * t51;
t30 = t32 * qJD(5);
t77 = t30 / 0.2e1;
t40 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t75 = -pkin(6) + t40;
t63 = qJD(1) * t51;
t64 = qJD(1) * t49;
t28 = -t53 * t64 + t54 * t63;
t74 = Ifges(6,4) * t28;
t19 = qJD(1) * t30;
t73 = t32 * t19;
t20 = qJD(1) * t29;
t72 = t32 * t20;
t71 = t81 * t19;
t70 = t81 * t20;
t27 = t32 * qJD(1);
t61 = t49 * mrSges(5,1) + t51 * mrSges(5,2);
t67 = t27 * mrSges(6,1) - t28 * mrSges(6,2) - qJD(1) * t61;
t31 = qJD(1) * t40 + qJD(3);
t15 = t51 * qJD(2) + t49 * t31;
t35 = qJD(4) + t82;
t62 = t20 * mrSges(6,1) + t19 * mrSges(6,2);
t14 = -t49 * qJD(2) + t51 * t31;
t12 = -pkin(6) * t63 + t14;
t13 = -pkin(6) * t64 + t15;
t3 = t54 * t12 - t53 * t13;
t4 = t53 * t12 + t54 * t13;
t60 = -t14 * t51 - t15 * t49;
t25 = t75 * t49;
t26 = t75 * t51;
t8 = t54 * t25 + t53 * t26;
t7 = -t53 * t25 + t54 * t26;
t58 = t32 * qJD(4);
t57 = t81 * qJD(4);
t1 = qJD(1) * t58 + t3 * qJD(5);
t2 = -qJD(1) * t57 - t4 * qJD(5);
t56 = -t1 * t32 + t2 * t81 + t4 * t29 + t3 * t30;
t55 = qJD(1) ^ 2;
t36 = t49 * pkin(4) + t41;
t24 = pkin(4) * t64 + t35;
t23 = Ifges(6,4) * t27;
t17 = qJD(5) * mrSges(6,1) - t28 * mrSges(6,3);
t16 = -qJD(5) * mrSges(6,2) + t27 * mrSges(6,3);
t10 = Ifges(6,1) * t28 + Ifges(6,5) * qJD(5) + t23;
t9 = Ifges(6,2) * t27 + Ifges(6,6) * qJD(5) + t74;
t6 = -t8 * qJD(5) - t57;
t5 = t7 * qJD(5) + t58;
t11 = [t36 * t62 + t9 * t78 + t24 * (t29 * mrSges(6,1) + t30 * mrSges(6,2)) + qJD(5) * (Ifges(6,5) * t30 - Ifges(6,6) * t29) / 0.2e1 + t10 * t77 + t5 * t16 + t6 * t17 + m(6) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5) + (t27 * t78 - t72) * Ifges(6,2) + (t28 * t77 + t71) * Ifges(6,1) + (m(5) * (-t40 * t80 + t60) + 0.2e1 * mrSges(5,3) * t80) * qJD(4) + ((-t32 * mrSges(6,1) + mrSges(6,2) * t81 + (2 * mrSges(4,3)) + t61) * qJD(1) + m(5) * (t35 + t82) + m(6) * (qJD(1) * t36 + t24) + 0.2e1 * t83 - t67) * qJD(3) + (-t7 * t19 - t8 * t20 - t56) * mrSges(6,3) + (t27 * t77 + t28 * t78 - t70 + t73) * Ifges(6,4); m(6) * (t1 * t81 + t2 * t32 - t3 * t29 + t4 * t30) + t30 * t16 - t29 * t17 + (-t70 - t73) * mrSges(6,3); m(6) * t56 + t29 * t16 + t30 * t17 - t55 * mrSges(4,3) + (-t71 + t72) * mrSges(6,3) + (-t83 - m(6) * t24 + (-t66 * qJD(4) - t35) * m(5) + t67) * qJD(1); -m(6) * (t4 * t27 - t3 * t28) - t27 * t16 + t28 * t17 - t66 * t55 * mrSges(5,3) + (m(6) * qJD(3) + (qJD(3) - t60) * m(5)) * qJD(1) + t62; Ifges(6,5) * t19 - Ifges(6,6) * t20 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t24 * (t28 * mrSges(6,1) + t27 * mrSges(6,2)) - t28 * (Ifges(6,1) * t27 - t74) / 0.2e1 + t28 * t9 / 0.2e1 - qJD(5) * (Ifges(6,5) * t27 - Ifges(6,6) * t28) / 0.2e1 - t3 * t16 + t4 * t17 + (t3 * t27 + t4 * t28) * mrSges(6,3) - (-Ifges(6,2) * t28 + t10 + t23) * t27 / 0.2e1;];
tauc = t11(:);
