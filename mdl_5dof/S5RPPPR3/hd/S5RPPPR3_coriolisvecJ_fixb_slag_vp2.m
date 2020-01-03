% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:53
% DurationCPUTime: 0.73s
% Computational Cost: add. (571->124), mult. (1411->183), div. (0->0), fcn. (870->6), ass. (0->61)
t48 = sin(pkin(8));
t50 = cos(pkin(8));
t85 = t48 ^ 2 + t50 ^ 2;
t84 = mrSges(5,2) + mrSges(4,3);
t83 = t85 * qJD(1) * qJD(3);
t52 = sin(qJ(5));
t53 = cos(qJ(5));
t68 = t48 * qJD(1);
t71 = qJD(1) * t50;
t29 = -t52 * t68 - t53 * t71;
t82 = -t29 / 0.2e1;
t59 = t48 * t52 + t50 * t53;
t27 = t59 * qJD(5);
t24 = qJD(1) * t27;
t81 = t59 * t24;
t36 = t48 * t53 - t50 * t52;
t28 = t36 * qJD(5);
t25 = qJD(1) * t28;
t80 = t36 * t25;
t58 = cos(pkin(7)) * pkin(1) + t48 * qJ(4) + pkin(2);
t78 = (t50 * pkin(3) + t58) * qJD(1);
t30 = -t52 * t71 + t53 * t68;
t76 = t30 / 0.2e1;
t44 = sin(pkin(7)) * pkin(1) + qJ(3);
t75 = -pkin(6) + t44;
t74 = Ifges(6,4) * t30;
t40 = t44 * qJD(1);
t23 = t48 * qJD(2) + t50 * t40;
t21 = (pkin(3) + pkin(4)) * t50 + t58;
t72 = qJD(1) * t21;
t70 = qJD(3) * t48;
t66 = t23 * t50 * qJD(3) + t83 * t44;
t63 = qJD(4) * t68;
t22 = t50 * qJD(2) - t48 * t40;
t60 = t25 * mrSges(6,1) - t24 * mrSges(6,2);
t20 = qJD(4) - t22;
t13 = -pkin(6) * t68 + t20;
t14 = -pkin(6) * t71 + t23;
t3 = t53 * t13 - t52 * t14;
t4 = t52 * t13 + t53 * t14;
t31 = t75 * t48;
t32 = t75 * t50;
t9 = t53 * t31 - t52 * t32;
t10 = t52 * t31 + t53 * t32;
t37 = (-t50 * mrSges(5,1) - t48 * mrSges(5,3)) * qJD(1);
t57 = t36 * qJD(3);
t56 = t59 * qJD(3);
t26 = Ifges(6,4) * t29;
t19 = qJD(5) * mrSges(6,1) - t30 * mrSges(6,3);
t18 = -qJD(5) * mrSges(6,2) + t29 * mrSges(6,3);
t17 = qJD(3) - t78;
t16 = t23 * t71;
t12 = -qJD(3) + t72;
t11 = -t29 * mrSges(6,1) + t30 * mrSges(6,2);
t8 = Ifges(6,1) * t30 + Ifges(6,5) * qJD(5) + t26;
t7 = Ifges(6,2) * t29 + Ifges(6,6) * qJD(5) + t74;
t6 = -t10 * qJD(5) + t57;
t5 = t9 * qJD(5) + t56;
t2 = qJD(1) * t57 - t4 * qJD(5);
t1 = qJD(1) * t56 + t3 * qJD(5);
t15 = [(mrSges(6,1) * t59 + t36 * mrSges(6,2)) * t63 - t27 * t8 / 0.2e1 + t12 * (t28 * mrSges(6,1) - t27 * mrSges(6,2)) + qJD(5) * (-Ifges(6,5) * t27 - Ifges(6,6) * t28) / 0.2e1 - t28 * t7 / 0.2e1 + t5 * t18 + t6 * t19 + t21 * t60 + m(4) * (-t22 * t70 + t66) + m(5) * (t20 * t70 + t66) + m(6) * (t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5) + (t11 - 0.2e1 * t37 + m(5) * (-t17 + t78) + m(6) * (t12 + t72)) * qJD(4) * t48 + (t25 * t59 + t28 * t82) * Ifges(6,2) + (-t24 * t36 - t27 * t76) * Ifges(6,1) + (-t1 * t59 - t10 * t25 - t2 * t36 + t9 * t24 + t3 * t27 - t4 * t28) * mrSges(6,3) + (t27 * t82 - t28 * t76 - t80 + t81) * Ifges(6,4) + 0.2e1 * t84 * t83; -t27 * t18 - t28 * t19 + m(6) * (t1 * t36 - t2 * t59 - t4 * t27 - t3 * t28) + (-t80 - t81) * mrSges(6,3); t29 * t18 - t30 * t19 + (-m(5) - m(6)) * t63 - m(4) * (-t22 * t68 + t16) - m(5) * (t20 * t68 + t16) - m(6) * (-t4 * t29 + t3 * t30) - t60 - t84 * qJD(1) ^ 2 * t85; m(6) * (t1 * t52 + t2 * t53) + (t53 * t24 - t52 * t25) * mrSges(6,3) + (-m(6) * t12 - t11 + t37 + (qJD(3) + t17) * m(5)) * t68 + (m(6) * (-t3 * t52 + t4 * t53) + t53 * t18 - t52 * t19) * qJD(5); -Ifges(6,5) * t24 - Ifges(6,6) * t25 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t12 * (t30 * mrSges(6,1) + t29 * mrSges(6,2)) - t30 * (Ifges(6,1) * t29 - t74) / 0.2e1 + t7 * t76 - qJD(5) * (Ifges(6,5) * t29 - Ifges(6,6) * t30) / 0.2e1 - t3 * t18 + t4 * t19 + (t3 * t29 + t4 * t30) * mrSges(6,3) + (-Ifges(6,2) * t30 + t26 + t8) * t82;];
tauc = t15(:);
