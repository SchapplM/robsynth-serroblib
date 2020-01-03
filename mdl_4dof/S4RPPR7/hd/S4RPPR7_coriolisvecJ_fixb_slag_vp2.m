% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:37
% DurationCPUTime: 0.46s
% Computational Cost: add. (454->103), mult. (1059->155), div. (0->0), fcn. (608->4), ass. (0->55)
t40 = sin(pkin(6));
t43 = sin(qJ(4));
t41 = cos(pkin(6));
t44 = cos(qJ(4));
t56 = t44 * t41;
t68 = -t43 * t40 + t56;
t22 = t68 * qJD(4);
t42 = -pkin(1) - qJ(3);
t30 = qJD(1) * t42 + qJD(2);
t55 = t40 ^ 2 + t41 ^ 2;
t69 = t30 * t55;
t67 = qJD(3) * t55;
t65 = -t22 / 0.2e1;
t24 = -t44 * t40 - t43 * t41;
t23 = t24 * qJD(4);
t64 = t23 / 0.2e1;
t63 = -pkin(5) + t42;
t20 = t24 * qJD(1);
t54 = qJD(1) * t40;
t21 = qJD(1) * t56 - t43 * t54;
t50 = t40 * mrSges(4,1) + t41 * mrSges(4,2);
t62 = -t20 * mrSges(5,1) + t21 * mrSges(5,2) + qJD(1) * t50;
t61 = (m(3) * qJ(2));
t60 = Ifges(5,4) * t21;
t16 = qJD(1) * t22;
t59 = t24 * t16;
t15 = qJD(1) * t23;
t58 = t68 * t15;
t36 = qJD(1) * qJ(2) + qJD(3);
t52 = t16 * mrSges(5,1) + t15 * mrSges(5,2);
t51 = -pkin(5) * qJD(1) + t30;
t17 = t51 * t40;
t18 = t51 * t41;
t6 = t44 * t17 + t43 * t18;
t5 = -t43 * t17 + t44 * t18;
t27 = t63 * t40;
t28 = t63 * t41;
t11 = t44 * t27 + t43 * t28;
t10 = -t43 * t27 + t44 * t28;
t48 = t24 * qJD(3);
t47 = t68 * qJD(3);
t1 = qJD(1) * t48 + qJD(4) * t5;
t2 = -qJD(1) * t47 - qJD(4) * t6;
t46 = -t1 * t24 + t2 * t68 + t6 * t22 + t5 * t23;
t45 = qJD(1) ^ 2;
t33 = t40 * pkin(3) + qJ(2);
t29 = pkin(3) * t54 + t36;
t19 = Ifges(5,4) * t20;
t13 = qJD(4) * mrSges(5,1) - t21 * mrSges(5,3);
t12 = -qJD(4) * mrSges(5,2) + t20 * mrSges(5,3);
t8 = Ifges(5,1) * t21 + Ifges(5,5) * qJD(4) + t19;
t7 = Ifges(5,2) * t20 + Ifges(5,6) * qJD(4) + t60;
t4 = -qJD(4) * t11 - t47;
t3 = qJD(4) * t10 + t48;
t9 = [t29 * (t22 * mrSges(5,1) + t23 * mrSges(5,2)) + t33 * t52 + t7 * t65 + qJD(4) * (Ifges(5,5) * t23 - Ifges(5,6) * t22) / 0.2e1 + t8 * t64 + t4 * t13 + t3 * t12 + t62 * qJD(2) + (t20 * t65 - t59) * Ifges(5,2) + (t21 * t64 + t58) * Ifges(5,1) + m(5) * (t29 * qJD(2) + t1 * t11 + t2 * t10 + t6 * t3 + t5 * t4) + m(4) * (t36 * qJD(2) - qJD(3) * t69) + ((m(4) * qJ(2) + m(5) * t33 - t24 * mrSges(5,1) + mrSges(5,2) * t68 + (2 * mrSges(3,3)) + t50 + (2 * t61)) * qJD(2) + (-m(4) * t42 + (2 * mrSges(4,3))) * t67) * qJD(1) + (-t10 * t15 - t11 * t16 - t46) * mrSges(5,3) + (t24 * t15 - t16 * t68 + t20 * t64 + t21 * t65) * Ifges(5,4); m(5) * t46 + t22 * t12 + t23 * t13 + (-mrSges(3,3) - t61) * t45 + (-t58 + t59) * mrSges(5,3) + (-m(5) * t29 + (-t36 - t67) * m(4) - t62) * qJD(1); -m(5) * (t6 * t20 - t5 * t21) - t20 * t12 + t21 * t13 - t55 * t45 * mrSges(4,3) + (m(5) * qJD(2) + (qJD(2) + t69) * m(4)) * qJD(1) + t52; Ifges(5,5) * t15 - Ifges(5,6) * t16 - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t29 * (t21 * mrSges(5,1) + t20 * mrSges(5,2)) - t21 * (Ifges(5,1) * t20 - t60) / 0.2e1 + t21 * t7 / 0.2e1 - qJD(4) * (Ifges(5,5) * t20 - Ifges(5,6) * t21) / 0.2e1 - t5 * t12 + t6 * t13 + (t5 * t20 + t6 * t21) * mrSges(5,3) - (-Ifges(5,2) * t21 + t19 + t8) * t20 / 0.2e1;];
tauc = t9(:);
