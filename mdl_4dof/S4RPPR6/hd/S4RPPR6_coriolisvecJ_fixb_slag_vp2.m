% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:39
% DurationCPUTime: 0.74s
% Computational Cost: add. (370->112), mult. (1049->176), div. (0->0), fcn. (624->4), ass. (0->59)
t40 = sin(pkin(6));
t38 = t40 ^ 2;
t41 = cos(pkin(6));
t39 = t41 ^ 2;
t71 = (mrSges(4,2) + mrSges(3,3)) * (t38 + t39);
t42 = sin(qJ(4));
t43 = cos(qJ(4));
t57 = t40 * qJD(1);
t59 = qJD(1) * t41;
t21 = -t42 * t57 - t43 * t59;
t68 = -t21 / 0.2e1;
t22 = -t42 * t59 + t43 * t57;
t65 = t22 / 0.2e1;
t64 = Ifges(5,4) * t22;
t63 = -pkin(5) + qJ(2);
t44 = qJD(1) ^ 2;
t61 = qJ(2) * t44;
t53 = t40 * qJ(3) + pkin(1);
t23 = (pkin(2) + pkin(3)) * t41 + t53;
t60 = qJD(1) * t23;
t58 = qJD(3) * t40;
t31 = qJ(2) * t57 + qJD(3);
t56 = qJD(1) * qJD(2);
t30 = t63 * t41;
t54 = qJD(3) * t57;
t52 = qJ(2) * t56;
t48 = t40 * t42 + t41 * t43;
t19 = t48 * qJD(4);
t15 = qJD(1) * t19;
t26 = t40 * t43 - t41 * t42;
t20 = t26 * qJD(4);
t16 = qJD(1) * t20;
t49 = t16 * mrSges(5,1) - t15 * mrSges(5,2);
t24 = -pkin(5) * t57 + t31;
t28 = qJD(1) * t30;
t8 = t43 * t24 - t42 * t28;
t9 = t42 * t24 + t43 * t28;
t29 = t63 * t40;
t10 = t43 * t29 - t42 * t30;
t11 = t42 * t29 + t43 * t30;
t47 = -t41 * pkin(2) - t53;
t27 = (-t41 * mrSges(4,1) - t40 * mrSges(4,3)) * qJD(1);
t46 = t26 * qJD(2);
t45 = t48 * qJD(2);
t37 = t39 * t61;
t32 = t38 * t52;
t18 = t47 * qJD(1) + qJD(2);
t17 = Ifges(5,4) * t21;
t14 = qJD(4) * mrSges(5,1) - t22 * mrSges(5,3);
t13 = -qJD(4) * mrSges(5,2) + t21 * mrSges(5,3);
t12 = -qJD(2) + t60;
t7 = -t21 * mrSges(5,1) + t22 * mrSges(5,2);
t6 = Ifges(5,1) * t22 + Ifges(5,5) * qJD(4) + t17;
t5 = Ifges(5,2) * t21 + Ifges(5,6) * qJD(4) + t64;
t4 = -t11 * qJD(4) + t46;
t3 = t10 * qJD(4) + t45;
t2 = qJD(1) * t46 - t9 * qJD(4);
t1 = qJD(1) * t45 + t8 * qJD(4);
t25 = [(mrSges(5,1) * t48 + t26 * mrSges(5,2)) * t54 + t7 * t58 - t19 * t6 / 0.2e1 + t12 * (t20 * mrSges(5,1) - t19 * mrSges(5,2)) + qJD(4) * (-Ifges(5,5) * t19 - Ifges(5,6) * t20) / 0.2e1 - t20 * t5 / 0.2e1 + t3 * t13 + t4 * t14 + t23 * t49 + 0.2e1 * m(3) * (t39 * t52 + t32) + m(4) * (t32 + (t31 * qJD(2) - t18 * qJD(3)) * t40 + (0.2e1 * t39 * qJD(2) * qJ(2) - t47 * t58) * qJD(1)) - 0.2e1 * t27 * t58 + m(5) * (t1 * t11 + t2 * t10 + t9 * t3 + t8 * t4 + (t12 + t60) * t58) + (t16 * t48 + t20 * t68) * Ifges(5,2) + (-t15 * t26 - t19 * t65) * Ifges(5,1) + (-t1 * t48 + t10 * t15 - t11 * t16 + t8 * t19 - t2 * t26 - t9 * t20) * mrSges(5,3) + (t15 * t48 - t16 * t26 + t19 * t68 - t20 * t65) * Ifges(5,4) + 0.2e1 * t56 * t71; t21 * t13 - t22 * t14 + (-m(4) - m(5)) * t54 - m(3) * (t38 * t61 + t37) - m(4) * (t31 * t57 + t37) - m(5) * (-t9 * t21 + t8 * t22) - t49 - t44 * t71; m(5) * (t1 * t42 + t2 * t43) + (t43 * t15 - t42 * t16) * mrSges(5,3) + (-m(5) * t12 + t27 - t7 + (qJD(2) + t18) * m(4)) * t57 + (m(5) * (-t42 * t8 + t43 * t9) + t43 * t13 - t42 * t14) * qJD(4); -Ifges(5,5) * t15 - Ifges(5,6) * t16 - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t12 * (t22 * mrSges(5,1) + t21 * mrSges(5,2)) - t22 * (Ifges(5,1) * t21 - t64) / 0.2e1 + t5 * t65 - qJD(4) * (Ifges(5,5) * t21 - Ifges(5,6) * t22) / 0.2e1 - t8 * t13 + t9 * t14 + (t8 * t21 + t9 * t22) * mrSges(5,3) + (-Ifges(5,2) * t22 + t17 + t6) * t68;];
tauc = t25(:);
