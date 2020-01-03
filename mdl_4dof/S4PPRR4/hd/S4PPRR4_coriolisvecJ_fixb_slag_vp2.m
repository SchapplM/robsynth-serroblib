% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:26
% EndTime: 2019-12-31 16:18:28
% DurationCPUTime: 0.33s
% Computational Cost: add. (271->83), mult. (807->135), div. (0->0), fcn. (512->6), ass. (0->41)
t43 = -Ifges(5,1) / 0.2e1;
t26 = cos(qJ(4));
t35 = qJD(3) * t26;
t21 = Ifges(5,4) * t35;
t42 = -t21 / 0.2e1;
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t16 = t25 * t22 - t27 * t23;
t9 = t16 * qJD(1);
t41 = 2 * m(5);
t17 = t27 * t22 + t25 * t23;
t12 = t17 * qJD(3);
t8 = qJD(1) * t12;
t40 = t8 * t16;
t24 = sin(qJ(4));
t39 = Ifges(5,4) * t24;
t38 = Ifges(5,5) * qJD(4);
t37 = Ifges(5,6) * qJD(4);
t36 = qJD(3) * t24;
t34 = t38 / 0.2e1;
t33 = -t37 / 0.2e1;
t32 = (-t26 * mrSges(5,1) + t24 * mrSges(5,2) - mrSges(4,1)) * qJD(3);
t10 = t17 * qJD(1);
t6 = qJD(3) * pkin(5) + t10;
t3 = t26 * qJD(2) - t24 * t6;
t4 = t24 * qJD(2) + t26 * t6;
t31 = -t24 * t3 + t26 * t4;
t19 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t36;
t20 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t35;
t30 = -t24 * t19 + t26 * t20;
t5 = -qJD(3) * pkin(3) + t9;
t29 = t3 * mrSges(5,3) + t36 * t43 + t42 - t38 / 0.2e1 - t5 * mrSges(5,2);
t28 = t4 * mrSges(5,3) + t37 / 0.2e1 + (Ifges(5,2) * t26 + t39) * qJD(3) / 0.2e1 - t5 * mrSges(5,1);
t11 = t16 * qJD(3);
t15 = (mrSges(5,1) * t24 + mrSges(5,2) * t26) * qJD(4) * qJD(3);
t7 = qJD(1) * t11;
t2 = -qJD(4) * t4 + t24 * t7;
t1 = qJD(4) * t3 - t26 * t7;
t13 = [t16 * t15 + t32 * t12 + (-t19 * t26 - t20 * t24) * t17 * qJD(4) - (-qJD(3) * mrSges(4,2) + t30) * t11 + m(4) * (-t10 * t11 + t9 * t12 - t7 * t17 + t40) + m(5) * (t5 * t12 + t40 - t31 * t11 + (t1 * t26 - t2 * t24 + (-t4 * t24 - t3 * t26) * qJD(4)) * t17); m(5) * (t1 * t24 + t2 * t26) + (m(5) * t31 + (-t24 ^ 2 - t26 ^ 2) * qJD(3) * mrSges(5,3) + t30) * qJD(4); -t8 * mrSges(4,1) + (-m(5) * t8 - t15) * pkin(3) + (-t9 * qJD(3) + t7) * mrSges(4,2) + (-m(5) * t5 - t32) * t10 + (-t8 * mrSges(5,1) + t1 * mrSges(5,3) + t9 * t20 + (pkin(5) * t1 / 0.2e1 + t4 * t9 / 0.2e1) * t41 + (0.3e1 / 0.2e1 * t21 + t34 + (-m(5) * t3 - t19) * pkin(5) - t29) * qJD(4)) * t26 + (t8 * mrSges(5,2) - t2 * mrSges(5,3) - t9 * t19 + (-pkin(5) * t2 / 0.2e1 - t3 * t9 / 0.2e1) * t41 + (t33 + (-m(5) * t4 - t20) * pkin(5) + (-0.3e1 / 0.2e1 * t39 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t26) * qJD(3) - t28) * qJD(4)) * t24; t2 * mrSges(5,1) - t1 * mrSges(5,2) + t4 * t19 - t3 * t20 + ((t34 + t42 + t29) * t26 + (t33 + (t39 / 0.2e1 + (t43 + Ifges(5,2) / 0.2e1) * t26) * qJD(3) + t28) * t24) * qJD(3);];
tauc = t13(:);
