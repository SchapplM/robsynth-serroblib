% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:42
% DurationCPUTime: 0.29s
% Computational Cost: add. (208->74), mult. (570->122), div. (0->0), fcn. (276->4), ass. (0->39)
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t12 = t21 * qJD(1) + t19 * qJD(2);
t7 = qJD(3) * pkin(5) + t12;
t29 = (t18 ^ 2 + t20 ^ 2) * t7;
t32 = qJD(3) * t18;
t13 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t32;
t31 = qJD(3) * t20;
t14 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t31;
t22 = qJD(3) ^ 2;
t26 = mrSges(5,1) * t18 + mrSges(5,2) * t20;
t5 = t26 * qJD(4) * qJD(3);
t44 = t22 * mrSges(4,2) + t5 + (t13 * t18 - t14 * t20) * qJD(3);
t43 = -t18 / 0.2e1;
t42 = t18 / 0.2e1;
t9 = qJD(3) * t12;
t40 = t9 * t19;
t39 = t9 * t21;
t38 = Ifges(5,4) * t18;
t37 = Ifges(5,2) * t20;
t35 = qJD(4) * t7;
t34 = Ifges(5,5) * qJD(4);
t33 = Ifges(5,6) * qJD(4);
t11 = -t19 * qJD(1) + t21 * qJD(2);
t8 = qJD(3) * t11;
t15 = Ifges(5,4) * t31;
t30 = m(5) * pkin(5) + mrSges(5,3);
t1 = -t18 * t35 + t20 * t8;
t2 = -t18 * t8 - t20 * t35;
t27 = t1 * t20 - t18 * t2;
t25 = t13 * t20 + t14 * t18;
t10 = (-t20 * mrSges(5,1) + t18 * mrSges(5,2)) * qJD(3);
t23 = -t22 * mrSges(4,1) + qJD(3) * t10 - t25 * qJD(4);
t6 = -qJD(3) * pkin(3) - t11;
t4 = Ifges(5,1) * t32 + t15 + t34;
t3 = t33 + (t37 + t38) * qJD(3);
t16 = [t23 * t21 + t44 * t19 + m(4) * (t40 + t8 * t21 + (-t11 * t21 - t12 * t19) * qJD(3)) + m(5) * (t40 + t27 * t21 + (-t19 * t29 + t21 * t6) * qJD(3)); -t44 * t21 + t23 * t19 + m(4) * (t8 * t19 - t39 + (-t11 * t19 + t12 * t21) * qJD(3)) + m(5) * (-t39 + t27 * t19 + (t19 * t6 + t21 * t29) * qJD(3)); -pkin(3) * t5 - t12 * t10 + (-t9 * mrSges(5,1) - t11 * t14 + t30 * t1 + (t4 / 0.2e1 - pkin(5) * t13 + t6 * mrSges(5,2) + 0.3e1 / 0.2e1 * t15 + t34 / 0.2e1) * qJD(4)) * t20 + (t9 * mrSges(5,2) + t11 * t13 - t30 * t2 + (-t3 / 0.2e1 - pkin(5) * t14 + t6 * mrSges(5,1) - t33 / 0.2e1 + (-0.3e1 / 0.2e1 * t38 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t20) * qJD(3)) * qJD(4)) * t18 + (-t9 * pkin(3) - t11 * t29 - t6 * t12) * m(5); t2 * mrSges(5,1) - t1 * mrSges(5,2) + t25 * t7 + (-t6 * t26 + t3 * t42 + ((Ifges(5,1) * t20 - t38) * t43 + t37 * t42) * qJD(3) + (Ifges(5,5) * t20 / 0.2e1 + Ifges(5,6) * t43) * qJD(4) - (t4 + t15) * t20 / 0.2e1) * qJD(3);];
tauc = t16(:);
