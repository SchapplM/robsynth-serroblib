% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:23
% DurationCPUTime: 0.30s
% Computational Cost: add. (176->66), mult. (505->107), div. (0->0), fcn. (241->4), ass. (0->36)
t19 = qJD(3) ^ 2;
t27 = Ifges(5,5) * qJD(4) / 0.2e1;
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t16 = sin(qJ(3));
t11 = qJD(3) * pkin(5) + t16 * qJD(2);
t22 = t15 * qJD(1) - t17 * t11;
t3 = -t17 * qJD(1) - t15 * t11;
t24 = t15 * t3 + t17 * t22;
t33 = qJD(3) * t17;
t10 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t33;
t34 = qJD(3) * t15;
t9 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t34;
t23 = t17 * t10 - t15 * t9;
t43 = -m(5) * t24 + t23;
t42 = t15 ^ 2;
t41 = t17 ^ 2;
t38 = Ifges(5,4) * t15;
t37 = qJD(3) * pkin(3);
t35 = Ifges(5,6) * qJD(4);
t32 = qJD(4) * t15;
t31 = qJD(4) * t17;
t18 = cos(qJ(3));
t30 = t18 * qJD(2);
t29 = qJD(4) * qJD(3);
t28 = qJD(3) * t30;
t26 = -t35 / 0.2e1;
t1 = qJD(4) * t3 + t17 * t28;
t2 = t22 * qJD(4) - t15 * t28;
t25 = t1 * t17 - t15 * t2;
t12 = -t30 - t37;
t13 = Ifges(5,4) * t33;
t21 = t12 * mrSges(5,2) + Ifges(5,1) * t34 / 0.2e1 + t13 / 0.2e1 + t27 - t3 * mrSges(5,3);
t20 = -t22 * mrSges(5,3) + t35 / 0.2e1 + (Ifges(5,2) * t17 + t38) * qJD(3) / 0.2e1 - t12 * mrSges(5,1);
t7 = (mrSges(5,1) * t15 + mrSges(5,2) * t17) * t29;
t4 = [m(5) * (-t1 * t15 - t2 * t17) + ((t41 + t42) * qJD(3) * mrSges(5,3) - t43) * qJD(4); (-t19 * mrSges(4,2) + t43 * qJD(3) - t7) * t18 + (-t10 * t32 - t9 * t31 + m(5) * (qJD(3) * t12 + t22 * t32 - t3 * t31 + t25 - t28) + (-t17 * mrSges(5,1) + t15 * mrSges(5,2) - mrSges(4,1)) * t19) * t16; -pkin(3) * t7 + (-t23 * t18 + (t24 * t18 + (-t12 - t37) * t16) * m(5)) * qJD(2) + ((t27 + t21) * t17 + (t26 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t33 - t20) * t15 + (-t15 * t10 - t17 * t9 + m(5) * (t15 * t22 - t17 * t3)) * pkin(5)) * qJD(4) + (0.3e1 / 0.2e1 * t41 - 0.3e1 / 0.2e1 * t42) * Ifges(5,4) * t29 + (m(5) * pkin(5) + mrSges(5,3)) * t25; t2 * mrSges(5,1) - t1 * mrSges(5,2) - t3 * t10 - t22 * t9 + ((t27 - t13 / 0.2e1 - t21) * t17 + (t26 + (t38 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t17) * qJD(3) + t20) * t15) * qJD(3);];
tauc = t4(:);
