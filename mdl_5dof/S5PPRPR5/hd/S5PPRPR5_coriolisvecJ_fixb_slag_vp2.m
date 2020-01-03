% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:24
% EndTime: 2019-12-31 17:33:25
% DurationCPUTime: 0.35s
% Computational Cost: add. (306->103), mult. (680->151), div. (0->0), fcn. (292->4), ass. (0->46)
t15 = sin(qJ(5));
t46 = t15 / 0.2e1;
t17 = cos(qJ(5));
t45 = -t17 / 0.2e1;
t19 = -pkin(3) - pkin(6);
t18 = cos(qJ(3));
t30 = t18 * qJD(2);
t27 = qJD(4) - t30;
t21 = t19 * qJD(3) + t27;
t3 = t15 * qJD(1) + t17 * t21;
t44 = t17 * t3;
t16 = sin(qJ(3));
t31 = t16 * qJD(2);
t28 = qJD(3) * t31;
t4 = t17 * qJD(1) - t15 * t21;
t2 = qJD(5) * t4 + t17 * t28;
t43 = t2 * t17;
t42 = Ifges(6,4) * t15;
t41 = Ifges(6,4) * t17;
t14 = qJD(3) * qJ(4) + t31;
t40 = t14 * t18;
t39 = qJD(3) * pkin(3);
t38 = Ifges(6,5) * qJD(5);
t37 = Ifges(6,6) * qJD(5);
t36 = qJD(3) * t15;
t35 = qJD(3) * t17;
t34 = qJD(5) * t15;
t33 = qJD(5) * t17;
t32 = t14 * qJD(3);
t29 = t4 * t33;
t1 = qJD(5) * t3 + t15 * t28;
t26 = t1 * t15 + t43;
t25 = -t15 * t3 - t17 * t4;
t24 = mrSges(6,1) * t17 - mrSges(6,2) * t15;
t12 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t36;
t13 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t35;
t23 = t15 * t12 + t17 * t13;
t10 = (qJD(4) + t30) * qJD(3);
t22 = t10 * qJ(4) + t14 * qJD(4);
t20 = qJD(3) ^ 2;
t11 = t27 - t39;
t9 = (t15 * mrSges(6,1) + t17 * mrSges(6,2)) * qJD(3);
t8 = t24 * qJD(5) * qJD(3);
t7 = t38 + (Ifges(6,1) * t17 - t42) * qJD(3);
t6 = t37 + (-Ifges(6,2) * t15 + t41) * qJD(3);
t5 = [m(6) * (-t1 * t17 + t2 * t15) + (m(6) * (-t15 * t4 + t44) + (t15 ^ 2 + t17 ^ 2) * qJD(3) * mrSges(6,3) + t23) * qJD(5); (t8 + (-mrSges(4,1) + mrSges(5,2)) * t20 + t23 * qJD(3) + m(5) * (qJD(3) * t11 + t10) + m(6) * (t3 * t35 - t4 * t36 + t10)) * t16 + (qJD(3) * t9 + (-mrSges(4,2) + mrSges(5,3)) * t20 + (-t17 * t12 + t15 * t13) * qJD(5) + m(5) * (-t28 + t32) + m(6) * (t3 * t34 - t26 + t29 + t32)) * t18; qJ(4) * t8 + t27 * t9 + (t27 * qJD(3) + t10) * mrSges(5,3) + (-t13 * t31 + t10 * mrSges(6,2) - t2 * mrSges(6,3) + (t14 * mrSges(6,1) + t4 * mrSges(6,3) + t19 * t12 - t6 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,4) * t35 - t37 / 0.2e1) * qJD(5)) * t17 + m(5) * t22 + m(6) * (t22 + (-t29 + t43) * t19) + 0.2e1 * (-m(6) * (t16 * t44 + t40) / 0.2e1 - (t40 + (t11 + t39) * t16) * m(5) / 0.2e1) * qJD(2) + (t10 * mrSges(6,1) + (-t14 * mrSges(6,2) + t3 * mrSges(6,3) - t7 / 0.2e1 - t38 / 0.2e1 + (-m(6) * t3 - t13) * t19 + (0.3e1 / 0.2e1 * t42 + (0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(6,1)) * t17) * qJD(3)) * qJD(5) + (m(6) * t4 - t12) * t31 + (t19 * m(6) - mrSges(6,3)) * t1) * t15; m(6) * (t25 * qJD(5) + t26) + t12 * t33 - t13 * t34 - t20 * mrSges(5,3) + (-m(6) * t14 - t9 + (-t14 + t31) * m(5)) * qJD(3); t2 * mrSges(6,1) - t1 * mrSges(6,2) - t3 * t12 - t4 * t13 + (-t14 * t24 + t7 * t46 + t17 * t6 / 0.2e1 + ((-Ifges(6,1) * t15 - t41) * t45 + (-Ifges(6,2) * t17 - t42) * t46) * qJD(3) + t25 * mrSges(6,3) + (-Ifges(6,5) * t15 / 0.2e1 + Ifges(6,6) * t45) * qJD(5)) * qJD(3);];
tauc = t5(:);
