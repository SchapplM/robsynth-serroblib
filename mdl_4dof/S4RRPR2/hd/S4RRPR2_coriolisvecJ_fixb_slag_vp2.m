% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-05-28 15:34
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-28 15:34:02
% EndTime: 2019-05-28 15:34:03
% DurationCPUTime: 0.27s
% Computational Cost: add. (454->80), mult. (716->112), div. (0->0), fcn. (275->4), ass. (0->44)
t45 = mrSges(3,1) + mrSges(4,1);
t22 = sin(qJ(4));
t23 = sin(qJ(2));
t24 = cos(qJ(4));
t25 = cos(qJ(2));
t26 = -pkin(2) - pkin(3);
t28 = qJ(3) * t24 + t22 * t26;
t44 = pkin(1) * qJD(1);
t51 = -t22 * qJD(3) - t28 * qJD(4) - (-t22 * t25 + t23 * t24) * t44;
t27 = -qJ(3) * t22 + t24 * t26;
t50 = -t24 * qJD(3) - t27 * qJD(4) + (t22 * t23 + t24 * t25) * t44;
t21 = qJD(1) + qJD(2);
t20 = qJD(4) - t21;
t49 = (mrSges(5,1) * t22 + mrSges(5,2) * t24) * t20;
t47 = mrSges(4,3) * t21;
t46 = t20 * mrSges(5,2);
t43 = pkin(1) * qJD(2);
t42 = t21 * qJD(3);
t41 = t23 * t43;
t40 = t25 * t43;
t39 = -pkin(1) * t25 - pkin(2);
t38 = qJD(1) * t43;
t34 = t25 * t38;
t14 = t34 + t42;
t35 = t23 * t38;
t32 = -t25 * t44 + qJD(3);
t10 = t21 * t26 + t32;
t16 = qJ(3) * t21 + t23 * t44;
t6 = t10 * t24 - t16 * t22;
t2 = qJD(4) * t6 + t24 * t14 + t22 * t35;
t1 = t2 * mrSges(5,2);
t7 = t10 * t22 + t16 * t24;
t3 = -qJD(4) * t7 - t22 * t14 + t24 * t35;
t36 = -t3 * mrSges(5,1) + t14 * mrSges(4,3) + t1;
t33 = -t22 * t6 + t24 * t7;
t18 = -pkin(3) + t39;
t19 = pkin(1) * t23 + qJ(3);
t30 = t18 * t24 - t19 * t22;
t29 = t18 * t22 + t19 * t24;
t17 = qJD(3) + t40;
t15 = -pkin(2) * t21 + t32;
t5 = -qJD(4) * t29 - t22 * t17 + t24 * t41;
t4 = qJD(4) * t30 + t24 * t17 + t22 * t41;
t8 = [-t4 * t46 + t5 * t20 * mrSges(5,1) + t17 * t47 + m(5) * (t2 * t29 + t3 * t30 + t7 * t4 + t6 * t5) + m(4) * (t14 * t19 + t16 * t17 + (qJD(1) * t39 + t15) * t41) + t36 + (-t21 * t40 - t34) * mrSges(3,2) + t45 * (-t21 * t41 - t35); mrSges(4,3) * t42 + m(4) * (t14 * qJ(3) + t16 * qJD(3)) + (mrSges(5,1) * t51 + mrSges(5,2) * t50) * t20 + (-m(4) * (t15 * t23 + t16 * t25) + ((mrSges(3,2) - mrSges(4,3)) * t25 + t45 * t23) * t21 + (-t25 * mrSges(3,2) + (-m(4) * pkin(2) - t45) * t23) * qJD(2)) * t44 + t36 + (t2 * t28 + t3 * t27 - t50 * t7 + t51 * t6) * m(5); m(4) * t35 + m(5) * (qJD(4) * t33 + t2 * t22 + t3 * t24) - qJD(4) * t49 + (-m(4) * t16 - m(5) * t33 - t47 + t49) * t21; t6 * t46 - t1 + (t20 * t7 + t3) * mrSges(5,1);];
tauc  = t8(:);
