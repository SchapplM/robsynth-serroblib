% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.18s
% Computational Cost: add. (161->42), mult. (430->70), div. (0->0), fcn. (198->4), ass. (0->30)
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t33 = mrSges(4,1) * t15 + mrSges(4,2) * t17;
t13 = qJD(2) + qJD(3);
t12 = qJD(4) + t13;
t30 = t12 * mrSges(5,1);
t14 = sin(qJ(4));
t29 = t14 * t15;
t16 = cos(qJ(4));
t28 = t15 * t16;
t27 = qJD(2) * pkin(1);
t26 = qJD(4) * t14;
t25 = qJD(4) * t16;
t22 = t15 * t27;
t21 = -t14 * t17 - t28;
t20 = t16 * t17 - t29;
t19 = (t20 * qJD(3) - t15 * t26) * pkin(1);
t18 = (t21 * qJD(3) - t15 * t25) * pkin(1);
t11 = t17 * pkin(1) + pkin(2);
t10 = t13 * pkin(2) + t17 * t27;
t9 = t20 * t27;
t8 = t21 * t27;
t7 = t14 * t10 + t16 * t22;
t6 = t16 * t10 - t14 * t22;
t5 = -t11 * t26 + t18;
t4 = t11 * t25 + t19;
t3 = qJD(2) * t18 - t10 * t26;
t2 = qJD(2) * t19 + t10 * t25;
t1 = t3 * mrSges(5,1);
t23 = [0; t1 + t5 * t30 + m(5) * (t2 * (pkin(1) * t28 + t14 * t11) + t7 * t4 + t3 * (-pkin(1) * t29 + t16 * t11) + t6 * t5) + (-t4 * t12 - t2) * mrSges(5,2) + t33 * pkin(1) * qJD(3) * (-qJD(2) - t13); t1 - t8 * t30 - m(5) * (t6 * t8 + t7 * t9) + (t9 * t12 - t2) * mrSges(5,2) + (m(5) * (t2 * t14 + t3 * t16 + (-t14 * t6 + t16 * t7) * qJD(4)) + (-mrSges(5,1) * t14 - mrSges(5,2) * t16) * qJD(4) * t12) * pkin(2) + t33 * t27 * (-qJD(3) + t13); t7 * t30 + t1 + (t12 * t6 - t2) * mrSges(5,2);];
tauc  = t23(:);
