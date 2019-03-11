% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:48
% EndTime: 2019-03-08 18:21:48
% DurationCPUTime: 0.21s
% Computational Cost: add. (279->57), mult. (728->93), div. (0->0), fcn. (528->6), ass. (0->36)
t27 = sin(pkin(6));
t28 = cos(pkin(6));
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t21 = t27 * t32 + t28 * t30;
t17 = t21 * qJD(1);
t22 = -t27 * t30 + t28 * t32;
t18 = t22 * qJD(1);
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t25 = t28 * pkin(2) + pkin(3);
t39 = pkin(2) * t27;
t34 = t31 * t25 - t29 * t39;
t41 = t34 * qJD(4) + t29 * t17 - t31 * t18;
t35 = t29 * t25 + t31 * t39;
t40 = -t35 * qJD(4) + t31 * t17 + t29 * t18;
t38 = qJD(1) * t30;
t24 = qJD(2) * pkin(2) + t32 * qJD(1);
t11 = t28 * t24 - t27 * t38;
t10 = qJD(2) * pkin(3) + t11;
t12 = t27 * t24 + t28 * t38;
t6 = t31 * t10 - t29 * t12;
t7 = t29 * t10 + t31 * t12;
t37 = t31 * t21 + t29 * t22;
t36 = -t29 * t21 + t31 * t22;
t19 = t21 * qJD(2);
t20 = t22 * qJD(2);
t26 = qJD(2) + qJD(4);
t16 = qJD(1) * t20;
t15 = qJD(1) * t19;
t5 = -qJD(4) * t37 - t31 * t19 - t29 * t20;
t4 = qJD(4) * t36 - t29 * t19 + t31 * t20;
t3 = -qJD(4) * t7 - t31 * t15 - t29 * t16;
t2 = qJD(4) * t6 - t29 * t15 + t31 * t16;
t1 = t3 * mrSges(5,1);
t8 = [(t5 * mrSges(5,1) - t4 * mrSges(5,2)) * t26 + m(4) * (-t11 * t19 + t12 * t20 - t15 * t22 + t16 * t21) + m(5) * (t2 * t37 + t3 * t36 + t7 * t4 + t6 * t5) + (-t19 * mrSges(4,1) - t20 * mrSges(4,2) + (-t30 * mrSges(3,1) - t32 * mrSges(3,2)) * qJD(2)) * qJD(2); -t2 * mrSges(5,2) + t1 + (t18 * qJD(2) - t16) * mrSges(4,2) + (t17 * qJD(2) - t15) * mrSges(4,1) + (t40 * mrSges(5,1) - t41 * mrSges(5,2)) * t26 + (t2 * t35 + t3 * t34 + t40 * t6 + t41 * t7) * m(5) + ((-t15 * t28 + t16 * t27) * pkin(2) + t11 * t17 - t12 * t18) * m(4); 0; t7 * t26 * mrSges(5,1) + t1 + (t26 * t6 - t2) * mrSges(5,2);];
tauc  = t8(:);
