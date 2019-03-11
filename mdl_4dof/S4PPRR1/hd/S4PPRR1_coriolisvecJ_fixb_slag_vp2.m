% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR1
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:43
% EndTime: 2019-03-08 18:15:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (105->33), mult. (284->59), div. (0->0), fcn. (176->4), ass. (0->24)
t11 = qJD(3) + qJD(4);
t24 = t11 * mrSges(5,1);
t13 = sin(qJ(3));
t23 = qJD(2) * t13;
t12 = sin(qJ(4));
t22 = qJD(4) * t12;
t14 = cos(qJ(4));
t21 = qJD(4) * t14;
t15 = cos(qJ(3));
t20 = t12 * t15 + t14 * t13;
t19 = -t12 * t13 + t14 * t15;
t18 = t20 * qJD(3);
t17 = t19 * qJD(3);
t10 = qJD(3) * pkin(3) + t15 * qJD(2);
t9 = t19 * qJD(2);
t8 = t20 * qJD(2);
t7 = t12 * t10 + t14 * t23;
t6 = t14 * t10 - t12 * t23;
t5 = -t20 * qJD(4) - t18;
t4 = t19 * qJD(4) + t17;
t3 = -t10 * t22 + (-t13 * t21 - t18) * qJD(2);
t2 = t10 * t21 + (-t13 * t22 + t17) * qJD(2);
t1 = t3 * mrSges(5,1);
t16 = [0; m(5) * (t3 * t19 + t2 * t20 + t7 * t4 + t6 * t5) + (-t13 * mrSges(4,1) - t15 * mrSges(4,2)) * qJD(3) ^ 2 + (t5 * mrSges(5,1) - t4 * mrSges(5,2)) * t11; t1 + t8 * t24 - m(5) * (-t6 * t8 + t7 * t9) + (t9 * t11 - t2) * mrSges(5,2) + (m(5) * (t2 * t12 + t3 * t14 + (-t12 * t6 + t14 * t7) * qJD(4)) + (-mrSges(5,1) * t12 - mrSges(5,2) * t14) * qJD(4) * t11) * pkin(3); t7 * t24 + t1 + (t11 * t6 - t2) * mrSges(5,2);];
tauc  = t16(:);
