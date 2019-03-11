% Calculate vector of centrifugal and Coriolis load on the joints for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3PRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:50
% EndTime: 2019-03-08 18:03:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (105->33), mult. (284->59), div. (0->0), fcn. (176->4), ass. (0->24)
t11 = qJD(2) + qJD(3);
t24 = t11 * mrSges(4,1);
t13 = sin(qJ(2));
t23 = qJD(1) * t13;
t12 = sin(qJ(3));
t22 = qJD(3) * t12;
t14 = cos(qJ(3));
t21 = qJD(3) * t14;
t15 = cos(qJ(2));
t20 = t12 * t15 + t14 * t13;
t19 = -t12 * t13 + t14 * t15;
t18 = t20 * qJD(2);
t17 = t19 * qJD(2);
t10 = qJD(2) * pkin(2) + t15 * qJD(1);
t9 = t19 * qJD(1);
t8 = t20 * qJD(1);
t7 = t12 * t10 + t14 * t23;
t6 = t14 * t10 - t12 * t23;
t5 = -t20 * qJD(3) - t18;
t4 = t19 * qJD(3) + t17;
t3 = -t10 * t22 + (-t13 * t21 - t18) * qJD(1);
t2 = t10 * t21 + (-t13 * t22 + t17) * qJD(1);
t1 = t3 * mrSges(4,1);
t16 = [m(4) * (t3 * t19 + t2 * t20 + t7 * t4 + t6 * t5) + (-t13 * mrSges(3,1) - t15 * mrSges(3,2)) * qJD(2) ^ 2 + (t5 * mrSges(4,1) - t4 * mrSges(4,2)) * t11; t1 + t8 * t24 - m(4) * (-t6 * t8 + t7 * t9) + (t9 * t11 - t2) * mrSges(4,2) + (m(4) * (t2 * t12 + t3 * t14 + (-t12 * t6 + t14 * t7) * qJD(3)) + (-mrSges(4,1) * t12 - mrSges(4,2) * t14) * qJD(3) * t11) * pkin(2); t7 * t24 + t1 + (t11 * t6 - t2) * mrSges(4,2);];
tauc  = t16(:);
