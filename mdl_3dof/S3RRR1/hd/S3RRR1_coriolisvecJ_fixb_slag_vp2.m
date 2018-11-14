% Calculate vector of centrifugal and coriolis load on the joints for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3RRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:53
% EndTime: 2018-11-14 10:15:54
% DurationCPUTime: 0.18s
% Computational Cost: add. (161->42), mult. (430->70), div. (0->0), fcn. (198->4), ass. (0->30)
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t33 = mrSges(3,1) * t15 + mrSges(3,2) * t17;
t13 = qJD(1) + qJD(2);
t12 = qJD(3) + t13;
t30 = t12 * mrSges(4,1);
t14 = sin(qJ(3));
t29 = t14 * t15;
t16 = cos(qJ(3));
t28 = t15 * t16;
t27 = pkin(1) * qJD(1);
t26 = qJD(3) * t14;
t25 = qJD(3) * t16;
t22 = t15 * t27;
t21 = -t14 * t17 - t28;
t20 = t16 * t17 - t29;
t19 = (t20 * qJD(2) - t15 * t26) * pkin(1);
t18 = (t21 * qJD(2) - t15 * t25) * pkin(1);
t11 = t17 * pkin(1) + pkin(2);
t10 = t13 * pkin(2) + t17 * t27;
t9 = t20 * t27;
t8 = t21 * t27;
t7 = t14 * t10 + t16 * t22;
t6 = t16 * t10 - t14 * t22;
t5 = -t11 * t26 + t18;
t4 = t11 * t25 + t19;
t3 = qJD(1) * t18 - t10 * t26;
t2 = qJD(1) * t19 + t10 * t25;
t1 = t3 * mrSges(4,1);
t23 = [t1 + t5 * t30 + m(4) * (t2 * (pkin(1) * t28 + t14 * t11) + t7 * t4 + t3 * (-pkin(1) * t29 + t16 * t11) + t6 * t5) + (-t4 * t12 - t2) * mrSges(4,2) + t33 * pkin(1) * qJD(2) * (-qJD(1) - t13); t1 - t8 * t30 - m(4) * (t6 * t8 + t7 * t9) + (t9 * t12 - t2) * mrSges(4,2) + (m(4) * (t2 * t14 + t3 * t16 + (-t14 * t6 + t16 * t7) * qJD(3)) + (-mrSges(4,1) * t14 - mrSges(4,2) * t16) * qJD(3) * t12) * pkin(2) + t33 * t27 * (-qJD(2) + t13); t7 * t30 + t1 + (t12 * t6 - t2) * mrSges(4,2);];
tauc  = t23(:);
