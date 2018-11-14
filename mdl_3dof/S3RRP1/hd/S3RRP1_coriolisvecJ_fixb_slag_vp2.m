% Calculate vector of centrifugal and coriolis load on the joints for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3RRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (63->31), mult. (142->45), div. (0->0), fcn. (33->2), ass. (0->18)
t18 = mrSges(3,1) + mrSges(4,1);
t17 = pkin(1) * qJD(1);
t16 = pkin(1) * qJD(2);
t6 = qJD(1) + qJD(2);
t15 = t6 * qJD(3);
t7 = sin(qJ(2));
t14 = t7 * t16;
t8 = cos(qJ(2));
t13 = t8 * t16;
t12 = qJD(1) * t16;
t10 = t8 * t12;
t9 = t7 * t12;
t5 = qJD(3) + t13;
t4 = qJ(3) * t6 + t7 * t17;
t3 = -t6 * pkin(2) - t8 * t17 + qJD(3);
t2 = t10 + t15;
t1 = t2 * mrSges(4,3);
t11 = [t1 + t5 * t6 * mrSges(4,3) + m(4) * (t2 * (pkin(1) * t7 + qJ(3)) + t4 * t5 + (t3 + (-pkin(1) * t8 - pkin(2)) * qJD(1)) * t14) + (-t6 * t13 - t10) * mrSges(3,2) + t18 * (-t6 * t14 - t9); mrSges(4,3) * t15 + t1 + m(4) * (t2 * qJ(3) + qJD(3) * t4) + (-m(4) * (t3 * t7 + t4 * t8) + ((mrSges(3,2) - mrSges(4,3)) * t8 + t18 * t7) * t6 + (-t8 * mrSges(3,2) + (-m(4) * pkin(2) - t18) * t7) * qJD(2)) * t17; -t6 ^ 2 * mrSges(4,3) + (-t4 * t6 + t9) * m(4);];
tauc  = t11(:);
