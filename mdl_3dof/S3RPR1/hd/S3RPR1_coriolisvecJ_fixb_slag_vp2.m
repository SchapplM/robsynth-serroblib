% Calculate vector of centrifugal and coriolis load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3RPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:22
% EndTime: 2018-11-14 10:14:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (117->33), mult. (228->51), div. (0->0), fcn. (74->2), ass. (0->22)
t10 = cos(qJ(3));
t8 = -qJD(1) + qJD(3);
t9 = sin(qJ(3));
t24 = (mrSges(4,1) * t9 + mrSges(4,2) * t10) * t8;
t23 = (m(3) * qJ(2) + mrSges(3,3)) * qJD(1);
t11 = -pkin(1) - pkin(2);
t7 = t11 * qJD(1) + qJD(2);
t22 = qJD(3) * t7;
t21 = t9 * qJD(2);
t20 = t10 * qJD(2);
t19 = qJ(2) * qJD(1);
t18 = qJ(2) * qJD(3);
t5 = t10 * t7 - t9 * t19;
t6 = t10 * t19 + t9 * t7;
t16 = t10 * t6 - t5 * t9;
t14 = -t9 * qJ(2) + t10 * t11;
t13 = t10 * qJ(2) + t9 * t11;
t4 = -t13 * qJD(3) - t21;
t3 = t14 * qJD(3) + t20;
t2 = -t9 * t22 + (-t10 * t18 - t21) * qJD(1);
t1 = t10 * t22 + (-t9 * t18 + t20) * qJD(1);
t12 = [m(4) * (t1 * t13 + t2 * t14 + t6 * t3 + t5 * t4) + 0.2e1 * qJD(2) * t23 + (-t3 * t8 + t1) * mrSges(4,2) + (t4 * t8 - t2) * mrSges(4,1); m(4) * (t16 * qJD(3) + t1 * t9 + t2 * t10) - qJD(3) * t24 + (-m(4) * t16 - t23 + t24) * qJD(1); (t5 * t8 - t1) * mrSges(4,2) + (t6 * t8 + t2) * mrSges(4,1);];
tauc  = t12(:);
