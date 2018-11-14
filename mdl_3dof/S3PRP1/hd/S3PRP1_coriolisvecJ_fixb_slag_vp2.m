% Calculate vector of centrifugal and coriolis load on the joints for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
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
% Datum: 2018-11-14 10:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3PRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:04:14
% EndTime: 2018-11-14 10:04:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (27->19), mult. (78->25), div. (0->0), fcn. (23->2), ass. (0->12)
t4 = sin(qJ(2));
t10 = qJD(1) * t4;
t3 = qJD(2) * qJ(3) + t10;
t5 = cos(qJ(2));
t11 = t3 * t5;
t9 = qJD(1) * t5;
t8 = qJD(2) * pkin(2);
t7 = qJD(3) - t9;
t6 = qJD(2) ^ 2;
t2 = t7 - t8;
t1 = (qJD(3) + t9) * qJD(2);
t12 = [m(4) * (t1 * t4 + (t11 + (t2 - t9) * t4) * qJD(2)) + ((-mrSges(3,2) + mrSges(4,3)) * t5 + (-mrSges(3,1) - mrSges(4,1)) * t4) * t6; (t7 * qJD(2) + t1) * mrSges(4,3) + (qJ(3) * t1 + qJD(3) * t3 - t8 * t10 - (t2 * t4 + t11) * qJD(1)) * m(4); -t6 * mrSges(4,3) + (-t3 + t10) * qJD(2) * m(4);];
tauc  = t12(:);
