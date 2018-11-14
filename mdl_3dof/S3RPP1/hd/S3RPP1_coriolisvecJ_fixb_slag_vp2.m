% Calculate vector of centrifugal and coriolis load on the joints for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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

function tauc = S3RPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:33
% EndTime: 2018-11-14 10:13:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (20->14), mult. (48->19), div. (0->0), fcn. (0->0), ass. (0->7)
t9 = (m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3));
t6 = (-pkin(1) - qJ(3));
t5 = (m(4) * qJD(1));
t3 = (qJD(1) ^ 2);
t2 = (qJD(1) * qJ(2) + qJD(3));
t1 = (qJD(1) * t6 + qJD(2));
t4 = [m(4) * (t2 * qJD(2) - t1 * qJD(3) + (qJD(2) * qJ(2) - qJD(3) * t6) * qJD(1)) + 2 * (qJD(3) * mrSges(4,3) + t9 * qJD(2)) * qJD(1); (-qJD(3) - t2) * t5 - t9 * t3; -t3 * mrSges(4,3) + (qJD(2) + t1) * t5;];
tauc  = t4(:);
