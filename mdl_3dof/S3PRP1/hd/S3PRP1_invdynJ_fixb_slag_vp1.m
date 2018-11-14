% Calculate vector of inverse dynamics joint torques for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S3PRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRP1_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:04:14
% EndTime: 2018-11-14 10:04:15
% DurationCPUTime: 0.24s
% Computational Cost: add. (171->43), mult. (343->45), div. (0->0), fcn. (175->2), ass. (0->23)
t30 = cos(qJ(2));
t44 = rSges(4,3) + qJ(3);
t39 = t44 * t30;
t29 = sin(qJ(2));
t45 = rSges(4,1) + pkin(2);
t40 = t45 * t29;
t36 = t39 - t40;
t41 = t44 * t29;
t9 = t45 * t30 + t41;
t31 = -qJD(2) * t45 + qJD(3);
t21 = qJD(3) * t29;
t32 = qJD(2) * t39 + t21;
t1 = -qJDD(3) * t30 + qJDD(1) + t9 * qJDD(2) + (t31 * t29 + t32) * qJD(2);
t43 = t1 - g(2);
t22 = qJD(3) * t30;
t42 = qJDD(3) * t29 + t36 * qJDD(2) + (-qJD(2) * t41 + t31 * t30 + t22) * qJD(2) - g(1);
t18 = t30 * rSges(3,1) - t29 * rSges(3,2);
t34 = qJD(2) * t18;
t15 = t29 * rSges(3,1) + t30 * rSges(3,2);
t11 = qJD(2) * t15;
t10 = qJD(1) + t34;
t7 = -qJD(2) * t11 + qJDD(2) * t18 + qJDD(1);
t2 = [m(2) * qJDD(1) + m(3) * t7 + m(4) * t1 + (-m(2) - m(3) - m(4)) * g(2); (Icges(3,3) + Icges(4,2)) * qJDD(2) + (-t10 * t11 + (t7 - g(2)) * t18 + (qJDD(2) * t15 + g(1) + (t10 + t34) * qJD(2)) * t15) * m(3) + (t43 * t9 + t42 * t36 + (-t21 + t32 + (-t36 - t40) * qJD(2)) * (t9 * qJD(2) + qJD(1) - t22)) * m(4); (t42 * t29 - t43 * t30) * m(4);];
tau  = t2;
