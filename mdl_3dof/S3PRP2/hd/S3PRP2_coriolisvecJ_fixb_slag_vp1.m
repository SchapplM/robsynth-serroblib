% Calculate vector of centrifugal and coriolis load on the joints for
% S3PRP2
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
% tauc [3x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3PRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP2_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:07:06
% EndTime: 2018-11-14 10:07:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (115->15), mult. (264->36), div. (0->0), fcn. (130->2), ass. (0->14)
t42 = (rSges(4,1) + pkin(2));
t43 = -qJD(2) * t42 + 2 * qJD(3);
t27 = cos(qJ(2));
t41 = rSges(4,3) + qJ(3);
t38 = t41 * t27;
t26 = sin(qJ(2));
t40 = t41 * t26;
t12 = t26 * rSges(3,1) + t27 * rSges(3,2);
t32 = qJD(2) * t12;
t15 = t27 * rSges(3,1) - t26 * rSges(3,2);
t31 = qJD(2) * t15;
t2 = (-qJD(2) * t40 + t43 * t27) * qJD(2);
t1 = (qJD(2) * t38 + t43 * t26) * qJD(2);
t3 = [-m(3) * qJD(2) * t32 + m(4) * t1; (-t31 * t32 + (0.2e1 * t12 * t31 - t15 * t32) * qJD(2)) * m(3) + (t1 * (t42 * t27 + t40) + t2 * (-t42 * t26 + t38)) * m(4); 0.2e1 * (-t1 * t27 / 0.2e1 + t2 * t26 / 0.2e1) * m(4);];
tauc  = t3(:);
