% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:04
% EndTime: 2018-11-14 13:38:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (63->17), mult. (176->33), div. (0->0), fcn. (140->4), ass. (0->14)
t11 = sin(pkin(5));
t12 = cos(pkin(5));
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t10 = t11 * t14 + t12 * t13;
t9 = -t11 * t13 + t12 * t14;
t17 = qJD(4) * (t10 * rSges(5,1) + t9 * rSges(5,2));
t18 = qJD(4) * (t9 * rSges(5,1) - t10 * rSges(5,2));
t16 = m(5) * qJD(4);
t8 = t10 * qJD(4);
t7 = t9 * qJD(4);
t4 = t8 * rSges(5,1) + t7 * rSges(5,2);
t3 = t7 * rSges(5,1) - t8 * rSges(5,2);
t1 = [0; (-t11 * t4 - t12 * t3) * t16; (t11 * t3 - t12 * t4) * t16; (t3 * t17 + (t3 - t18) * (-qJD(2) * t12 + qJD(3) * t11 + t17) + (-t4 + t17) * (qJD(2) * t11 + qJD(3) * t12) + (-0.2e1 * t4 + t17) * t18) * m(5);];
tauc  = t1(:);
