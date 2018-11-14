% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:03
% EndTime: 2018-11-14 14:06:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (63->27), mult. (172->43), div. (0->0), fcn. (70->2), ass. (0->22)
t17 = sin(qJ(3));
t18 = cos(qJ(3));
t26 = rSges(5,2) * t18;
t20 = rSges(5,1) * t17 + t26;
t21 = rSges(4,1) * t17 + rSges(4,2) * t18;
t7 = t21 * qJD(3);
t29 = -rSges(5,1) - pkin(3);
t32 = t29 * t18;
t19 = qJD(3) ^ 2;
t31 = t19 * pkin(3);
t27 = rSges(5,1) * t18;
t25 = m(4) * qJD(3);
t12 = rSges(4,1) * t18 - t17 * rSges(4,2);
t24 = qJD(3) * t12;
t23 = qJD(3) * t17;
t15 = t17 * rSges(5,2);
t13 = rSges(5,2) * t23;
t6 = qJD(1) - t7;
t5 = -qJD(2) - t24;
t2 = -t18 * t31 + qJD(3) * (-qJD(3) * t27 + t13);
t1 = t17 * t31 + t20 * t19;
t3 = [m(5) * t2 - t24 * t25; -m(5) * t1 - t7 * t25; -(-t12 * t6 + t21 * t5) * t25 + m(4) * (t5 * t7 - t6 * t24 + (-t12 * t7 + t21 * t24) * qJD(3)) + (t2 * (t29 * t17 - t26) + t1 * (t15 + t32) + (t13 + (pkin(3) * t18 - t15 + t27 + t32) * qJD(3)) * (-pkin(3) * t23 - qJD(3) * t20 + qJD(1))) * m(5); 0;];
tauc  = t3(:);
