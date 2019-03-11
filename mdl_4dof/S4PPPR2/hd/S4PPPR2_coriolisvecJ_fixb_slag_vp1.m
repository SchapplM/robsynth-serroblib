% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:07
% EndTime: 2019-03-08 18:10:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (58->16), mult. (156->29), div. (0->0), fcn. (126->4), ass. (0->14)
t12 = sin(pkin(5));
t13 = cos(pkin(5));
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t10 = -t12 * t15 + t13 * t14;
t9 = -t12 * t14 - t13 * t15;
t18 = qJD(4) * (-rSges(5,1) * t10 + rSges(5,2) * t9);
t19 = qJD(4) * (rSges(5,1) * t9 + rSges(5,2) * t10);
t17 = m(5) * qJD(4);
t8 = t9 * qJD(4);
t7 = t10 * qJD(4);
t4 = rSges(5,1) * t8 + rSges(5,2) * t7;
t3 = -rSges(5,1) * t7 + rSges(5,2) * t8;
t1 = [-t3 * t17; 0; (t12 * t4 + t13 * t3) * t17; (t4 * t18 + (t4 - t19) * (qJD(3) * t12 + t18) + (-t3 + t18) * (-qJD(3) * t13 + qJD(1)) + (-0.2e1 * t3 + t18) * t19) * m(5);];
tauc  = t1(:);
