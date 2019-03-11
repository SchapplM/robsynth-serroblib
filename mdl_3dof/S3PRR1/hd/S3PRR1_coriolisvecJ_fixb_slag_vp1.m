% Calculate vector of centrifugal and Coriolis load on the joints for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3PRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:50
% EndTime: 2019-03-08 18:03:50
% DurationCPUTime: 0.17s
% Computational Cost: add. (165->31), mult. (224->58), div. (0->0), fcn. (95->4), ass. (0->31)
t18 = qJD(2) + qJD(3);
t21 = cos(qJ(2));
t25 = pkin(2) * qJD(2);
t19 = qJ(2) + qJ(3);
t17 = cos(t19);
t15 = t17 * rSges(4,1);
t16 = sin(t19);
t27 = t16 * rSges(4,2);
t9 = t15 - t27;
t3 = t18 * t9 + t21 * t25 + qJD(1);
t8 = t16 * rSges(4,1) + t17 * rSges(4,2);
t32 = t3 * t8;
t20 = sin(qJ(2));
t13 = t20 * rSges(3,1) + t21 * rSges(3,2);
t10 = qJD(2) * t13;
t31 = -t10 / 0.2e1;
t30 = t18 * t8;
t29 = t20 * pkin(2);
t28 = t21 * pkin(2);
t26 = m(3) * qJD(2);
t14 = t21 * rSges(3,1) - t20 * rSges(3,2);
t24 = qJD(2) * t14;
t4 = -t20 * t25 - t30;
t23 = m(4) * (-t4 * t9 - t32);
t22 = qJD(2) ^ 2;
t12 = t18 * t27;
t7 = qJD(1) + t24;
t6 = t18 * t15 - t12;
t2 = -t18 * t6 - t22 * t28;
t1 = -t18 * t30 - t22 * t29;
t5 = [m(4) * t1 - t10 * t26; m(4) * (t2 * (-t8 - t29) + t4 * t12 + t1 * (t9 + t28)) + 0.2e1 * (m(4) * (-t4 * t15 - t32) / 0.2e1 - t23 / 0.2e1) * t18 + 0.2e1 * (t7 * t31 - t22 * t13 * t14 / 0.2e1) * m(3) + 0.2e1 * (t14 * t31 + (t24 + t7 / 0.2e1) * t13) * t26; m(4) * (t1 * t9 - t2 * t8 - t3 * t30 - t4 * t6) - t18 * t23;];
tauc  = t5(:);
