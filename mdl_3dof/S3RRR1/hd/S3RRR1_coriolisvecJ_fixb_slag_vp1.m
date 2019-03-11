% Calculate vector of centrifugal and Coriolis load on the joints for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3RRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:57
% DurationCPUTime: 0.27s
% Computational Cost: add. (550->45), mult. (479->61), div. (0->0), fcn. (212->6), ass. (0->40)
t37 = cos(qJ(1));
t49 = pkin(1) * qJD(1);
t45 = t37 * t49;
t35 = qJ(1) + qJ(2);
t30 = cos(t35);
t34 = qJD(1) + qJD(2);
t51 = t30 * t34;
t31 = qJ(3) + t35;
t26 = sin(t31);
t27 = cos(t31);
t15 = t27 * rSges(4,1) - t26 * rSges(4,2);
t28 = qJD(3) + t34;
t60 = t28 * t15;
t4 = pkin(2) * t51 + t45 + t60;
t36 = sin(qJ(1));
t46 = t36 * t49;
t29 = sin(t35);
t16 = t29 * rSges(3,1) + t30 * rSges(3,2);
t50 = t34 * t16;
t9 = -t46 - t50;
t56 = pkin(2) * t34 ^ 2;
t55 = t36 * pkin(1);
t32 = t37 * pkin(1);
t14 = t26 * rSges(4,1) + t27 * rSges(4,2);
t11 = t28 * t14;
t52 = t29 * t34;
t38 = qJD(1) ^ 2;
t48 = t38 * t55;
t47 = t38 * t32;
t17 = t30 * rSges(3,1) - t29 * rSges(3,2);
t13 = rSges(3,1) * t51 - rSges(3,2) * t52;
t44 = pkin(2) * t30 + t15;
t39 = -pkin(2) * t29 - t14;
t10 = t34 * t17 + t45;
t6 = -t34 * t13 - t47;
t5 = -t34 * t50 - t48;
t3 = -pkin(2) * t52 - t11 - t46;
t2 = -t28 * t60 - t30 * t56 - t47;
t1 = -t11 * t28 - t29 * t56 - t48;
t7 = [m(3) * (t6 * (-t16 - t55) + t5 * (t17 + t32) + (-t13 - t45 + t10) * t9) + (t2 * (t39 - t55) + t1 * (t32 + t44)) * m(4); (t1 * t44 + t2 * t39) * m(4) + (-t10 * t50 - t9 * t13 - t6 * t16 + t5 * t17 - (-t10 * t16 - t17 * t9) * t34) * m(3); (t1 * t15 - t2 * t14 - (-t14 * t4 - t15 * t3) * t28 - t3 * t60 - t4 * t11) * m(4);];
tauc  = t7(:);
