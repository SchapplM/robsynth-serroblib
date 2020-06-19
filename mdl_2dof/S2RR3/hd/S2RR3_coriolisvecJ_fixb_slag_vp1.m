% Calculate vector of centrifugal and Coriolis load on the joints for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S2RR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_coriolisvecJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:24
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (147->16), mult. (197->27), div. (0->0), fcn. (86->4), ass. (0->19)
t19 = cos(qJ(1));
t25 = pkin(1) * qJD(1);
t16 = qJD(1) + qJD(2);
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t15 = cos(t17);
t8 = t15 * rSges(3,1) - t14 * rSges(3,2);
t6 = t16 * t8;
t4 = t19 * t25 + t6;
t18 = sin(qJ(1));
t7 = t14 * rSges(3,1) + t15 * rSges(3,2);
t29 = t16 * t7;
t3 = -t18 * t25 - t29;
t28 = t18 * pkin(1);
t27 = t19 * pkin(1);
t20 = qJD(1) ^ 2;
t2 = -t16 * t6 - t20 * t27;
t1 = -t16 * t29 - t20 * t28;
t5 = [m(3) * (t2 * (-t7 - t28) + t1 * (t8 + t27)); (t1 * t8 - t2 * t7 - t3 * t6 - t4 * t29 - (-t3 * t8 - t4 * t7) * t16) * m(3);];
tauc = t5(:);
