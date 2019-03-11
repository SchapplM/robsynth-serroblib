% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:14
% EndTime: 2019-03-08 18:12:14
% DurationCPUTime: 0.22s
% Computational Cost: add. (267->41), mult. (628->63), div. (0->0), fcn. (568->4), ass. (0->32)
t44 = rSges(5,1) + pkin(3);
t39 = -rSges(5,3) - qJ(4);
t37 = sin(pkin(5));
t38 = cos(pkin(5));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t24 = -t37 * t42 + t38 * t41;
t18 = qJD(4) * t24;
t23 = -t37 * t41 - t38 * t42;
t43 = t39 * t23 - t44 * t24;
t46 = t43 * qJD(3) + t18;
t19 = t23 * qJD(4);
t32 = t44 * t23 + t39 * t24;
t45 = -t32 * qJD(3) + t19;
t40 = m(4) * qJD(3);
t36 = qJD(2) * t38;
t33 = t23 * rSges(4,1) + t24 * rSges(4,2);
t16 = t24 * qJD(3);
t17 = qJD(3) * t23;
t31 = t39 * t16 + t44 * t17 - t19;
t30 = t44 * t16 - t39 * t17 - t18;
t29 = qJD(2) * t37;
t11 = -t24 * rSges(4,1) + t23 * rSges(4,2);
t8 = t17 * rSges(4,1) + t16 * rSges(4,2);
t7 = -t16 * rSges(4,1) + t17 * rSges(4,2);
t6 = qJD(3) * t33 - t36;
t5 = qJD(3) * t11 + t29;
t4 = -t36 - t45;
t3 = t29 + t46;
t2 = t31 * qJD(3) - qJD(4) * t17;
t1 = t30 * qJD(3) - qJD(4) * t16;
t9 = [0; m(5) * (-t1 * t38 + t2 * t37) + (t37 * t8 + t38 * t7) * t40; -(-t6 * t11 + t5 * t33) * t40 + m(4) * (t5 * t8 - t6 * t7 + (t11 * t8 - t33 * t7) * qJD(3)) + (t1 * t32 + t2 * t43 + (t30 + t46) * t4 + (t31 + t45) * t3) * m(5); (-t1 * t23 - t4 * t16 - t3 * t17 + t2 * t24 - (-t23 * t3 - t24 * t4) * qJD(3)) * m(5);];
tauc  = t9(:);
