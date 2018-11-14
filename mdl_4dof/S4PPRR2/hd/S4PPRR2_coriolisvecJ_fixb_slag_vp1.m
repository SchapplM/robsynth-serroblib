% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 13:59
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:22
% EndTime: 2018-11-14 13:59:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (260->32), mult. (224->59), div. (0->0), fcn. (95->4), ass. (0->31)
t21 = pkin(6) + qJ(3);
t18 = sin(t21);
t19 = cos(t21);
t12 = t18 * rSges(4,1) + t19 * rSges(4,2);
t8 = qJD(3) * t12;
t32 = -t8 / 0.2e1;
t23 = qJD(3) ^ 2;
t31 = pkin(3) * t23;
t20 = qJ(4) + t21;
t16 = sin(t20);
t17 = cos(t20);
t10 = t16 * rSges(5,1) + t17 * rSges(5,2);
t15 = t17 * rSges(5,1);
t29 = t16 * rSges(5,2);
t11 = t15 - t29;
t22 = qJD(3) + qJD(4);
t26 = pkin(3) * qJD(3);
t3 = t22 * t11 + t19 * t26 + qJD(1);
t30 = t3 * t10;
t28 = t22 * t10;
t27 = m(4) * qJD(3);
t13 = t19 * rSges(4,1) - t18 * rSges(4,2);
t25 = qJD(3) * t13;
t4 = -t18 * t26 - t28;
t24 = m(5) * (-t4 * t11 - t30);
t14 = t22 * t29;
t7 = qJD(1) + t25;
t6 = t22 * t15 - t14;
t2 = -t19 * t31 - t22 * t6;
t1 = -t18 * t31 - t22 * t28;
t5 = [m(5) * t1 - t8 * t27; 0; m(5) * (t2 * (-pkin(3) * t18 - t10) + t4 * t14 + t1 * (pkin(3) * t19 + t11)) + 0.2e1 * (m(5) * (-t4 * t15 - t30) / 0.2e1 - t24 / 0.2e1) * t22 + 0.2e1 * (t7 * t32 - t23 * t12 * t13 / 0.2e1) * m(4) + 0.2e1 * (t13 * t32 + (t25 + t7 / 0.2e1) * t12) * t27; m(5) * (t1 * t11 - t2 * t10 - t28 * t3 - t4 * t6) - t22 * t24;];
tauc  = t5(:);
