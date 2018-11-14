% Calculate vector of centrifugal and coriolis load on the joints for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:11
% EndTime: 2018-11-14 14:02:12
% DurationCPUTime: 0.40s
% Computational Cost: add. (314->32), mult. (358->50), div. (0->0), fcn. (153->6), ass. (0->27)
t30 = qJ(2) + pkin(6);
t26 = cos(t30);
t32 = cos(qJ(2));
t49 = t32 * pkin(2) + pkin(3) * t26;
t27 = qJ(4) + t30;
t22 = sin(t27);
t23 = cos(t27);
t11 = t23 * rSges(5,1) - t22 * rSges(5,2);
t29 = qJD(2) + qJD(4);
t8 = t29 * t11;
t25 = sin(t30);
t31 = sin(qJ(2));
t44 = t31 * pkin(2);
t38 = -pkin(3) * t25 - t44;
t10 = t22 * rSges(5,1) + t23 * rSges(5,2);
t41 = t29 * t10;
t4 = t38 * qJD(2) - t41;
t33 = qJD(2) ^ 2;
t19 = t32 * rSges(3,1) - t31 * rSges(3,2);
t40 = qJD(2) * t19;
t18 = t31 * rSges(3,1) + t32 * rSges(3,2);
t15 = qJD(2) * t18;
t14 = qJD(1) + t40;
t3 = t49 * qJD(2) + qJD(1) + t8;
t2 = -t29 * t8 - t33 * t49;
t1 = -t29 * t41 + t38 * t33;
t5 = [-m(3) * qJD(2) * t15 + m(4) * t33 * (-t25 * rSges(4,1) - t26 * rSges(4,2) - t44) + m(5) * t1; (-(-qJD(2) * t14 + t19 * t33) * t18 - t14 * t15 + (-t15 * t19 + 0.2e1 * t18 * t40) * qJD(2)) * m(3) + m(5) * (t2 * (-t10 + t38) + t1 * (t11 + t49)); 0; (t1 * t11 - t2 * t10 - t3 * t41 - t4 * t8 - (-t3 * t10 - t4 * t11) * t29) * m(5);];
tauc  = t5(:);
