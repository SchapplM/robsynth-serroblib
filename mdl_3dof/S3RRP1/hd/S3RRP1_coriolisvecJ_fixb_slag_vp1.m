% Calculate vector of centrifugal and coriolis load on the joints for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3RRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (562->55), mult. (543->70), div. (0->0), fcn. (268->4), ass. (0->44)
t40 = qJ(1) + qJ(2);
t37 = cos(t40);
t30 = t37 * qJ(3);
t36 = sin(t40);
t67 = rSges(4,1) + pkin(2);
t72 = -t37 * rSges(4,3) + t36 * t67 - t30;
t27 = qJD(3) * t36;
t39 = qJD(1) + qJD(2);
t65 = t72 * t39;
t71 = t27 - t65;
t70 = t67 * t37;
t41 = sin(qJ(1));
t58 = pkin(1) * qJD(1);
t53 = t41 * t58;
t17 = t36 * rSges(3,1) + t37 * rSges(3,2);
t59 = t39 * t17;
t9 = -t53 - t59;
t66 = rSges(4,3) + qJ(3);
t60 = t37 * t39;
t64 = rSges(4,3) * t60 + t39 * t30;
t63 = t41 * pkin(1);
t42 = cos(qJ(1));
t38 = t42 * pkin(1);
t61 = t36 * t39;
t51 = t66 * t36 + t70;
t43 = qJD(1) ^ 2;
t56 = t43 * t63;
t55 = t43 * t38;
t54 = t27 + t64;
t52 = t42 * t58;
t20 = t37 * rSges(3,1) - t36 * rSges(3,2);
t12 = rSges(3,1) * t60 - rSges(3,2) * t61;
t28 = qJD(3) * t37;
t49 = t28 - t52;
t47 = -t39 * t67 + qJD(3);
t5 = -t53 + t71;
t6 = t51 * t39 - t49;
t44 = (-t5 * t70 + (-t5 * t66 - t6 * t67) * t36) * t39;
t10 = t39 * t20 + t52;
t8 = -t39 * t12 - t55;
t7 = -t39 * t59 - t56;
t2 = -t55 + (t47 * t37 - t61 * t66 + t28) * t39;
t1 = -t56 + (t47 * t36 + t54) * t39;
t3 = [m(3) * (t8 * (-t17 - t63) + t7 * (t20 + t38) + (-t12 - t52 + t10) * t9) + (t2 * (-t72 - t63) + t5 * t49 + t1 * (t38 + t51) + t44 + (t5 + t64 + t65) * t6) * m(4); (-t2 * t72 + t44 + (t54 - t71) * t6 + (t5 * t39 + t1) * t51) * m(4) + (-t10 * t59 - t9 * t12 - t8 * t17 + t7 * t20 - (-t10 * t17 - t20 * t9) * t39) * m(3); 0.2e1 * (-t1 * t37 / 0.2e1 + t2 * t36 / 0.2e1) * m(4);];
tauc  = t3(:);
