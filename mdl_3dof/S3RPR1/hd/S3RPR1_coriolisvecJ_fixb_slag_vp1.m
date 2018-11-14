% Calculate vector of centrifugal and coriolis load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3RPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:22
% EndTime: 2018-11-14 10:14:23
% DurationCPUTime: 0.37s
% Computational Cost: add. (368->76), mult. (775->101), div. (0->0), fcn. (544->4), ass. (0->50)
t47 = qJD(1) - qJD(3);
t48 = sin(qJ(3));
t49 = sin(qJ(1));
t50 = cos(qJ(1));
t66 = cos(qJ(3));
t53 = t49 * t48 + t50 * t66;
t57 = t49 * t66;
t73 = -t50 * t48 + t57;
t55 = -rSges(4,1) * t73 + rSges(4,2) * t53;
t76 = t47 * t55;
t28 = t50 * pkin(1) + t49 * qJ(2);
t61 = t50 * rSges(3,1) + t49 * rSges(3,3);
t75 = t28 + t61;
t68 = t50 * pkin(2);
t74 = t28 + t68;
t72 = t49 / 0.2e1;
t71 = -t50 / 0.2e1;
t70 = -pkin(1) - pkin(2);
t69 = t49 * pkin(2);
t67 = -rSges(3,1) - pkin(1);
t58 = qJD(1) * qJD(2);
t60 = qJD(1) * t49;
t39 = qJD(2) * t49;
t59 = qJD(1) * t50;
t63 = qJ(2) * t59 + t39;
t64 = qJD(1) * (-pkin(1) * t60 + t63) + t49 * t58;
t42 = t50 * qJ(2);
t25 = t49 * pkin(1) - t42;
t62 = -qJD(1) * t25 + t39;
t13 = t47 * t53;
t14 = -qJD(1) * t57 + qJD(3) * t73 + t48 * t59;
t4 = t14 * rSges(4,1) + t13 * rSges(4,2);
t3 = t13 * rSges(4,1) - t14 * rSges(4,2);
t54 = -rSges(4,1) * t53 - rSges(4,2) * t73;
t51 = qJD(1) ^ 2;
t44 = t50 * rSges(3,3);
t40 = qJD(2) * t50;
t38 = rSges(3,3) * t59;
t33 = t50 * t58;
t26 = t49 * rSges(3,1) - t44;
t20 = qJD(1) * t28 - t40;
t18 = qJD(1) * t75 - t40;
t17 = t39 + (-t25 - t26) * qJD(1);
t10 = t33 + (-qJD(1) * t61 - t20) * qJD(1);
t9 = qJD(1) * (-rSges(3,1) * t60 + t38) + t64;
t8 = qJD(1) * t74 - t54 * t47 - t40;
t7 = t76 + t39 + (-t25 - t69) * qJD(1);
t2 = -qJD(1) * t20 - t47 * t3 - t51 * t68 + t33;
t1 = t47 * t4 - t51 * t69 + t64;
t5 = [(-(-pkin(2) * t60 + t62 - t7 + t76) * t8 + t2 * (t70 * t49 + t42 + t55) + t7 * (-t3 + t40) + t1 * (-t54 + t74) + t8 * (t4 + t63) + (t7 * t70 * t50 + (-t7 * qJ(2) + t8 * t70) * t49) * qJD(1)) * m(4) + (-(-qJD(1) * t26 - t17 + t62) * t18 + t10 * (t67 * t49 + t42 + t44) + t17 * t40 + t9 * t75 + t18 * (t38 + t63) + (t17 * t67 * t50 + (t17 * (-rSges(3,3) - qJ(2)) + t18 * t67) * t49) * qJD(1)) * m(3); 0.2e1 * (t1 * t71 + t2 * t72) * m(4) + 0.2e1 * (t10 * t72 + t9 * t71) * m(3); (t1 * t54 - t2 * t55 + t7 * t3 - t8 * t4 - (-t7 * t54 - t55 * t8) * t47) * m(4);];
tauc  = t5(:);
