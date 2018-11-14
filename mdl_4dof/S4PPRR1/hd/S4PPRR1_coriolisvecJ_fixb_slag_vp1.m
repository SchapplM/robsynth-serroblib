% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:08
% EndTime: 2018-11-14 13:40:09
% DurationCPUTime: 0.31s
% Computational Cost: add. (373->50), mult. (526->88), div. (0->0), fcn. (426->6), ass. (0->42)
t37 = qJD(3) + qJD(4);
t49 = qJ(3) + qJ(4);
t36 = sin(t49);
t38 = sin(pkin(6));
t39 = cos(pkin(6));
t47 = cos(t49);
t22 = -t38 * t36 - t39 * t47;
t23 = t39 * t36 - t38 * t47;
t9 = t22 * rSges(5,1) + t23 * rSges(5,2);
t59 = t37 * t9;
t41 = cos(qJ(3));
t57 = t41 * pkin(3);
t35 = qJD(2) * t38;
t48 = t57 * t38;
t40 = sin(qJ(3));
t51 = t39 * t40;
t10 = -t23 * rSges(5,1) + t22 * rSges(5,2);
t54 = t37 * t10;
t3 = t35 + qJD(3) * (-pkin(3) * t51 + t48) + t54;
t53 = t38 * t40;
t4 = -qJD(3) * pkin(3) * t53 + t59 + (-qJD(3) * t57 - qJD(2)) * t39;
t18 = t23 * t37;
t19 = t22 * t37;
t44 = t18 * rSges(5,1) - t19 * rSges(5,2);
t6 = t19 * rSges(5,1) + t18 * rSges(5,2);
t58 = t3 * t6 + t4 * t44;
t56 = qJD(3) ^ 2;
t50 = m(4) * qJD(3);
t26 = -t39 * t41 - t53;
t43 = t38 * t41 - t51;
t45 = t26 * rSges(4,1) - rSges(4,2) * t43;
t42 = t43 * pkin(3);
t25 = qJD(3) * t26;
t24 = t43 * qJD(3);
t14 = rSges(4,1) * t43 + t26 * rSges(4,2);
t13 = t25 * rSges(4,1) - t24 * rSges(4,2);
t12 = t24 * rSges(4,1) + t25 * rSges(4,2);
t8 = -qJD(2) * t39 + qJD(3) * t45;
t7 = qJD(3) * t14 + t35;
t2 = t37 * t44 - t56 * t42;
t1 = t56 * t26 * pkin(3) + t37 * t6;
t5 = [0; m(5) * (t1 * t38 - t2 * t39) + (t12 * t39 + t13 * t38) * t50; -(-t8 * t14 + t7 * t45) * t50 + m(4) * (-t8 * t12 + t7 * t13 + (-t12 * t45 + t13 * t14) * qJD(3)) + (-t3 * (pkin(3) * t25 + t59) - t4 * (-qJD(3) * t42 - t54) + t1 * (t48 + t10) + t2 * (-t57 * t39 + t9) + ((-t1 * t39 - t2 * t38) * t40 + (t3 * t26 - t4 * t43) * qJD(3)) * pkin(3) + t58) * m(5); (t1 * t10 + t2 * t9 - (-t4 * t10 + t3 * t9) * t37 + t58) * m(5);];
tauc  = t5(:);
