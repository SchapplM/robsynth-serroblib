% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:27
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.52s
% Computational Cost: add. (906->59), mult. (713->79), div. (0->0), fcn. (346->6), ass. (0->43)
t48 = qJ(1) + pkin(6);
t45 = qJ(3) + t48;
t42 = cos(t45);
t33 = t42 * qJ(4);
t41 = sin(t45);
t86 = rSges(5,1) + pkin(3);
t60 = t42 * rSges(5,3) - t41 * t86 + t33;
t30 = qJD(4) * t41;
t47 = qJD(1) + qJD(3);
t89 = t47 * t60 + t30;
t88 = t86 * t42;
t65 = -pkin(2) * sin(t48) - sin(qJ(1)) * pkin(1);
t57 = t65 * qJD(1);
t5 = t57 + t89;
t84 = rSges(5,3) + qJ(4);
t81 = pkin(2) * cos(t48) + cos(qJ(1)) * pkin(1);
t31 = qJD(4) * t42;
t66 = t84 * t41 + t88;
t78 = -t66 * t47 + t31;
t51 = qJD(1) ^ 2;
t37 = t42 * rSges(4,1);
t22 = -rSges(4,2) * t41 + t37;
t56 = t81 * qJD(1);
t10 = t22 * t47 + t56;
t19 = rSges(4,1) * t41 + rSges(4,2) * t42;
t74 = t10 * t19;
t73 = t41 * t47;
t72 = t42 * t47;
t71 = t47 * t19;
t68 = rSges(5,3) * t72 + t47 * t33 + t30;
t61 = -t47 * t86 + qJD(4);
t59 = t65 * t51;
t58 = t81 * t51;
t9 = t57 - t71;
t6 = t56 - t78;
t52 = t5 * t31 + t6 * t68 + (-t5 * t88 + (-t5 * t84 - t6 * t86) * t41) * t47;
t26 = rSges(4,2) * t73;
t12 = rSges(4,1) * t72 - t26;
t8 = -t12 * t47 - t58;
t7 = -t47 * t71 + t59;
t2 = -t58 + (t61 * t42 - t73 * t84 + t31) * t47;
t1 = t59 + (t61 * t41 + t68) * t47;
t3 = [(t2 * (t60 + t65) + t1 * (t66 + t81) + (-t5 * t81 + t6 * t65) * qJD(1) + t52) * m(5) + (t8 * (-t19 + t65) + t9 * t26 + t7 * (t22 + t81) + (-t9 * t37 - t74) * t47 + (t10 * t65 - t81 * t9) * qJD(1)) * m(4); 0; (t1 * t66 + t2 * t60 - t5 * t78 - t6 * t89 + t52) * m(5) + (-t10 * t71 - t9 * t12 - t8 * t19 + t7 * t22 - (-t22 * t9 - t74) * t47) * m(4); 0.2e1 * (-t1 * t42 / 0.2e1 + t2 * t41 / 0.2e1) * m(5);];
tauc  = t3(:);
