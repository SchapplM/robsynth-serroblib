% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:34
% EndTime: 2018-11-14 13:46:35
% DurationCPUTime: 0.69s
% Computational Cost: add. (828->94), mult. (945->103), div. (0->0), fcn. (622->6), ass. (0->62)
t55 = qJ(1) + pkin(6);
t51 = sin(t55);
t102 = (-rSges(4,1) - pkin(2)) * t51;
t52 = cos(t55);
t45 = t52 * qJ(3);
t89 = pkin(2) * t51;
t26 = -t45 + t89;
t42 = qJD(3) * t51;
t78 = qJD(1) * t52;
t82 = qJ(3) * t78 + t42;
t101 = qJD(1) * t26 - t42 + t82;
t43 = qJD(3) * t52;
t53 = cos(qJ(1)) * pkin(1);
t88 = pkin(3) * t52;
t68 = t53 + t88;
t50 = t52 * pkin(2);
t96 = t51 * qJ(3) + t50;
t100 = -(t96 + t68) * qJD(1) + t43;
t59 = qJD(1) ^ 2;
t54 = qJD(1) - qJD(4);
t56 = sin(qJ(4));
t86 = cos(qJ(4));
t62 = t51 * t56 + t52 * t86;
t73 = t51 * t86;
t94 = -t52 * t56 + t73;
t65 = -rSges(5,1) * t94 + rSges(5,2) * t62;
t99 = t54 * t65;
t74 = t53 + t96;
t48 = t52 * rSges(4,1);
t80 = t51 * rSges(4,3) + t48;
t98 = t74 + t80;
t70 = t52 * rSges(3,1) - rSges(3,2) * t51;
t95 = t53 + t70;
t92 = t51 / 0.2e1;
t91 = -t52 / 0.2e1;
t90 = sin(qJ(1)) * pkin(1);
t77 = qJD(1) * qJD(3);
t79 = qJD(1) * t51;
t83 = qJD(1) * (-pkin(2) * t79 + t82) + t51 * t77;
t76 = t59 * t90;
t75 = t59 * t53;
t47 = t52 * rSges(4,3);
t72 = -rSges(4,1) * t51 + t47 - t90;
t71 = t45 - t90;
t69 = -pkin(3) * t51 - t90;
t28 = rSges(3,1) * t51 + rSges(3,2) * t52;
t15 = t54 * t62;
t16 = -qJD(1) * t73 + qJD(4) * t94 + t56 * t78;
t6 = rSges(5,1) * t16 + rSges(5,2) * t15;
t5 = rSges(5,1) * t15 - rSges(5,2) * t16;
t64 = -rSges(5,1) * t62 - rSges(5,2) * t94;
t41 = rSges(4,3) * t78;
t36 = t52 * t77;
t20 = qJD(1) * t96 - t43;
t13 = t42 + (-t26 + t72) * qJD(1);
t10 = -t75 + t36 + (-qJD(1) * t80 - t20) * qJD(1);
t9 = -t76 + qJD(1) * (-rSges(4,1) * t79 + t41) + t83;
t8 = -t64 * t54 - t100;
t7 = t99 + t42 + (-t26 + t69) * qJD(1);
t2 = -qJD(1) * t20 - t5 * t54 - t68 * t59 + t36;
t1 = t54 * t6 + t69 * t59 + t83;
t3 = [m(3) * ((-t28 * t59 - t76) * t95 + (-t75 + (-0.2e1 * t70 - t53 + t95) * t59) * (-t28 - t90)) + (t2 * ((-pkin(2) - pkin(3)) * t51 + t65 + t71) + t1 * (-t64 + t74 + t88) + (t100 - t5) * t7 + (-t89 * qJD(1) + t101 + t6 + t7 - t99) * t8) * m(5) + (t10 * (t47 + t71 + t102) + t9 * t98 + (t43 + (-t48 - t53 - t50 + (-rSges(4,3) - qJ(3)) * t51) * qJD(1)) * t13 + (t41 + t13 + (-t72 - t90 + t102) * qJD(1) + t101) * (t98 * qJD(1) - t43)) * m(4); 0; 0.2e1 * (t1 * t91 + t2 * t92) * m(5) + 0.2e1 * (t10 * t92 + t9 * t91) * m(4); (t1 * t64 - t65 * t2 + t5 * t7 - t6 * t8 - (-t7 * t64 - t65 * t8) * t54) * m(5);];
tauc  = t3(:);
