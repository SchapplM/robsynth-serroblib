% Calculate matrix of centrifugal and coriolis load on the joints for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
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
% Cq [2x2]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S2RR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_coriolismatJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_coriolismatJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_coriolismatJ_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_coriolismatJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR2_coriolismatJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR2_coriolismatJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:53
% EndTime: 2019-03-08 18:00:53
% DurationCPUTime: 0.50s
% Computational Cost: add. (688->86), mult. (1660->138), div. (0->0), fcn. (1556->4), ass. (0->63)
t55 = sin(qJ(1));
t56 = cos(qJ(2));
t57 = cos(qJ(1));
t76 = t56 * t57;
t54 = sin(qJ(2));
t80 = t54 * t57;
t63 = -rSges(3,1) * t76 + rSges(3,2) * t80;
t84 = rSges(3,3) + pkin(1);
t18 = -t84 * t55 + t63;
t81 = t54 * t55;
t83 = rSges(3,1) * t56;
t62 = rSges(3,2) * t81 - t55 * t83;
t19 = t84 * t57 + t62;
t43 = rSges(3,1) * t54 + rSges(3,2) * t56;
t36 = t43 * t55;
t37 = t43 * t57;
t67 = Icges(3,2) * t56;
t71 = Icges(3,4) * t54;
t72 = Icges(3,1) * t56;
t96 = -(-t71 + t72 / 0.2e1 - t67 / 0.2e1) * t54 - m(3) * (t18 * t37 + t19 * t36);
t66 = Icges(3,2) * t57;
t70 = Icges(3,4) * t57;
t26 = Icges(3,6) * t55 - t54 * t66 + t56 * t70;
t50 = t54 * t70;
t28 = Icges(3,5) * t55 + t57 * t72 - t50;
t77 = t56 * t28;
t95 = (t54 * t26 - t77) * t57;
t53 = Icges(3,4) * t56;
t94 = Icges(3,1) * t54 + t53;
t93 = Icges(3,2) * t54 - t53;
t91 = t55 ^ 2;
t90 = t57 ^ 2;
t88 = t55 / 0.2e1;
t85 = t57 / 0.2e1;
t65 = Icges(3,6) * t57;
t25 = t93 * t55 + t65;
t82 = t54 * t25;
t79 = t55 * t57;
t49 = t55 * t71;
t69 = Icges(3,5) * t57;
t27 = -t55 * t72 + t49 + t69;
t78 = t56 * t27;
t59 = Icges(3,5) * t56 - Icges(3,6) * t54;
t23 = Icges(3,3) * t57 - t59 * t55;
t75 = t57 * t23 + t25 * t81;
t74 = t55 * t23 + t27 * t76;
t24 = Icges(3,3) * t55 - t54 * t65 + t56 * t69;
t64 = -t24 - t78;
t58 = Icges(3,5) * t54 + Icges(3,6) * t56;
t9 = t57 * t24 + t26 * t81 - t55 * t77;
t45 = -rSges(3,2) * t54 + t83;
t31 = t58 * t57;
t30 = t55 * t58;
t11 = t55 * t24 - t95;
t10 = -t25 * t80 + t74;
t8 = -t55 * t78 + t75;
t6 = (t94 / 0.2e1 - t93 / 0.2e1) * t56 - t96;
t5 = t10 * t57 + t11 * t55;
t4 = t9 * t55 + t57 * t8;
t3 = (t9 + (-t24 + t82) * t57 - t74) * t57 + (t64 * t55 + t75 - t8) * t55;
t2 = (t11 + t75 + t95) * t57 + (-t10 + (t64 - t82) * t57 + t9 + t74) * t55;
t1 = (t5 / 0.2e1 + t3 / 0.2e1) * t57 + (t2 / 0.2e1 - t4 / 0.2e1) * t55;
t7 = [t6 * qJD(2), t6 * qJD(1) + (m(3) * ((t18 * t55 - t19 * t57) * t45 + (-t36 * t57 + t37 * t55) * t43) + ((t55 * t67 + t27 + t49) * t56 + (t55 * t94 - t25) * t54) * t85 - t55 * t2 / 0.2e1 + (t91 / 0.2e1 + t90 / 0.2e1) * t59 + ((-t56 * t66 + t28 - t50) * t56 + (-t57 * t94 - t26) * t54 + t4) * t88 - (t5 + t3) * t57 / 0.2e1) * qJD(2); (-(-t93 + t94) * t56 / 0.2e1 + t96) * qJD(1) + t1 * qJD(2), t1 * qJD(1) + (m(3) * ((-t55 * t62 - t57 * t63) * (-t55 * t36 - t37 * t57) + (t90 + t91) * t45 * t43) + (t30 * t79 - t91 * t31) * t88 + (t90 * t30 - t31 * t79) * t85) * qJD(2);];
Cq  = t7;
