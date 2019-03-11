% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:47
% EndTime: 2019-03-08 18:30:47
% DurationCPUTime: 0.73s
% Computational Cost: add. (788->126), mult. (1567->151), div. (0->0), fcn. (1150->4), ass. (0->79)
t73 = sin(qJ(3));
t75 = cos(qJ(1));
t110 = t73 * t75;
t111 = cos(qJ(3));
t74 = sin(qJ(1));
t97 = t74 * t111;
t38 = -t97 + t110;
t109 = t74 * t73;
t96 = t75 * t111;
t84 = t96 + t109;
t129 = -t38 * rSges(5,1) - rSges(5,2) * t84 - pkin(3) * t110;
t72 = qJD(1) - qJD(3);
t88 = rSges(4,1) * t38 + rSges(4,2) * t84;
t128 = t72 * t88;
t61 = t111 * pkin(3) + pkin(2);
t127 = -t74 * t61 - t129;
t115 = pkin(2) * t74;
t108 = t115 + t127;
t126 = t108 * t72;
t102 = qJD(1) * t75;
t78 = t84 * qJD(3);
t23 = t84 * qJD(1) - t78;
t124 = t38 * qJD(3) - t73 * t102;
t24 = -qJD(1) * t97 - t124;
t54 = pkin(3) * t109;
t95 = -t23 * rSges(5,1) + t24 * rSges(5,2) - qJD(1) * t54 - t61 * t102;
t125 = -pkin(2) * t102 - pkin(3) * t78 - t95;
t71 = t75 * pkin(1);
t45 = t74 * qJ(2) + t71;
t70 = t75 * pkin(2);
t122 = t45 + t70;
t104 = t75 * rSges(3,1) + t74 * rSges(3,3);
t121 = t45 + t104;
t120 = rSges(5,1) * t84 - t38 * rSges(5,2) + t54;
t103 = qJD(1) * t74;
t77 = -t24 * rSges(5,1) - t23 * rSges(5,2) + t124 * pkin(3);
t119 = (pkin(2) - t61) * t103 - t77;
t118 = t74 / 0.2e1;
t117 = -t75 / 0.2e1;
t116 = -pkin(1) - pkin(2);
t114 = qJD(1) ^ 2 * pkin(2);
t113 = -rSges(3,1) - pkin(1);
t101 = qJD(1) * qJD(2);
t63 = qJD(2) * t74;
t106 = qJ(2) * t102 + t63;
t107 = qJD(1) * (-pkin(1) * t103 + t106) + t74 * t101;
t66 = t75 * qJ(2);
t42 = pkin(1) * t74 - t66;
t105 = -qJD(1) * t42 + t63;
t94 = t75 * t61 + t120;
t93 = pkin(3) * t96;
t90 = t70 - t94;
t10 = rSges(4,1) * t24 + rSges(4,2) * t23;
t9 = rSges(4,1) * t23 - rSges(4,2) * t24;
t87 = -rSges(4,1) * t84 + rSges(4,2) * t38;
t86 = -t74 * t114 + t107;
t85 = -pkin(2) * t103 + t105;
t83 = t63 + (-t42 - t115) * qJD(1);
t64 = qJD(2) * t75;
t82 = t122 * qJD(1) - t64;
t32 = qJD(1) * t45 - t64;
t56 = t75 * t101;
t80 = -qJD(1) * t32 - t75 * t114 + t56;
t68 = t75 * rSges(3,3);
t62 = rSges(3,3) * t102;
t43 = rSges(3,1) * t74 - t68;
t29 = t121 * qJD(1) - t64;
t28 = t63 + (-t42 - t43) * qJD(1);
t16 = t56 + (-t104 * qJD(1) - t32) * qJD(1);
t15 = qJD(1) * (-rSges(3,1) * t103 + t62) + t107;
t14 = -t87 * t72 + t82;
t13 = t83 + t128;
t8 = -t90 * t72 + t82;
t7 = t83 + t126;
t4 = -t72 * t9 + t80;
t3 = t10 * t72 + t86;
t2 = -t125 * t72 + t80;
t1 = t119 * t72 + t86;
t5 = [(t2 * (-t42 + t127) + t7 * (t64 + t95 + (t54 + t93) * qJD(3)) + t1 * (t94 + t45) + t8 * (-t77 + t106) + (-t7 * t71 + (-t7 * qJ(2) + t8 * (-pkin(1) - t61)) * t74) * qJD(1) - (-t7 + t85 + t126) * t8) * m(5) + (t4 * (t116 * t74 + t66 + t88) + t13 * (t64 - t9) + t3 * (-t87 + t122) + t14 * (t10 + t106) + (t13 * t116 * t75 + (-t13 * qJ(2) + t14 * t116) * t74) * qJD(1) - (-t13 + t85 + t128) * t14) * m(4) + (t16 * (t113 * t74 + t66 + t68) + t28 * t64 + t15 * t121 + t29 * (t62 + t106) + (t28 * t113 * t75 + (t28 * (-rSges(3,3) - qJ(2)) + t29 * t113) * t74) * qJD(1) - (-qJD(1) * t43 + t105 - t28) * t29) * m(3); 0.2e1 * (t1 * t117 + t2 * t118) * m(5) + 0.2e1 * (t3 * t117 + t4 * t118) * m(4) + 0.2e1 * (t15 * t117 + t16 * t118) * m(3); (-(t7 * (t93 + t120) + t8 * (pkin(3) * t97 + t129)) * t72 - t2 * t108 + t7 * t125 + t1 * t90 - t8 * t119) * m(5) + (-(-t13 * t87 - t88 * t14) * t72 - t10 * t14 + t13 * t9 - t88 * t4 + t87 * t3) * m(4); 0;];
tauc  = t5(:);
