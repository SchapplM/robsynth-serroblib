% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:26
% EndTime: 2019-03-08 18:28:27
% DurationCPUTime: 0.62s
% Computational Cost: add. (756->121), mult. (1233->148), div. (0->0), fcn. (880->6), ass. (0->79)
t71 = qJD(1) - qJD(4);
t95 = pkin(6) + qJ(4);
t62 = sin(t95);
t74 = cos(qJ(1));
t73 = sin(qJ(1));
t89 = cos(t95);
t80 = t73 * t89;
t30 = t62 * t74 - t80;
t76 = t73 * t62 + t74 * t89;
t83 = rSges(5,1) * t30 + rSges(5,2) * t76;
t120 = t71 * t83;
t112 = pkin(2) * t74;
t72 = sin(pkin(6));
t107 = t72 * t74;
t99 = cos(pkin(6));
t90 = t73 * t99;
t37 = -t90 + t107;
t45 = t74 * pkin(1) + t73 * qJ(2);
t106 = t73 * t72;
t79 = t74 * t99 + t106;
t119 = rSges(4,1) * t79 - t37 * rSges(4,2) + t112 + t45;
t57 = t99 * pkin(3) + pkin(2);
t103 = pkin(3) * t106 + t74 * t57;
t100 = t74 * rSges(3,1) + t73 * rSges(3,3);
t118 = t45 + t100;
t117 = t45 + t103;
t116 = t73 / 0.2e1;
t115 = -t74 / 0.2e1;
t114 = -pkin(1) - pkin(2);
t113 = pkin(2) * t73;
t111 = qJD(1) ^ 2 * pkin(2);
t110 = -rSges(3,1) - pkin(1);
t109 = -pkin(1) - t57;
t63 = qJD(2) * t73;
t97 = qJD(1) * t74;
t102 = qJ(2) * t97 + t63;
t96 = qJD(1) * qJD(2);
t98 = qJD(1) * t73;
t105 = qJD(1) * (-pkin(1) * t98 + t102) + t73 * t96;
t32 = t79 * qJD(1);
t93 = t72 * t97;
t33 = -qJD(1) * t90 + t93;
t104 = t33 * rSges(4,1) + t32 * rSges(4,2);
t66 = t74 * qJ(2);
t42 = pkin(1) * t73 - t66;
t101 = -qJD(1) * t42 + t63;
t94 = pkin(3) * t107;
t92 = (pkin(2) - t57) * t73;
t91 = t92 + t94 - t113;
t59 = t74 * t96;
t88 = -t74 * t111 + t59;
t85 = -rSges(4,1) * t32 + rSges(4,2) * t33;
t84 = rSges(4,1) * t37 + rSges(4,2) * t79;
t17 = t71 * t76;
t18 = -qJD(1) * t80 - t30 * qJD(4) + t62 * t97;
t8 = rSges(5,1) * t18 + rSges(5,2) * t17;
t7 = rSges(5,1) * t17 - rSges(5,2) * t18;
t82 = -rSges(5,1) * t76 + rSges(5,2) * t30;
t81 = -t73 * t111 + t105;
t77 = t84 - t113;
t68 = t74 * rSges(3,3);
t64 = qJD(2) * t74;
t61 = rSges(3,3) * t97;
t48 = pkin(3) * t93;
t43 = rSges(3,1) * t73 - t68;
t31 = qJD(1) * t45 - t64;
t24 = t118 * qJD(1) - t64;
t23 = t63 + (-t42 - t43) * qJD(1);
t16 = t59 + (-t100 * qJD(1) - t31) * qJD(1);
t15 = qJD(1) * (-rSges(3,1) * t98 + t61) + t105;
t14 = t119 * qJD(1) - t64;
t13 = t63 + (-t42 + t77) * qJD(1);
t10 = (-t31 + t85) * qJD(1) + t88;
t9 = qJD(1) * t104 + t81;
t6 = t117 * qJD(1) - t82 * t71 - t64;
t5 = t120 + t63 + (-t42 + t91) * qJD(1);
t2 = -t71 * t7 + (-t31 + (t112 - t103) * qJD(1)) * qJD(1) + t88;
t1 = t71 * t8 + (qJD(1) * t92 + t48) * qJD(1) + t81;
t3 = [(t2 * (t109 * t73 + t66 + t83 + t94) + t5 * (t64 - t7) + t1 * (-t82 + t117) + t6 * (t48 + t8 + t102) + (t5 * t109 * t74 + (t5 * (-pkin(3) * t72 - qJ(2)) + t6 * t109) * t73) * qJD(1) - (t91 * qJD(1) + t101 + t120 - t5) * t6) * m(5) + (t10 * (t114 * t73 + t66 + t84) + t13 * (t64 + t85) + t9 * t119 + t14 * (t102 + t104) + (t13 * t114 * t74 + (-t13 * qJ(2) + t14 * t114) * t73) * qJD(1) - (t77 * qJD(1) + t101 - t13) * t14) * m(4) + (t16 * (t110 * t73 + t66 + t68) + t23 * t64 + t15 * t118 + t24 * (t61 + t102) + (t23 * t110 * t74 + (t23 * (-rSges(3,3) - qJ(2)) + t24 * t110) * t73) * qJD(1) - (-qJD(1) * t43 + t101 - t23) * t24) * m(3); 0.2e1 * (t1 * t115 + t2 * t116) * m(5) + 0.2e1 * (t10 * t116 + t9 * t115) * m(4) + 0.2e1 * (t15 * t115 + t16 * t116) * m(3); 0; (t1 * t82 - t2 * t83 + t5 * t7 - t6 * t8 - (-t5 * t82 - t6 * t83) * t71) * m(5);];
tauc  = t3(:);
