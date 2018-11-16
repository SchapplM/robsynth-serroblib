% Calculate vector of inverse dynamics joint torques for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S2RR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR2_invdynJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:48:45
% EndTime: 2018-11-16 16:48:47
% DurationCPUTime: 1.80s
% Computational Cost: add. (1143->219), mult. (3010->318), div. (0->0), fcn. (2394->4), ass. (0->129)
t90 = cos(qJ(2));
t85 = Icges(3,4) * t90;
t88 = sin(qJ(2));
t102 = -Icges(3,2) * t88 + t85;
t91 = cos(qJ(1));
t131 = Icges(3,6) * t91;
t89 = sin(qJ(1));
t43 = -t102 * t89 + t131;
t132 = Icges(3,6) * t89;
t133 = Icges(3,2) * t91;
t138 = Icges(3,4) * t91;
t44 = -t133 * t88 + t138 * t90 + t132;
t170 = Icges(3,1) * t88 + t85;
t53 = t170 * t89;
t54 = t170 * t91;
t165 = (t44 + t54) * t89 + (t43 - t53) * t91;
t134 = Icges(3,2) * t90;
t136 = Icges(3,5) * t91;
t140 = Icges(3,1) * t90;
t139 = Icges(3,4) * t88;
t81 = t89 * t139;
t45 = -t140 * t89 + t136 + t81;
t137 = Icges(3,5) * t89;
t82 = t88 * t138;
t46 = t140 * t91 + t137 - t82;
t96 = (-t133 * t90 + t46 - t82) * t89 + (t134 * t89 + t45 + t81) * t91;
t175 = t165 * t90 + t96 * t88;
t125 = t89 * qJD(2);
t72 = rSges(3,1) * t88 + rSges(3,2) * t90;
t174 = t72 * t125;
t142 = t170 + t102;
t68 = t134 + t139;
t71 = -t139 + t140;
t143 = t68 - t71;
t173 = (t142 * t88 + t143 * t90) * qJD(1);
t92 = qJD(1) ^ 2;
t107 = t44 * t88 - t46 * t90;
t171 = t107 * t91;
t67 = Icges(3,5) * t90 - Icges(3,6) * t88;
t41 = Icges(3,3) * t91 - t67 * t89;
t129 = qJD(1) * t41;
t19 = t43 * t90 + t45 * t88;
t28 = t68 * t125 + (-t102 * t91 - t132) * qJD(1);
t30 = qJD(2) * t53 + (-t71 * t91 - t137) * qJD(1);
t169 = qJD(2) * t19 + t28 * t88 - t30 * t90 - t129;
t20 = t44 * t90 + t46 * t88;
t128 = qJD(2) * t91;
t27 = qJD(1) * t43 - t128 * t68;
t29 = -qJD(2) * t54 + (-t71 * t89 + t136) * qJD(1);
t130 = Icges(3,3) * t89;
t42 = -t131 * t88 + t136 * t90 + t130;
t168 = -qJD(1) * t42 + qJD(2) * t20 + t27 * t88 - t29 * t90;
t105 = t170 * t88 + t90 * t68;
t59 = t102 * qJD(2);
t60 = t71 * qJD(2);
t66 = Icges(3,5) * t88 + Icges(3,6) * t90;
t167 = -qJD(1) * t66 + qJD(2) * t105 + t59 * t88 - t60 * t90;
t124 = qJD(1) * qJD(2);
t64 = qJDD(2) * t89 + t124 * t91;
t164 = t64 / 0.2e1;
t65 = qJDD(2) * t91 - t124 * t89;
t163 = t65 / 0.2e1;
t162 = -rSges(3,3) - pkin(1);
t160 = rSges(3,1) * t90;
t158 = t43 * t88;
t157 = t45 * t90;
t49 = t66 * t89;
t156 = t66 * t91;
t155 = t88 * t89;
t154 = t88 * t91;
t153 = t89 * t90;
t152 = t90 * t91;
t87 = t91 * rSges(3,3);
t149 = t43 * t155 + t91 * t41;
t148 = t45 * t152 + t89 * t41;
t104 = -t170 * t90 + t88 * t68;
t22 = -t104 * t91 + t49;
t127 = t22 * qJD(1);
t126 = t67 * qJD(1);
t123 = rSges(3,1) * t152;
t84 = rSges(3,2) * t154;
t122 = qJD(1) * t84 + t174;
t121 = t72 * t128;
t47 = -rSges(3,1) * t153 + rSges(3,2) * t155 + t87;
t37 = pkin(1) * t91 + t47;
t120 = -t128 / 0.2e1;
t119 = t128 / 0.2e1;
t118 = -t125 / 0.2e1;
t117 = t125 / 0.2e1;
t116 = -t42 - t157;
t108 = -t157 + t158;
t94 = t49 * qJD(2) + (-t67 * t91 + t108 - t130) * qJD(1);
t95 = qJD(1) * t107 - qJD(2) * t156 + t129;
t115 = (-t169 * t91 + t94 * t89) * t91 + (-t168 * t91 + t95 * t89) * t89;
t114 = (t169 * t89 + t94 * t91) * t91 + (t168 * t89 + t95 * t91) * t89;
t75 = rSges(2,1) * t91 - t89 * rSges(2,2);
t113 = rSges(2,1) * t89 + rSges(2,2) * t91;
t74 = -rSges(3,2) * t88 + t160;
t13 = -t153 * t45 + t149;
t14 = -t153 * t46 + t44 * t155 + t91 * t42;
t112 = t13 * t91 + t14 * t89;
t15 = -t154 * t43 + t148;
t16 = t89 * t42 - t171;
t111 = t15 * t91 + t16 * t89;
t101 = -rSges(3,3) * t89 - t123;
t48 = -t101 - t84;
t23 = t174 + (-pkin(1) * t89 - t48) * qJD(1);
t24 = qJD(1) * t37 - t121;
t110 = t23 * t89 - t24 * t91;
t56 = t72 * t91;
t31 = -qJD(2) * t56 + (-t74 * t89 + t87) * qJD(1);
t32 = qJD(1) * t101 + t122;
t109 = t31 * t91 - t32 * t89;
t106 = -t47 * t89 + t48 * t91;
t93 = t104 * qJD(1) + t67 * qJD(2);
t61 = t74 * qJD(2);
t55 = t72 * t89;
t21 = t104 * t89 + t156;
t18 = t21 * qJD(1);
t17 = t106 * qJD(2);
t12 = -t61 * t128 + qJD(1) * t32 + qJDD(1) * t47 - t65 * t72 + (qJDD(1) * t91 - t89 * t92) * pkin(1);
t11 = t61 * t125 - qJD(1) * t31 - qJDD(1) * t48 + t64 * t72 + (-qJDD(1) * t89 - t91 * t92) * pkin(1);
t10 = t167 * t89 + t93 * t91;
t9 = -t167 * t91 + t93 * t89;
t8 = -t107 * qJD(2) + t27 * t90 + t29 * t88;
t7 = -qJD(2) * t108 + t28 * t90 + t30 * t88;
t6 = qJD(2) * t111 + t127;
t5 = qJD(2) * t112 + t18;
t1 = [(t18 + ((t149 + t16 + t171) * t91 + (-t15 + (t116 - t158) * t91 + t14 + t148) * t89) * qJD(2)) * t118 + (-t104 * qJD(2) + t59 * t90 + t60 * t88) * qJD(1) + (t22 + t20) * t164 + (t19 + t21) * t163 + (t6 - t127 + ((t14 + (-t42 + t158) * t91 - t148) * t91 + (t116 * t89 - t13 + t149) * t89) * qJD(2)) * t120 + (t7 + t10) * t119 + (Icges(2,3) + t105) * qJDD(1) + (t24 * t122 + t23 * t121 + ((-t160 * t24 + t23 * t162) * t91 + (t24 * t162 + t23 * t74) * t89) * qJD(1) + (-g(1) + t12) * t37 + (-g(3) + t11) * (t162 * t89 - t123 + t84)) * m(3) + ((qJDD(1) * t113 + g(1)) * t113 + (qJDD(1) * t75 + g(3)) * t75) * m(2) + (t9 + t8 + t5) * t117; t89 * (t9 * qJD(1) + t115 * qJD(2) + t22 * qJDD(1) + t15 * t65 + t16 * t64) / 0.2e1 + t111 * t164 + ((-t15 * t89 + t16 * t91) * qJD(1) + t115) * t117 + qJDD(1) * (t19 * t91 + t20 * t89) / 0.2e1 + t91 * (t10 * qJD(1) + t114 * qJD(2) + t21 * qJDD(1) + t13 * t65 + t14 * t64) / 0.2e1 + t112 * t163 + ((-t13 * t89 + t14 * t91) * qJD(1) + t114) * t119 + ((-t125 * t156 + t126) * t89 + (-t173 + (t89 * t49 - t175) * qJD(2)) * t91) * t118 + ((t49 * t128 + t126) * t91 + (t173 + (-t156 * t91 + t175) * qJD(2)) * t89) * t120 - (t89 * t5 + (t142 * t90 - t143 * t88) * qJD(1) + (-t165 * t88 + t90 * t96) * qJD(2)) * qJD(1) / 0.2e1 + ((-qJD(1) * t19 + t8) * t89 + (qJD(1) * t20 + t6 + t7) * t91) * qJD(1) / 0.2e1 + ((qJD(2) * t109 - t47 * t64 + t48 * t65) * t106 + t17 * ((-t47 * t91 - t48 * t89) * qJD(1) + t109) + t110 * t61 + (t11 * t89 - t12 * t91 + (t23 * t91 + t24 * t89) * qJD(1)) * t72 - (t23 * t56 + t24 * t55) * qJD(1) - (t17 * (-t55 * t89 - t56 * t91) + t110 * t74) * qJD(2) + g(1) * t56 - g(2) * t74 - g(3) * t55) * m(3);];
tau  = t1;
