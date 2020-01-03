% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR3
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:21
% EndTime: 2019-12-31 16:17:25
% DurationCPUTime: 2.49s
% Computational Cost: add. (2479->214), mult. (5965->321), div. (0->0), fcn. (6382->6), ass. (0->122)
t95 = cos(qJ(4));
t150 = Icges(5,4) * t95;
t94 = sin(qJ(4));
t114 = -Icges(5,2) * t94 + t150;
t145 = sin(pkin(6));
t146 = cos(pkin(6));
t173 = sin(qJ(3));
t174 = cos(qJ(3));
t75 = -t145 * t173 - t146 * t174;
t76 = -t145 * t174 + t146 * t173;
t43 = -Icges(5,6) * t76 + t114 * t75;
t151 = Icges(5,4) * t94;
t116 = Icges(5,1) * t95 - t151;
t46 = -Icges(5,5) * t76 + t116 * t75;
t119 = -t43 * t94 + t46 * t95;
t112 = Icges(5,5) * t95 - Icges(5,6) * t94;
t40 = -Icges(5,3) * t76 + t112 * t75;
t15 = t119 * t75 - t76 * t40;
t41 = -Icges(5,3) * t75 - t112 * t76;
t191 = t76 * t41;
t124 = t94 * rSges(5,1) + rSges(5,2) * t95;
t186 = qJD(4) * t124;
t189 = t75 * t186;
t136 = t76 * t186;
t113 = Icges(5,2) * t95 + t151;
t115 = Icges(5,1) * t94 + t150;
t117 = -t113 * t94 + t115 * t95;
t111 = Icges(5,5) * t94 + Icges(5,6) * t95;
t54 = t111 * t76;
t32 = t117 * t75 - t54;
t188 = qJD(3) * t32;
t22 = -t43 * t95 - t94 * t46;
t167 = t94 * rSges(5,2);
t125 = rSges(5,1) * t95 - t167;
t49 = -t76 * rSges(5,3) + t125 * t75;
t187 = t136 + (t75 * pkin(3) - t76 * pkin(5) + t49) * qJD(3);
t55 = t111 * t75;
t185 = 0.2e1 * qJD(4);
t71 = t76 * qJD(3);
t72 = qJD(3) * t75;
t182 = t72 * pkin(3) - t71 * pkin(5);
t168 = t76 * t95;
t155 = -rSges(5,1) * t168 - t75 * rSges(5,3);
t50 = t167 * t76 + t155;
t74 = t76 * pkin(3);
t181 = qJD(3) * (-pkin(5) * t75 + t50 - t74) + t189;
t47 = -Icges(5,5) * t75 - t116 * t76;
t159 = t113 * t76 + t47;
t44 = -Icges(5,6) * t75 - t114 * t76;
t161 = -t115 * t76 + t44;
t180 = t159 * t94 + t161 * t95;
t175 = -rSges(5,3) - pkin(5);
t172 = t47 * t95;
t171 = t72 * t95;
t166 = t94 * t44;
t163 = t47 * t168 + t75 * t41;
t162 = -t75 * t172 + t191;
t160 = -t115 * t75 - t43;
t158 = t113 * t75 - t46;
t156 = rSges(5,1) * t171 - t71 * rSges(5,3);
t154 = -t113 + t116;
t153 = -t114 - t115;
t152 = m(4) * qJD(3);
t144 = qJD(4) * t75;
t143 = qJD(4) * t76;
t141 = qJD(4) * t94;
t140 = qJD(4) * t95;
t31 = t117 * t76 + t55;
t139 = t31 * qJD(3);
t138 = t112 * qJD(3);
t135 = t76 * t141;
t132 = t144 / 0.2e1;
t131 = -t143 / 0.2e1;
t130 = qJD(2) * t146;
t126 = rSges(4,1) * t75 + rSges(4,2) * t76;
t93 = qJD(2) * t145;
t19 = t181 + t93;
t20 = -t130 + t187;
t123 = -t19 * t75 - t20 * t76;
t121 = t166 - t172;
t118 = -t49 * t75 + t50 * t76;
t13 = t119 * t76 + t75 * t40;
t110 = pkin(3) + t125;
t109 = t140 * t75 + t71 * t94;
t108 = t141 * t75 - t71 * t95;
t107 = t140 * t76 - t72 * t94;
t106 = t135 + t171;
t12 = t166 * t76 - t163;
t105 = (-t12 * t75 + t13 * t76) * qJD(4);
t14 = t166 * t75 + t162;
t104 = (-t14 * t75 + t15 * t76) * qJD(4);
t103 = t158 * t94 + t160 * t95;
t102 = (t153 * t94 + t154 * t95) * qJD(3);
t29 = rSges(5,1) * t108 + rSges(5,2) * t109 - t72 * rSges(5,3);
t30 = rSges(5,1) * t135 + rSges(5,2) * t107 + t156;
t101 = t29 * t75 + t30 * t76 - t49 * t71 - t50 * t72;
t78 = t114 * qJD(4);
t79 = t116 * qJD(4);
t96 = -t78 * t94 + t79 * t95 + (-t113 * t95 - t115 * t94) * qJD(4);
t80 = t125 * qJD(4);
t77 = t112 * qJD(4);
t64 = -rSges(4,1) * t76 + rSges(4,2) * t75;
t63 = t124 * t75;
t62 = t124 * t76;
t61 = rSges(4,1) * t72 + rSges(4,2) * t71;
t60 = -rSges(4,1) * t71 + rSges(4,2) * t72;
t53 = qJD(3) * t126 - t130;
t52 = qJD(3) * t64 + t93;
t24 = Icges(5,5) * t106 + Icges(5,6) * t107 - Icges(5,3) * t71;
t23 = Icges(5,5) * t108 + Icges(5,6) * t109 - Icges(5,3) * t72;
t21 = t44 * t95 + t94 * t47;
t18 = qJD(4) * t118 + qJD(1);
t17 = (t124 * t71 + t75 * t80) * qJD(4) + (t182 + t30) * qJD(3);
t16 = (-t124 * t72 + t76 * t80) * qJD(4) + (t71 * pkin(3) + t72 * pkin(5) - t29) * qJD(3);
t11 = t111 * t71 - t117 * t72 + t75 * t77 + t76 * t96;
t10 = t111 * t72 + t117 * t71 + t75 * t96 - t76 * t77;
t9 = qJD(4) * t119 - (Icges(5,4) * t108 + Icges(5,2) * t109 - Icges(5,6) * t72) * t95 - t94 * (Icges(5,1) * t108 + Icges(5,4) * t109 - Icges(5,5) * t72);
t8 = qJD(4) * t121 - (Icges(5,4) * t106 + Icges(5,2) * t107 - Icges(5,6) * t71) * t95 - t94 * (Icges(5,1) * t106 + Icges(5,4) * t107 - Icges(5,5) * t71);
t7 = t101 * qJD(4);
t6 = t104 - t188;
t5 = t105 - t139;
t1 = [m(5) * t7; m(5) * (t145 * t17 - t146 * t16) + (t145 * t61 + t146 * t60) * t152; t5 * t143 / 0.2e1 - qJD(3) * (-t117 * qJD(4) - t78 * t95 - t94 * t79) + m(4) * (t52 * t61 - t53 * t60 + (-t126 * t60 + t61 * t64) * qJD(3)) - (t126 * t52 - t53 * t64) * t152 + (t17 * (-t74 + t155) + (t16 * t175 + t167 * t17) * t76 + (-t17 * pkin(5) + t110 * t16) * t75 + (t110 * t71 - t175 * t72 + t181 - t189) * t20 + (-t167 * t72 + t136 + t156 + t182 - t187) * t19) * m(5) + (-t139 + ((t13 - t14 + t162) * t76 + t163 * t75) * qJD(4) + t9 + t10) * t131 + (t6 + ((-t12 + (t40 + t166) * t76 - t163) * t76 + (-t13 + (t121 + t40) * t75 + t191) * t75) * qJD(4) + t8 + t11 + t188) * t132 + ((-t21 + t31) * t71 + (-t22 + t32) * t72) * qJD(4) / 0.2e1; -qJD(3) * (t21 * t71 + t22 * t72 - t75 * t8 + t76 * t9) / 0.2e1 + ((t143 * t55 + t138) * t76 + (-t102 + (-t180 * t75 + (-t54 + t103) * t76) * qJD(4)) * t75) * t131 + ((t144 * t54 - t138) * t75 + (-t102 + (t103 * t76 + (-t180 - t55) * t75) * qJD(4)) * t76) * t132 + qJD(3) * (-(-t153 * t95 + t154 * t94) * qJD(3) + ((-t158 * t76 + t159 * t75) * t95 + (t160 * t76 - t161 * t75) * t94) * qJD(4)) / 0.2e1 - (-qJD(3) * t11 + (-t12 * t71 - t13 * t72 - (-t121 * t72 - t75 * t24 - t71 * t41) * t75 + (-t119 * t72 - t75 * t23 + t40 * t71) * t76) * t185) * t75 / 0.2e1 + (-qJD(3) * t10 + (-(t121 * t71 + t76 * t24 - t72 * t41) * t75 - t14 * t71 - t15 * t72 + (t119 * t71 + t76 * t23 + t40 * t72) * t76) * t185) * t76 / 0.2e1 - (t5 + t105) * t71 / 0.2e1 - (t6 + t104) * t72 / 0.2e1 + (t7 * t118 + t18 * t101 - t123 * t80 - (-t16 * t76 - t17 * t75 - t19 * t71 + t20 * t72) * t124 - (t19 * t62 - t20 * t63) * qJD(3) - (t18 * (t62 * t76 + t63 * t75) - t123 * t125) * qJD(4)) * m(5);];
tauc = t1(:);
