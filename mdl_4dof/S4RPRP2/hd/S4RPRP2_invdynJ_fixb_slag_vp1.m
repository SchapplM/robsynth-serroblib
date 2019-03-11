% Calculate vector of inverse dynamics joint torques for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:43
% DurationCPUTime: 1.12s
% Computational Cost: add. (1034->163), mult. (1923->180), div. (0->0), fcn. (1460->4), ass. (0->89)
t144 = cos(qJ(3));
t96 = cos(qJ(1));
t117 = t96 * t144;
t94 = sin(qJ(3));
t95 = sin(qJ(1));
t138 = t95 * t94;
t107 = t117 + t138;
t139 = t94 * t96;
t118 = t95 * t144;
t51 = -t118 + t139;
t168 = -rSges(5,1) * t51 - rSges(5,2) * t107 - pkin(3) * t139;
t101 = t107 * qJD(3);
t25 = qJD(1) * t107 - t101;
t124 = qJD(1) * t96;
t157 = qJD(3) * t51 - t94 * t124;
t26 = -qJD(1) * t118 - t157;
t11 = rSges(4,1) * t25 - rSges(4,2) * t26;
t132 = rSges(4,1) * t51 + rSges(4,2) * t107;
t122 = qJD(1) * qJD(2);
t128 = qJDD(2) * t95 + t122 * t96;
t91 = t96 * pkin(1);
t60 = t95 * qJ(2) + t91;
t84 = qJD(2) * t96;
t39 = qJD(1) * t60 - t84;
t86 = t96 * qJ(2);
t57 = pkin(1) * t95 - t86;
t97 = qJD(1) ^ 2;
t151 = -qJD(1) * t39 - qJDD(1) * t57 + t128 + (-qJDD(1) * t95 - t96 * t97) * pkin(2);
t92 = qJDD(1) - qJDD(3);
t93 = qJD(1) - qJD(3);
t3 = -t11 * t93 + t132 * t92 + t151;
t167 = -g(1) + t3;
t79 = pkin(3) * t144 + pkin(2);
t166 = -t79 * t95 - t168;
t12 = rSges(4,1) * t26 + rSges(4,2) * t25;
t27 = -rSges(4,1) * t107 + rSges(4,2) * t51;
t125 = qJD(1) * t95;
t83 = qJD(2) * t95;
t127 = qJ(2) * t124 + t83;
t108 = -qJDD(2) * t96 + qJD(1) * (-pkin(1) * t125 + t127) + qJDD(1) * t60 + t95 * t122;
t147 = pkin(2) * t95;
t90 = t96 * pkin(2);
t99 = qJDD(1) * t90 - t147 * t97 + t108;
t4 = t12 * t93 - t27 * t92 + t99;
t165 = -g(2) + t4;
t135 = t147 + t166;
t164 = t135 * t93;
t163 = t93 * t132;
t160 = t27 * t93;
t30 = rSges(5,1) * t107 - rSges(5,2) * t51;
t71 = pkin(3) * t138;
t116 = -t79 * t96 - t30 - t71;
t159 = t90 + t116;
t158 = -rSges(5,1) * t25 + rSges(5,2) * t26 - qJD(1) * t71 - t124 * t79;
t156 = t60 + t90;
t61 = rSges(3,1) * t96 + rSges(3,3) * t95;
t155 = pkin(3) * t118 + t168;
t154 = pkin(3) * t117 + t71;
t153 = -pkin(2) * t124 - pkin(3) * t101 - t158;
t100 = -t26 * rSges(5,1) - t25 * rSges(5,2) + pkin(3) * t157;
t152 = (pkin(2) - t79) * t125 - t100;
t150 = t95 / 0.2e1;
t149 = -t96 / 0.2e1;
t148 = -pkin(1) - pkin(2);
t146 = -rSges(3,1) - pkin(1);
t88 = t96 * rSges(3,3);
t58 = rSges(3,1) * t95 - t88;
t130 = -t57 - t58;
t38 = t60 + t61;
t126 = -qJD(1) * t57 + t83;
t62 = rSges(2,1) * t96 - rSges(2,2) * t95;
t59 = rSges(2,1) * t95 + rSges(2,2) * t96;
t110 = (-Icges(4,3) - Icges(5,3)) * t92;
t109 = -pkin(2) * t125 + t126;
t106 = t83 + (-t57 - t147) * qJD(1);
t105 = qJD(1) * t156 - t84;
t104 = t30 + t154;
t81 = rSges(3,3) * t124;
t33 = qJD(1) * t38 - t84;
t32 = qJD(1) * t130 + t83;
t16 = t105 - t160;
t15 = t106 + t163;
t10 = qJDD(1) * t61 + qJD(1) * (-rSges(3,1) * t125 + t81) + t108;
t9 = t130 * qJDD(1) + (-qJD(1) * t61 - t39) * qJD(1) + t128;
t8 = -t159 * t93 + t105;
t7 = t106 + t164;
t2 = t152 * t93 - t159 * t92 + t99;
t1 = t135 * t92 - t153 * t93 + t151;
t5 = [-m(2) * (-g(1) * t59 + g(2) * t62) - t110 + (Icges(2,3) + Icges(3,2) + m(2) * (t59 ^ 2 + t62 ^ 2)) * qJDD(1) + (-(-t7 + t109 + t164) * t8 + t7 * (qJD(3) * t154 + t158 + t84) + t8 * (-t100 + t127) + (-t7 * t91 + (-t7 * qJ(2) + t8 * (-pkin(1) - t79)) * t95) * qJD(1) + (-g(2) + t2) * (-t116 + t60) + (-g(1) + t1) * (-t57 + t166)) * m(5) + (-(t109 - t15 + t163) * t16 + t15 * (-t11 + t84) + t16 * (t12 + t127) + (t15 * t148 * t96 + (-qJ(2) * t15 + t148 * t16) * t95) * qJD(1) + t165 * (-t27 + t156) + t167 * (t148 * t95 + t132 + t86)) * m(4) + (-(-qJD(1) * t58 + t126 - t32) * t33 + t32 * t84 + t33 * (t81 + t127) + (t32 * t146 * t96 + (t32 * (-rSges(3,3) - qJ(2)) + t33 * t146) * t95) * qJD(1) + (-g(2) + t10) * t38 + (-g(1) + t9) * (t146 * t95 + t86 + t88)) * m(3); (-m(3) - m(4) - m(5)) * (g(1) * t95 - g(2) * t96) + 0.2e1 * (t1 * t150 + t149 * t2) * m(5) + 0.2e1 * (t149 * t4 + t150 * t3) * m(4) + 0.2e1 * (t10 * t149 + t150 * t9) * m(3); t110 + (-g(1) * t155 + g(2) * t104 - t1 * t135 + t7 * t153 + t2 * t159 - t8 * t152 - (t7 * t104 + t155 * t8) * t93) * m(5) + (-t12 * t16 + (t11 + t160) * t15 - (-t16 * t93 + t167) * t132 + t165 * t27) * m(4); (g(3) + qJDD(4)) * m(5);];
tau  = t5;
