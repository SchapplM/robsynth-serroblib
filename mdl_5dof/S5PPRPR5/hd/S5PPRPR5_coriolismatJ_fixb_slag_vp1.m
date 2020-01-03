% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:27
% DurationCPUTime: 1.12s
% Computational Cost: add. (2393->147), mult. (5272->240), div. (0->0), fcn. (6236->6), ass. (0->104)
t101 = sin(qJ(5));
t102 = cos(qJ(5));
t108 = Icges(6,5) * t101 + Icges(6,6) * t102;
t122 = sin(pkin(7));
t123 = cos(pkin(7));
t136 = sin(qJ(3));
t137 = cos(qJ(3));
t84 = -t122 * t136 - t123 * t137;
t85 = -t122 * t137 + t123 * t136;
t52 = -Icges(6,3) * t84 + t108 * t85;
t134 = t84 * t52;
t100 = Icges(6,4) * t101;
t154 = Icges(6,2) * t102 + t100;
t55 = -Icges(6,6) * t84 + t154 * t85;
t159 = t102 * t55;
t57 = -Icges(6,6) * t85 - t154 * t84;
t158 = t102 * t57;
t95 = t101 * rSges(6,1) + rSges(6,2) * t102;
t157 = t95 * t84;
t156 = t95 * t85;
t155 = rSges(6,3) + pkin(3) + pkin(6);
t94 = -Icges(6,1) * t102 + t100;
t119 = Icges(6,4) * t102;
t121 = Icges(6,1) * t101;
t93 = t119 + t121;
t153 = -t93 / 0.2e1 + Icges(6,2) * t101 / 0.2e1 - t119 / 0.2e1;
t83 = t84 ^ 2;
t152 = t85 ^ 2;
t151 = 2 * qJD(3);
t150 = m(5) / 0.2e1;
t149 = m(6) / 0.2e1;
t78 = t84 * qJ(4);
t29 = -t155 * t85 - t157 - t78;
t79 = t85 * qJ(4);
t31 = t155 * t84 - t156 - t79;
t124 = t102 * t85;
t126 = t101 * t85;
t68 = rSges(6,1) * t124 - rSges(6,2) * t126;
t96 = -rSges(6,1) * t102 + t101 * rSges(6,2);
t69 = t96 * t84;
t148 = m(6) * (t29 * t69 - t31 * t68);
t80 = t84 * rSges(6,3);
t114 = rSges(6,1) * t126 + rSges(6,2) * t124 - t80;
t61 = -t80 + t156;
t21 = (-t61 + t114) * t84;
t147 = -t21 / 0.2e1;
t146 = -t84 / 0.2e1;
t139 = rSges(5,2) - pkin(3);
t48 = -rSges(5,3) * t84 + t139 * t85 - t78;
t37 = t84 * t48;
t49 = -rSges(5,3) * t85 - t139 * t84 - t79;
t17 = -t49 * t85 - t37;
t143 = m(5) * t17;
t26 = t84 * t29;
t10 = -t31 * t85 - t26;
t142 = m(6) * t10;
t141 = m(6) * t96;
t125 = t102 * t84;
t127 = t101 * t84;
t130 = Icges(6,5) * t84;
t73 = t85 * t119;
t58 = t85 * t121 - t130 + t73;
t13 = -t55 * t125 - t58 * t127 - t85 * t52;
t135 = t13 * t84;
t133 = t85 * t84;
t81 = t84 * pkin(3);
t30 = t84 * pkin(6) - t114 - t79 + t81;
t132 = t30 - t31;
t131 = m(6) * qJD(5);
t116 = t83 + t152;
t54 = -Icges(6,3) * t85 - t108 * t84;
t74 = t84 * t119;
t60 = -Icges(6,5) * t85 - t84 * t121 - t74;
t115 = t57 * t125 + t60 * t127 + t85 * t54;
t113 = t95 * t131;
t112 = t21 * t149;
t109 = Icges(6,5) * t102 - Icges(6,6) * t101;
t106 = m(6) * (t68 * t84 + t69 * t85);
t15 = -t106 / 0.2e1;
t105 = t84 * t122 - t85 * t123;
t103 = t105 * t141 / 0.2e1;
t104 = m(6) * (t69 * t122 + t68 * t123);
t18 = -t104 / 0.2e1 + t103;
t107 = -t18 * qJD(2) - t15 * qJD(4);
t63 = t109 * t84;
t62 = t109 * t85;
t59 = -t93 * t85 + t130;
t47 = -rSges(5,2) * t84 + t81 + (-rSges(5,3) - qJ(4)) * t85;
t34 = -0.2e1 * (t150 + t149) * t105;
t28 = t68 * t85 - t69 * t84;
t22 = t85 * t61 + (t85 * rSges(6,3) + t157) * t84;
t20 = qJD(3) * t112;
t19 = t104 / 0.2e1 + t103;
t16 = t106 / 0.2e1;
t11 = -t134 + (t101 * t58 + t159) * t85;
t8 = t148 + t153 * t102 + (t94 / 0.2e1 + t154 / 0.2e1) * t101;
t7 = t22 * t21;
t6 = t142 + t143;
t5 = t115 * t85 - t135;
t4 = -t11 * t84 - (t57 * t124 + t60 * t126 - t84 * t54) * t85;
t3 = -t135 + (-t134 - t11 + (-t101 * t59 + t159) * t85 + t115) * t85;
t2 = -t134 * t84 + (-t13 + (t101 * t60 + t158 - t52) * t85 + (-t54 + (t58 + t59) * t101) * t84) * t85;
t1 = m(6) * t7 + (-t4 / 0.2e1 - t2 / 0.2e1) * t85 + (-t3 / 0.2e1 + t5 / 0.2e1) * t84;
t9 = [0, 0, qJD(5) * t112, 0, t28 * t131 + t20; 0, 0, t34 * qJD(4) + t19 * qJD(5) + (m(4) * ((rSges(4,1) * t84 + rSges(4,2) * t85) * t122 + (-rSges(4,1) * t85 + rSges(4,2) * t84) * t123) / 0.2e1 + (t47 * t122 + t48 * t123) * t150 + (t30 * t122 + t29 * t123) * t149) * t151, t34 * qJD(3), t19 * qJD(3) + (-t122 * t85 - t123 * t84) * t113; t131 * t147, -t18 * qJD(5), (m(5) * (t47 - t49) * t48 + m(6) * t132 * t29) * qJD(3) + t6 * qJD(4) + t8 * qJD(5), t6 * qJD(3) - t15 * qJD(5), qJD(1) * t147 * m(6) + t8 * qJD(3) + t107 + (t5 * t146 + ((-t29 * t95 - t69 * t96) * t85 + (t31 * t95 - t68 * t96) * t84 - t7) * m(6) + (-t83 / 0.2e1 - t152 / 0.2e1) * t108 + (t159 + (t58 + t73) * t101 + t3) * t84 / 0.2e1 + (t158 + (t60 - t74) * t101 + t4 + t2) * t85 / 0.2e1) * qJD(5); 0, 0, t16 * qJD(5) + 0.4e1 * (-t143 / 0.4e1 - t142 / 0.4e1) * qJD(3) + ((t47 * t85 + t17 + t37) * t150 + (t30 * t85 + t10 + t26) * t149) * t151, 0, t16 * qJD(3) - t116 * t113; t20, t18 * qJD(3), qJD(1) * t112 + t1 * qJD(5) - t107 + (-t148 - t132 * t141 * t85 - (t154 + t94) * t101 / 0.2e1 + ((-t58 / 0.2e1 - t59 / 0.2e1) * t85 - t153) * t102) * qJD(3), t15 * qJD(3), t1 * qJD(3) + (m(6) * (t116 * t96 * t95 + t22 * t28) + (-t63 * t133 + t83 * t62) * t146 - t85 * (t62 * t133 - t152 * t63) / 0.2e1) * qJD(5);];
Cq = t9;
