% Calculate time derivative of joint inertia matrix for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:35
% DurationCPUTime: 1.48s
% Computational Cost: add. (1986->184), mult. (4626->278), div. (0->0), fcn. (4530->6), ass. (0->94)
t84 = sin(qJ(4));
t127 = qJD(4) * t84;
t140 = sin(qJ(1));
t142 = cos(qJ(1));
t130 = t142 * pkin(1) + t140 * qJ(2);
t125 = t142 * pkin(2) + t130;
t52 = t142 * rSges(3,1) + t140 * rSges(3,3) + t130;
t155 = qJD(1) - qJD(3);
t139 = sin(qJ(3));
t141 = cos(qJ(3));
t61 = -t140 * t139 - t142 * t141;
t45 = t155 * t61;
t85 = cos(qJ(4));
t136 = t45 * t85;
t62 = t142 * t139 - t140 * t141;
t101 = t62 * t127 - t136;
t46 = t155 * t62;
t135 = t46 * t85;
t99 = t61 * t127 + t135;
t154 = -t45 * t62 + t46 * t61;
t153 = t45 * t84;
t152 = t46 * t84;
t128 = Icges(5,4) * t85;
t104 = Icges(5,2) * t84 - t128;
t34 = Icges(5,6) * t62 + t104 * t61;
t129 = Icges(5,4) * t84;
t105 = -Icges(5,1) * t85 + t129;
t36 = Icges(5,5) * t62 + t105 * t61;
t109 = t34 * t84 - t36 * t85;
t151 = t109 * t62;
t33 = -Icges(5,6) * t61 + t104 * t62;
t35 = -Icges(5,5) * t61 + t105 * t62;
t110 = t33 * t84 - t35 * t85;
t150 = t110 * t61;
t126 = qJD(4) * t85;
t100 = t61 * t126 - t152;
t97 = rSges(5,1) * t99 + t45 * rSges(5,3);
t8 = -t100 * rSges(5,2) - t46 * pkin(3) - t45 * pkin(6) - t97;
t102 = t62 * t126 + t153;
t98 = rSges(5,1) * t101 + t46 * rSges(5,3);
t7 = t102 * rSges(5,2) - t45 * pkin(3) + t46 * pkin(6) + t98;
t137 = rSges(5,2) * t84;
t121 = -pkin(3) + t137;
t138 = rSges(5,1) * t85;
t132 = t62 * rSges(5,3) - t61 * t138;
t26 = -pkin(6) * t62 - t121 * t61 - t132;
t133 = t61 * rSges(5,3) + t62 * t138;
t25 = -t61 * pkin(6) + t121 * t62 - t133;
t120 = (Icges(5,2) * t85 + t105 + t129) * t127;
t59 = t104 * qJD(4);
t68 = -Icges(5,1) * t84 - t128;
t149 = -(qJD(4) * t68 + t59) * t85 - t120;
t148 = 2 * m(4);
t147 = 2 * m(5);
t81 = t142 * qJ(2);
t131 = qJD(1) * t81 + qJD(2) * t140;
t124 = t140 * pkin(1);
t30 = -t46 * rSges(4,1) + t45 * rSges(4,2);
t29 = -t45 * rSges(4,1) - t46 * rSges(4,2);
t47 = -t62 * rSges(4,1) + t61 * rSges(4,2);
t48 = rSges(4,1) * t61 + rSges(4,2) * t62;
t103 = -Icges(5,5) * t85 + Icges(5,6) * t84;
t96 = -t140 * pkin(2) - t124;
t93 = t81 + t96;
t90 = -t140 * rSges(3,1) + t142 * rSges(3,3) - t124;
t89 = t96 * qJD(1) + t131;
t79 = qJD(2) * t142;
t87 = -qJD(1) * t125 + t79;
t58 = t103 * qJD(4);
t86 = -(-t34 * t85 - t84 * t36) * t45 / 0.2e1 - (-t33 * t85 - t84 * t35) * t46 / 0.2e1 + (t110 * qJD(4) - (t101 * Icges(5,4) + t102 * Icges(5,2) + Icges(5,6) * t46) * t85 - t84 * (t101 * Icges(5,1) + t102 * Icges(5,4) + Icges(5,5) * t46) - t61 * t58) * t61 / 0.2e1 - (t109 * qJD(4) - (t99 * Icges(5,4) + t100 * Icges(5,2) + Icges(5,6) * t45) * t85 - t84 * (t99 * Icges(5,1) + t100 * Icges(5,4) + Icges(5,5) * t45) + t62 * t58) * t62 / 0.2e1 + t154 * (-Icges(5,5) * t84 - Icges(5,6) * t85);
t69 = -t84 * rSges(5,1) - rSges(5,2) * t85;
t64 = (t137 - t138) * qJD(4);
t51 = t81 + t90;
t50 = -qJD(1) * t52 + t79;
t49 = t90 * qJD(1) + t131;
t40 = -t48 + t125;
t39 = -t47 + t93;
t38 = t61 * t137 + t132;
t37 = t62 * t137 - t133;
t32 = Icges(5,3) * t62 + t103 * t61;
t31 = -Icges(5,3) * t61 + t103 * t62;
t24 = -t29 + t87;
t23 = -t30 + t89;
t22 = t125 - t26;
t21 = -t25 + t93;
t14 = t99 * Icges(5,5) + t100 * Icges(5,6) + Icges(5,3) * t45;
t13 = t101 * Icges(5,5) + t102 * Icges(5,6) + Icges(5,3) * t46;
t12 = t109 * t61 + t62 * t32;
t11 = t62 * t31 + t150;
t10 = -t61 * t32 + t151;
t9 = t110 * t62 - t61 * t31;
t6 = -t7 + t87;
t5 = -t8 + t89;
t1 = [0.2e1 * m(3) * (t49 * t52 + t50 * t51) + (t23 * t40 + t24 * t39) * t148 - t68 * t126 - t85 * t59 + (t21 * t6 + t22 * t5) * t147 - t120; m(3) * (t140 * t50 - t142 * t49 + (t140 * t52 + t142 * t51) * qJD(1)) + m(4) * (t140 * t24 - t142 * t23 + (t140 * t40 + t142 * t39) * qJD(1)) + m(5) * (t140 * t6 - t142 * t5 + (t140 * t22 + t142 * t21) * qJD(1)); 0; m(4) * (t23 * t48 + t24 * t47 + t29 * t39 + t30 * t40) + m(5) * (t21 * t7 + t22 * t8 + t25 * t6 + t26 * t5) - t149; m(4) * (t29 * t140 - t30 * t142 + (t140 * t48 + t142 * t47) * qJD(1)) + m(5) * (t7 * t140 - t8 * t142 + (t140 * t26 + t142 * t25) * qJD(1)); (t29 * t47 + t30 * t48) * t148 + (t25 * t7 + t26 * t8) * t147 + t149; m(5) * ((-t21 * t61 - t22 * t62) * t64 + (t21 * t46 - t22 * t45 - t5 * t62 - t6 * t61) * t69) - t86; m(5) * ((-t140 * t61 + t142 * t62) * t64 + (t140 * t46 + t142 * t45 + (-t140 * t62 - t142 * t61) * qJD(1)) * t69); m(5) * ((-t25 * t61 - t26 * t62) * t64 + (t25 * t46 - t26 * t45 - t61 * t7 - t62 * t8) * t69) + t86; ((t37 * t62 + t38 * t61) * (t45 * t37 + t62 * t98 - t46 * t38 + t61 * t97 + (t61 * t100 + t62 * t102) * rSges(5,2)) + ((t61 ^ 2 + t62 ^ 2) * t64 - t154 * t69) * t69) * t147 + t45 * (-t11 * t61 + t12 * t62) + t62 * ((t62 * t14 + t45 * t32) * t62 + t12 * t45 + (t11 - t151) * t46 + (-t62 * t13 - t35 * t135 + t33 * t152 - t45 * t31) * t61) + t46 * (t10 * t62 - t9 * t61) - t61 * (-(-t61 * t13 + t46 * t31) * t61 + t9 * t46 - (-t10 + t150) * t45 + (-t36 * t136 - t61 * t14 + t34 * t153 + t46 * t32) * t62);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
