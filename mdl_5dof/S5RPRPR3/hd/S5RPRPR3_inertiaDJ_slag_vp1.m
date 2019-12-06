% Calculate time derivative of joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:21
% DurationCPUTime: 2.14s
% Computational Cost: add. (5523->249), mult. (4724->369), div. (0->0), fcn. (4266->10), ass. (0->121)
t111 = qJ(1) + pkin(8);
t166 = pkin(2) * sin(t111) + sin(qJ(1)) * pkin(1);
t114 = sin(qJ(5));
t116 = cos(qJ(5));
t112 = sin(pkin(9));
t113 = cos(pkin(9));
t150 = Icges(6,4) * t116;
t79 = -Icges(6,6) * t113 + (-Icges(6,2) * t114 + t150) * t112;
t151 = Icges(6,4) * t114;
t80 = -Icges(6,5) * t113 + (Icges(6,1) * t116 - t151) * t112;
t165 = -t114 * t80 - t116 * t79;
t164 = 2 * m(4);
t163 = 2 * m(5);
t162 = 2 * m(6);
t110 = qJD(1) + qJD(3);
t109 = qJ(3) + t111;
t104 = sin(t109);
t105 = cos(t109);
t143 = t113 * t116;
t124 = t104 * t143 - t105 * t114;
t144 = t113 * t114;
t125 = -t104 * t116 + t105 * t144;
t61 = t124 * qJD(5) + t125 * t110;
t74 = t104 * t144 + t105 * t116;
t77 = t104 * t114 + t105 * t143;
t62 = t74 * qJD(5) - t77 * t110;
t158 = t62 * rSges(6,1) + t61 * rSges(6,2);
t157 = -rSges(6,1) * t124 + t74 * rSges(6,2);
t146 = t110 * t105;
t147 = t110 * t104;
t72 = rSges(4,1) * t147 + rSges(4,2) * t146;
t141 = qJD(5) * t112;
t83 = (-Icges(6,2) * t116 - t151) * t141;
t155 = t114 * t83;
t82 = (-Icges(6,5) * t114 - Icges(6,6) * t116) * t141;
t84 = (-Icges(6,1) * t114 - t150) * t141;
t133 = t112 * t116 * t84 - t113 * t82;
t26 = (t165 * qJD(5) - t155) * t112 + t133;
t153 = t26 * t113;
t152 = -rSges(5,3) - qJ(4);
t149 = t104 * t112;
t148 = t105 * t112;
t145 = t110 * t112;
t142 = t166 * qJD(1);
t97 = rSges(5,2) * t149;
t140 = t104 * t145;
t139 = t113 * t147;
t138 = t105 * t145;
t134 = -rSges(5,1) * t113 - pkin(3);
t87 = -t105 * rSges(4,1) + t104 * rSges(4,2);
t132 = pkin(3) * t147 - qJD(4) * t104;
t73 = -rSges(4,1) * t146 + rSges(4,2) * t147;
t129 = -t77 * rSges(6,1) + rSges(6,2) * t125;
t127 = -pkin(2) * cos(t111) - cos(qJ(1)) * pkin(1);
t86 = -rSges(4,1) * t104 - rSges(4,2) * t105;
t59 = -t77 * qJD(5) + t74 * t110;
t60 = -t125 * qJD(5) - t124 * t110;
t33 = rSges(6,1) * t60 + t59 * rSges(6,2) - rSges(6,3) * t140;
t123 = t127 * qJD(1);
t122 = -pkin(4) * t113 - pkin(3) + (-rSges(6,3) - pkin(7)) * t112;
t100 = t105 * qJ(4);
t63 = t105 * rSges(5,3) + t134 * t104 + t100 + t97;
t47 = -Icges(6,5) * t124 + Icges(6,6) * t74 - Icges(6,3) * t149;
t49 = -Icges(6,4) * t124 + Icges(6,2) * t74 - Icges(6,6) * t149;
t51 = -Icges(6,1) * t124 + Icges(6,4) * t74 - Icges(6,5) * t149;
t18 = -t113 * t47 + (-t114 * t49 + t116 * t51) * t112;
t48 = Icges(6,5) * t77 - Icges(6,6) * t125 + Icges(6,3) * t148;
t50 = Icges(6,4) * t77 - Icges(6,2) * t125 + Icges(6,6) * t148;
t52 = Icges(6,1) * t77 - Icges(6,4) * t125 + Icges(6,5) * t148;
t19 = -t113 * t48 + (-t114 * t50 + t116 * t52) * t112;
t28 = Icges(6,5) * t62 + Icges(6,6) * t61 - Icges(6,3) * t138;
t30 = Icges(6,4) * t62 + Icges(6,2) * t61 - Icges(6,6) * t138;
t32 = Icges(6,1) * t62 + Icges(6,4) * t61 - Icges(6,5) * t138;
t3 = -t113 * t28 + (-t114 * t30 + t116 * t32 + (-t114 * t51 - t116 * t49) * qJD(5)) * t112;
t78 = -Icges(6,3) * t113 + (Icges(6,5) * t116 - Icges(6,6) * t114) * t112;
t35 = -t124 * t80 - t78 * t149 + t74 * t79;
t36 = -t125 * t79 + t78 * t148 + t77 * t80;
t27 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t140;
t29 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t140;
t31 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t140;
t4 = -t113 * t27 + (-t114 * t29 + t116 * t31 + (-t114 * t52 - t116 * t50) * qJD(5)) * t112;
t8 = t59 * t79 + t60 * t80 - t125 * t83 + t77 * t84 + (t105 * t82 - t78 * t147) * t112;
t9 = t61 * t79 + t62 * t80 + t74 * t83 - t124 * t84 + (-t104 * t82 - t78 * t146) * t112;
t120 = -t153 - (t3 + t9) * t149 / 0.2e1 + (t4 + t8) * t148 / 0.2e1 - ((t19 + t36) * t104 + (t18 + t35) * t105) * t145 / 0.2e1;
t119 = t152 * t104 + t134 * t105;
t64 = rSges(5,2) * t148 + t119;
t118 = -t104 * qJ(4) + t122 * t105;
t39 = t122 * t104 + t100 + t157;
t99 = qJD(4) * t105;
t46 = rSges(5,2) * t138 + t119 * t110 + t99;
t45 = rSges(5,1) * t139 + (t152 * t105 - t97) * t110 + t132;
t22 = pkin(4) * t139 + pkin(7) * t140 - qJ(4) * t146 + t132 - t33;
t40 = t118 + t129;
t23 = t118 * t110 + t158 + t99;
t85 = (-rSges(6,1) * t114 - rSges(6,2) * t116) * t141;
t81 = -rSges(6,3) * t113 + (rSges(6,1) * t116 - rSges(6,2) * t114) * t112;
t68 = t127 + t87;
t67 = -t166 + t86;
t66 = t123 + t73;
t65 = t142 + t72;
t56 = t127 + t64;
t55 = -t166 + t63;
t54 = rSges(6,3) * t148 - t129;
t53 = -rSges(6,3) * t149 + t157;
t44 = t123 + t46;
t43 = t45 + t142;
t42 = t113 * t54 + t81 * t148;
t41 = -t113 * t53 + t81 * t149;
t38 = t127 + t40;
t37 = -t166 + t39;
t34 = -rSges(6,3) * t138 + t158;
t21 = -t113 * t34 + (t104 * t85 + t81 * t146) * t112;
t20 = t113 * t33 + (t105 * t85 - t81 * t147) * t112;
t17 = t123 + t23;
t16 = t22 + t142;
t13 = -t125 * t50 + t48 * t148 + t52 * t77;
t12 = -t125 * t49 + t47 * t148 + t51 * t77;
t11 = -t124 * t52 - t48 * t149 + t50 * t74;
t10 = -t124 * t51 - t47 * t149 + t49 * t74;
t5 = ((-t110 * t54 - t34) * t105 + (t110 * t53 - t33) * t104) * t112;
t1 = [(t65 * t68 + t66 * t67) * t164 + (t43 * t56 + t44 * t55) * t163 - t112 * t155 + (t16 * t38 + t17 * t37) * t162 + t133 + t165 * t141; 0; 0; m(4) * (t65 * t87 + t66 * t86 + t67 * t73 + t68 * t72) + m(5) * (t43 * t64 + t44 * t63 + t45 * t56 + t46 * t55) + m(6) * (t16 * t40 + t17 * t39 + t22 * t38 + t23 * t37) + t26; 0; (t22 * t40 + t23 * t39) * t162 + (t45 * t64 + t46 * t63) * t163 + (t72 * t87 + t73 * t86) * t164 + t26; m(5) * ((t110 * t55 + t43) * t105 + (-t110 * t56 + t44) * t104) + m(6) * ((t110 * t37 + t16) * t105 + (-t110 * t38 + t17) * t104); 0; m(6) * ((t110 * t39 + t22) * t105 + (-t110 * t40 + t23) * t104) + m(5) * ((t110 * t63 + t45) * t105 + (-t110 * t64 + t46) * t104); 0; m(6) * (t16 * t42 + t17 * t41 + t20 * t38 + t21 * t37) + t120; m(6) * t5; m(6) * (t20 * t40 + t21 * t39 + t22 * t42 + t23 * t41) + t120; m(6) * ((t110 * t41 + t20) * t105 + (-t110 * t42 + t21) * t104); (t42 * t20 + t41 * t21 + (-t104 * t54 - t105 * t53) * t5 * t112) * t162 - t113 * (-t153 + ((-t110 * t18 + t4) * t105 + (-t110 * t19 - t3) * t104) * t112) - (-t113 * t35 + (-t10 * t104 + t105 * t11) * t112) * t138 - (-t9 * t113 + (-(-t124 * t32 + t74 * t30 + t61 * t49 + t62 * t51) * t104 - t10 * t146 + (-t124 * t31 + t74 * t29 + t61 * t50 + t62 * t52) * t105 - t11 * t147 + (-(-t104 * t28 - t47 * t146) * t104 + (-t104 * t27 - t48 * t146) * t105) * t112) * t112) * t149 - (-t113 * t36 + (-t104 * t12 + t105 * t13) * t112) * t140 + (-t8 * t113 + (-(-t125 * t30 + t77 * t32 + t59 * t49 + t60 * t51) * t104 - t12 * t146 + (-t125 * t29 + t77 * t31 + t59 * t50 + t60 * t52) * t105 - t13 * t147 + (-(t105 * t28 - t47 * t147) * t104 + (t105 * t27 - t48 * t147) * t105) * t112) * t112) * t148;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
