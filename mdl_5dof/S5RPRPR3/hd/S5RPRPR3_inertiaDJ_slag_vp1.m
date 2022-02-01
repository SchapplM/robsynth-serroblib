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
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:40
% DurationCPUTime: 1.52s
% Computational Cost: add. (5523->250), mult. (4724->370), div. (0->0), fcn. (4266->10), ass. (0->121)
t117 = sin(qJ(5));
t119 = cos(qJ(5));
t115 = sin(pkin(9));
t116 = cos(pkin(9));
t152 = Icges(6,4) * t119;
t79 = -Icges(6,6) * t116 + (-Icges(6,2) * t117 + t152) * t115;
t153 = Icges(6,4) * t117;
t80 = -Icges(6,5) * t116 + (Icges(6,1) * t119 - t153) * t115;
t164 = -t117 * t80 - t119 * t79;
t114 = qJ(1) + pkin(8);
t142 = pkin(2) * cos(t114) + cos(qJ(1)) * pkin(1);
t163 = 2 * m(4);
t162 = 2 * m(5);
t161 = 2 * m(6);
t113 = qJD(1) + qJD(3);
t111 = qJ(3) + t114;
t107 = sin(t111);
t108 = cos(t111);
t145 = t116 * t117;
t127 = t107 * t145 + t108 * t119;
t144 = t116 * t119;
t77 = t107 * t117 + t108 * t144;
t59 = -t77 * qJD(5) + t127 * t113;
t75 = t107 * t144 - t108 * t117;
t76 = t107 * t119 - t108 * t145;
t60 = t76 * qJD(5) - t75 * t113;
t159 = t60 * rSges(6,1) + t59 * rSges(6,2);
t147 = t113 * t108;
t158 = qJ(4) * t147 + qJD(4) * t107;
t141 = qJD(5) * t115;
t83 = (-Icges(6,2) * t119 - t153) * t141;
t156 = t117 * t83;
t82 = (-Icges(6,5) * t117 - Icges(6,6) * t119) * t141;
t84 = (-Icges(6,1) * t117 - t152) * t141;
t135 = t115 * t119 * t84 - t116 * t82;
t26 = (qJD(5) * t164 - t156) * t115 + t135;
t154 = t26 * t116;
t100 = t107 * qJ(4);
t151 = t107 * t115;
t150 = t108 * t115;
t149 = t108 * t116;
t148 = t113 * t107;
t146 = t113 * t115;
t143 = t108 * pkin(3) + t100;
t54 = t77 * rSges(6,1) + t76 * rSges(6,2) + rSges(6,3) * t150;
t140 = t107 * t146;
t139 = t108 * t146;
t136 = -rSges(5,1) * t116 - pkin(3);
t87 = t108 * rSges(4,1) - rSges(4,2) * t107;
t73 = -rSges(4,1) * t147 + rSges(4,2) * t148;
t61 = -t75 * qJD(5) + t76 * t113;
t62 = -t127 * qJD(5) + t77 * t113;
t132 = -t62 * rSges(6,1) - t61 * rSges(6,2);
t131 = -rSges(6,1) * t75 + rSges(6,2) * t127;
t130 = -pkin(2) * sin(t114) - sin(qJ(1)) * pkin(1);
t129 = t136 * t107;
t86 = -rSges(4,1) * t107 - rSges(4,2) * t108;
t126 = t130 * qJD(1);
t125 = t142 * qJD(1);
t40 = pkin(4) * t149 + pkin(7) * t150 + t143 + t54;
t72 = t86 * t113;
t124 = -pkin(4) * t116 - pkin(3) + (-rSges(6,3) - pkin(7)) * t115;
t64 = rSges(5,1) * t149 - rSges(5,2) * t150 + t107 * rSges(5,3) + t143;
t101 = t108 * qJ(4);
t63 = rSges(5,2) * t151 + t108 * rSges(5,3) + t101 + t129;
t122 = t124 * t107;
t47 = Icges(6,5) * t75 - Icges(6,6) * t127 + Icges(6,3) * t151;
t49 = Icges(6,4) * t75 - Icges(6,2) * t127 + Icges(6,6) * t151;
t51 = Icges(6,1) * t75 - Icges(6,4) * t127 + Icges(6,5) * t151;
t18 = -t116 * t47 + (-t117 * t49 + t119 * t51) * t115;
t48 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t150;
t50 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t150;
t52 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t150;
t19 = -t116 * t48 + (-t117 * t50 + t119 * t52) * t115;
t28 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t139;
t30 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t139;
t32 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t139;
t3 = -t116 * t28 + (-t117 * t30 + t119 * t32 + (-t117 * t51 - t119 * t49) * qJD(5)) * t115;
t78 = -Icges(6,3) * t116 + (Icges(6,5) * t119 - Icges(6,6) * t117) * t115;
t35 = -t127 * t79 + t78 * t151 + t75 * t80;
t36 = t78 * t150 + t76 * t79 + t77 * t80;
t27 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t140;
t29 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t140;
t31 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t140;
t4 = -t116 * t27 + (-t117 * t29 + t119 * t31 + (-t117 * t52 - t119 * t50) * qJD(5)) * t115;
t8 = t59 * t79 + t60 * t80 + t76 * t83 + t77 * t84 + (t108 * t82 - t78 * t148) * t115;
t9 = t61 * t79 + t62 * t80 - t127 * t83 + t75 * t84 + (t107 * t82 + t78 * t147) * t115;
t121 = -t154 + (t3 + t9) * t151 / 0.2e1 - (t36 + t19) * t140 / 0.2e1 + (t4 + t8 + (t18 + t35) * t113) * t150 / 0.2e1;
t45 = rSges(5,2) * t140 + rSges(5,3) * t147 + t113 * t129 + t158;
t99 = qJD(4) * t108;
t46 = rSges(5,2) * t139 + t99 + (t136 * t108 + (-rSges(5,3) - qJ(4)) * t107) * t113;
t22 = t113 * t122 + t158 + t159;
t39 = t101 + t122 + t131;
t23 = t99 + (t124 * t108 - t100) * t113 + t132;
t85 = (-rSges(6,1) * t117 - rSges(6,2) * t119) * t141;
t81 = -rSges(6,3) * t116 + (rSges(6,1) * t119 - rSges(6,2) * t117) * t115;
t68 = t87 + t142;
t67 = t130 + t86;
t66 = -t125 + t73;
t65 = t72 + t126;
t56 = t64 + t142;
t55 = t130 + t63;
t53 = rSges(6,3) * t151 - t131;
t44 = -t125 + t46;
t43 = t126 + t45;
t42 = -t116 * t54 - t81 * t150;
t41 = t116 * t53 + t81 * t151;
t38 = t40 + t142;
t37 = t130 + t39;
t34 = rSges(6,3) * t139 - t132;
t33 = -rSges(6,3) * t140 + t159;
t21 = t116 * t34 + (t107 * t85 + t81 * t147) * t115;
t20 = -t116 * t33 + (-t108 * t85 + t81 * t148) * t115;
t17 = -t125 + t23;
t16 = t126 + t22;
t13 = t48 * t150 + t50 * t76 + t52 * t77;
t12 = t47 * t150 + t49 * t76 + t51 * t77;
t11 = -t127 * t50 + t48 * t151 + t52 * t75;
t10 = -t127 * t49 + t47 * t151 + t51 * t75;
t5 = ((-t113 * t54 + t34) * t108 + (-t113 * t53 - t33) * t107) * t115;
t1 = [(t16 * t38 + t17 * t37) * t161 + (t43 * t56 + t44 * t55) * t162 - t115 * t156 + (t65 * t68 + t66 * t67) * t163 + t135 + t164 * t141; 0; 0; m(6) * (t16 * t40 + t17 * t39 + t22 * t38 + t23 * t37) + m(5) * (t43 * t64 + t44 * t63 + t45 * t56 + t46 * t55) + m(4) * (t65 * t87 + t66 * t86 + t67 * t73 + t68 * t72) + t26; 0; (t22 * t40 + t23 * t39) * t161 + (t45 * t64 + t46 * t63) * t162 + (t72 * t87 + t73 * t86) * t163 + t26; m(6) * ((t113 * t37 - t16) * t108 + (t113 * t38 + t17) * t107) + m(5) * ((t113 * t55 - t43) * t108 + (t113 * t56 + t44) * t107); 0; m(6) * ((t113 * t39 - t22) * t108 + (t113 * t40 + t23) * t107) + m(5) * ((t113 * t63 - t45) * t108 + (t113 * t64 + t46) * t107); 0; m(6) * (t16 * t42 + t17 * t41 + t20 * t38 + t21 * t37) + t121; m(6) * t5; m(6) * (t20 * t40 + t21 * t39 + t22 * t42 + t23 * t41) + t121; m(6) * ((t113 * t41 - t20) * t108 + (t113 * t42 + t21) * t107); (t42 * t20 + t41 * t21 + (-t107 * t54 + t108 * t53) * t5 * t115) * t161 - (-t36 * t116 + (t107 * t12 + t108 * t13) * t115) * t140 + (-t8 * t116 + ((t76 * t29 + t77 * t31 + t59 * t50 + t60 * t52) * t108 - t13 * t148 + (t76 * t30 + t77 * t32 + t59 * t49 + t60 * t51) * t107 + t12 * t147 + ((t108 * t27 - t48 * t148) * t108 + (t108 * t28 - t47 * t148) * t107) * t115) * t115) * t150 + (-t116 * t35 + (t10 * t107 + t108 * t11) * t115) * t139 + (-t9 * t116 + ((-t127 * t29 + t75 * t31 + t61 * t50 + t62 * t52) * t108 - t11 * t148 + (-t127 * t30 + t75 * t32 + t61 * t49 + t62 * t51) * t107 + t10 * t147 + ((t107 * t27 + t48 * t147) * t108 + (t107 * t28 + t47 * t147) * t107) * t115) * t115) * t151 - t116 * (-t154 + ((t113 * t18 + t4) * t108 + (-t113 * t19 + t3) * t107) * t115);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
