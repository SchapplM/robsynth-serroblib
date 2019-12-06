% Calculate time derivative of joint inertia matrix for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:38
% DurationCPUTime: 1.28s
% Computational Cost: add. (3598->173), mult. (2214->250), div. (0->0), fcn. (1596->8), ass. (0->109)
t98 = pkin(9) + qJ(5);
t93 = sin(t98);
t135 = qJD(5) * t93;
t95 = cos(t98);
t134 = qJD(5) * t95;
t100 = qJD(2) + qJD(3);
t165 = t100 * t93;
t164 = t100 * t95;
t99 = pkin(8) + qJ(2);
t97 = qJ(3) + t99;
t90 = sin(t97);
t163 = (rSges(5,3) + qJ(4)) * t90;
t102 = cos(pkin(9));
t126 = -rSges(5,1) * t102 - pkin(3);
t91 = cos(t97);
t161 = t91 * qJ(4) + t126 * t90;
t149 = rSges(6,2) * t93;
t150 = rSges(6,1) * t95;
t160 = -t149 + t150;
t144 = Icges(6,4) * t95;
t113 = -Icges(6,2) * t93 + t144;
t39 = Icges(6,6) * t90 + t113 * t91;
t145 = Icges(6,4) * t93;
t114 = Icges(6,1) * t95 - t145;
t41 = Icges(6,5) * t90 + t114 * t91;
t116 = t39 * t93 - t41 * t95;
t159 = t116 * t90;
t38 = -Icges(6,6) * t91 + t113 * t90;
t40 = -Icges(6,5) * t91 + t114 * t90;
t118 = t38 * t93 - t40 * t95;
t158 = t118 * t91;
t112 = Icges(6,5) * t95 - Icges(6,6) * t93;
t62 = Icges(6,2) * t95 + t145;
t63 = Icges(6,1) * t93 + t144;
t115 = t62 * t93 - t63 * t95;
t157 = t112 * qJD(5) + t100 * t115;
t36 = -Icges(6,3) * t91 + t112 * t90;
t156 = 2 * m(4);
t155 = 2 * m(5);
t154 = 2 * m(6);
t94 = sin(t99);
t151 = pkin(2) * t94;
t83 = t90 * rSges(6,3);
t148 = t91 * rSges(6,3) + t90 * t149;
t146 = rSges(5,2) * sin(pkin(9));
t140 = pkin(2) * qJD(2);
t37 = Icges(6,3) * t90 + t112 * t91;
t139 = t100 * t37;
t138 = t100 * t90;
t137 = t100 * t91;
t103 = -pkin(7) - qJ(4);
t136 = t103 * t90;
t132 = t100 * t149;
t133 = -t91 * t132 + (-rSges(6,1) * t135 - rSges(6,2) * t134) * t90;
t131 = t94 * t140;
t96 = cos(t99);
t130 = t96 * t140;
t129 = t100 * t146;
t92 = pkin(4) * t102 + pkin(3);
t125 = -t92 - t150;
t58 = t91 * rSges(4,1) - rSges(4,2) * t90;
t47 = -rSges(4,1) * t137 + rSges(4,2) * t138;
t57 = -rSges(4,1) * t90 - rSges(4,2) * t91;
t68 = rSges(6,1) * t93 + rSges(6,2) * t95;
t106 = -t103 * t91 + t125 * t90;
t30 = t106 + t148;
t28 = t30 - t151;
t43 = t160 * t91 + t83;
t31 = t91 * t92 - t136 + t43;
t89 = pkin(2) * t96;
t29 = t31 + t89;
t121 = t28 * t91 + t29 * t90;
t120 = t30 * t91 + t31 * t90;
t119 = t38 * t95 + t40 * t93;
t117 = t39 * t95 + t41 * t93;
t61 = Icges(6,5) * t93 + Icges(6,6) * t95;
t111 = (t114 - t62) * t135 + (t113 + t63) * t134;
t46 = t57 * t100;
t110 = (-qJD(5) * t116 + t157 * t90 - t38 * t164 - t40 * t165) * t90 / 0.2e1 - (-qJD(5) * t118 - t157 * t91 + t39 * t164 + t41 * t165) * t91 / 0.2e1 + (-t115 * t90 - t61 * t91 + t119) * t138 / 0.2e1 + (-t115 * t91 + t61 * t90 + t117) * t137 / 0.2e1;
t107 = qJD(5) * t61;
t35 = t163 + (-t146 - t126) * t91;
t34 = t91 * rSges(5,3) + t90 * t146 + t161;
t105 = -qJD(5) * t68 * t91 + rSges(6,3) * t137 + t90 * t132;
t79 = qJD(4) * t90;
t26 = rSges(5,3) * t137 + t100 * t161 + t90 * t129 + t79;
t80 = qJD(4) * t91;
t27 = t91 * t129 + t80 + (t126 * t91 - t163) * t100;
t13 = t80 - t133 + (t125 * t91 + t136 - t83) * t100;
t12 = t100 * t106 + t105 + t79;
t56 = t160 * qJD(5);
t49 = t58 + t89;
t48 = t57 - t151;
t45 = t47 - t130;
t44 = t46 - t131;
t42 = t150 * t90 - t148;
t33 = t35 + t89;
t32 = t34 - t151;
t21 = -t107 * t90 + t139;
t20 = -t100 * t36 - t91 * t107;
t19 = t27 - t130;
t18 = t26 - t131;
t11 = t13 - t130;
t10 = t12 - t131;
t9 = -t116 * t91 + t37 * t90;
t8 = t36 * t90 - t158;
t7 = -t37 * t91 - t159;
t6 = -t118 * t90 - t36 * t91;
t3 = ((-t43 + t83) * t100 + t133) * t90 + (t100 * t42 + t105) * t91;
t1 = [0; 0; (t18 * t33 + t19 * t32) * t155 + (t10 * t29 + t11 * t28) * t154 + (t44 * t49 + t45 * t48) * t156 + t111; 0; m(5) * (t18 * t35 + t19 * t34 + t26 * t33 + t27 * t32) + m(6) * (t10 * t31 + t11 * t30 + t12 * t29 + t13 * t28) + m(4) * (t44 * t58 + t45 * t57 + t46 * t49 + t47 * t48) + t111; (t12 * t31 + t13 * t30) * t154 + (t26 * t35 + t27 * t34) * t155 + (t46 * t58 + t47 * t57) * t156 + t111; 0; m(5) * (-t18 * t91 + t19 * t90 + (t32 * t91 + t33 * t90) * t100) + m(6) * (-t10 * t91 + t100 * t121 + t11 * t90); m(6) * (t100 * t120 - t12 * t91 + t13 * t90) + m(5) * (-t26 * t91 + t27 * t90 + (t34 * t91 + t35 * t90) * t100); 0; m(6) * t3; m(6) * (-t121 * t56 + (-t10 * t90 - t11 * t91 + (t28 * t90 - t29 * t91) * t100) * t68) + t110; m(6) * (-t120 * t56 + (-t12 * t90 - t13 * t91 + (t30 * t90 - t31 * t91) * t100) * t68) + t110; 0; ((t42 * t90 + t43 * t91) * t3 + (t90 ^ 2 + t91 ^ 2) * t68 * t56) * t154 + (-t8 * t91 + t9 * t90) * t137 + t90 * ((t90 * t20 + (t8 + t159) * t100) * t90 + (t9 * t100 + (t134 * t38 + t135 * t40) * t91 + (-t21 - t117 * qJD(5) + (-t118 + t37) * t100) * t90) * t91) + (-t6 * t91 + t7 * t90) * t138 - t91 * ((t91 * t21 + (t7 + t158) * t100) * t91 + (t6 * t100 + (-t134 * t39 - t41 * t135 + t139) * t90 + (t119 * qJD(5) - t100 * t116 - t20) * t91) * t90);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
