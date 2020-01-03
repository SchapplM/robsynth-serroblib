% Calculate time derivative of joint inertia matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:30
% DurationCPUTime: 1.18s
% Computational Cost: add. (2004->180), mult. (2112->264), div. (0->0), fcn. (1498->6), ass. (0->107)
t91 = sin(qJ(4));
t93 = cos(qJ(4));
t116 = rSges(5,1) * t91 + rSges(5,2) * t93;
t90 = qJ(1) + qJ(2);
t87 = cos(t90);
t135 = t116 * t87;
t86 = sin(t90);
t154 = t86 / 0.2e1;
t153 = rSges(4,2) - pkin(2);
t133 = Icges(5,4) * t91;
t106 = Icges(5,2) * t93 + t133;
t39 = Icges(5,6) * t87 + t106 * t86;
t132 = Icges(5,4) * t93;
t107 = Icges(5,1) * t91 + t132;
t41 = Icges(5,5) * t87 + t107 * t86;
t113 = t39 * t93 + t41 * t91;
t152 = t113 * t87;
t53 = t106 * qJD(4);
t54 = t107 * qJD(4);
t61 = Icges(5,5) * t93 - Icges(5,6) * t91;
t62 = -Icges(5,2) * t91 + t132;
t63 = Icges(5,1) * t93 - t133;
t89 = qJD(1) + qJD(2);
t151 = (t62 * t91 - t63 * t93) * qJD(4) + t53 * t93 + t54 * t91 + t61 * t89;
t150 = -Icges(5,3) * t89 + qJD(4) * t61;
t149 = -Icges(5,6) * t89 + qJD(4) * t62;
t148 = -Icges(5,5) * t89 + qJD(4) * t63;
t147 = 2 * m(3);
t146 = 2 * m(4);
t145 = 2 * m(5);
t92 = sin(qJ(1));
t142 = pkin(1) * t92;
t83 = t87 * pkin(2);
t139 = t86 * rSges(5,3);
t138 = t86 * t89;
t137 = t87 * t89;
t76 = t87 * qJ(3);
t136 = qJD(3) * t86 + t89 * t76;
t134 = t86 * qJ(3) + t83;
t128 = pkin(1) * qJD(1);
t127 = qJD(4) * t91;
t126 = qJD(4) * t93;
t125 = -rSges(5,3) - pkin(2) - pkin(6);
t121 = rSges(5,1) * t126;
t124 = -t86 * t121 - t135 * t89;
t43 = t87 * rSges(5,3) + t116 * t86;
t123 = t92 * t128;
t94 = cos(qJ(1));
t122 = t94 * t128;
t120 = rSges(5,2) * t127;
t57 = t116 * qJD(4);
t118 = (t86 ^ 2 + t87 ^ 2) * t57;
t51 = t87 * rSges(3,1) - rSges(3,2) * t86;
t117 = t91 * t53 - t93 * t54;
t46 = -rSges(3,1) * t137 + rSges(3,2) * t138;
t50 = -rSges(3,1) * t86 - rSges(3,2) * t87;
t112 = -t39 * t91 + t41 * t93;
t40 = Icges(5,6) * t86 - t106 * t87;
t42 = Icges(5,5) * t86 - t107 * t87;
t111 = t40 * t93 + t42 * t91;
t110 = t40 * t91 - t42 * t93;
t109 = t93 * t62 + t91 * t63;
t33 = t87 * rSges(4,3) + t153 * t86 + t76;
t105 = Icges(5,5) * t91 + Icges(5,6) * t93;
t34 = -rSges(4,2) * t87 + t86 * rSges(4,3) + t134;
t30 = t87 * pkin(6) + t134 + t43;
t45 = t50 * t89;
t104 = t111 * t86;
t103 = t107 * t89;
t102 = t106 * t89;
t101 = t105 * t89;
t98 = -t105 * qJD(4) + t109 * t89;
t100 = (-t110 / 0.2e1 - t109 * t87 / 0.2e1 + t61 * t154) * t137 + (-qJD(4) * t111 + t151 * t87 - (t86 * t102 - t149 * t87) * t91 + (t86 * t103 - t148 * t87) * t93 + t98 * t86) * t154 + (-qJD(4) * t113 - t151 * t86 - (t87 * t102 + t149 * t86) * t91 + (t87 * t103 + t148 * t86) * t93 + t98 * t87) * t87 / 0.2e1 - (t109 * t86 + t61 * t87 + t112) * t138 / 0.2e1;
t27 = rSges(4,3) * t137 + t153 * t138 + t136;
t29 = t125 * t86 + t135 + t76;
t74 = qJD(3) * t87;
t28 = rSges(4,2) * t137 + t74 + (-t83 + (-rSges(4,3) - qJ(3)) * t86) * t89;
t97 = -qJD(4) * t109 + t117;
t59 = t87 * t121;
t12 = -t87 * t120 + t59 + t74 + (t125 * t87 + (-qJ(3) - t116) * t86) * t89;
t10 = t12 - t122;
t25 = t29 - t142;
t88 = t94 * pkin(1);
t26 = t88 + t30;
t11 = (t125 * t89 - t120) * t86 - t124 + t136;
t9 = t11 - t123;
t96 = (t25 * t89 - t9) * t87 + (t26 * t89 + t10) * t86;
t95 = (t29 * t89 - t11) * t87 + (t30 * t89 + t12) * t86;
t68 = rSges(5,1) * t93 - rSges(5,2) * t91;
t48 = t51 + t88;
t47 = t50 - t142;
t44 = -t135 + t139;
t38 = Icges(5,3) * t86 - t105 * t87;
t37 = Icges(5,3) * t87 + t105 * t86;
t36 = t46 - t122;
t35 = t45 - t123;
t32 = t34 + t88;
t31 = t33 - t142;
t24 = t28 - t122;
t23 = t27 - t123;
t18 = t87 * t101 + t150 * t86;
t17 = t86 * t101 - t150 * t87;
t8 = -t111 * t87 + t86 * t38;
t7 = t37 * t86 - t152;
t6 = t38 * t87 + t104;
t5 = t113 * t86 + t37 * t87;
t1 = [(t35 * t48 + t36 * t47) * t147 + (t23 * t32 + t24 * t31) * t146 - t63 * t127 - t62 * t126 + (t10 * t25 + t26 * t9) * t145 + t117; m(3) * (t35 * t51 + t36 * t50 + t45 * t48 + t46 * t47) + m(4) * (t23 * t34 + t24 * t33 + t27 * t32 + t28 * t31) + m(5) * (t10 * t29 + t11 * t26 + t12 * t25 + t30 * t9) + t97; (t45 * t51 + t46 * t50) * t147 + (t27 * t34 + t28 * t33) * t146 + (t11 * t30 + t12 * t29) * t145 + t97; m(4) * ((t31 * t89 - t23) * t87 + (t32 * t89 + t24) * t86) + m(5) * t96; m(4) * ((t33 * t89 - t27) * t87 + (t34 * t89 + t28) * t86) + m(5) * t95; 0; m(5) * (-(t25 * t86 - t26 * t87) * t57 + t96 * t68) + t100; m(5) * (-(t29 * t86 - t30 * t87) * t57 + t95 * t68) + t100; -m(5) * t118; ((-t43 * t86 + t44 * t87) * ((-t89 * t43 - t59 + (rSges(5,3) * t89 + t120) * t87) * t87 + (t86 * t120 + (t139 - t44 + t135) * t89 + t124) * t86) - t68 * t118) * t145 - (t5 * t87 + t6 * t86) * t138 + t87 * ((t87 * t18 + (t6 + t152) * t89) * t87 + (-t5 * t89 + (t126 * t42 - t127 * t40) * t86 + (t112 * qJD(4) + t111 * t89 + t17) * t87) * t86) + (t7 * t87 + t8 * t86) * t137 + t86 * ((t86 * t17 + (-t7 + t104) * t89) * t86 + (t8 * t89 + (-t126 * t41 + t127 * t39) * t87 + (t110 * qJD(4) + t113 * t89 + t18) * t86) * t87);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
