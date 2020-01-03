% Calculate time derivative of joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:24
% DurationCPUTime: 1.61s
% Computational Cost: add. (1786->181), mult. (3218->272), div. (0->0), fcn. (3152->8), ass. (0->95)
t114 = sin(pkin(7));
t115 = cos(pkin(7));
t130 = sin(qJ(1));
t131 = cos(qJ(1));
t54 = -t130 * t114 - t131 * t115;
t55 = t131 * t114 - t130 * t115;
t71 = pkin(8) + qJ(5);
t63 = sin(t71);
t64 = cos(t71);
t92 = -Icges(6,5) * t64 + Icges(6,6) * t63;
t25 = Icges(6,3) * t55 + t92 * t54;
t147 = t25 * t54;
t24 = -Icges(6,3) * t54 + t92 * t55;
t148 = t24 * t55;
t149 = t147 - t148;
t146 = t54 ^ 2;
t145 = t55 ^ 2;
t142 = t55 * t54;
t119 = t131 * pkin(1) + t130 * qJ(2);
t110 = t131 * pkin(2) + t119;
t39 = t131 * rSges(3,1) + t130 * rSges(3,3) + t119;
t113 = qJD(5) * t63;
t108 = t54 * t113;
t46 = t55 * qJD(1);
t86 = t46 * t64 + t108;
t106 = t55 * t113;
t45 = t54 * qJD(1);
t88 = -t45 * t64 + t106;
t136 = 2 * m(6);
t129 = rSges(6,2) * t63;
t124 = t54 * t63;
t123 = t54 * t64;
t74 = -pkin(6) - qJ(4);
t122 = rSges(6,3) - t74;
t121 = -rSges(6,1) * t123 + t55 * rSges(6,3);
t68 = t131 * qJ(2);
t120 = qJD(1) * t68 + qJD(2) * t130;
t118 = Icges(6,4) * t63;
t117 = Icges(6,4) * t64;
t116 = rSges(5,3) + qJ(4);
t112 = qJD(5) * t64;
t111 = m(5) / 0.2e1 + m(6) / 0.2e1;
t109 = t130 * pkin(1);
t107 = t54 * t112;
t105 = t55 * t112;
t100 = -rSges(6,1) * t64 + t129;
t53 = -rSges(6,1) * t63 - rSges(6,2) * t64;
t97 = t45 * t55 - t46 * t54;
t94 = -Icges(6,1) * t64 + t118;
t93 = Icges(6,2) * t63 - t117;
t73 = cos(pkin(8));
t91 = rSges(5,1) * t73 - rSges(5,2) * sin(pkin(8)) + pkin(3);
t61 = pkin(4) * t73 + pkin(3);
t90 = -t100 + t61;
t89 = t45 * t63 + t105;
t87 = -t46 * t63 + t107;
t85 = rSges(6,1) * t86 + t45 * rSges(6,3);
t84 = -t130 * pkin(2) - t109;
t81 = -t130 * t54 + t131 * t55;
t80 = t68 + t84;
t79 = -t130 * rSges(3,1) + t131 * rSges(3,3) - t109;
t78 = t84 * qJD(1) + t120;
t77 = t55 * qJD(4) + t78;
t66 = qJD(2) * t131;
t76 = -qJD(1) * t110 + t66;
t75 = t54 * qJD(4) + t76;
t52 = -Icges(6,1) * t63 - t117;
t51 = -Icges(6,2) * t64 - t118;
t44 = t100 * qJD(5);
t41 = t92 * qJD(5);
t38 = t68 + t79;
t35 = -qJD(1) * t39 + t66;
t34 = t79 * qJD(1) + t120;
t33 = -rSges(4,1) * t54 - rSges(4,2) * t55 + t110;
t32 = t55 * rSges(4,1) - t54 * rSges(4,2) + t80;
t31 = rSges(6,2) * t124 + t121;
t30 = -t54 * rSges(6,3) + t100 * t55;
t29 = Icges(6,5) * t55 + t94 * t54;
t28 = -Icges(6,5) * t54 + t94 * t55;
t27 = Icges(6,6) * t55 + t93 * t54;
t26 = -Icges(6,6) * t54 + t93 * t55;
t23 = t45 * rSges(4,1) + t46 * rSges(4,2) + t76;
t22 = t46 * rSges(4,1) - t45 * rSges(4,2) + t78;
t19 = t116 * t55 - t91 * t54 + t110;
t18 = t116 * t54 + t91 * t55 + t80;
t17 = -t55 * t74 + (-t61 + t129) * t54 + t110 + t121;
t16 = t122 * t54 + t90 * t55 + t80;
t11 = t86 * Icges(6,5) + t87 * Icges(6,6) + Icges(6,3) * t45;
t10 = t88 * Icges(6,5) + t89 * Icges(6,6) + Icges(6,3) * t46;
t9 = -t116 * t46 + t91 * t45 + t75;
t8 = t116 * t45 + t91 * t46 + t77;
t7 = t53 * t55 * qJD(5) - t122 * t46 + t90 * t45 + t75;
t6 = t87 * rSges(6,2) - t45 * t74 + t46 * t61 + t77 + t85;
t1 = t45 * t30 + t55 * (t88 * rSges(6,1) + t46 * rSges(6,3)) - t46 * t31 + t54 * t85 + (t54 * t87 + t55 * t89) * rSges(6,2);
t2 = [0.2e1 * m(3) * (t34 * t39 + t35 * t38) + 0.2e1 * m(4) * (t22 * t33 + t23 * t32) + 0.2e1 * m(5) * (t18 * t9 + t19 * t8) + (t16 * t7 + t17 * t6) * t136 + (-t94 + t51) * t113 + (-t52 - t93) * t112; m(3) * (t130 * t35 - t131 * t34 + (t130 * t39 + t131 * t38) * qJD(1)) + m(4) * (t130 * t23 - t131 * t22 + (t130 * t33 + t131 * t32) * qJD(1)) + m(5) * (t130 * t9 - t131 * t8 + (t130 * t19 + t131 * t18) * qJD(1)) + m(6) * (t130 * t7 - t131 * t6 + (t130 * t17 + t131 * t16) * qJD(1)); 0; 0; 0; 0; m(5) * (t45 * t18 + t46 * t19 - t54 * t8 + t55 * t9) + m(6) * (t45 * t16 + t46 * t17 - t54 * t6 + t55 * t7); 0.2e1 * t111 * (t81 * qJD(1) + t45 * t130 - t46 * t131); 0; 0.4e1 * t111 * t97; m(6) * ((-t16 * t54 - t17 * t55) * t44 + (t16 * t46 - t17 * t45 - t54 * t7 - t55 * t6) * t53) + t97 * (-Icges(6,5) * t63 - Icges(6,6) * t64) + (-t26 * t64 - t28 * t63) * t46 / 0.2e1 - ((t26 * t63 - t28 * t64) * qJD(5) + t51 * t105 + t52 * t106 - (t88 * Icges(6,4) + t89 * Icges(6,2) + Icges(6,6) * t46) * t64 - (t88 * Icges(6,1) + t89 * Icges(6,4) + Icges(6,5) * t46) * t63 - t41 * t54) * t54 / 0.2e1 + ((t27 * t63 - t29 * t64) * qJD(5) + t51 * t107 + t52 * t108 - (t86 * Icges(6,4) + t87 * Icges(6,2) + Icges(6,6) * t45) * t64 - (t86 * Icges(6,1) + t87 * Icges(6,4) + Icges(6,5) * t45) * t63 + t41 * t55) * t55 / 0.2e1 + (-t52 * t123 + t51 * t124 - t27 * t64 - t29 * t63 + t54 * (-t51 * t63 + t52 * t64)) * t45 / 0.2e1; m(6) * (t81 * t44 + (t130 * t46 + t131 * t45 + (-t130 * t55 - t131 * t54) * qJD(1)) * t53); -m(6) * t1; 0; ((t30 * t55 + t31 * t54) * t1 + ((t145 + t146) * t44 + t97 * t53) * t53) * t136 + t55 * (-t10 * t142 + t11 * t145) - t54 * (t10 * t146 - t11 * t142) + (0.3e1 * t24 * t146 + (-t147 - t149) * t55) * t46 + (0.3e1 * t25 * t145 + (-t148 + t149) * t54) * t45;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
