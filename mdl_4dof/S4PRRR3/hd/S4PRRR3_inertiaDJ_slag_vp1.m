% Calculate time derivative of joint inertia matrix for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:36
% DurationCPUTime: 0.95s
% Computational Cost: add. (2212->135), mult. (1728->206), div. (0->0), fcn. (1260->6), ass. (0->87)
t78 = sin(qJ(4));
t108 = qJD(4) * t78;
t79 = cos(qJ(4));
t107 = qJD(4) * t79;
t76 = pkin(7) + qJ(2);
t75 = qJ(3) + t76;
t71 = sin(t75);
t72 = cos(t75);
t110 = Icges(5,4) * t79;
t90 = -Icges(5,2) * t78 + t110;
t86 = t90 * t72;
t33 = Icges(5,6) * t71 + t86;
t111 = Icges(5,4) * t78;
t91 = Icges(5,1) * t79 - t111;
t87 = t91 * t72;
t35 = Icges(5,5) * t71 + t87;
t93 = t33 * t78 - t35 * t79;
t126 = t93 * t71;
t32 = -Icges(5,6) * t72 + t90 * t71;
t34 = -Icges(5,5) * t72 + t91 * t71;
t95 = t32 * t78 - t34 * t79;
t125 = t95 * t72;
t77 = qJD(2) + qJD(3);
t89 = Icges(5,5) * t79 - Icges(5,6) * t78;
t61 = Icges(5,2) * t79 + t111;
t62 = Icges(5,1) * t78 + t110;
t92 = t61 * t78 - t62 * t79;
t124 = t89 * qJD(4) + t92 * t77;
t123 = 2 * m(4);
t122 = 2 * m(5);
t73 = sin(t76);
t119 = pkin(2) * t73;
t118 = rSges(5,1) * t79;
t117 = rSges(5,2) * t78;
t64 = t71 * rSges(5,3);
t116 = t71 * t77;
t115 = t71 * t79;
t114 = t72 * t77;
t58 = t71 * t117;
t113 = rSges(5,3) * t114 + t77 * t58;
t112 = t72 * rSges(5,3) + t58;
t109 = pkin(2) * qJD(2);
t106 = t72 * t117;
t102 = rSges(5,2) * t107;
t105 = -t77 * t106 + (-t108 * rSges(5,1) - t102) * t71;
t104 = t73 * t109;
t74 = cos(t76);
t103 = t74 * t109;
t99 = -pkin(3) - t118;
t43 = t72 * rSges(4,1) - rSges(4,2) * t71;
t39 = -rSges(4,1) * t114 + rSges(4,2) * t116;
t42 = -rSges(4,1) * t71 - rSges(4,2) * t72;
t63 = t78 * rSges(5,1) + rSges(5,2) * t79;
t96 = t32 * t79 + t78 * t34;
t94 = t33 * t79 + t78 * t35;
t60 = Icges(5,5) * t78 + Icges(5,6) * t79;
t37 = t72 * t118 - t106 + t64;
t38 = t42 * t77;
t88 = (-t61 + t91) * t108 + (t62 + t90) * t107;
t85 = t89 * t72;
t84 = (-t93 * qJD(4) + t124 * t71 + (-t78 * t91 - t79 * t90) * t116) * t71 / 0.2e1 - (-t95 * qJD(4) - t124 * t72 + (t78 * t87 + t79 * t86) * t77) * t72 / 0.2e1 + (-t72 * t60 - t92 * t71 + t96) * t116 / 0.2e1 + (t71 * t60 - t92 * t72 + t94) * t114 / 0.2e1;
t27 = t72 * pkin(3) + t71 * pkin(6) + t37;
t26 = t72 * pkin(6) + t99 * t71 + t112;
t81 = Icges(5,3) * t77 - t60 * qJD(4);
t13 = (t99 * t72 + (-rSges(5,3) - pkin(6)) * t71) * t77 - t105;
t12 = -t72 * t102 - pkin(3) * t116 + pkin(6) * t114 + (-t72 * t108 - t77 * t115) * rSges(5,1) + t113;
t70 = pkin(2) * t74;
t52 = (-t117 + t118) * qJD(4);
t41 = t43 + t70;
t40 = t42 - t119;
t36 = rSges(5,1) * t115 - t112;
t31 = Icges(5,3) * t71 + t85;
t30 = -Icges(5,3) * t72 + t89 * t71;
t29 = t39 - t103;
t28 = t38 - t104;
t25 = t27 + t70;
t24 = t26 - t119;
t19 = t81 * t71 + t77 * t85;
t18 = -t89 * t116 + t81 * t72;
t11 = t13 - t103;
t10 = t12 - t104;
t9 = t71 * t31 - t93 * t72;
t8 = t71 * t30 - t125;
t7 = -t72 * t31 - t126;
t6 = -t72 * t30 - t95 * t71;
t1 = ((-t37 + t64) * t77 + t105) * t71 + (-t63 * t72 * qJD(4) + t77 * t36 + t113) * t72;
t2 = [0; 0; (t28 * t41 + t29 * t40) * t123 + (t10 * t25 + t11 * t24) * t122 + t88; 0; m(4) * (t28 * t43 + t29 * t42 + t38 * t41 + t39 * t40) + m(5) * (t10 * t27 + t11 * t26 + t12 * t25 + t13 * t24) + t88; (t38 * t43 + t39 * t42) * t123 + (t12 * t27 + t13 * t26) * t122 + t88; m(5) * t1; m(5) * ((-t24 * t72 - t25 * t71) * t52 + ((-t25 * t77 - t11) * t72 + (t24 * t77 - t10) * t71) * t63) + t84; m(5) * ((-t26 * t72 - t27 * t71) * t52 + ((-t27 * t77 - t13) * t72 + (t26 * t77 - t12) * t71) * t63) + t84; ((t36 * t71 + t37 * t72) * t1 + (t71 ^ 2 + t72 ^ 2) * t63 * t52) * t122 + (t9 * t71 - t72 * t8) * t114 + t71 * ((t71 * t18 + (t8 + t126) * t77) * t71 + (t9 * t77 + (t32 * t107 + t34 * t108) * t72 + (-t94 * qJD(4) - t95 * t77 - t19) * t71) * t72) + (-t6 * t72 + t7 * t71) * t116 - t72 * ((t72 * t19 + (t7 + t125) * t77) * t72 + (t6 * t77 + (-t33 * t107 - t35 * t108) * t71 + (t96 * qJD(4) - t93 * t77 - t18) * t72) * t71);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
