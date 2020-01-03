% Calculate time derivative of joint inertia matrix for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:09
% DurationCPUTime: 0.97s
% Computational Cost: add. (2244->137), mult. (1784->207), div. (0->0), fcn. (1294->8), ass. (0->84)
t79 = sin(qJ(4));
t114 = qJD(4) * t79;
t81 = cos(qJ(4));
t113 = qJD(4) * t81;
t78 = qJ(1) + pkin(7);
t75 = qJ(3) + t78;
t71 = sin(t75);
t72 = cos(t75);
t115 = Icges(5,4) * t81;
t96 = -Icges(5,2) * t79 + t115;
t90 = t96 * t72;
t33 = Icges(5,6) * t71 + t90;
t116 = Icges(5,4) * t79;
t97 = Icges(5,1) * t81 - t116;
t91 = t97 * t72;
t35 = Icges(5,5) * t71 + t91;
t99 = t33 * t79 - t35 * t81;
t132 = t99 * t71;
t117 = pkin(2) * cos(t78) + cos(qJ(1)) * pkin(1);
t32 = -Icges(5,6) * t72 + t96 * t71;
t34 = -Icges(5,5) * t72 + t97 * t71;
t101 = t32 * t79 - t34 * t81;
t131 = t101 * t72;
t77 = qJD(1) + qJD(3);
t95 = Icges(5,5) * t81 - Icges(5,6) * t79;
t61 = Icges(5,2) * t81 + t116;
t62 = Icges(5,1) * t79 + t115;
t98 = t61 * t79 - t62 * t81;
t130 = t95 * qJD(4) + t98 * t77;
t129 = 2 * m(4);
t128 = 2 * m(5);
t124 = rSges(5,1) * t81;
t123 = rSges(5,2) * t79;
t64 = t71 * rSges(5,3);
t122 = t71 * t77;
t121 = t71 * t81;
t120 = t72 * t77;
t58 = t71 * t123;
t119 = rSges(5,3) * t120 + t77 * t58;
t118 = t72 * rSges(5,3) + t58;
t112 = t72 * t123;
t110 = rSges(5,2) * t113;
t111 = -t77 * t112 + (-t114 * rSges(5,1) - t110) * t71;
t107 = -pkin(3) - t124;
t43 = t72 * rSges(4,1) - rSges(4,2) * t71;
t41 = -rSges(4,1) * t120 + rSges(4,2) * t122;
t106 = -pkin(2) * sin(t78) - sin(qJ(1)) * pkin(1);
t42 = -rSges(4,1) * t71 - rSges(4,2) * t72;
t63 = rSges(5,1) * t79 + rSges(5,2) * t81;
t102 = t32 * t81 + t34 * t79;
t100 = t33 * t81 + t35 * t79;
t60 = Icges(5,5) * t79 + Icges(5,6) * t81;
t37 = t72 * t124 - t112 + t64;
t40 = t42 * t77;
t94 = t106 * qJD(1);
t93 = t117 * qJD(1);
t92 = (-t61 + t97) * t114 + (t62 + t96) * t113;
t89 = t95 * t72;
t88 = (-t99 * qJD(4) + t130 * t71 + (-t79 * t97 - t81 * t96) * t122) * t71 / 0.2e1 - (-t101 * qJD(4) - t130 * t72 + (t79 * t91 + t81 * t90) * t77) * t72 / 0.2e1 + (-t72 * t60 - t98 * t71 + t102) * t122 / 0.2e1 + (t71 * t60 - t98 * t72 + t100) * t120 / 0.2e1;
t27 = t72 * pkin(3) + t71 * pkin(6) + t37;
t26 = t72 * pkin(6) + t107 * t71 + t118;
t84 = Icges(5,3) * t77 - t60 * qJD(4);
t13 = (t107 * t72 + (-rSges(5,3) - pkin(6)) * t71) * t77 - t111;
t12 = -t72 * t110 - pkin(3) * t122 + pkin(6) * t120 + (-t72 * t114 - t77 * t121) * rSges(5,1) + t119;
t52 = (-t123 + t124) * qJD(4);
t39 = t43 + t117;
t38 = t106 + t42;
t36 = rSges(5,1) * t121 - t118;
t31 = Icges(5,3) * t71 + t89;
t30 = -Icges(5,3) * t72 + t95 * t71;
t29 = t41 - t93;
t28 = t40 + t94;
t25 = t27 + t117;
t24 = t106 + t26;
t19 = t84 * t71 + t77 * t89;
t18 = -t95 * t122 + t84 * t72;
t11 = -t93 + t13;
t10 = t94 + t12;
t9 = t71 * t31 - t99 * t72;
t8 = t71 * t30 - t131;
t7 = -t72 * t31 - t132;
t6 = -t101 * t71 - t72 * t30;
t1 = ((-t37 + t64) * t77 + t111) * t71 + (-t63 * t72 * qJD(4) + t77 * t36 + t119) * t72;
t2 = [(t28 * t39 + t29 * t38) * t129 + (t10 * t25 + t11 * t24) * t128 + t92; 0; 0; m(4) * (t28 * t43 + t29 * t42 + t38 * t41 + t39 * t40) + m(5) * (t10 * t27 + t11 * t26 + t12 * t25 + t13 * t24) + t92; 0; (t40 * t43 + t41 * t42) * t129 + (t12 * t27 + t13 * t26) * t128 + t92; m(5) * ((-t24 * t72 - t25 * t71) * t52 + ((-t25 * t77 - t11) * t72 + (t24 * t77 - t10) * t71) * t63) + t88; m(5) * t1; m(5) * ((-t26 * t72 - t27 * t71) * t52 + ((-t27 * t77 - t13) * t72 + (t26 * t77 - t12) * t71) * t63) + t88; ((t36 * t71 + t37 * t72) * t1 + (t71 ^ 2 + t72 ^ 2) * t63 * t52) * t128 + (t9 * t71 - t72 * t8) * t120 + t71 * ((t71 * t18 + (t8 + t132) * t77) * t71 + (t9 * t77 + (t32 * t113 + t34 * t114) * t72 + (-t100 * qJD(4) - t101 * t77 - t19) * t71) * t72) + (-t6 * t72 + t7 * t71) * t122 - t72 * ((t72 * t19 + (t7 + t131) * t77) * t72 + (t6 * t77 + (-t33 * t113 - t35 * t114) * t71 + (t102 * qJD(4) - t99 * t77 - t18) * t72) * t71);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
