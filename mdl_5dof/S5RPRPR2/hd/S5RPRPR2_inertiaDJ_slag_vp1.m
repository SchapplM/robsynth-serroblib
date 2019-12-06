% Calculate time derivative of joint inertia matrix for
% S5RPRPR2
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:21
% EndTime: 2019-12-05 17:49:25
% DurationCPUTime: 1.49s
% Computational Cost: add. (3646->176), mult. (2294->249), div. (0->0), fcn. (1646->10), ass. (0->105)
t98 = pkin(9) + qJ(5);
t93 = sin(t98);
t142 = qJD(5) * t93;
t95 = cos(t98);
t141 = qJD(5) * t95;
t118 = Icges(6,5) * t95 - Icges(6,6) * t93;
t148 = Icges(6,4) * t95;
t119 = -Icges(6,2) * t93 + t148;
t149 = Icges(6,4) * t93;
t120 = Icges(6,1) * t95 - t149;
t63 = Icges(6,2) * t95 + t149;
t64 = Icges(6,1) * t93 + t148;
t121 = t63 * t93 - t64 * t95;
t99 = qJD(1) + qJD(3);
t183 = t118 * qJD(5) + (-t119 * t95 - t120 * t93 + t121) * t99;
t150 = rSges(5,2) * sin(pkin(9));
t176 = rSges(5,3) + qJ(4);
t100 = qJ(1) + pkin(8);
t97 = qJ(3) + t100;
t89 = sin(t97);
t90 = cos(t97);
t182 = -t89 * t150 - t176 * t90;
t179 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t100);
t102 = cos(pkin(9));
t177 = rSges(5,1) * t102 + pkin(3);
t35 = -t176 * t89 + (-t177 + t150) * t90;
t178 = t35 * t99;
t103 = -pkin(7) - qJ(4);
t162 = rSges(6,1) * t95;
t91 = pkin(4) * t102 + pkin(3);
t132 = -t91 - t162;
t156 = t89 * rSges(6,3);
t175 = t89 * t103 + t132 * t90 - t156;
t174 = rSges(6,1) * t142 + rSges(6,2) * t141;
t40 = Icges(6,6) * t90 - t119 * t89;
t42 = Icges(6,5) * t90 - t120 * t89;
t124 = t40 * t93 - t42 * t95;
t173 = t124 * t90;
t62 = Icges(6,5) * t93 + Icges(6,6) * t95;
t172 = -Icges(6,3) * t99 + t62 * qJD(5);
t168 = 2 * m(4);
t167 = 2 * m(5);
t166 = 2 * m(6);
t161 = rSges(6,2) * t93;
t155 = t89 * t99;
t87 = t90 * rSges(6,3);
t154 = t90 * t99;
t153 = t89 * t161 + t87;
t48 = rSges(4,1) * t155 + rSges(4,2) * t154;
t152 = t179 * qJD(1);
t144 = t103 * t90;
t140 = t89 * t162;
t75 = t90 * t161;
t139 = -t99 * t140 - t174 * t90;
t138 = t174 * t89 + t99 * t75;
t58 = -t90 * rSges(4,1) + t89 * rSges(4,2);
t12 = -rSges(6,3) * t154 + t91 * t155 + t99 * t144 + (-t99 * t161 - qJD(4)) * t89 - t139;
t10 = t12 + t152;
t32 = t132 * t89 - t144 + t153;
t28 = -t179 + t32;
t131 = t28 * t99 + t10;
t126 = -pkin(2) * cos(t100) - cos(qJ(1)) * pkin(1);
t115 = t126 * qJD(1);
t83 = qJD(4) * t90;
t13 = t175 * t99 + t138 + t83;
t11 = t115 + t13;
t33 = t75 + t175;
t29 = t126 + t33;
t130 = t29 * t99 - t11;
t129 = t32 * t99 + t12;
t128 = t33 * t99 - t13;
t49 = -rSges(4,1) * t154 + rSges(4,2) * t155;
t57 = -rSges(4,1) * t89 - rSges(4,2) * t90;
t41 = Icges(6,6) * t89 + t119 * t90;
t43 = Icges(6,5) * t89 + t120 * t90;
t123 = t41 * t93 - t43 * t95;
t117 = t123 * t89;
t116 = (t120 - t63) * t142 + (t119 + t64) * t141;
t112 = t118 * t99;
t111 = (-t123 * qJD(5) + t183 * t89) * t89 / 0.2e1 + (-t124 * qJD(5) + t183 * t90) * t90 / 0.2e1 - (t121 * t89 + t40 * t95 + t42 * t93 + t62 * t90) * t155 / 0.2e1 + (-t121 * t90 + t41 * t95 + t43 * t93 + t62 * t89) * t154 / 0.2e1;
t34 = -t177 * t89 - t182;
t27 = t83 + t178;
t26 = -qJD(4) * t89 + t177 * t155 + t182 * t99;
t71 = rSges(6,1) * t93 + rSges(6,2) * t95;
t56 = (-t161 + t162) * qJD(5);
t47 = t126 + t58;
t46 = -t179 + t57;
t45 = t90 * t162 + t156 - t75;
t44 = -t140 + t153;
t39 = Icges(6,3) * t89 + t118 * t90;
t38 = Icges(6,3) * t90 - t118 * t89;
t37 = t115 + t49;
t36 = t152 + t48;
t31 = t126 + t35;
t30 = -t179 + t34;
t21 = -t90 * t112 + t172 * t89;
t20 = -t89 * t112 - t172 * t90;
t19 = t115 + t27;
t18 = t26 + t152;
t9 = -t123 * t90 + t39 * t89;
t8 = t38 * t89 - t173;
t7 = t39 * t90 + t117;
t6 = t124 * t89 + t38 * t90;
t3 = t90 * t139 - t89 * t138 + ((-t44 + t87) * t90 + (t156 - t45 + (t161 + t162) * t90) * t89) * t99;
t1 = [(t36 * t47 + t37 * t46) * t168 + (t18 * t31 + t19 * t30) * t167 + (t10 * t29 + t11 * t28) * t166 + t116; 0; 0; m(4) * (t36 * t58 + t37 * t57 + t49 * t46 + t48 * t47) + m(5) * (t18 * t35 + t19 * t34 + t26 * t31 + t27 * t30) + m(6) * (t10 * t33 + t11 * t32 + t12 * t29 + t13 * t28) + t116; 0; (t26 * t35 + t27 * t34) * t167 + (t12 * t33 + t13 * t32) * t166 + (t48 * t58 + t49 * t57) * t168 + t116; m(5) * ((t30 * t99 + t18) * t90 + (-t31 * t99 + t19) * t89) + m(6) * (-t130 * t89 + t131 * t90); 0; m(5) * ((t34 * t99 + t26) * t90 + (t27 - t178) * t89) + m(6) * (-t128 * t89 + t129 * t90); 0; m(6) * ((-t28 * t90 + t29 * t89) * t56 + (t130 * t90 + t131 * t89) * t71) + t111; m(6) * t3; m(6) * ((-t32 * t90 + t33 * t89) * t56 + (t128 * t90 + t129 * t89) * t71) + t111; 0; ((-t44 * t89 + t45 * t90) * t3 + (t89 ^ 2 + t90 ^ 2) * t71 * t56) * t166 - (t6 * t90 + t7 * t89) * t155 + t90 * ((t90 * t21 + (t7 + t173) * t99) * t90 + (-t6 * t99 + (t41 * t141 + t43 * t142) * t89 + (t20 + (qJD(5) * t40 - t43 * t99) * t95 + (qJD(5) * t42 + t41 * t99) * t93) * t90) * t89) + (t8 * t90 + t9 * t89) * t154 + t89 * ((t89 * t20 + (-t8 + t117) * t99) * t89 + (t9 * t99 + (-t40 * t141 - t42 * t142) * t90 + (t21 + (-qJD(5) * t41 - t42 * t99) * t95 + (-qJD(5) * t43 + t40 * t99) * t93) * t89) * t90);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
