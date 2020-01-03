% Calculate time derivative of joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:51
% DurationCPUTime: 3.04s
% Computational Cost: add. (2131->216), mult. (2932->317), div. (0->0), fcn. (2232->6), ass. (0->121)
t163 = rSges(5,3) + qJ(4);
t181 = rSges(5,1) + pkin(3);
t95 = pkin(6) + qJ(3);
t89 = sin(t95);
t90 = cos(t95);
t190 = t163 * t89 + t181 * t90;
t189 = t163 * t90;
t175 = t181 * t89 - t189;
t102 = cos(qJ(1));
t166 = t102 * t89;
t101 = sin(qJ(1));
t93 = t101 * rSges(4,3);
t188 = -rSges(4,2) * t166 + t93;
t171 = Icges(4,4) * t90;
t128 = -Icges(4,2) * t89 + t171;
t46 = Icges(4,6) * t101 + t102 * t128;
t172 = Icges(4,4) * t89;
t132 = Icges(4,1) * t90 - t172;
t50 = Icges(4,5) * t101 + t102 * t132;
t133 = t46 * t89 - t50 * t90;
t114 = t133 * t101;
t45 = -Icges(4,6) * t102 + t101 * t128;
t49 = -Icges(4,5) * t102 + t101 * t132;
t134 = t45 * t89 - t49 * t90;
t115 = t134 * t102;
t169 = Icges(5,5) * t90;
t124 = Icges(5,3) * t89 + t169;
t40 = Icges(5,6) * t101 + t102 * t124;
t170 = Icges(5,5) * t89;
t130 = Icges(5,1) * t90 + t170;
t48 = Icges(5,4) * t101 + t102 * t130;
t135 = t40 * t89 + t48 * t90;
t116 = t135 * t101;
t39 = -Icges(5,6) * t102 + t101 * t124;
t47 = -Icges(5,4) * t102 + t101 * t130;
t136 = t39 * t89 + t47 * t90;
t117 = t136 * t102;
t100 = -pkin(5) - qJ(2);
t179 = rSges(4,2) * t89;
t180 = rSges(4,1) * t90;
t143 = -t179 + t180;
t99 = cos(pkin(6));
t86 = pkin(2) * t99 + pkin(1);
t119 = -t143 - t86;
t33 = (rSges(4,3) - t100) * t102 + t119 * t101;
t144 = -t101 * t100 + t102 * t86;
t165 = t102 * t90;
t55 = rSges(4,1) * t165 + t188;
t34 = t144 + t55;
t187 = t101 * t34 + t102 * t33;
t125 = Icges(4,5) * t90 - Icges(4,6) * t89;
t41 = -Icges(4,3) * t102 + t101 * t125;
t126 = Icges(5,4) * t90 + Icges(5,6) * t89;
t43 = -Icges(5,2) * t102 + t101 * t126;
t120 = rSges(3,1) * t99 - rSges(3,2) * sin(pkin(6)) + pkin(1);
t164 = rSges(3,3) + qJ(2);
t38 = t101 * t164 + t102 * t120;
t186 = 2 * m(4);
t185 = 2 * m(5);
t96 = t101 ^ 2;
t97 = t102 ^ 2;
t178 = -t190 * qJD(3) + qJD(4) * t90;
t177 = -rSges(5,2) * t102 + t190 * t101;
t94 = t101 * rSges(5,2);
t176 = t163 * t166 + t181 * t165 + t94;
t150 = qJD(1) * t101;
t92 = qJD(2) * t102;
t174 = t100 * t150 + t92;
t173 = t96 + t97;
t42 = Icges(4,3) * t101 + t102 * t125;
t156 = qJD(1) * t42;
t44 = Icges(5,2) * t101 + t102 * t126;
t155 = qJD(1) * t44;
t154 = qJD(3) * t89;
t153 = qJD(3) * t90;
t152 = qJD(4) * t89;
t151 = t100 * t102;
t149 = qJD(1) * t102;
t148 = qJD(3) * t101;
t145 = t181 * qJD(3);
t36 = t175 * t102;
t74 = rSges(4,1) * t89 + rSges(4,2) * t90;
t106 = -t190 - t86;
t103 = t106 * t101;
t17 = (rSges(5,2) - t100) * t102 + t103;
t18 = t144 + t176;
t122 = t101 * t18 + t102 * t17;
t35 = t175 * t101;
t121 = -t101 * t35 - t102 * t36;
t113 = qJD(3) * t74;
t109 = qJD(3) * (-Icges(5,4) * t89 + Icges(5,6) * t90);
t108 = qJD(3) * (-Icges(4,5) * t89 - Icges(4,6) * t90);
t105 = rSges(5,2) * t149 - t145 * t166 + (t153 * t163 + t152) * t102;
t104 = rSges(4,3) * t149 - t102 * t113 + t150 * t179;
t37 = -t101 * t120 + t102 * t164;
t91 = qJD(2) * t101;
t65 = t143 * qJD(3);
t53 = -rSges(4,3) * t102 + t101 * t143;
t32 = -t38 * qJD(1) + t92;
t31 = qJD(1) * t37 + t91;
t24 = t101 * t109 + t155;
t23 = -qJD(1) * t43 + t102 * t109;
t22 = t101 * t108 + t156;
t21 = -qJD(1) * t41 + t102 * t108;
t16 = t74 * t148 + (t102 * t119 - t93) * qJD(1) + t174;
t15 = t91 + (-t151 + (-t86 - t180) * t101) * qJD(1) + t104;
t14 = -qJD(1) * t36 + t101 * t178;
t13 = t102 * t178 + t150 * t175;
t12 = t101 * t42 - t133 * t102;
t11 = t101 * t41 - t115;
t10 = t101 * t44 + t135 * t102;
t9 = t101 * t43 + t117;
t8 = -t102 * t42 - t114;
t7 = -t101 * t134 - t102 * t41;
t6 = -t102 * t44 + t116;
t5 = t101 * t136 - t102 * t43;
t4 = t101 * t177 + t102 * t176;
t3 = (t175 * qJD(3) - t152) * t101 + (t102 * t106 - t94) * qJD(1) + t174;
t2 = t91 + (t103 - t151) * qJD(1) + t105;
t1 = (qJD(1) * t177 + t105) * t102 + (t148 * t189 + (t94 - t176) * qJD(1) + (qJD(4) - t145) * t101 * t89) * t101;
t19 = [0.2e1 * m(3) * (t31 * t38 + t32 * t37) + (t15 * t34 + t16 * t33) * t186 + (t17 * t3 + t18 * t2) * t185 + (t130 + t132 + t170 - t172 + (-Icges(4,2) - Icges(5,3)) * t90) * t154 + (-t124 + t128 - t169 + t171 + (Icges(4,1) + Icges(5,1)) * t89) * t153; m(3) * (t101 * t32 - t102 * t31 + (t101 * t38 + t102 * t37) * qJD(1)) + m(4) * (t187 * qJD(1) + t101 * t16 - t102 * t15) + m(5) * (qJD(1) * t122 + t101 * t3 - t102 * t2); 0; m(5) * (t13 * t17 + t14 * t18 - t2 * t35 - t3 * t36) + (-t114 / 0.2e1 + t116 / 0.2e1 + t115 / 0.2e1 - t117 / 0.2e1 + (t125 + t126) * (t96 / 0.2e1 + t97 / 0.2e1)) * qJD(3) + (-t187 * t65 + (-t101 * t15 - t102 * t16 + (t33 * t101 - t34 * t102) * qJD(1)) * t74) * m(4); m(5) * (qJD(1) * t121 + t13 * t101 - t102 * t14); ((t101 * t53 + t102 * t55) * ((qJD(1) * t53 + t104) * t102 + (-t101 * t113 + (-t55 + t188) * qJD(1)) * t101) + t173 * t74 * t65) * t186 + t101 * ((t101 * t21 + (t11 + t114) * qJD(1)) * t101 + (t12 * qJD(1) + (t153 * t45 + t154 * t49) * t102 + (-t22 + (-t46 * t90 - t50 * t89) * qJD(3) + (-t134 + t42) * qJD(1)) * t101) * t102) - t102 * ((t102 * t22 + (t8 + t115) * qJD(1)) * t102 + (t7 * qJD(1) + (-t153 * t46 - t154 * t50 + t156) * t101 + (-t21 + (t45 * t90 + t49 * t89) * qJD(3) - t133 * qJD(1)) * t102) * t101) + (t1 * t4 - t13 * t36 - t14 * t35) * t185 + t101 * ((t101 * t23 + (t9 - t116) * qJD(1)) * t101 + (t10 * qJD(1) + (-t153 * t39 + t154 * t47) * t102 + (-t24 + (t40 * t90 - t48 * t89) * qJD(3) + (t136 + t44) * qJD(1)) * t101) * t102) - t102 * ((t102 * t24 + (t6 - t117) * qJD(1)) * t102 + (t5 * qJD(1) + (t153 * t40 - t154 * t48 + t155) * t101 + (-t23 + (-t39 * t90 + t47 * t89) * qJD(3) + t135 * qJD(1)) * t102) * t101) + ((-t5 - t7) * t102 + (t6 + t8) * t101) * t150 + ((-t11 - t9) * t102 + (t10 + t12) * t101) * t149; m(5) * (t122 * t153 + (t101 * t2 + t102 * t3 + (-t101 * t17 + t102 * t18) * qJD(1)) * t89); 0; m(5) * ((qJD(3) * t121 - t1) * t90 + (qJD(3) * t4 + t101 * t14 + t102 * t13 + (t101 * t36 - t102 * t35) * qJD(1)) * t89); (-0.1e1 + t173) * t89 * t153 * t185;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1), t19(2), t19(4), t19(7); t19(2), t19(3), t19(5), t19(8); t19(4), t19(5), t19(6), t19(9); t19(7), t19(8), t19(9), t19(10);];
Mq = res;
