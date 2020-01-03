% Calculate time derivative of joint inertia matrix for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:50
% DurationCPUTime: 2.22s
% Computational Cost: add. (2120->179), mult. (2714->258), div. (0->0), fcn. (2069->4), ass. (0->113)
t209 = Icges(5,4) + Icges(4,5);
t208 = -Icges(4,6) + Icges(5,6);
t207 = Icges(5,2) + Icges(4,3);
t94 = sin(qJ(3));
t95 = cos(qJ(3));
t205 = t208 * t94 + t209 * t95;
t168 = rSges(5,1) + pkin(3);
t206 = t168 * t95;
t93 = pkin(6) + qJ(2);
t91 = sin(t93);
t92 = cos(t93);
t203 = t205 * t92 + t207 * t91;
t204 = t205 * t91 - t207 * t92;
t150 = Icges(5,5) * t95;
t113 = Icges(5,3) * t94 + t150;
t36 = -Icges(5,6) * t92 + t113 * t91;
t154 = Icges(4,4) * t95;
t117 = -Icges(4,2) * t94 + t154;
t42 = -Icges(4,6) * t92 + t117 * t91;
t202 = t42 - t36;
t37 = Icges(5,6) * t91 + t113 * t92;
t43 = Icges(4,6) * t91 + t117 * t92;
t201 = t43 - t37;
t151 = Icges(5,5) * t94;
t119 = Icges(5,1) * t95 + t151;
t44 = -Icges(5,4) * t92 + t119 * t91;
t155 = Icges(4,4) * t94;
t121 = Icges(4,1) * t95 - t155;
t46 = -Icges(4,5) * t92 + t121 * t91;
t200 = t46 + t44;
t45 = Icges(5,4) * t91 + t119 * t92;
t47 = Icges(4,5) * t91 + t121 * t92;
t199 = t47 + t45;
t122 = t43 * t94 - t47 * t95;
t124 = t37 * t94 + t45 * t95;
t175 = t122 - t124;
t198 = t175 * t92;
t197 = (t208 * t95 - t209 * t94) * qJD(3);
t123 = t42 * t94 - t46 * t95;
t125 = t36 * t94 + t44 * t95;
t196 = t123 - t125;
t144 = rSges(5,3) + qJ(4);
t158 = -t144 * t95 + t168 * t94;
t145 = qJ(4) * t94;
t193 = rSges(5,3) * t94 + t145 + t206;
t192 = t203 * qJD(2);
t106 = t122 * t91;
t107 = t123 * t92;
t108 = t124 * t91;
t109 = t125 * t92;
t191 = t203 * t92 - t204 * t91 + t106 + t107 - t108 - t109;
t190 = t199 * t91 + t200 * t92;
t189 = t201 * t91 + t202 * t92;
t89 = t91 ^ 2;
t90 = t92 ^ 2;
t188 = t196 * t91 + t204 * t92;
t173 = 2 * m(4);
t75 = t94 * rSges(4,1) + rSges(4,2) * t95;
t105 = t75 * qJD(3);
t140 = qJD(2) * t92;
t141 = qJD(2) * t91;
t167 = rSges(4,2) * t94;
t159 = rSges(4,3) * t140 + t141 * t167;
t164 = t92 * t94;
t176 = -rSges(4,2) * t164 + t91 * rSges(4,3);
t177 = t91 * t105;
t132 = rSges(4,1) * t95 - t167;
t165 = t92 * rSges(4,3);
t49 = t132 * t91 - t165;
t163 = t92 * t95;
t51 = rSges(4,1) * t163 + t176;
t2 = (qJD(2) * t49 - t105 * t92 + t159) * t92 + (-t177 + (-t51 + t176) * qJD(2)) * t91;
t186 = t173 * t2;
t184 = t203 * t91 - t198;
t183 = -t204 * qJD(2) + t197 * t92;
t182 = -t197 * t91 - t192;
t172 = 2 * m(5);
t169 = m(4) * t75;
t85 = t91 * rSges(5,2);
t166 = t92 * rSges(5,2);
t162 = t193 * t91 - t166;
t161 = rSges(5,3) * t164 + t92 * t145 + t168 * t163 + t85;
t160 = -t193 * qJD(3) + qJD(4) * t95;
t157 = t92 * pkin(2) + t91 * pkin(5);
t156 = t89 + t90;
t139 = qJD(3) * t94;
t138 = qJD(3) * t95;
t137 = qJD(4) * t94;
t135 = m(5) * t139;
t134 = t92 * t138;
t35 = t158 * t92;
t133 = t168 * qJD(3);
t111 = -pkin(2) - t132;
t98 = -t144 * t94 - pkin(2) - t206;
t97 = rSges(5,2) * t140 - t133 * t164 + t144 * t134 + t92 * t137;
t96 = t98 * t91;
t87 = t92 * pkin(5);
t83 = pkin(5) * t140;
t62 = t132 * qJD(3);
t34 = t158 * t91;
t33 = t51 + t157;
t32 = t111 * t91 + t165 + t87;
t19 = t157 + t161;
t18 = t87 + t96 + t166;
t17 = t177 + ((-rSges(4,3) - pkin(5)) * t91 + t111 * t92) * qJD(2);
t16 = -rSges(4,2) * t134 - pkin(2) * t141 + t83 + (-t139 * t92 - t141 * t95) * rSges(4,1) + t159;
t15 = -qJD(2) * t35 + t160 * t91;
t14 = t141 * t158 + t160 * t92;
t5 = t161 * t92 + t162 * t91;
t4 = (t158 * qJD(3) - t137) * t91 + ((-rSges(5,2) - pkin(5)) * t91 + t98 * t92) * qJD(2);
t3 = qJD(2) * t96 + t83 + t97;
t1 = (qJD(2) * t162 + t97) * t92 + ((-t161 + t85) * qJD(2) + (t144 * t138 + (qJD(4) - t133) * t94) * t91) * t91;
t6 = [0; 0; (t16 * t33 + t17 * t32) * t173 + (t18 * t4 + t19 * t3) * t172 + (t119 + t121 + t151 - t155 + (-Icges(4,2) - Icges(5,3)) * t95) * t139 + (-t113 + t117 - t150 + t154 + (Icges(4,1) + Icges(5,1)) * t94) * t138; m(4) * t2 + m(5) * t1; m(4) * ((-t16 * t91 - t17 * t92) * t75 + (-t32 * t92 - t33 * t91) * t62) + m(5) * (t14 * t18 + t15 * t19 - t3 * t34 - t35 * t4) + (-t106 / 0.2e1 + t108 / 0.2e1 + t107 / 0.2e1 - t109 / 0.2e1 + t205 * (t89 / 0.2e1 + t90 / 0.2e1)) * qJD(3) + ((-t33 * t169 + (t43 / 0.2e1 - t37 / 0.2e1) * t95 + (t47 / 0.2e1 + t45 / 0.2e1) * t94) * t92 + (t32 * t169 + (t42 / 0.2e1 - t36 / 0.2e1) * t95 + (t46 / 0.2e1 + t44 / 0.2e1) * t94) * t91 + (-t200 * t94 - t202 * t95) * t91 / 0.2e1 - (t199 * t94 + t201 * t95) * t92 / 0.2e1) * qJD(2); t156 * t75 * t62 * t173 + (t1 * t5 - t14 * t35 - t15 * t34) * t172 + (t188 * t141 + t182 * t90 + t51 * t186 + (-t196 * t92 + t191) * t140) * t92 + (t49 * t186 + t183 * t89 + t184 * t140 + ((t182 - t192) * t91 + t183 * t92 + t190 * t139 + t189 * t138 + (-t189 * t95 - t190 * t94) * qJD(3) + ((-t196 + t203) * t91 + t198 + t184 + t188) * qJD(2)) * t92 + (t175 * t91 - t191) * t141) * t91; t135; m(5) * ((t18 * t92 + t19 * t91) * t138 + (t3 * t91 + t4 * t92 + (-t18 * t91 + t19 * t92) * qJD(2)) * t94); m(5) * ((-t1 + (-t34 * t91 - t35 * t92) * qJD(3)) * t95 + (qJD(3) * t5 + t14 * t92 + t15 * t91 + (-t34 * t92 + t35 * t91) * qJD(2)) * t94); 0.2e1 * (-0.1e1 + t156) * t95 * t135;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
