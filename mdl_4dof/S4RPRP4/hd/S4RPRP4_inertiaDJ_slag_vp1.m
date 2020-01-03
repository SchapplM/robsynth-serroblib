% Calculate time derivative of joint inertia matrix for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:42
% DurationCPUTime: 2.28s
% Computational Cost: add. (2154->183), mult. (2772->255), div. (0->0), fcn. (2105->6), ass. (0->116)
t217 = Icges(5,4) + Icges(4,5);
t216 = -Icges(4,6) + Icges(5,6);
t215 = Icges(5,2) + Icges(4,3);
t95 = sin(qJ(3));
t97 = cos(qJ(3));
t213 = t216 * t95 + t217 * t97;
t174 = rSges(5,1) + pkin(3);
t214 = t174 * t97;
t94 = qJ(1) + pkin(6);
t91 = sin(t94);
t92 = cos(t94);
t211 = t213 * t92 + t215 * t91;
t212 = t213 * t91 - t215 * t92;
t157 = Icges(5,5) * t97;
t118 = Icges(5,3) * t95 + t157;
t36 = -Icges(5,6) * t92 + t118 * t91;
t161 = Icges(4,4) * t97;
t122 = -Icges(4,2) * t95 + t161;
t42 = -Icges(4,6) * t92 + t122 * t91;
t210 = t42 - t36;
t37 = Icges(5,6) * t91 + t118 * t92;
t43 = Icges(4,6) * t91 + t122 * t92;
t209 = t43 - t37;
t158 = Icges(5,5) * t95;
t124 = Icges(5,1) * t97 + t158;
t44 = -Icges(5,4) * t92 + t124 * t91;
t162 = Icges(4,4) * t95;
t126 = Icges(4,1) * t97 - t162;
t46 = -Icges(4,5) * t92 + t126 * t91;
t208 = t46 + t44;
t45 = Icges(5,4) * t91 + t124 * t92;
t47 = Icges(4,5) * t91 + t126 * t92;
t207 = t47 + t45;
t127 = t43 * t95 - t47 * t97;
t129 = t37 * t95 + t45 * t97;
t182 = t127 - t129;
t206 = t182 * t92;
t145 = qJD(3) * t97;
t151 = rSges(5,3) + qJ(4);
t205 = t151 * t145;
t204 = (t216 * t97 - t217 * t95) * qJD(3);
t128 = t42 * t95 - t46 * t97;
t130 = t36 * t95 + t44 * t97;
t203 = t128 - t130;
t164 = -t151 * t97 + t174 * t95;
t152 = qJ(4) * t95;
t200 = rSges(5,3) * t95 + t152 + t214;
t199 = t211 * qJD(1);
t111 = t127 * t91;
t112 = t128 * t92;
t113 = t129 * t91;
t114 = t130 * t92;
t198 = t211 * t92 - t212 * t91 + t111 + t112 - t113 - t114;
t197 = t207 * t91 + t208 * t92;
t196 = t209 * t91 + t210 * t92;
t89 = t91 ^ 2;
t90 = t92 ^ 2;
t195 = t203 * t91 + t212 * t92;
t180 = 2 * m(4);
t75 = rSges(4,1) * t95 + rSges(4,2) * t97;
t110 = qJD(3) * t75;
t147 = qJD(1) * t92;
t148 = qJD(1) * t91;
t172 = rSges(4,2) * t95;
t100 = rSges(4,3) * t147 - t110 * t92 + t148 * t172;
t169 = t92 * t95;
t183 = -rSges(4,2) * t169 + t91 * rSges(4,3);
t184 = t91 * t110;
t173 = rSges(4,1) * t97;
t137 = -t172 + t173;
t170 = rSges(4,3) * t92;
t49 = t137 * t91 - t170;
t168 = t92 * t97;
t51 = rSges(4,1) * t168 + t183;
t2 = (qJD(1) * t49 + t100) * t92 + (-t184 + (-t51 + t183) * qJD(1)) * t91;
t193 = t180 * t2;
t191 = t211 * t91 - t206;
t190 = -t212 * qJD(1) + t204 * t92;
t189 = -t204 * t91 - t199;
t179 = 2 * m(5);
t176 = m(4) * t75;
t175 = sin(qJ(1)) * pkin(1);
t93 = cos(qJ(1)) * pkin(1);
t171 = rSges(5,2) * t92;
t85 = t91 * rSges(5,2);
t167 = t200 * t91 - t171;
t166 = rSges(5,3) * t169 + t92 * t152 + t174 * t168 + t85;
t165 = -t200 * qJD(3) + qJD(4) * t97;
t163 = t89 + t90;
t146 = qJD(3) * t95;
t144 = qJD(4) * t95;
t142 = t92 * pkin(2) + t91 * pkin(5) + t93;
t141 = m(5) * t146;
t35 = t164 * t92;
t139 = t174 * qJD(3);
t116 = -pkin(2) - t137;
t102 = -t151 * t95 - pkin(2) - t214;
t101 = rSges(5,2) * t147 - t139 * t169 + (t144 + t205) * t92;
t99 = t102 * t91 - t175;
t87 = t92 * pkin(5);
t83 = pkin(5) * t147;
t62 = t137 * qJD(3);
t34 = t164 * t91;
t33 = t51 + t142;
t32 = t116 * t91 + t170 - t175 + t87;
t19 = t142 + t166;
t18 = t87 + t99 + t171;
t17 = t184 + (-t93 + (-rSges(4,3) - pkin(5)) * t91 + t116 * t92) * qJD(1);
t16 = t83 + (-t175 + (-pkin(2) - t173) * t91) * qJD(1) + t100;
t15 = -qJD(1) * t35 + t165 * t91;
t14 = t148 * t164 + t165 * t92;
t5 = t166 * t92 + t167 * t91;
t4 = (t164 * qJD(3) - t144) * t91 + (-t93 + (-rSges(5,2) - pkin(5)) * t91 + t102 * t92) * qJD(1);
t3 = qJD(1) * t99 + t101 + t83;
t1 = (qJD(1) * t167 + t101) * t92 + ((-t166 + t85) * qJD(1) + (t205 + (qJD(4) - t139) * t95) * t91) * t91;
t6 = [(t16 * t33 + t17 * t32) * t180 + (t18 * t4 + t19 * t3) * t179 + (t124 + t126 + t158 - t162 + (-Icges(4,2) - Icges(5,3)) * t97) * t146 + (-t118 + t122 - t157 + t161 + (Icges(4,1) + Icges(5,1)) * t95) * t145; 0; 0; m(4) * ((-t16 * t91 - t17 * t92) * t75 + (-t32 * t92 - t33 * t91) * t62) + m(5) * (t14 * t18 + t15 * t19 - t3 * t34 - t35 * t4) + (-t111 / 0.2e1 + t113 / 0.2e1 + t112 / 0.2e1 - t114 / 0.2e1 + t213 * (t89 / 0.2e1 + t90 / 0.2e1)) * qJD(3) + ((-t33 * t176 + (t43 / 0.2e1 - t37 / 0.2e1) * t97 + (t47 / 0.2e1 + t45 / 0.2e1) * t95) * t92 + (t32 * t176 + (t42 / 0.2e1 - t36 / 0.2e1) * t97 + (t46 / 0.2e1 + t44 / 0.2e1) * t95) * t91 + (-t208 * t95 - t210 * t97) * t91 / 0.2e1 - (t207 * t95 + t209 * t97) * t92 / 0.2e1) * qJD(1); m(4) * t2 + m(5) * t1; t163 * t75 * t62 * t180 + (t1 * t5 - t14 * t35 - t15 * t34) * t179 + (t195 * t148 + t189 * t90 + t51 * t193 + (-t203 * t92 + t198) * t147) * t92 + (t49 * t193 + t190 * t89 + t191 * t147 + ((t189 - t199) * t91 + t190 * t92 + t197 * t146 + t196 * t145 + (-t196 * t97 - t197 * t95) * qJD(3) + ((-t203 + t211) * t91 + t206 + t191 + t195) * qJD(1)) * t92 + (t182 * t91 - t198) * t148) * t91; m(5) * ((t18 * t92 + t19 * t91) * t145 + (t3 * t91 + t4 * t92 + (-t18 * t91 + t19 * t92) * qJD(1)) * t95); t141; m(5) * ((-t1 + (-t34 * t91 - t35 * t92) * qJD(3)) * t97 + (qJD(3) * t5 + t14 * t92 + t15 * t91 + (-t34 * t92 + t35 * t91) * qJD(1)) * t95); 0.2e1 * (-0.1e1 + t163) * t97 * t141;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
