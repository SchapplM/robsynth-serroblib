% Calculate time derivative of joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:10
% DurationCPUTime: 3.93s
% Computational Cost: add. (2338->229), mult. (5544->317), div. (0->0), fcn. (5568->6), ass. (0->131)
t111 = cos(qJ(4));
t110 = sin(qJ(4));
t239 = -Icges(5,4) - Icges(6,4);
t249 = t239 * t110;
t251 = Icges(5,2) + Icges(6,2);
t254 = -t251 * t111 + t249;
t250 = t239 * t111;
t252 = Icges(5,1) + Icges(6,1);
t253 = -t252 * t110 + t250;
t238 = Icges(5,5) + Icges(6,5);
t237 = Icges(5,6) + Icges(6,6);
t248 = t251 * t110 + t250;
t247 = t252 * t111 + t249;
t170 = sin(pkin(7));
t171 = cos(pkin(7));
t192 = sin(qJ(1));
t193 = cos(qJ(1));
t77 = -t170 * t192 - t171 * t193;
t78 = t193 * t170 - t192 * t171;
t219 = -t237 * t77 + t248 * t78;
t218 = t237 * t78 + t248 * t77;
t217 = t238 * t77 + t247 * t78;
t216 = -t238 * t78 + t247 * t77;
t236 = Icges(5,3) + Icges(6,3);
t246 = t237 * t110 - t238 * t111;
t230 = -qJ(5) - pkin(6) - rSges(6,3);
t243 = t230 * t77;
t221 = -t236 * t77 + t246 * t78;
t231 = t236 * t78 + t246 * t77;
t102 = pkin(4) * t111 + pkin(3);
t183 = rSges(6,2) * t110;
t139 = -rSges(6,1) * t111 + t183;
t127 = t102 - t139;
t240 = (-pkin(3) + t127) * t78;
t210 = t218 * t110 + t216 * t111;
t209 = t219 * t110 + t217 * t111;
t233 = -t110 * t238 - t111 * t237;
t163 = t193 * pkin(1) + t192 * qJ(2);
t153 = t193 * pkin(2) + t163;
t162 = qJD(4) * t110;
t73 = t78 * qJD(1);
t176 = t111 * t73;
t123 = t77 * t162 + t176;
t229 = -rSges(3,1) * t193 - rSges(3,3) * t192 - t163;
t226 = -0.2e1 * t78;
t182 = rSges(6,2) * t111;
t207 = t110 * (rSges(6,1) + pkin(4)) + t182;
t117 = qJD(4) * t207;
t72 = t77 * qJD(1);
t113 = rSges(6,1) * t176 + t73 * t102 + t117 * t77 - t230 * t72;
t196 = pkin(6) * t78;
t175 = t111 * t77;
t180 = t110 * t77;
t212 = -rSges(6,1) * t175 + rSges(6,2) * t180 - t77 * t102 - t230 * t78;
t187 = pkin(3) * t77 - t196 + t212;
t76 = t77 * pkin(6);
t188 = -t240 + t76 + t243;
t195 = t72 * pkin(6);
t214 = t230 * t73;
t70 = t73 * pkin(6);
t1 = -t187 * t73 + (t117 * t78 - t214 - t70) * t78 - (-t188 + t240) * t72 + (-t195 + (-pkin(3) - t183) * t73 + t113) * t77;
t203 = 2 * m(6);
t225 = t1 * t203;
t122 = t123 * rSges(5,1) + t72 * rSges(5,3);
t161 = qJD(4) * t111;
t181 = t110 * t73;
t124 = t161 * t77 - t181;
t125 = -t111 * t72 + t78 * t162;
t126 = t110 * t72 + t161 * t78;
t191 = t73 * rSges(5,3);
t184 = rSges(5,2) * t110;
t140 = -rSges(5,1) * t111 + t184;
t190 = t77 * rSges(5,3);
t50 = t140 * t78 - t190;
t185 = -rSges(5,1) * t175 + t78 * rSges(5,3);
t52 = rSges(5,2) * t180 + t185;
t2 = t72 * t50 + t78 * (t125 * rSges(5,1) + t191) - t73 * t52 + t77 * t122 + (t124 * t77 + t126 * t78) * rSges(5,2);
t204 = 2 * m(5);
t224 = t2 * t204;
t223 = -t125 * t238 - t126 * t237 - t236 * t73;
t222 = t123 * t238 + t124 * t237 + t236 * t72;
t215 = t246 * qJD(4);
t96 = -t110 * rSges(5,1) - rSges(5,2) * t111;
t213 = t96 ^ 2 * t204;
t87 = t140 * qJD(4);
t211 = t87 * t96 * t204;
t208 = (t210 + t221) * t78 + (t209 - t231) * t77;
t206 = t161 * t254 + t162 * t253;
t205 = -t110 * t254 + t111 * t253;
t202 = t72 / 0.2e1;
t201 = t73 / 0.2e1;
t200 = -t77 / 0.2e1;
t199 = t78 / 0.2e1;
t198 = m(5) * t87;
t197 = m(5) * t96;
t165 = qJD(4) * t78;
t106 = t193 * qJ(2);
t164 = qJD(1) * t106 + qJD(2) * t192;
t160 = t192 * pkin(1);
t95 = -t110 * rSges(6,1) - t182;
t152 = pkin(4) * t110 - t95;
t145 = t72 * t78 - t73 * t77;
t128 = pkin(3) - t140;
t121 = -pkin(2) * t192 - t160;
t118 = -t192 * t77 + t193 * t78;
t116 = t106 + t121;
t115 = -rSges(3,1) * t192 + rSges(3,3) * t193 - t160;
t114 = qJD(1) * t121 + t164;
t104 = qJD(2) * t193;
t112 = -t153 * qJD(1) + t104;
t86 = t139 * qJD(4);
t64 = t106 + t115;
t58 = t229 * qJD(1) + t104;
t57 = qJD(1) * t115 + t164;
t56 = t152 * t77;
t55 = t152 * t78;
t54 = -rSges(4,1) * t77 - rSges(4,2) * t78 + t153;
t53 = t78 * rSges(4,1) - t77 * rSges(4,2) + t116;
t36 = t72 * rSges(4,1) + t73 * rSges(4,2) + t112;
t35 = t73 * rSges(4,1) - t72 * rSges(4,2) + t114;
t32 = pkin(4) * t124 + t73 * t95 - t77 * t86;
t31 = pkin(4) * t126 - t72 * t95 - t78 * t86;
t30 = t196 + (-pkin(3) + t184) * t77 + t153 + t185;
t29 = t128 * t78 + t116 + t190 + t76;
t16 = t153 + t212;
t15 = t127 * t78 + t116 - t243;
t14 = t128 * t72 + t165 * t96 + t112 - t191 - t70;
t13 = rSges(5,2) * t124 + t73 * pkin(3) + t114 + t122 + t195;
t12 = t77 * qJD(5) + t127 * t72 - t207 * t165 + t112 + t214;
t11 = -rSges(6,2) * t181 + t78 * qJD(5) + t113 + t114;
t3 = [0.2e1 * m(3) * (-t229 * t57 + t58 * t64) + 0.2e1 * m(4) * (t35 * t54 + t36 * t53) + (t13 * t30 + t14 * t29) * t204 + (t11 * t16 + t12 * t15) * t203 + (t254 + t247) * t162 + (-t253 - t248) * t161; m(3) * (t192 * t58 - t193 * t57 + (-t192 * t229 + t193 * t64) * qJD(1)) + m(4) * (t192 * t36 - t193 * t35 + (t192 * t54 + t193 * t53) * qJD(1)) + m(5) * (t192 * t14 - t193 * t13 + (t192 * t30 + t193 * t29) * qJD(1)) + m(6) * (t192 * t12 - t193 * t11 + (t15 * t193 + t16 * t192) * qJD(1)); 0; 0; 0; 0; m(6) * (t11 * t55 + t12 * t56 + t15 * t32 + t16 * t31) + (t29 * t73 - t30 * t72) * t197 + (t216 * t110 - t218 * t111) * t202 + (t217 * t110 - t219 * t111) * t201 + t233 * t145 + (t209 * qJD(4) + t125 * t253 + t126 * t254 - t205 * t72 + t233 * t73) * t200 + (t210 * qJD(4) + t123 * t253 + t124 * t254 + t205 * t73 + t233 * t72) * t199 + (-t13 * t197 - t30 * t198 + t215 * t199 + t206 * t200 - t205 * t201) * t78 + (-t14 * t197 - t29 * t198 + t206 * t199 - t215 * t200 - t205 * t202) * t77; m(6) * (t32 * t192 - t31 * t193 + (t192 * t55 + t193 * t56) * qJD(1)) + t118 * t198 + (t192 * t73 + t193 * t72 + (-t192 * t78 - t193 * t77) * qJD(1)) * t197; -m(5) * t2 - m(6) * t1; (t55 * t31 + t56 * t32) * t203 + (t50 * t224 + t188 * t225 + (t222 * t78 + t211) * t78 - (-0.2e1 * t210 * t77 - t213 + (t226 - t78) * t231) * t72 + (-t210 * t78 + t208) * t73) * t78 + (t52 * t224 + t187 * t225 + (t223 * t77 + t211) * t77 + (t209 * t226 + 0.3e1 * t221 * t77 - t213) * t73 - (-t209 * t77 + t208) * t72 + (t222 * t77 + t223 * t78 - t231 * t73 - t221 * t72 + (-t216 * t72 + t217 * t73) * t111 + (-t218 * t72 + t219 * t73) * t110) * t78) * t77; m(6) * (-t11 * t77 + t12 * t78 + t15 * t72 + t16 * t73); m(6) * (qJD(1) * t118 + t192 * t72 - t193 * t73); 0; m(6) * (-t31 * t77 + t32 * t78 + t55 * t73 + t56 * t72); t145 * t203;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
