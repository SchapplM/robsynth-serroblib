% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:47
% EndTime: 2019-12-05 17:28:53
% DurationCPUTime: 2.85s
% Computational Cost: add. (9030->199), mult. (8175->319), div. (0->0), fcn. (8642->10), ass. (0->116)
t145 = pkin(9) + qJ(5);
t141 = sin(t145);
t143 = cos(t145);
t147 = sin(pkin(8));
t148 = cos(pkin(8));
t187 = Icges(6,4) * t141;
t110 = -Icges(6,5) * t148 + (Icges(6,1) * t143 - t187) * t147;
t124 = (-Icges(6,2) * t143 - t187) * t147;
t175 = t110 + t124;
t186 = Icges(6,4) * t143;
t109 = -Icges(6,6) * t148 + (-Icges(6,2) * t141 + t186) * t147;
t125 = (-Icges(6,1) * t141 - t186) * t147;
t176 = t109 - t125;
t123 = (-Icges(6,5) * t141 - Icges(6,6) * t143) * t147;
t181 = t148 * t123;
t235 = (-t181 + (-t141 * t175 - t143 * t176) * t147) * t148;
t146 = qJ(1) + pkin(7);
t142 = sin(t146);
t144 = cos(t146);
t189 = sin(pkin(9));
t163 = t144 * t189;
t164 = sin(qJ(1)) * pkin(1) - t144 * qJ(3);
t184 = t142 * t148;
t111 = t141 * t184 + t143 * t144;
t112 = -t144 * t141 + t143 * t184;
t178 = t112 * rSges(6,1) - t111 * rSges(6,2);
t190 = cos(pkin(9));
t229 = -(pkin(4) * t190 + pkin(3)) * t148 - pkin(2) + (-pkin(6) - qJ(4) - rSges(6,3)) * t147;
t54 = -pkin(4) * t163 - t229 * t142 + t164 + t178;
t203 = cos(qJ(1)) * pkin(1);
t182 = t144 * t148;
t113 = t141 * t182 - t142 * t143;
t114 = t141 * t142 + t143 * t182;
t226 = -t114 * rSges(6,1) + t113 * rSges(6,2);
t55 = (-pkin(4) * t189 - qJ(3)) * t142 - t203 + t229 * t144 + t226;
t94 = rSges(6,1) * t111 + rSges(6,2) * t112;
t95 = -rSges(6,1) * t113 - rSges(6,2) * t114;
t16 = (-(t110 / 0.2e1 + t124 / 0.2e1) * t141 + (t125 / 0.2e1 - t109 / 0.2e1) * t143) * t147 + m(6) * (-t54 * t94 - t55 * t95) - t181 / 0.2e1;
t234 = t16 * qJD(1);
t183 = t144 * t147;
t77 = Icges(6,5) * t114 - Icges(6,6) * t113 + Icges(6,3) * t183;
t232 = t183 * t77;
t153 = t147 * (-rSges(6,3) * t148 + (rSges(6,1) * t143 - rSges(6,2) * t141) * t147);
t86 = rSges(6,3) * t183 - t226;
t58 = t144 * t153 + t148 * t86;
t103 = Icges(6,4) * t114;
t80 = -Icges(6,2) * t113 + Icges(6,6) * t183 + t103;
t102 = Icges(6,4) * t113;
t84 = -Icges(6,1) * t114 - Icges(6,5) * t183 + t102;
t199 = t111 * t80 + t112 * t84;
t228 = -t147 / 0.2e1;
t204 = -t148 / 0.2e1;
t227 = m(6) * t147;
t48 = (t142 * t95 + t144 * t94) * t147;
t223 = pkin(3) * t148 + pkin(2) + (rSges(5,3) + qJ(4)) * t147;
t221 = t142 ^ 2;
t129 = (t144 ^ 2 + t221) * t147;
t220 = 0.2e1 * t129;
t219 = 0.4e1 * qJD(1);
t185 = t142 * t147;
t216 = m(4) * (-t142 * (-t203 + (-rSges(4,3) - qJ(3)) * t142 + rSges(4,2) * t183) - t144 * (-rSges(4,2) * t185 - t144 * rSges(4,3) + t164));
t161 = t148 * t189;
t162 = t148 * t190;
t151 = -(t142 * t189 + t144 * t162) * rSges(5,1) + (-t142 * t190 + t144 * t161) * rSges(5,2) - t203 - t142 * qJ(3) - t223 * t144;
t63 = -(t142 * t161 + t144 * t190) * rSges(5,2) - (-t142 * t162 + t163) * rSges(5,1) + t223 * t142 + t164;
t215 = m(5) * (t142 * t63 - t144 * t151) * t147;
t214 = m(5) * (-t142 * t151 - t144 * t63);
t211 = (t142 * t54 - t144 * t55) * t227;
t210 = m(6) * (-t142 * t55 - t144 * t54);
t85 = -rSges(6,3) * t185 - t178;
t57 = -t142 * t153 + t148 * t85;
t209 = (t142 * t57 - t144 * t58) * t227;
t208 = m(6) * (-t142 * t58 - t57 * t144);
t207 = m(6) * t48;
t206 = m(6) * (t142 * t94 - t144 * t95);
t188 = Icges(6,4) * t112;
t79 = Icges(6,2) * t111 - Icges(6,6) * t185 - t188;
t101 = Icges(6,4) * t111;
t82 = -Icges(6,1) * t112 - Icges(6,5) * t185 + t101;
t198 = -t111 * t79 + t112 * t82;
t197 = -Icges(6,1) * t111 - t188 + t79;
t196 = Icges(6,1) * t113 + t103 + t80;
t195 = Icges(6,2) * t112 + t101 + t82;
t194 = -Icges(6,2) * t114 - t102 - t84;
t193 = m(6) * qJD(5);
t71 = (m(5) / 0.2e1 + m(6) / 0.2e1) * t220;
t179 = t71 * qJD(1);
t76 = -Icges(6,5) * t112 + Icges(6,6) * t111 - Icges(6,3) * t185;
t157 = t185 * t76 + t198;
t25 = -t185 * t77 + t199;
t172 = (t142 * t157 + t144 * t25) * t228 + (t199 * t144 + (t157 - t232) * t142) * t147 / 0.2e1;
t171 = (-t77 * t221 * t147 + (t199 - t25) * t142 + (-t198 + t157 + (-t142 * t76 - t144 * t77) * t147 + t232) * t144) * t228 + (t204 + t148 / 0.2e1) * ((-Icges(6,3) * t148 + (Icges(6,5) * t143 - Icges(6,6) * t141) * t147) * t183 - t109 * t113 + t110 * t114);
t126 = (-rSges(6,1) * t141 - rSges(6,2) * t143) * t147;
t89 = -Icges(6,5) * t113 - Icges(6,6) * t114;
t88 = Icges(6,5) * t111 + Icges(6,6) * t112;
t74 = t126 * t183 + t148 * t95;
t73 = t126 * t185 - t148 * t94;
t72 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t220 - (m(5) + m(6)) * t129 / 0.2e1;
t50 = t206 / 0.2e1;
t47 = t207 / 0.2e1;
t34 = t208 / 0.2e1;
t32 = -t113 * t175 - t114 * t176 + t123 * t183;
t31 = t111 * t175 + t112 * t176 - t123 * t185;
t30 = t209 / 0.2e1;
t23 = -t148 * t89 + (-t194 * t141 - t196 * t143) * t147;
t22 = -t148 * t88 + (-t195 * t141 - t197 * t143) * t147;
t17 = t211 + t215;
t15 = t210 + t214 + t216;
t9 = t34 - t206 / 0.2e1;
t8 = t50 + t34;
t7 = t50 - t208 / 0.2e1;
t6 = t30 - t207 / 0.2e1;
t5 = t47 + t30;
t4 = t47 - t209 / 0.2e1;
t1 = (t142 * t171 + t144 * t172) * t147;
t2 = [t15 * qJD(3) + t17 * qJD(4) + t16 * qJD(5), 0, qJD(1) * t15 + qJD(4) * t72 + qJD(5) * t8, qJD(1) * t17 + qJD(3) * t72 + qJD(5) * t5, t234 + t8 * qJD(3) + t5 * qJD(4) + (-t235 + (-t73 * t54 + t74 * t55 - t57 * t94 - t58 * t95) * m(6) + ((t23 / 0.2e1 + t32 / 0.2e1 - t172) * t144 + (-t22 / 0.2e1 - t31 / 0.2e1 - t171) * t142) * t147) * qJD(5); 0, 0, 0, 0, -t48 * t193; -t71 * qJD(4) + t7 * qJD(5) + (-t216 / 0.4e1 - t214 / 0.4e1 - t210 / 0.4e1) * t219, 0, 0, -t179, t7 * qJD(1) + (t142 * t73 + t144 * t74) * t193; t71 * qJD(3) + t4 * qJD(5) + (-t215 / 0.4e1 - t211 / 0.4e1) * t219, 0, t179, 0, t4 * qJD(1) + (t148 * t48 + (-t142 * t74 + t144 * t73) * t147) * t193; t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) - t234, 0, t9 * qJD(1), t6 * qJD(1), t1 * qJD(1) + (m(6) * (-(-t142 * t86 - t144 * t85) * t147 * t48 - t57 * t73 + t58 * t74) + (-t235 + (-t142 * t22 + t144 * t23) * t147) * t204 - (-t31 * t148 - (t195 * t111 + t197 * t112 - t88 * t185) * t185 + (t194 * t111 + t196 * t112 - t89 * t185) * t183) * t185 / 0.2e1 + (-t32 * t148 - (-t195 * t113 - t197 * t114 + t88 * t183) * t185 + (-t194 * t113 - t196 * t114 + t89 * t183) * t183) * t183 / 0.2e1) * qJD(5);];
Cq = t2;
