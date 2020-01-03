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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:06
% EndTime: 2020-01-03 11:20:15
% DurationCPUTime: 2.72s
% Computational Cost: add. (9030->213), mult. (8175->336), div. (0->0), fcn. (8642->10), ass. (0->122)
t154 = pkin(9) + qJ(5);
t149 = sin(t154);
t151 = cos(t154);
t157 = sin(pkin(8));
t159 = cos(pkin(8));
t189 = Icges(6,4) * t149;
t112 = -Icges(6,5) * t159 + (Icges(6,1) * t151 - t189) * t157;
t125 = (-Icges(6,2) * t151 - t189) * t157;
t178 = t112 + t125;
t188 = Icges(6,4) * t151;
t111 = -Icges(6,6) * t159 + (-Icges(6,2) * t149 + t188) * t157;
t126 = (-Icges(6,1) * t149 - t188) * t157;
t179 = t111 - t126;
t124 = (-Icges(6,5) * t149 - Icges(6,6) * t151) * t157;
t182 = t159 * t124;
t238 = (-t182 + (-t178 * t149 - t179 * t151) * t157) * t159;
t155 = qJ(1) + pkin(7);
t152 = cos(t155);
t150 = sin(t155);
t185 = t150 * t159;
t114 = -t149 * t185 - t151 * t152;
t115 = -t152 * t149 + t151 * t185;
t101 = rSges(6,1) * t114 - rSges(6,2) * t115;
t183 = t152 * t159;
t116 = t149 * t183 - t150 * t151;
t117 = t149 * t150 + t151 * t183;
t102 = rSges(6,1) * t116 + rSges(6,2) * t117;
t158 = cos(pkin(9));
t148 = pkin(4) * t158 + pkin(3);
t170 = t157 * (-pkin(6) - qJ(4));
t173 = t152 * pkin(2) + t150 * qJ(3) + cos(qJ(1)) * pkin(1);
t156 = sin(pkin(9));
t187 = t150 * t156;
t184 = t152 * t157;
t93 = t117 * rSges(6,1) - t116 * rSges(6,2) + rSges(6,3) * t184;
t225 = pkin(4) * t187 + t148 * t183 - t152 * t170 + t173 + t93;
t168 = rSges(6,1) * t115 + rSges(6,2) * t114;
t169 = sin(qJ(1)) * pkin(1) - t152 * qJ(3);
t56 = -t152 * t156 * pkin(4) + t168 + t169 + (rSges(6,3) * t157 + t148 * t159 + pkin(2) - t170) * t150;
t16 = (-(t112 / 0.2e1 + t125 / 0.2e1) * t149 + (t126 / 0.2e1 - t111 / 0.2e1) * t151) * t157 + m(6) * (t101 * t56 - t102 * t225) - t182 / 0.2e1;
t237 = t16 * qJD(1);
t84 = Icges(6,5) * t117 - Icges(6,6) * t116 + Icges(6,3) * t184;
t235 = t184 * t84;
t113 = -t159 * rSges(6,3) + (rSges(6,1) * t151 - rSges(6,2) * t149) * t157;
t234 = -t113 * t184 - t159 * t93;
t107 = Icges(6,4) * t117;
t88 = Icges(6,2) * t116 - Icges(6,6) * t184 - t107;
t106 = Icges(6,4) * t116;
t90 = Icges(6,1) * t117 + Icges(6,5) * t184 - t106;
t233 = -t114 * t88 + t115 * t90;
t230 = rSges(5,3) + qJ(4);
t186 = t150 * t157;
t83 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t186;
t229 = t83 * t184;
t228 = t157 / 0.2e1;
t203 = -t159 / 0.2e1;
t48 = (t101 * t152 + t102 * t150) * t157;
t223 = t152 ^ 2;
t130 = (t150 ^ 2 + t223) * t157;
t222 = 0.2e1 * t130;
t220 = 0.4e1 * qJD(1);
t219 = m(5) / 0.2e1;
t217 = m(6) / 0.2e1;
t215 = m(4) * (-(-rSges(4,2) * t186 - t152 * rSges(4,3) + t150 * pkin(2) + t169) * t152 + (-rSges(4,2) * t184 + rSges(4,3) * t150 + t173) * t150);
t162 = rSges(5,1) * (t158 * t183 + t187) - rSges(5,2) * (-t150 * t158 + t156 * t183) + pkin(3) * t183 + t173 + t230 * t184;
t69 = (-rSges(5,1) * t156 - rSges(5,2) * t158) * t152 + (pkin(2) + t230 * t157 + (rSges(5,1) * t158 - rSges(5,2) * t156 + pkin(3)) * t159) * t150 + t169;
t38 = t162 * t184 + t69 * t186;
t214 = m(5) * t38;
t213 = m(5) * (t162 * t150 - t152 * t69);
t92 = rSges(6,3) * t186 + t168;
t61 = t113 * t186 + t159 * t92;
t33 = t184 * t234 - t61 * t186;
t212 = m(6) * ((t150 * t61 - t152 * t234) * t157 + t33);
t28 = t184 * t225 + t56 * t186;
t210 = m(6) * t28;
t209 = m(6) * (t225 * t150 - t152 * t56);
t208 = m(6) * t33;
t207 = m(6) * (t234 * t150 + t152 * t61);
t206 = m(6) * t48;
t205 = m(6) * (-t101 * t150 + t102 * t152);
t190 = Icges(6,4) * t115;
t86 = Icges(6,2) * t114 + Icges(6,6) * t186 + t190;
t199 = -Icges(6,1) * t114 + t190 + t86;
t105 = Icges(6,4) * t114;
t89 = Icges(6,1) * t115 + Icges(6,5) * t186 + t105;
t198 = -Icges(6,2) * t115 + t105 + t89;
t197 = Icges(6,2) * t117 + t106 - t90;
t196 = m(6) * qJD(5);
t194 = t150 * t84;
t193 = t152 * t83;
t191 = Icges(6,1) * t116 + t107 - t88;
t77 = (t219 + t217) * t222;
t180 = t77 * qJD(1);
t24 = t114 * t86 + t115 * t89 + t83 * t186;
t25 = -t186 * t84 - t233;
t176 = (t150 * t24 - t152 * t25) * t228 - ((t24 + t235) * t150 + ((t193 - t194) * t157 - t25 - t229) * t152) * t157 / 0.2e1;
t175 = (t223 * t84 * t157 - t235 * t152 + (t25 + (t193 + t194) * t157 - t229 + t233) * t150) * t228 + (t203 + t159 / 0.2e1) * ((-Icges(6,3) * t159 + (Icges(6,5) * t151 - Icges(6,6) * t149) * t157) * t184 - t111 * t116 + t112 * t117);
t127 = (-rSges(6,1) * t149 - rSges(6,2) * t151) * t157;
t96 = Icges(6,5) * t116 + Icges(6,6) * t117;
t95 = Icges(6,5) * t114 - Icges(6,6) * t115;
t80 = t102 * t159 - t127 * t184;
t79 = -t101 * t159 - t127 * t186;
t78 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t222 - (m(5) + m(6)) * t130 / 0.2e1;
t52 = t205 / 0.2e1;
t47 = -t206 / 0.2e1;
t34 = t207 / 0.2e1;
t32 = t178 * t116 + t179 * t117 - t124 * t184;
t31 = t178 * t114 - t179 * t115 + t124 * t186;
t30 = t208 / 0.2e1;
t23 = -t159 * t96 + (-t197 * t149 + t191 * t151) * t157;
t22 = -t159 * t95 + (-t198 * t149 - t199 * t151) * t157;
t17 = t210 + t214;
t15 = t209 + t213 + t215;
t13 = t212 / 0.2e1;
t9 = t34 - t205 / 0.2e1;
t8 = t52 + t34;
t7 = t52 - t207 / 0.2e1;
t6 = t30 + t13 + t206 / 0.2e1;
t5 = t47 + t30 - t212 / 0.2e1;
t4 = t47 + t13 - t208 / 0.2e1;
t1 = (t175 * t150 + t176 * t152) * t157;
t2 = [t15 * qJD(3) + t17 * qJD(4) + t16 * qJD(5), 0, qJD(1) * t15 + qJD(4) * t78 + qJD(5) * t8, qJD(1) * t17 + qJD(3) * t78 + qJD(5) * t5, t237 + t8 * qJD(3) + t5 * qJD(4) + (-t238 + (-t61 * t101 - t102 * t234 + t225 * t80 + t79 * t56) * m(6) + ((-t23 / 0.2e1 - t32 / 0.2e1 - t176) * t152 + (t22 / 0.2e1 + t31 / 0.2e1 - t175) * t150) * t157) * qJD(5); 0, 0, 0, 0, t48 * t196; -t77 * qJD(4) + t7 * qJD(5) + (-t215 / 0.4e1 - t213 / 0.4e1 - t209 / 0.4e1) * t220, 0, 0, -t180, t7 * qJD(1) + (-t150 * t79 - t152 * t80) * t196; t77 * qJD(3) + t4 * qJD(5) + (-t214 / 0.4e1 - t210 / 0.4e1) * t220 + 0.2e1 * (((-t150 * t69 - t152 * t162) * t157 + t38) * t219 + ((-t150 * t56 - t152 * t225) * t157 + t28) * t217) * qJD(1), 0, t180, 0, t4 * qJD(1) + (-t159 * t48 + (t150 * t80 - t152 * t79) * t157) * t196; t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) - t237, 0, t9 * qJD(1), t6 * qJD(1), t1 * qJD(1) + (m(6) * (t234 * t80 + (-t150 * t93 + t152 * t92) * t157 * t48 - t61 * t79) + (-t238 + (t150 * t22 - t152 * t23) * t157) * t203 + (-t31 * t159 + (t198 * t114 - t199 * t115 + t95 * t186) * t186 - (t197 * t114 + t191 * t115 + t96 * t186) * t184) * t186 / 0.2e1 - (-t32 * t159 + (t198 * t116 + t199 * t117 - t95 * t184) * t186 - (t197 * t116 - t191 * t117 - t96 * t184) * t184) * t184 / 0.2e1) * qJD(5);];
Cq = t2;
