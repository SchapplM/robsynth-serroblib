% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP3
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:47
% DurationCPUTime: 2.11s
% Computational Cost: add. (4200->141), mult. (4307->204), div. (0->0), fcn. (3910->4), ass. (0->97)
t139 = pkin(6) + qJ(2);
t135 = sin(t139);
t136 = cos(t139);
t140 = sin(qJ(3));
t141 = cos(qJ(3));
t222 = rSges(5,2) * t141 + (rSges(5,1) + pkin(3)) * t140;
t226 = t222 * t136;
t227 = t222 * t135;
t236 = m(5) * (t135 * t227 + t136 * t226);
t35 = -t236 / 0.2e1;
t30 = t236 / 0.2e1;
t238 = -Icges(4,4) - Icges(5,4);
t237 = Icges(4,5) + Icges(5,5);
t234 = Icges(4,2) + Icges(5,2);
t232 = Icges(4,6) + Icges(5,6);
t233 = Icges(4,3) + Icges(5,3);
t175 = Icges(5,4) * t140;
t104 = Icges(5,1) * t141 - t175;
t176 = Icges(4,4) * t140;
t106 = Icges(4,1) * t141 - t176;
t228 = (t104 + t106) * t136 + t237 * t135;
t170 = t135 * t140;
t231 = t238 * t170;
t214 = t238 * t141;
t97 = Icges(5,5) * t141 - Icges(5,6) * t140;
t98 = Icges(4,5) * t141 - Icges(4,6) * t140;
t229 = (t97 + t98) * t136 + t233 * t135;
t169 = t135 * t141;
t221 = t136 * t233 - t169 * t237 + t170 * t232;
t224 = t232 * t136 + t169 * t238 + t234 * t170;
t220 = t231 + (Icges(4,1) + Icges(5,1)) * t169 - t237 * t136;
t225 = t228 * t169;
t173 = Icges(5,2) * t140;
t174 = Icges(4,2) * t140;
t223 = (-t173 - t174 - t214) * t136 + t232 * t135;
t219 = t136 * t229 - t225;
t218 = t224 * t140;
t167 = t136 * t141;
t211 = t135 * t229 + t167 * t228;
t217 = t135 * t221 - t167 * t220;
t193 = pkin(3) * t141;
t132 = pkin(2) + t193;
t192 = -qJ(4) - pkin(5);
t44 = -rSges(5,1) * t169 + rSges(5,2) * t170 - t135 * t132 + (rSges(5,3) - t192) * t136;
t184 = rSges(5,1) * t141;
t152 = t132 + t184;
t168 = t136 * t140;
t165 = -rSges(5,2) * t168 - t135 * t192;
t45 = t135 * rSges(5,3) + t152 * t136 + t165;
t216 = m(5) * (t135 * t45 + t136 * t44);
t215 = -t168 * t223 + t211;
t213 = -t140 * t237 - t232 * t141;
t212 = t140 * t223 + t221;
t210 = t168 * t224 - t170 * t223 - t217 - t219;
t209 = -t135 / 0.2e1;
t197 = t135 / 0.2e1;
t196 = -t136 / 0.2e1;
t208 = t136 / 0.2e1;
t207 = qJD(2) * t216;
t101 = Icges(4,2) * t141 + t176;
t177 = Icges(5,1) * t140;
t178 = Icges(4,1) * t140;
t131 = t136 * pkin(5);
t185 = rSges(4,1) * t141;
t159 = pkin(2) + t185;
t166 = rSges(4,2) * t170 + rSges(4,3) * t136;
t46 = -t135 * t159 + t131 + t166;
t120 = rSges(4,2) * t168;
t47 = -t120 + t159 * t136 + (rSges(4,3) + pkin(5)) * t135;
t108 = rSges(4,1) * t140 + rSges(4,2) * t141;
t92 = t108 * t135;
t93 = t108 * t136;
t99 = Icges(5,2) * t141 + t175;
t6 = (t106 / 0.2e1 - t101 / 0.2e1 + t104 / 0.2e1 - t99 / 0.2e1) * t140 + (t178 / 0.2e1 - t174 / 0.2e1 + t177 / 0.2e1 - t173 / 0.2e1 - t214) * t141 + m(5) * (-t226 * t45 + t227 * t44) + m(4) * (t46 * t92 - t47 * t93);
t206 = t6 * qJD(2);
t133 = t135 ^ 2;
t134 = t136 ^ 2;
t164 = t133 + t134;
t205 = -t177 - t178 + t214;
t186 = m(5) * qJD(3);
t163 = qJD(3) * t135;
t162 = t213 * t209;
t161 = t213 * t208;
t160 = t98 / 0.2e1 + t97 / 0.2e1;
t157 = rSges(5,2) * t140 - t184 - t193;
t145 = (t136 * t212 - t211 + t215) * t208 + (t221 * t136 + (t141 * t220 + t218) * t135) * t196 + (t135 * t212 + t210 + t219) * t197;
t144 = t211 * t209 + t215 * t197 + ((-t218 + t229) * t136 + t210 + t217 - t225) * t196;
t110 = -rSges(4,2) * t140 + t185;
t79 = t157 * t136;
t77 = t157 * t135;
t43 = -t135 * t92 - t136 * t93;
t31 = t222 * t164;
t16 = 0.2e1 * t35;
t15 = t35 + t30;
t14 = 0.2e1 * t30;
t1 = t135 * t145 + t136 * t144;
t2 = [0, 0, 0.2e1 * (m(4) * t43 / 0.2e1 - m(5) * t31 / 0.2e1) * qJD(3), 0; 0, qJD(3) * t6 + qJD(4) * t216, t206 + (t44 * t79 + t45 * t77) * t186 + t15 * qJD(4) + (m(4) * (t108 * t93 - t110 * t47) + t160 * t135 - t145) * t163 + ((m(4) * (-t108 * t92 - t110 * t46) + t160 * t136 - t144) * t136 + (((-t101 - t99) * t136 + t228) * t197 + (-t169 * t234 + t220 + t231) * t196) * t141 + ((t136 * t205 - t223) * t197 + (t135 * t205 + t224) * t196) * t140) * qJD(3), qJD(3) * t15 + t207; 0, t1 * qJD(3) + t16 * qJD(4) - t206, t1 * qJD(2) + (t135 * t162 + t136 * t161) * t163 * t136 + (m(4) * (t164 * t110 * t108 + (t135 * (rSges(4,1) * t169 - t166) + t136 * (rSges(4,1) * t167 + t135 * rSges(4,3) - t120)) * t43) + m(5) * (-((-t135 * pkin(2) + t131 - t44) * t135 + ((-pkin(2) + t152) * t136 + (rSges(5,3) - pkin(5)) * t135 + t165) * t136) * t31 - t227 * t77 - t226 * t79) + t161 * t135 * t133 + t162 * t134 * t136) * qJD(3), t16 * qJD(2); 0, t14 * qJD(3) - t207, t14 * qJD(2) + (t135 * t79 - t136 * t77) * t186, 0;];
Cq = t2;
