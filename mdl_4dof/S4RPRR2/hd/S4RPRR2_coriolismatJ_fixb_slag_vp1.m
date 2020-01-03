% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:09
% DurationCPUTime: 0.97s
% Computational Cost: add. (7440->137), mult. (5024->196), div. (0->0), fcn. (4440->8), ass. (0->93)
t132 = cos(qJ(4));
t128 = Icges(5,4) * t132;
t130 = sin(qJ(4));
t110 = -Icges(5,2) * t130 + t128;
t111 = Icges(5,1) * t130 + t128;
t197 = t110 + t111;
t129 = qJ(1) + pkin(7);
t127 = qJ(3) + t129;
t123 = sin(t127);
t124 = cos(t127);
t144 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t129);
t145 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t129);
t196 = m(4) * (t144 * (-rSges(4,1) * t123 - rSges(4,2) * t124) - (t124 * rSges(4,1) - t123 * rSges(4,2)) * t145);
t167 = rSges(5,1) * t132;
t150 = pkin(3) + t167;
t159 = t123 * t130;
t151 = rSges(5,2) * t159 + t124 * rSges(5,3);
t63 = t124 * pkin(6) - t150 * t123 + t151;
t61 = t145 + t63;
t157 = t124 * t130;
t107 = rSges(5,2) * t157;
t64 = -t107 + t150 * t124 + (rSges(5,3) + pkin(6)) * t123;
t62 = t144 + t64;
t195 = m(5) * (-t64 * t61 + t62 * t63);
t10 = t195 + t196;
t193 = t10 * qJD(1);
t162 = Icges(5,4) * t130;
t109 = Icges(5,2) * t132 + t162;
t112 = Icges(5,1) * t132 - t162;
t147 = t197 * t132 / 0.2e1 + (-t109 / 0.2e1 + t112 / 0.2e1) * t130;
t192 = t123 ^ 2;
t191 = t124 ^ 2;
t113 = rSges(5,1) * t130 + rSges(5,2) * t132;
t89 = t113 * t123;
t90 = t113 * t124;
t25 = t61 * t89 - t62 * t90;
t26 = t63 * t89 - t64 * t90;
t186 = m(5) * (t26 + t25);
t185 = m(5) * ((-t62 + t64) * t124 + (t61 - t63) * t123) * t113;
t182 = m(5) * t25;
t181 = m(5) * t26;
t180 = -t123 / 0.2e1;
t179 = t123 / 0.2e1;
t178 = -t124 / 0.2e1;
t158 = t123 * t132;
t76 = Icges(5,4) * t158 - Icges(5,2) * t159 - Icges(5,6) * t124;
t164 = t130 * t76;
t79 = Icges(5,5) * t123 + t112 * t124;
t67 = t79 * t158;
t108 = Icges(5,5) * t132 - Icges(5,6) * t130;
t160 = t108 * t124;
t75 = Icges(5,3) * t123 + t160;
t149 = t124 * t75 - t67;
t77 = Icges(5,6) * t123 + t110 * t124;
t32 = -t77 * t159 - t149;
t74 = Icges(5,5) * t158 - Icges(5,6) * t159 - Icges(5,3) * t124;
t105 = Icges(5,4) * t159;
t78 = Icges(5,1) * t158 - Icges(5,5) * t124 - t105;
t19 = t123 * t32 - t124 * (-(-t132 * t78 + t164) * t123 - t124 * t74);
t156 = t124 * t132;
t173 = -t123 * t74 - t78 * t156;
t33 = -t76 * t157 - t173;
t172 = t123 * t75 + t79 * t156;
t34 = -t77 * t157 + t172;
t20 = t123 * t34 - t124 * t33;
t148 = t130 * t77 - t74;
t8 = (t148 * t124 - t172 + t34) * t124 + (t148 * t123 + t149 + t33) * t123;
t9 = (t32 - t67 + (t75 + t164) * t124 + t173) * t124 + t172 * t123;
t2 = (t20 / 0.2e1 - t9 / 0.2e1) * t124 + (t8 / 0.2e1 + t19 / 0.2e1) * t123;
t177 = t2 * qJD(4);
t171 = t111 * t123 + t76;
t170 = -t111 * t124 - t77;
t169 = -Icges(5,2) * t158 - t105 + t78;
t168 = -t109 * t124 + t79;
t146 = t186 / 0.2e1 + t147;
t140 = Icges(5,5) * t130 + Icges(5,6) * t132;
t139 = (t123 * t90 - t124 * t89) * t113;
t134 = (-t109 + t112) * t132 - t197 * t130;
t138 = t124 * t9 / 0.2e1 + (t8 + t19) * t180 + (t108 * t123 + t134 * t124 + t170 * t130 + t168 * t132) * t179 + (t134 * t123 - t171 * t130 + t169 * t132 - t160 + t20) * t178;
t137 = -t147 + (t179 + t180) * (t130 * t78 + t132 * t76);
t136 = t169 * t130 + t171 * t132;
t135 = -t168 * t130 + t170 * t132;
t115 = -rSges(5,2) * t130 + t167;
t84 = t124 * t140;
t83 = t140 * t123;
t59 = -t123 * t89 - t124 * t90;
t23 = t147 + t181;
t21 = t147 + t182;
t13 = t185 / 0.2e1;
t5 = -t185 / 0.2e1 + t146;
t4 = t13 + t146;
t3 = t13 - t186 / 0.2e1 + t137;
t1 = [qJD(3) * t10 + qJD(4) * t21, 0, t193 + t4 * qJD(4) + 0.2e1 * (t195 / 0.2e1 + t196 / 0.2e1) * qJD(3), t21 * qJD(1) + t4 * qJD(3) + (((-t123 * t62 - t124 * t61) * t115 + t139) * m(5) + t138) * qJD(4); 0, 0, 0, m(5) * t59 * qJD(4); t5 * qJD(4) - t193, 0, qJD(4) * t23, t5 * qJD(1) + t23 * qJD(3) + (((-t123 * t64 - t124 * t63) * t115 + t139) * m(5) + t138) * qJD(4); (t137 - t182) * qJD(1) + t3 * qJD(3) + t177, 0, t3 * qJD(1) + (t137 - t181) * qJD(3) + t177, (m(5) * ((t123 * (rSges(5,1) * t158 - t151) + t124 * (rSges(5,1) * t156 + t123 * rSges(5,3) - t107)) * t59 + (t191 + t192) * t115 * t113) + (-t192 * t84 + (t136 * t124 + (t83 + t135) * t123) * t124) * t179 + (-t191 * t83 + (t135 * t123 + (t84 + t136) * t124) * t123) * t178) * qJD(4) + (qJD(1) + qJD(3)) * t2;];
Cq = t1;
