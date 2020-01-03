% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:36
% DurationCPUTime: 0.97s
% Computational Cost: add. (7364->135), mult. (4940->194), div. (0->0), fcn. (4362->6), ass. (0->93)
t129 = cos(qJ(4));
t126 = Icges(5,4) * t129;
t128 = sin(qJ(4));
t110 = -Icges(5,2) * t128 + t126;
t111 = Icges(5,1) * t128 + t126;
t193 = t110 + t111;
t127 = pkin(7) + qJ(2);
t125 = qJ(3) + t127;
t121 = sin(t125);
t122 = cos(t125);
t171 = pkin(2) * cos(t127);
t172 = pkin(2) * sin(t127);
t192 = m(4) * (t171 * (-rSges(4,1) * t121 - rSges(4,2) * t122) + (rSges(4,1) * t122 - t121 * rSges(4,2)) * t172);
t161 = rSges(5,1) * t129;
t144 = pkin(3) + t161;
t154 = t121 * t128;
t145 = rSges(5,2) * t154 + t122 * rSges(5,3);
t63 = t122 * pkin(6) - t144 * t121 + t145;
t61 = t63 - t172;
t151 = t122 * t128;
t107 = rSges(5,2) * t151;
t64 = -t107 + t144 * t122 + (rSges(5,3) + pkin(6)) * t121;
t62 = t64 + t171;
t191 = m(5) * (-t64 * t61 + t62 * t63);
t11 = t191 + t192;
t189 = t11 * qJD(2);
t156 = Icges(5,4) * t128;
t109 = Icges(5,2) * t129 + t156;
t112 = Icges(5,1) * t129 - t156;
t141 = t193 * t129 / 0.2e1 + (-t109 / 0.2e1 + t112 / 0.2e1) * t128;
t188 = t121 ^ 2;
t187 = t122 ^ 2;
t113 = t128 * rSges(5,1) + rSges(5,2) * t129;
t89 = t113 * t121;
t90 = t113 * t122;
t25 = t61 * t89 - t62 * t90;
t26 = t63 * t89 - t64 * t90;
t182 = m(5) * (t26 + t25);
t181 = m(5) * ((-t62 + t64) * t122 + (t61 - t63) * t121) * t113;
t178 = m(5) * t25;
t177 = m(5) * t26;
t176 = -t121 / 0.2e1;
t175 = t121 / 0.2e1;
t174 = -t122 / 0.2e1;
t153 = t121 * t129;
t76 = Icges(5,4) * t153 - Icges(5,2) * t154 - Icges(5,6) * t122;
t158 = t128 * t76;
t79 = Icges(5,5) * t121 + t112 * t122;
t67 = t79 * t153;
t108 = Icges(5,5) * t129 - Icges(5,6) * t128;
t152 = t122 * t108;
t75 = Icges(5,3) * t121 + t152;
t143 = t122 * t75 - t67;
t77 = Icges(5,6) * t121 + t110 * t122;
t32 = -t77 * t154 - t143;
t74 = Icges(5,5) * t153 - Icges(5,6) * t154 - Icges(5,3) * t122;
t105 = Icges(5,4) * t154;
t78 = Icges(5,1) * t153 - Icges(5,5) * t122 - t105;
t19 = t121 * t32 - t122 * (-(-t129 * t78 + t158) * t121 - t122 * t74);
t150 = t122 * t129;
t167 = -t121 * t74 - t78 * t150;
t33 = -t76 * t151 - t167;
t166 = t121 * t75 + t79 * t150;
t34 = -t77 * t151 + t166;
t20 = t121 * t34 - t122 * t33;
t142 = t128 * t77 - t74;
t8 = (t142 * t122 - t166 + t34) * t122 + (t142 * t121 + t143 + t33) * t121;
t9 = (t32 - t67 + (t75 + t158) * t122 + t167) * t122 + t166 * t121;
t2 = (t20 / 0.2e1 - t9 / 0.2e1) * t122 + (t8 / 0.2e1 + t19 / 0.2e1) * t121;
t173 = t2 * qJD(4);
t165 = t111 * t121 + t76;
t164 = -t111 * t122 - t77;
t163 = -Icges(5,2) * t153 - t105 + t78;
t162 = -t109 * t122 + t79;
t140 = t182 / 0.2e1 + t141;
t136 = Icges(5,5) * t128 + Icges(5,6) * t129;
t135 = (t121 * t90 - t122 * t89) * t113;
t130 = (-t109 + t112) * t129 - t193 * t128;
t134 = t122 * t9 / 0.2e1 + (t8 + t19) * t176 + (t121 * t108 + t130 * t122 + t164 * t128 + t162 * t129) * t175 + (t130 * t121 - t165 * t128 + t163 * t129 - t152 + t20) * t174;
t133 = -t141 + (t175 + t176) * (t128 * t78 + t129 * t76);
t132 = t163 * t128 + t165 * t129;
t131 = -t162 * t128 + t164 * t129;
t114 = -t128 * rSges(5,2) + t161;
t84 = t122 * t136;
t83 = t136 * t121;
t59 = -t121 * t89 - t122 * t90;
t22 = t141 + t177;
t21 = t141 + t178;
t13 = t181 / 0.2e1;
t7 = -t181 / 0.2e1 + t140;
t6 = t13 + t140;
t3 = t13 - t182 / 0.2e1 + t133;
t1 = [0, 0, 0, m(5) * t59 * qJD(4); 0, qJD(3) * t11 + qJD(4) * t21, t189 + t6 * qJD(4) + 0.2e1 * (t191 / 0.2e1 + t192 / 0.2e1) * qJD(3), t21 * qJD(2) + t6 * qJD(3) + (((-t121 * t62 - t122 * t61) * t114 + t135) * m(5) + t134) * qJD(4); 0, t7 * qJD(4) - t189, qJD(4) * t22, t7 * qJD(2) + t22 * qJD(3) + (((-t121 * t64 - t122 * t63) * t114 + t135) * m(5) + t134) * qJD(4); 0, (t133 - t178) * qJD(2) + t3 * qJD(3) + t173, t3 * qJD(2) + (t133 - t177) * qJD(3) + t173, (m(5) * ((t121 * (rSges(5,1) * t153 - t145) + t122 * (rSges(5,1) * t150 + t121 * rSges(5,3) - t107)) * t59 + (t187 + t188) * t114 * t113) + (-t188 * t84 + (t132 * t122 + (t83 + t131) * t121) * t122) * t175 + (-t187 * t83 + (t131 * t121 + (t84 + t132) * t122) * t121) * t174) * qJD(4) + (qJD(2) + qJD(3)) * t2;];
Cq = t1;
