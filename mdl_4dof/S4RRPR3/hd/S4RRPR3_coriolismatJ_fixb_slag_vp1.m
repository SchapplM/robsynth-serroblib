% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:31
% DurationCPUTime: 1.07s
% Computational Cost: add. (7826->142), mult. (5296->199), div. (0->0), fcn. (4650->8), ass. (0->96)
t141 = cos(qJ(4));
t137 = Icges(5,4) * t141;
t139 = sin(qJ(4));
t118 = -Icges(5,2) * t139 + t137;
t119 = Icges(5,1) * t139 + t137;
t210 = t118 + t119;
t138 = qJ(1) + qJ(2);
t134 = pkin(7) + t138;
t132 = sin(t134);
t133 = cos(t134);
t135 = sin(t138);
t136 = cos(t138);
t185 = cos(qJ(1)) * pkin(1);
t186 = sin(qJ(1)) * pkin(1);
t187 = pkin(2) * t136;
t188 = pkin(2) * t135;
t172 = t141 * rSges(5,1);
t157 = pkin(3) + t172;
t167 = t132 * t139;
t158 = rSges(5,2) * t167 + t133 * rSges(5,3);
t65 = t133 * pkin(6) - t157 * t132 + t158 - t188;
t62 = t65 - t186;
t164 = t133 * t139;
t115 = rSges(5,2) * t164;
t66 = t187 - t115 + t157 * t133 + (rSges(5,3) + pkin(6)) * t132;
t63 = t66 + t185;
t10 = m(4) * (t185 * (-t132 * rSges(4,1) - t133 * rSges(4,2) - t188) + (t133 * rSges(4,1) - t132 * rSges(4,2) + t187) * t186) + m(5) * (-t66 * t62 + t63 * t65) + m(3) * (t185 * (-t135 * rSges(3,1) - t136 * rSges(3,2)) + (t136 * rSges(3,1) - t135 * rSges(3,2)) * t186);
t208 = t10 * qJD(1);
t207 = t10 * qJD(2);
t170 = Icges(5,4) * t139;
t117 = Icges(5,2) * t141 + t170;
t120 = Icges(5,1) * t141 - t170;
t154 = t210 * t141 / 0.2e1 + (-t117 / 0.2e1 + t120 / 0.2e1) * t139;
t206 = t132 ^ 2;
t205 = t133 ^ 2;
t121 = t139 * rSges(5,1) + t141 * rSges(5,2);
t97 = t121 * t132;
t98 = t121 * t133;
t25 = t62 * t97 - t63 * t98;
t26 = t65 * t97 - t66 * t98;
t198 = m(5) * (t26 + t25);
t197 = m(5) * ((-t63 + t66) * t133 + (t62 - t65) * t132) * t121;
t194 = m(5) * t25;
t193 = m(5) * t26;
t192 = -t132 / 0.2e1;
t191 = t132 / 0.2e1;
t190 = -t133 / 0.2e1;
t166 = t132 * t141;
t80 = Icges(5,4) * t166 - Icges(5,2) * t167 - Icges(5,6) * t133;
t173 = t139 * t80;
t83 = Icges(5,5) * t132 + t120 * t133;
t71 = t83 * t166;
t116 = Icges(5,5) * t141 - Icges(5,6) * t139;
t165 = t133 * t116;
t79 = Icges(5,3) * t132 + t165;
t156 = t133 * t79 - t71;
t81 = Icges(5,6) * t132 + t118 * t133;
t32 = -t81 * t167 - t156;
t78 = Icges(5,5) * t166 - Icges(5,6) * t167 - Icges(5,3) * t133;
t113 = Icges(5,4) * t167;
t82 = Icges(5,1) * t166 - Icges(5,5) * t133 - t113;
t19 = t32 * t132 - (-(-t141 * t82 + t173) * t132 - t133 * t78) * t133;
t163 = t133 * t141;
t181 = -t132 * t78 - t82 * t163;
t33 = -t80 * t164 - t181;
t180 = t132 * t79 + t83 * t163;
t34 = -t81 * t164 + t180;
t20 = t34 * t132 - t33 * t133;
t155 = t139 * t81 - t78;
t8 = (t155 * t133 - t180 + t34) * t133 + (t155 * t132 + t156 + t33) * t132;
t9 = (t32 - t71 + (t79 + t173) * t133 + t181) * t133 + t180 * t132;
t2 = (t20 / 0.2e1 - t9 / 0.2e1) * t133 + (t8 / 0.2e1 + t19 / 0.2e1) * t132;
t189 = t2 * qJD(4);
t179 = t119 * t132 + t80;
t178 = -t119 * t133 - t81;
t177 = -Icges(5,2) * t166 - t113 + t82;
t176 = -t117 * t133 + t83;
t153 = t198 / 0.2e1 + t154;
t149 = Icges(5,5) * t139 + Icges(5,6) * t141;
t148 = (t132 * t98 - t133 * t97) * t121;
t143 = (-t117 + t120) * t141 - t210 * t139;
t147 = t133 * t9 / 0.2e1 + (t8 + t19) * t192 + (t132 * t116 + t143 * t133 + t178 * t139 + t176 * t141) * t191 + (t143 * t132 - t179 * t139 + t177 * t141 - t165 + t20) * t190;
t146 = -t154 + (t191 + t192) * (t139 * t82 + t141 * t80);
t145 = t177 * t139 + t179 * t141;
t144 = -t176 * t139 + t178 * t141;
t123 = -t139 * rSges(5,2) + t172;
t92 = t133 * t149;
t91 = t149 * t132;
t60 = -t132 * t97 - t133 * t98;
t23 = t154 + t193;
t21 = t154 + t194;
t15 = t197 / 0.2e1;
t5 = -t197 / 0.2e1 + t153;
t4 = t15 + t153;
t3 = t15 - t198 / 0.2e1 + t146;
t1 = [t21 * qJD(4) + t207, t4 * qJD(4) + t207 + t208, 0, t21 * qJD(1) + t4 * qJD(2) + (((-t132 * t63 - t133 * t62) * t123 + t148) * m(5) + t147) * qJD(4); t5 * qJD(4) - t208, qJD(4) * t23, 0, t5 * qJD(1) + t23 * qJD(2) + (((-t132 * t66 - t133 * t65) * t123 + t148) * m(5) + t147) * qJD(4); 0, 0, 0, m(5) * t60 * qJD(4); (t146 - t194) * qJD(1) + t3 * qJD(2) + t189, t3 * qJD(1) + (t146 - t193) * qJD(2) + t189, 0, (m(5) * ((t132 * (rSges(5,1) * t166 - t158) + t133 * (rSges(5,1) * t163 + t132 * rSges(5,3) - t115)) * t60 + (t205 + t206) * t123 * t121) + (-t206 * t92 + (t145 * t133 + (t91 + t144) * t132) * t133) * t191 + (-t205 * t91 + (t144 * t132 + (t92 + t145) * t133) * t132) * t190) * qJD(4) + (qJD(1) + qJD(2)) * t2;];
Cq = t1;
