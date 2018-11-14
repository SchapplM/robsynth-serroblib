% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:27
% DurationCPUTime: 1.17s
% Computational Cost: add. (2305->158), mult. (3333->181), div. (0->0), fcn. (2644->9), ass. (0->90)
t93 = sin(pkin(4));
t96 = cos(qJ(1));
t147 = t93 * t96;
t159 = rSges(5,1) + pkin(3);
t92 = sin(pkin(6));
t95 = sin(qJ(1));
t128 = pkin(4) + pkin(6);
t129 = pkin(4) - pkin(6);
t101 = cos(t129) / 0.2e1 + cos(t128) / 0.2e1;
t98 = t96 * t101;
t57 = t92 * t95 - t98;
t113 = -t57 * rSges(5,2) + t159 * t147;
t155 = rSges(5,3) + qJ(4);
t100 = sin(t128) / 0.2e1 - sin(t129) / 0.2e1;
t94 = cos(pkin(6));
t58 = t100 * t96 + t95 * t94;
t160 = -t155 * t58 + t113;
t59 = t101 * t95 + t96 * t92;
t99 = t95 * t100;
t60 = t94 * t96 - t99;
t117 = t60 * pkin(2) + qJ(3) * t59;
t130 = qJD(3) * t57;
t135 = qJ(2) * t93;
t91 = t96 * pkin(1);
t136 = t95 * t135 + t91;
t131 = qJD(2) * t93;
t80 = t96 * t131;
t157 = qJD(1) * t136 - t80;
t158 = qJD(1) * t117 + t130 + t157;
t156 = pkin(1) * t95;
t148 = t93 * t95;
t107 = -rSges(3,1) * t58 + rSges(3,2) * t57 + rSges(3,3) * t147;
t67 = -t96 * t135 + t156;
t154 = t107 - t67;
t48 = t57 * qJ(3);
t153 = -pkin(2) * t58 - t48;
t133 = qJD(1) * t95;
t38 = -qJD(1) * t98 + t92 * t133;
t39 = t58 * qJD(1);
t46 = qJD(3) * t59;
t105 = -t39 * pkin(2) - qJ(3) * t38 + t46;
t79 = t95 * t131;
t137 = -qJD(1) * t67 + t79;
t132 = qJD(1) * t96;
t121 = t93 * t132;
t138 = qJ(2) * t121 + t79;
t152 = -qJD(1) * t153 + t105 - t137 + t138 - t46;
t22 = rSges(4,1) * t148 - t60 * rSges(4,2) + t59 * rSges(4,3);
t12 = qJD(1) * t22 + t158;
t43 = qJD(4) * t58;
t151 = t43 + t158;
t40 = t59 * qJD(1);
t143 = t40 * qJ(3) + t130;
t41 = -qJD(1) * t99 + t94 * t132;
t20 = pkin(2) * t41 + t143;
t146 = -t20 - t157;
t145 = t59 * rSges(5,2) + t159 * t148 + t155 * t60;
t144 = -t67 + t153;
t142 = t46 + t79;
t120 = qJD(1) * t131;
t141 = qJD(1) * (-pkin(1) * t133 + t138) + t95 * t120;
t73 = t96 * t120;
t139 = -qJD(3) * t38 + t73;
t44 = qJD(4) * t60;
t127 = -pkin(2) - t155;
t125 = -t39 * rSges(3,1) + t38 * rSges(3,2) + rSges(3,3) * t121;
t124 = t60 * rSges(3,1) - t59 * rSges(3,2) + rSges(3,3) * t148;
t123 = rSges(4,1) * t121 + t39 * rSges(4,2) - t38 * rSges(4,3);
t122 = t93 * t133;
t119 = -t40 * rSges(5,2) - t43;
t118 = rSges(4,1) * t147 - rSges(4,3) * t57;
t116 = qJD(1) * t105 + qJD(3) * t40 + t141;
t114 = -t48 - t67;
t111 = -rSges(3,1) * t41 + rSges(3,2) * t40;
t110 = rSges(4,2) * t41 - rSges(4,3) * t40;
t109 = t117 + t136;
t106 = rSges(4,2) * t58 + t118;
t102 = -t38 * rSges(5,2) + t159 * t121 - t155 * t39 + t44;
t19 = qJD(1) * t124 + t157;
t18 = t154 * qJD(1) + t79;
t14 = t73 + (-rSges(3,3) * t122 + t111 - t157) * qJD(1);
t13 = qJD(1) * t125 + t141;
t11 = (t106 + t144) * qJD(1) + t142;
t8 = (-rSges(4,1) * t122 + t110 + t146) * qJD(1) + t139;
t7 = qJD(1) * t123 + t116;
t6 = t145 * qJD(1) + t151;
t5 = t44 + (t144 + t160) * qJD(1) + t142;
t2 = -qJD(4) * t39 + (-t122 * t159 - t155 * t41 + t119 + t146) * qJD(1) + t139;
t1 = qJD(1) * t102 + qJD(4) * t41 + t116;
t3 = [(t1 * (t109 + t145) + (t127 * t58 + t113 + t114) * t2 + (t102 - t44 + (-t156 - t160) * qJD(1) + t152) * t6 + (t127 * t41 + t119 - t143 + t80 + (-t91 + (-qJ(2) - t159) * t148 + t145) * qJD(1) + t151) * t5) * m(5) + (t8 * ((rSges(4,2) - pkin(2)) * t58 + t114 + t118) + t7 * (t109 + t22) + (t123 + (-t106 - t156) * qJD(1) + t152) * t12 + (t110 - t20 + t80 + (-t91 + (-rSges(4,1) - qJ(2)) * t148) * qJD(1) + t12) * t11) * m(4) + (t14 * t154 + t18 * (t111 + t80) + t13 * (t124 + t136) + t19 * (t125 + t138) + (-t18 * t91 + (-t19 * pkin(1) + t18 * (-rSges(3,3) - qJ(2)) * t93) * t95) * qJD(1) - (qJD(1) * t107 + t137 - t18) * t19) * m(3); 0.2e1 * (m(3) * (-t13 * t96 + t14 * t95) / 0.2e1 + m(4) * (-t7 * t96 + t8 * t95) / 0.2e1 + m(5) * (-t1 * t96 + t2 * t95) / 0.2e1) * t93; m(4) * (-t11 * t38 + t12 * t40 + t7 * t57 + t8 * t59) + m(5) * (t1 * t57 + t2 * t59 - t5 * t38 + t6 * t40) + 0.2e1 * (-m(4) * (-t11 * t57 + t12 * t59) / 0.2e1 - m(5) * (-t5 * t57 + t6 * t59) / 0.2e1) * qJD(1); (t1 * t58 + t2 * t60 - t5 * t39 + t6 * t41 - (-t5 * t58 + t6 * t60) * qJD(1)) * m(5);];
tauc  = t3(:);
