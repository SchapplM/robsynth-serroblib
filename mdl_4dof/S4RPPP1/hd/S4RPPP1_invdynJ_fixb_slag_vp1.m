% Calculate vector of inverse dynamics joint torques for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:20
% EndTime: 2019-03-08 18:26:22
% DurationCPUTime: 1.92s
% Computational Cost: add. (1514->218), mult. (3814->249), div. (0->0), fcn. (3553->6), ass. (0->106)
t199 = rSges(5,1) + pkin(3);
t197 = rSges(5,3) + qJ(4);
t119 = sin(pkin(4));
t123 = cos(qJ(1));
t164 = t119 * t123;
t118 = sin(pkin(6));
t120 = cos(pkin(6));
t122 = sin(qJ(1));
t121 = cos(pkin(4));
t162 = t121 * t123;
t78 = t118 * t122 - t120 * t162;
t201 = -t78 * rSges(5,2) + t199 * t164;
t79 = t118 * t162 + t120 * t122;
t130 = -t197 * t79 + t201;
t63 = t78 * qJ(3);
t140 = -t79 * pkin(2) - t63;
t198 = pkin(1) * t122;
t86 = -qJ(2) * t164 + t198;
t196 = t140 - t86;
t139 = t196 + t130;
t141 = rSges(4,1) * t164 - t78 * rSges(4,3);
t151 = t79 * rSges(4,2) + t141;
t153 = t196 + t151;
t159 = qJD(2) * t119;
t143 = qJD(1) * t159;
t158 = qJDD(2) * t119;
t104 = t122 * t159;
t160 = qJD(1) * t123;
t146 = t119 * t160;
t170 = qJ(2) * t146 + t104;
t161 = qJD(1) * t122;
t192 = pkin(1) * t161;
t117 = t123 * pkin(1);
t165 = t119 * t122;
t87 = qJ(2) * t165 + t117;
t129 = -t123 * t158 + qJD(1) * (t170 - t192) + qJDD(1) * t87 + t122 * t143;
t145 = t120 * t160;
t148 = t118 * t161;
t52 = -t121 * t145 + t148;
t53 = t79 * qJD(1);
t163 = t121 * t122;
t80 = t118 * t123 + t120 * t163;
t62 = qJD(3) * t80;
t131 = -t53 * pkin(2) - qJ(3) * t52 + t62;
t81 = -t118 * t163 + t120 * t123;
t35 = t81 * pkin(2) + qJ(3) * t80;
t54 = t80 * qJD(1);
t125 = qJD(1) * t131 + qJD(3) * t54 + qJDD(1) * t35 + qJDD(3) * t78 + t129;
t60 = qJD(4) * t81;
t126 = -t52 * rSges(5,2) + t199 * t146 - t197 * t53 + t60;
t177 = t80 * rSges(5,2) + t199 * t165 + t197 * t81;
t55 = -t121 * t148 + t145;
t2 = qJD(1) * t126 + qJD(4) * t55 + qJDD(1) * t177 + qJDD(4) * t79 + t125;
t200 = -g(2) + t2;
t105 = t123 * t159;
t195 = qJD(1) * t87 - t105;
t171 = -qJD(1) * t86 + t104;
t194 = -qJD(1) * t140 + t131 + t170 - t171 - t62;
t168 = qJD(3) * t78;
t193 = -qJD(1) * t35 - t168 - t195;
t149 = -t79 * rSges(3,1) + t78 * rSges(3,2) + rSges(3,3) * t164;
t175 = -t86 + t149;
t25 = rSges(4,1) * t165 - t81 * rSges(4,2) + t80 * rSges(4,3);
t189 = -qJD(1) * t25 + t193;
t59 = qJD(4) * t79;
t188 = -t59 + t193;
t186 = t121 ^ 2;
t184 = -m(4) - m(5);
t183 = t121 / 0.2e1;
t182 = t122 / 0.2e1;
t181 = -t123 / 0.2e1;
t176 = t54 * qJ(3) + t168;
t22 = pkin(2) * t55 + t176;
t179 = -t22 - t195;
t173 = t122 * t158 + t123 * t143;
t172 = t104 + t62;
t166 = t119 * t120;
t157 = qJDD(3) * t120;
t156 = 0.2e1 * t121;
t155 = -pkin(2) - t197;
t154 = -t53 * rSges(3,1) + t52 * rSges(3,2) + rSges(3,3) * t146;
t152 = rSges(4,1) * t146 + t53 * rSges(4,2) - t52 * rSges(4,3);
t29 = t81 * rSges(3,1) - t80 * rSges(3,2) + rSges(3,3) * t165;
t147 = t119 * t161;
t144 = -t54 * rSges(5,2) - t59;
t138 = -qJD(3) * t52 + qJDD(3) * t80 + t173;
t136 = -t63 - t86;
t135 = -rSges(3,1) * t55 + rSges(3,2) * t54;
t134 = rSges(4,2) * t55 - rSges(4,3) * t54;
t93 = rSges(2,1) * t123 - rSges(2,2) * t122;
t92 = rSges(2,1) * t122 + rSges(2,2) * t123;
t133 = t35 + t87;
t116 = qJDD(2) * t121;
t85 = -t119 * t157 + t116;
t40 = t116 + (qJDD(4) * t118 - t157) * t119;
t20 = qJD(1) * t29 + t195;
t19 = qJD(1) * t175 + t104;
t13 = qJD(1) * t153 + t172;
t10 = qJD(1) * t154 + qJDD(1) * t29 + t129;
t9 = t175 * qJDD(1) + (-rSges(3,3) * t147 + t135 - t195) * qJD(1) + t173;
t8 = t177 * qJD(1) - t188;
t7 = qJD(1) * t139 + t172 + t60;
t4 = qJD(1) * t152 + qJDD(1) * t25 + t125;
t3 = t153 * qJDD(1) + (-rSges(4,1) * t147 + t134 + t179) * qJD(1) + t138;
t1 = -qJD(4) * t53 + qJDD(4) * t81 + t139 * qJDD(1) + (-t147 * t199 - t197 * t55 + t144 + t179) * qJD(1) + t138;
t5 = [-m(2) * (-g(1) * t92 + g(2) * t93) + (-g(1) * t139 + (t155 * t79 + t136 + t201) * t1 + (-t130 * qJD(1) + t126 - t192 + t194 - t60) * t8 + (t155 * t55 + t105 + t144 - t176 - t188 + (t177 - t117 + (-qJ(2) - t199) * t165) * qJD(1)) * t7 + t200 * (t133 + t177)) * m(5) + (-g(1) * t153 + t3 * ((rSges(4,2) - pkin(2)) * t79 + t136 + t141) + (-g(2) + t4) * (t133 + t25) - (t152 + (-t151 - t198) * qJD(1) + t194) * t189 + (-t189 + t105 + t134 - t22 + (-t117 + (-rSges(4,1) - qJ(2)) * t165) * qJD(1)) * t13) * m(4) + (-(qJD(1) * t149 + t171 - t19) * t20 + t19 * (t105 + t135) + t20 * (t154 + t170) + (-t19 * t117 + (-t20 * pkin(1) + t19 * (-rSges(3,3) - qJ(2)) * t119) * t122) * qJD(1) + (-g(2) + t10) * (t29 + t87) + (-g(1) + t9) * t175) * m(3) + (Icges(2,3) + m(2) * (t92 ^ 2 + t93 ^ 2) + (Icges(3,3) + Icges(4,1) + Icges(5,1)) * t186 + ((Icges(3,2) + Icges(4,3) + Icges(5,2)) * t119 * t120 ^ 2 + (-Icges(5,4) - Icges(4,5) + Icges(3,6)) * t156 * t120 + ((-0.2e1 * Icges(5,6) * t120 + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t118 + 0.2e1 * (Icges(3,4) + Icges(4,6)) * t120) * t119 + (Icges(3,5) - Icges(4,4) + Icges(5,5)) * t156) * t118) * t119) * qJDD(1); (-m(3) + t184) * (g(3) * t121 + (g(1) * t122 - g(2) * t123) * t119) + 0.2e1 * (t40 * t183 + (t1 * t182 + t181 * t2) * t119) * m(5) + 0.2e1 * (t85 * t183 + (t181 * t4 + t182 * t3) * t119) * m(4) + 0.2e1 * (qJDD(2) * t186 / 0.2e1 + (t10 * t181 + t182 * t9) * t119) * m(3); t184 * (g(1) * t80 + g(2) * t78 - g(3) * t166) + m(4) * (-t13 * t52 - t85 * t166 - t189 * t54 + t3 * t80 + t4 * t78) + m(5) * (t1 * t80 - t40 * t166 + t2 * t78 - t7 * t52 + t8 * t54) + 0.2e1 * (-m(4) * (-t13 * t78 - t189 * t80) / 0.2e1 - m(5) * (-t7 * t78 + t8 * t80) / 0.2e1) * qJD(1); (-t53 * t7 + t55 * t8 + (-t8 * qJD(1) - g(1) + t1) * t81 + (qJD(1) * t7 + t200) * t79 + (t40 - g(3)) * t118 * t119) * m(5);];
tau  = t5;
