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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:27
% DurationCPUTime: 2.06s
% Computational Cost: add. (3194->223), mult. (4374->249), div. (0->0), fcn. (3553->10), ass. (0->111)
t208 = rSges(5,1) + pkin(3);
t206 = rSges(5,3) + qJ(4);
t125 = sin(pkin(4));
t129 = cos(qJ(1));
t176 = t125 * t129;
t124 = sin(pkin(6));
t128 = sin(qJ(1));
t169 = pkin(4) - pkin(6);
t145 = cos(t169) / 0.2e1;
t168 = pkin(4) + pkin(6);
t150 = cos(t168);
t133 = t145 + t150 / 0.2e1;
t131 = t129 * t133;
t77 = t124 * t128 - t131;
t210 = -t77 * rSges(5,2) + t176 * t208;
t126 = cos(pkin(6));
t144 = sin(t168) / 0.2e1;
t149 = sin(t169);
t141 = t144 - t149 / 0.2e1;
t78 = t126 * t128 + t129 * t141;
t139 = -t206 * t78 + t210;
t61 = t77 * qJ(3);
t154 = -pkin(2) * t78 - t61;
t207 = pkin(1) * t128;
t93 = -qJ(2) * t176 + t207;
t205 = t154 - t93;
t153 = t205 + t139;
t155 = rSges(4,1) * t176 - t77 * rSges(4,3);
t164 = rSges(4,2) * t78 + t155;
t166 = t205 + t164;
t171 = qJD(2) * t125;
t157 = qJD(1) * t171;
t170 = qJDD(2) * t125;
t107 = t128 * t171;
t172 = qJD(1) * t129;
t159 = t125 * t172;
t174 = qJ(2) * t159 + t107;
t173 = qJD(1) * t128;
t201 = pkin(1) * t173;
t123 = t129 * pkin(1);
t177 = t125 * t128;
t94 = qJ(2) * t177 + t123;
t138 = -t129 * t170 + qJD(1) * (t174 - t201) + qJDD(1) * t94 + t128 * t157;
t49 = -qJD(1) * t131 + t124 * t173;
t50 = t78 * qJD(1);
t79 = t124 * t129 + t128 * t133;
t59 = qJD(3) * t79;
t140 = -pkin(2) * t50 - qJ(3) * t49 + t59;
t136 = t128 * t141;
t80 = t126 * t129 - t136;
t35 = pkin(2) * t80 + qJ(3) * t79;
t51 = t79 * qJD(1);
t132 = qJD(1) * t140 + qJD(3) * t51 + qJDD(1) * t35 + qJDD(3) * t77 + t138;
t57 = qJD(4) * t80;
t134 = -t49 * rSges(5,2) + t159 * t208 - t206 * t50 + t57;
t185 = t79 * rSges(5,2) + t177 * t208 + t206 * t80;
t52 = -qJD(1) * t136 + t126 * t172;
t2 = qJD(1) * t134 + qJD(4) * t52 + qJDD(1) * t185 + qJDD(4) * t78 + t132;
t209 = -g(2) + t2;
t108 = t129 * t171;
t204 = qJD(1) * t94 - t108;
t180 = -qJD(1) * t93 + t107;
t203 = -qJD(1) * t154 + t140 + t174 - t180 - t59;
t178 = qJD(3) * t77;
t202 = -qJD(1) * t35 - t178 - t204;
t161 = -rSges(3,1) * t78 + rSges(3,2) * t77 + rSges(3,3) * t176;
t183 = -t93 + t161;
t26 = rSges(4,1) * t177 - rSges(4,2) * t80 + rSges(4,3) * t79;
t198 = -qJD(1) * t26 + t202;
t56 = qJD(4) * t78;
t197 = -t56 + t202;
t92 = t145 - t150 / 0.2e1;
t194 = -0.2e1 * t92;
t193 = 0.2e1 * t92;
t192 = -m(4) - m(5);
t127 = cos(pkin(4));
t191 = t127 / 0.2e1;
t190 = t128 / 0.2e1;
t189 = -t129 / 0.2e1;
t184 = qJ(3) * t51 + t178;
t22 = pkin(2) * t52 + t184;
t187 = -t22 - t204;
t181 = t107 + t59;
t91 = t144 + t149 / 0.2e1;
t81 = qJDD(2) * t127 - qJDD(3) * t91;
t175 = t128 * t170 + t129 * t157;
t167 = -pkin(2) - t206;
t165 = rSges(4,1) * t159 + rSges(4,2) * t50 - rSges(4,3) * t49;
t162 = -rSges(3,1) * t50 + rSges(3,2) * t49 + rSges(3,3) * t159;
t30 = rSges(3,1) * t80 - rSges(3,2) * t79 + rSges(3,3) * t177;
t160 = t125 * t173;
t158 = -t51 * rSges(5,2) - t56;
t151 = -qJD(3) * t49 + qJDD(3) * t79 + t175;
t148 = -t61 - t93;
t147 = -rSges(3,1) * t52 + rSges(3,2) * t51;
t146 = rSges(4,2) * t52 - rSges(4,3) * t51;
t98 = rSges(2,1) * t129 - rSges(2,2) * t128;
t97 = rSges(2,1) * t128 + rSges(2,2) * t129;
t143 = t35 + t94;
t36 = qJDD(4) * t92 + t81;
t21 = qJD(1) * t30 + t204;
t20 = qJD(1) * t183 + t107;
t13 = qJD(1) * t166 + t181;
t10 = qJD(1) * t162 + qJDD(1) * t30 + t138;
t9 = t183 * qJDD(1) + (-rSges(3,3) * t160 + t147 - t204) * qJD(1) + t175;
t8 = t185 * qJD(1) - t197;
t7 = qJD(1) * t153 + t181 + t57;
t4 = qJD(1) * t165 + qJDD(1) * t26 + t132;
t3 = t166 * qJDD(1) + (-rSges(4,1) * t160 + t146 + t187) * qJD(1) + t151;
t1 = -qJD(4) * t50 + qJDD(4) * t80 + t153 * qJDD(1) + (-t160 * t208 - t206 * t52 + t158 + t187) * qJD(1) + t151;
t5 = [-m(2) * (-g(1) * t97 + g(2) * t98) + (-g(1) * t153 + (t167 * t78 + t148 + t210) * t1 + (-qJD(1) * t139 + t134 - t201 + t203 - t57) * t8 + (t167 * t52 + t108 + t158 - t184 - t197 + (t185 - t123 + (-qJ(2) - t208) * t177) * qJD(1)) * t7 + t209 * (t143 + t185)) * m(5) + (-g(1) * t166 + t3 * ((rSges(4,2) - pkin(2)) * t78 + t148 + t155) + (-g(2) + t4) * (t143 + t26) - (t165 + (-t164 - t207) * qJD(1) + t203) * t198 + (-t198 + t108 + t146 - t22 + (-t123 + (-rSges(4,1) - qJ(2)) * t177) * qJD(1)) * t13) * m(4) + (-(qJD(1) * t161 + t180 - t20) * t21 + t20 * (t108 + t147) + t21 * (t162 + t174) + (-t20 * t123 + (-t21 * pkin(1) + t20 * (-rSges(3,3) - qJ(2)) * t125) * t128) * qJD(1) + (-g(2) + t10) * (t30 + t94) + (-g(1) + t9) * t183) * m(3) + (Icges(2,3) + m(2) * (t97 ^ 2 + t98 ^ 2) + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t92 ^ 2 + (Icges(5,6) * t194 + (Icges(3,2) + Icges(4,3) + Icges(5,2)) * t91 + (Icges(3,4) + Icges(4,6)) * t193) * t91 + (Icges(4,4) * t194 + 0.2e1 * Icges(3,6) * t91 + (Icges(3,3) + Icges(4,1) + Icges(5,1)) * t127 - 0.2e1 * (Icges(5,4) + Icges(4,5)) * t91 + (Icges(3,5) + Icges(5,5)) * t193) * t127) * qJDD(1); (-m(3) + t192) * (g(3) * t127 + (g(1) * t128 - g(2) * t129) * t125) + 0.2e1 * (t36 * t191 + (t1 * t190 + t189 * t2) * t125) * m(5) + 0.2e1 * (t81 * t191 + (t189 * t4 + t190 * t3) * t125) * m(4) + 0.2e1 * (qJDD(2) * t127 ^ 2 / 0.2e1 + (t10 * t189 + t190 * t9) * t125) * m(3); t192 * (g(1) * t79 + g(2) * t77 - g(3) * t91) + m(4) * (-t13 * t49 - t198 * t51 + t3 * t79 + t4 * t77 - t81 * t91) + m(5) * (t1 * t79 + t2 * t77 - t36 * t91 - t7 * t49 + t8 * t51) + 0.2e1 * (-m(4) * (-t13 * t77 - t198 * t79) / 0.2e1 - m(5) * (-t7 * t77 + t8 * t79) / 0.2e1) * qJD(1); (-t50 * t7 + t52 * t8 + (t36 - g(3)) * t92 + (-qJD(1) * t8 - g(1) + t1) * t80 + (qJD(1) * t7 + t209) * t78) * m(5);];
tau  = t5;
