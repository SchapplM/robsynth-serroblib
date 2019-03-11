% Calculate potential energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:36
% EndTime: 2019-03-09 04:21:36
% DurationCPUTime: 0.49s
% Computational Cost: add. (477->127), mult. (1189->175), div. (0->0), fcn. (1483->14), ass. (0->61)
t135 = sin(pkin(12));
t140 = cos(pkin(6));
t147 = cos(qJ(1));
t138 = cos(pkin(12));
t144 = sin(qJ(1));
t163 = t144 * t138;
t122 = -t135 * t147 - t140 * t163;
t136 = sin(pkin(7));
t139 = cos(pkin(7));
t137 = sin(pkin(6));
t169 = t137 * t144;
t113 = -t122 * t136 + t139 * t169;
t170 = t137 * t138;
t119 = -t136 * t170 + t139 * t140;
t178 = rSges(6,3) + pkin(10);
t177 = pkin(11) + rSges(7,3);
t176 = cos(qJ(3));
t175 = t140 * qJ(2) + pkin(8);
t174 = rSges(5,3) + qJ(4);
t172 = t135 * t137;
t143 = sin(qJ(3));
t171 = t136 * t143;
t168 = t137 * t147;
t166 = t139 * t143;
t165 = t140 * t147;
t164 = t144 * t135;
t162 = t147 * pkin(1) + qJ(2) * t169;
t159 = t136 * t176;
t158 = t139 * t176;
t157 = t137 * t159;
t120 = t138 * t165 - t164;
t112 = -t120 * t136 - t139 * t168;
t123 = t138 * t147 - t140 * t164;
t156 = t123 * pkin(2) + t113 * pkin(9) + t162;
t155 = pkin(2) * t172 + t119 * pkin(9) + t175;
t105 = t123 * t176 + (t122 * t139 + t136 * t169) * t143;
t154 = t105 * pkin(3) + t156;
t109 = t140 * t171 + (t176 * t135 + t138 * t166) * t137;
t153 = t109 * pkin(3) + t155;
t104 = -t122 * t158 + t123 * t143 - t144 * t157;
t152 = t113 * pkin(4) + t104 * qJ(4) + t154;
t108 = -t140 * t159 + t143 * t172 - t158 * t170;
t151 = t119 * pkin(4) + t108 * qJ(4) + t153;
t121 = t135 * t165 + t163;
t133 = t144 * pkin(1);
t150 = t121 * pkin(2) + t112 * pkin(9) - qJ(2) * t168 + t133;
t103 = t120 * t166 + t121 * t176 - t168 * t171;
t149 = t103 * pkin(3) + t150;
t102 = -t120 * t158 + t121 * t143 + t147 * t157;
t148 = t112 * pkin(4) + t102 * qJ(4) + t149;
t146 = cos(qJ(5));
t145 = cos(qJ(6));
t142 = sin(qJ(5));
t141 = sin(qJ(6));
t101 = t108 * t142 + t119 * t146;
t100 = -t108 * t146 + t119 * t142;
t95 = t104 * t142 + t113 * t146;
t94 = -t104 * t146 + t113 * t142;
t93 = t102 * t142 + t112 * t146;
t92 = -t102 * t146 + t112 * t142;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t147 - t144 * rSges(2,2)) + g(2) * (t144 * rSges(2,1) + rSges(2,2) * t147) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t123 + rSges(3,2) * t122 + t162) + g(2) * (t121 * rSges(3,1) + t120 * rSges(3,2) + t133) + g(3) * (rSges(3,3) * t140 + t175) + (g(1) * rSges(3,3) * t144 + g(3) * (rSges(3,1) * t135 + rSges(3,2) * t138) + g(2) * (-rSges(3,3) - qJ(2)) * t147) * t137) - m(4) * (g(1) * (rSges(4,1) * t105 - rSges(4,2) * t104 + rSges(4,3) * t113 + t156) + g(2) * (t103 * rSges(4,1) - t102 * rSges(4,2) + t112 * rSges(4,3) + t150) + g(3) * (rSges(4,1) * t109 - rSges(4,2) * t108 + rSges(4,3) * t119 + t155)) - m(5) * (g(1) * (rSges(5,1) * t113 - rSges(5,2) * t105 + t174 * t104 + t154) + g(2) * (t112 * rSges(5,1) - t103 * rSges(5,2) + t174 * t102 + t149) + g(3) * (rSges(5,1) * t119 - rSges(5,2) * t109 + t174 * t108 + t153)) - m(6) * (g(1) * (rSges(6,1) * t95 - rSges(6,2) * t94 + t178 * t105 + t152) + g(2) * (t93 * rSges(6,1) - t92 * rSges(6,2) + t178 * t103 + t148) + g(3) * (rSges(6,1) * t101 - rSges(6,2) * t100 + t178 * t109 + t151)) - m(7) * (g(1) * (t95 * pkin(5) + t105 * pkin(10) + (t105 * t141 + t145 * t95) * rSges(7,1) + (t105 * t145 - t141 * t95) * rSges(7,2) + t177 * t94 + t152) + g(2) * (t93 * pkin(5) + t103 * pkin(10) + (t103 * t141 + t145 * t93) * rSges(7,1) + (t103 * t145 - t141 * t93) * rSges(7,2) + t177 * t92 + t148) + g(3) * (t101 * pkin(5) + t109 * pkin(10) + (t101 * t145 + t109 * t141) * rSges(7,1) + (-t101 * t141 + t109 * t145) * rSges(7,2) + t177 * t100 + t151));
U  = t1;
