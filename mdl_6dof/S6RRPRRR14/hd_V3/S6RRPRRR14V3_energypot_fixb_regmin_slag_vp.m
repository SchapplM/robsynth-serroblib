% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR14V3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energypot_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:05
% EndTime: 2019-04-12 15:10:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (59->36), mult. (152->64), div. (0->0), fcn. (183->10), ass. (0->29)
t140 = sin(qJ(2));
t145 = cos(qJ(2));
t141 = sin(qJ(1));
t146 = cos(qJ(1));
t147 = g(1) * t146 + g(2) * t141;
t129 = -g(3) * t145 + t147 * t140;
t139 = sin(qJ(4));
t154 = t139 * t140;
t153 = t140 * t141;
t144 = cos(qJ(4));
t152 = t140 * t144;
t151 = t140 * t146;
t150 = t141 * t145;
t149 = t146 * t139;
t148 = t146 * t144;
t143 = cos(qJ(5));
t142 = cos(qJ(6));
t138 = sin(qJ(5));
t137 = sin(qJ(6));
t136 = g(1) * t141 - g(2) * t146;
t135 = t141 * t139 + t145 * t148;
t134 = -t141 * t144 + t145 * t149;
t133 = t144 * t150 - t149;
t132 = t139 * t150 + t148;
t131 = -t145 * t138 + t143 * t152;
t130 = -g(3) * t140 - t147 * t145;
t128 = t135 * t143 + t138 * t151;
t127 = t133 * t143 + t138 * t153;
t1 = [0, -t147, t136, 0, 0, 0, 0, 0, t130, t129, t130, -t136, -t129, -t129 * qJ(3), 0, 0, 0, 0, 0, -g(1) * t135 - g(2) * t133 - g(3) * t152, g(1) * t134 + g(2) * t132 + g(3) * t154, 0, 0, 0, 0, 0, -g(1) * t128 - g(2) * t127 - g(3) * t131, -g(1) * (-t135 * t138 + t143 * t151) - g(2) * (-t133 * t138 + t143 * t153) - g(3) * (-t138 * t152 - t145 * t143) 0, 0, 0, 0, 0, -g(1) * (t128 * t142 + t134 * t137) - g(2) * (t127 * t142 + t132 * t137) - g(3) * (t131 * t142 + t137 * t154) -g(1) * (-t128 * t137 + t134 * t142) - g(2) * (-t127 * t137 + t132 * t142) - g(3) * (-t131 * t137 + t142 * t154);];
U_reg  = t1;
