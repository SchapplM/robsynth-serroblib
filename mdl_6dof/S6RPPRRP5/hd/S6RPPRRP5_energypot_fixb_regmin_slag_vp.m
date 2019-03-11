% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:46
% EndTime: 2019-03-09 02:08:46
% DurationCPUTime: 0.06s
% Computational Cost: add. (55->41), mult. (91->48), div. (0->0), fcn. (86->6), ass. (0->25)
t150 = pkin(2) + pkin(6);
t138 = cos(qJ(4));
t149 = g(3) * t138;
t134 = sin(qJ(5));
t136 = sin(qJ(1));
t148 = t136 * t134;
t137 = cos(qJ(5));
t147 = t136 * t137;
t139 = cos(qJ(1));
t146 = t139 * t134;
t145 = t139 * t137;
t144 = t139 * pkin(1) + t136 * qJ(2);
t143 = -pkin(5) * t134 - pkin(7);
t130 = t136 * pkin(1);
t142 = -t139 * qJ(2) + t130;
t141 = g(1) * (t139 * qJ(3) + t144);
t125 = g(1) * t139 + g(2) * t136;
t126 = t137 * pkin(5) + pkin(4);
t133 = -qJ(6) - pkin(8);
t135 = sin(qJ(4));
t140 = t126 * t135 + t133 * t138;
t127 = t136 * qJ(3);
t124 = g(1) * t136 - g(2) * t139;
t123 = -g(3) * t135 + t125 * t138;
t1 = [0, -t125, t124, t125, -t124, -g(3) * pkin(6) - g(1) * t144 - g(2) * t142, -t124, -t125, -t141 - g(2) * (t127 + t142) - g(3) * t150, 0, 0, 0, 0, 0, -t125 * t135 - t149, -t123, 0, 0, 0, 0, 0, -g(1) * (t135 * t145 - t148) - g(2) * (t135 * t147 + t146) - t137 * t149, -g(1) * (-t135 * t146 - t147) - g(2) * (-t135 * t148 + t145) + t134 * t149, t123, -t141 - g(2) * (t127 + t130) - g(3) * (t138 * t126 - t135 * t133 + pkin(3) + t150) + (-g(1) * t143 - g(2) * t140) * t136 + (-g(1) * t140 - g(2) * (-qJ(2) - t143)) * t139;];
U_reg  = t1;
