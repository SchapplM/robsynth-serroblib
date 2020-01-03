% Calculate minimal parameter regressor of potential energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (81->28), mult. (67->38), div. (0->0), fcn. (58->8), ass. (0->21)
t149 = qJ(2) + pkin(5);
t139 = sin(qJ(3));
t148 = t139 * pkin(3) + t149;
t141 = cos(qJ(3));
t128 = t141 * pkin(3) + pkin(2);
t137 = qJ(1) + pkin(7);
t130 = sin(t137);
t132 = cos(t137);
t138 = -qJ(4) - pkin(6);
t140 = sin(qJ(1));
t147 = t140 * pkin(1) + t130 * t128 + t132 * t138;
t146 = g(1) * t132 + g(2) * t130;
t142 = cos(qJ(1));
t145 = -g(1) * t142 - g(2) * t140;
t144 = t142 * pkin(1) + t132 * t128 - t130 * t138;
t136 = qJ(3) + pkin(8);
t129 = sin(t136);
t131 = cos(t136);
t143 = pkin(4) * t131 + qJ(5) * t129;
t124 = -g(1) * t130 + g(2) * t132;
t1 = [0, t145, g(1) * t140 - g(2) * t142, t145 * pkin(1) - g(3) * t149, 0, 0, 0, 0, 0, -g(3) * t139 - t146 * t141, -g(3) * t141 + t146 * t139, t124, -g(1) * t144 - g(2) * t147 - g(3) * t148, -g(3) * t129 - t146 * t131, t124, g(3) * t131 - t146 * t129, -g(1) * (t143 * t132 + t144) - g(2) * (t143 * t130 + t147) - g(3) * (t129 * pkin(4) - t131 * qJ(5) + t148);];
U_reg = t1;
