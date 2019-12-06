% Calculate minimal parameter regressor of potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:32
% EndTime: 2019-12-05 18:25:32
% DurationCPUTime: 0.04s
% Computational Cost: add. (35->21), mult. (56->34), div. (0->0), fcn. (57->8), ass. (0->19)
t134 = cos(qJ(2));
t143 = pkin(1) * t134;
t129 = qJ(2) + qJ(4);
t127 = sin(t129);
t142 = g(3) * t127;
t131 = sin(qJ(2));
t141 = g(3) * t131;
t130 = sin(qJ(5));
t132 = sin(qJ(1));
t140 = t132 * t130;
t133 = cos(qJ(5));
t139 = t132 * t133;
t135 = cos(qJ(1));
t138 = t135 * t130;
t137 = t135 * t133;
t136 = g(1) * t135 + g(2) * t132;
t128 = cos(t129);
t126 = g(1) * t132 - g(2) * t135;
t1 = [0, -t136, t126, 0, 0, 0, 0, 0, -t136 * t134 - t141, -g(3) * t134 + t136 * t131, -t126, -g(1) * (t132 * qJ(3) + t135 * t143) - g(2) * (-t135 * qJ(3) + t132 * t143) - pkin(1) * t141, 0, 0, 0, 0, 0, -t136 * t128 - t142, -g(3) * t128 + t136 * t127, 0, 0, 0, 0, 0, -g(1) * (t128 * t137 + t140) - g(2) * (t128 * t139 - t138) - t133 * t142, -g(1) * (-t128 * t138 + t139) - g(2) * (-t128 * t140 - t137) + t130 * t142;];
U_reg = t1;
