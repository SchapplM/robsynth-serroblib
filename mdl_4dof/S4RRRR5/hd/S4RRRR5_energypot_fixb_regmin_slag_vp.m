% Calculate minimal parameter regressor of potential energy for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (32->22), mult. (54->37), div. (0->0), fcn. (62->8), ass. (0->17)
t129 = sin(qJ(2));
t140 = g(3) * t129;
t130 = sin(qJ(1));
t132 = cos(qJ(2));
t139 = t130 * t132;
t127 = qJ(3) + qJ(4);
t125 = sin(t127);
t133 = cos(qJ(1));
t138 = t133 * t125;
t126 = cos(t127);
t137 = t133 * t126;
t128 = sin(qJ(3));
t136 = t133 * t128;
t131 = cos(qJ(3));
t135 = t133 * t131;
t134 = g(1) * t133 + g(2) * t130;
t1 = [0, -t134, g(1) * t130 - g(2) * t133, 0, 0, 0, 0, 0, -t132 * t134 - t140, -g(3) * t132 + t129 * t134, 0, 0, 0, 0, 0, -g(1) * (t130 * t128 + t132 * t135) - g(2) * (t131 * t139 - t136) - t131 * t140, -g(1) * (t130 * t131 - t132 * t136) - g(2) * (-t128 * t139 - t135) + t128 * t140, 0, 0, 0, 0, 0, -g(1) * (t130 * t125 + t132 * t137) - g(2) * (t126 * t139 - t138) - t126 * t140, -g(1) * (t130 * t126 - t132 * t138) - g(2) * (-t125 * t139 - t137) + t125 * t140;];
U_reg = t1;
