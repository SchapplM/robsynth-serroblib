% Calculate minimal parameter regressor of potential energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:53
% EndTime: 2019-12-31 18:12:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->39), mult. (104->43), div. (0->0), fcn. (95->6), ass. (0->19)
t133 = pkin(7) + qJ(3);
t129 = sin(t133);
t135 = cos(pkin(7));
t146 = t135 * pkin(2) + qJ(4) * t129 + pkin(1);
t130 = cos(t133);
t137 = sin(qJ(1));
t144 = t130 * t137;
t138 = cos(qJ(1));
t143 = t130 * t138;
t142 = pkin(3) * t143 + t146 * t138;
t136 = -pkin(6) - qJ(2);
t141 = pkin(3) * t144 + t138 * t136 + t146 * t137;
t140 = g(1) * t138 + g(2) * t137;
t134 = sin(pkin(7));
t139 = t134 * pkin(2) + t129 * pkin(3) - t130 * qJ(4) + pkin(5);
t122 = g(1) * t137 - g(2) * t138;
t117 = g(3) * t129 + t140 * t130;
t116 = -g(3) * t130 + t140 * t129;
t1 = [0, -t140, t122, -g(3) * t134 - t140 * t135, -g(3) * t135 + t140 * t134, -t122, -g(1) * (t138 * pkin(1) + t137 * qJ(2)) - g(2) * (t137 * pkin(1) - t138 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t117, t116, -t122, t117, -t116, -g(1) * (-t137 * t136 + t142) - g(2) * t141 - g(3) * t139, -t122, -t116, -t117, -g(1) * (qJ(5) * t143 + (pkin(4) - t136) * t137 + t142) - g(2) * (-t138 * pkin(4) + qJ(5) * t144 + t141) - g(3) * (t129 * qJ(5) + t139);];
U_reg = t1;
