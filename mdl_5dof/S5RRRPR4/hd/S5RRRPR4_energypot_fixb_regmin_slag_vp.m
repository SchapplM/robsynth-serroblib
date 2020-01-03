% Calculate minimal parameter regressor of potential energy for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->24), mult. (75->32), div. (0->0), fcn. (78->8), ass. (0->18)
t135 = qJ(1) + qJ(2);
t133 = sin(t135);
t134 = cos(t135);
t146 = g(1) * t134 + g(2) * t133;
t138 = sin(qJ(1));
t141 = cos(qJ(1));
t145 = -g(1) * t141 - g(2) * t138;
t136 = sin(qJ(5));
t137 = sin(qJ(3));
t139 = cos(qJ(5));
t140 = cos(qJ(3));
t144 = t140 * t136 - t137 * t139;
t143 = t137 * t136 + t140 * t139;
t142 = pkin(3) * t140 + qJ(4) * t137 + pkin(2);
t132 = g(1) * t133 - g(2) * t134;
t131 = -g(3) * t137 - t146 * t140;
t130 = -g(3) * t140 + t146 * t137;
t1 = [0, t145, g(1) * t138 - g(2) * t141, 0, -t146, t132, 0, 0, 0, 0, 0, t131, t130, t131, -t132, -t130, -g(3) * (t137 * pkin(3) - t140 * qJ(4) + pkin(5) + pkin(6)) + t145 * pkin(1) + (g(2) * pkin(7) - g(1) * t142) * t134 + (-g(1) * pkin(7) - g(2) * t142) * t133, 0, 0, 0, 0, 0, g(3) * t144 - t146 * t143, g(3) * t143 + t146 * t144;];
U_reg = t1;
