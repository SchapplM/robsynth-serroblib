% Calculate minimal parameter regressor of potential energy for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:57
% EndTime: 2020-01-03 12:03:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (64->22), mult. (49->29), div. (0->0), fcn. (46->10), ass. (0->16)
t138 = pkin(9) + qJ(4);
t139 = qJ(1) + qJ(2);
t136 = sin(t139);
t137 = cos(t139);
t144 = g(2) * t136 - g(3) * t137;
t143 = cos(qJ(1));
t142 = sin(qJ(1));
t141 = cos(pkin(9));
t140 = sin(pkin(9));
t135 = qJ(5) + t138;
t134 = cos(t138);
t133 = sin(t138);
t132 = cos(t135);
t131 = sin(t135);
t130 = g(2) * t137 + g(3) * t136;
t1 = [0, -g(2) * t142 + g(3) * t143, -g(2) * t143 - g(3) * t142, 0, -t144, -t130, -g(1) * t140 - t144 * t141, -g(1) * t141 + t144 * t140, t130, -g(1) * (pkin(6) + pkin(5)) - g(2) * (t142 * pkin(1) + t136 * pkin(2) - t137 * qJ(3)) - g(3) * (-t143 * pkin(1) - t137 * pkin(2) - t136 * qJ(3)), 0, 0, 0, 0, 0, -g(1) * t133 - t144 * t134, -g(1) * t134 + t144 * t133, 0, 0, 0, 0, 0, -g(1) * t131 - t144 * t132, -g(1) * t132 + t144 * t131;];
U_reg = t1;
