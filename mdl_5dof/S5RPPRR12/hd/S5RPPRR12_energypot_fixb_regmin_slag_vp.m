% Calculate minimal parameter regressor of potential energy for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (44->27), mult. (66->38), div. (0->0), fcn. (64->8), ass. (0->19)
t129 = pkin(8) + qJ(4);
t125 = cos(t129);
t142 = g(3) * t125;
t132 = sin(qJ(5));
t133 = sin(qJ(1));
t141 = t133 * t132;
t134 = cos(qJ(5));
t140 = t133 * t134;
t135 = cos(qJ(1));
t139 = t135 * t132;
t138 = t135 * t134;
t137 = t135 * pkin(1) + t133 * qJ(2);
t136 = t133 * pkin(1) - t135 * qJ(2);
t122 = g(1) * t133 - g(2) * t135;
t131 = cos(pkin(8));
t130 = sin(pkin(8));
t124 = sin(t129);
t123 = g(1) * t135 + g(2) * t133;
t1 = [0, -t123, t122, t123, -t122, -g(3) * pkin(5) - g(1) * t137 - g(2) * t136, -g(3) * t131 - t122 * t130, g(3) * t130 - t122 * t131, -t123, -g(1) * (t135 * qJ(3) + t137) - g(2) * (t133 * qJ(3) + t136) - g(3) * (pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -t122 * t124 - t142, g(3) * t124 - t122 * t125, 0, 0, 0, 0, 0, -g(1) * (t124 * t140 + t139) - g(2) * (-t124 * t138 + t141) - t134 * t142, -g(1) * (-t124 * t141 + t138) - g(2) * (t124 * t139 + t140) + t132 * t142;];
U_reg = t1;
