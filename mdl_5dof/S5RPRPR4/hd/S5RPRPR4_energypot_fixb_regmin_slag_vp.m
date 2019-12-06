% Calculate minimal parameter regressor of potential energy for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:56
% EndTime: 2019-12-05 17:53:56
% DurationCPUTime: 0.07s
% Computational Cost: add. (52->21), mult. (43->29), div. (0->0), fcn. (37->8), ass. (0->16)
t149 = qJ(2) + pkin(5);
t141 = qJ(1) + pkin(8);
t138 = sin(t141);
t139 = cos(t141);
t148 = g(2) * t138 - g(3) * t139;
t144 = sin(qJ(1));
t146 = cos(qJ(1));
t147 = g(2) * t144 - g(3) * t146;
t145 = cos(qJ(3));
t143 = sin(qJ(3));
t142 = -qJ(4) - pkin(6);
t140 = qJ(3) + pkin(9) + qJ(5);
t137 = t145 * pkin(3) + pkin(2);
t136 = cos(t140);
t135 = sin(t140);
t1 = [0, t147, g(2) * t146 + g(3) * t144, t147 * pkin(1) - g(1) * t149, 0, 0, 0, 0, 0, -g(1) * t143 + t148 * t145, -g(1) * t145 - t148 * t143, -g(2) * t139 - g(3) * t138, -g(1) * (t143 * pkin(3) + t149) - g(2) * (-t144 * pkin(1) - t138 * t137 - t139 * t142) - g(3) * (t146 * pkin(1) + t139 * t137 - t138 * t142), 0, 0, 0, 0, 0, -g(1) * t135 + t148 * t136, -g(1) * t136 - t148 * t135;];
U_reg = t1;
