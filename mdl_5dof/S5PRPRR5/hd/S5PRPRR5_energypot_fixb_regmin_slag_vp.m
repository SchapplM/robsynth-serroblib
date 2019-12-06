% Calculate minimal parameter regressor of potential energy for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:40
% EndTime: 2019-12-05 15:54:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (71->39), mult. (93->62), div. (0->0), fcn. (101->10), ass. (0->19)
t142 = sin(qJ(2));
t148 = g(3) * t142;
t139 = sin(pkin(8));
t143 = cos(qJ(2));
t147 = t139 * t143;
t141 = cos(pkin(8));
t146 = t141 * t143;
t137 = pkin(9) + qJ(4);
t145 = g(1) * t141 + g(2) * t139;
t144 = pkin(2) * t143 + qJ(3) * t142 + pkin(1);
t140 = cos(pkin(9));
t138 = sin(pkin(9));
t136 = qJ(5) + t137;
t135 = cos(t137);
t134 = sin(t137);
t133 = cos(t136);
t132 = sin(t136);
t131 = -g(3) * t143 + t145 * t142;
t1 = [-g(3) * qJ(1), 0, -t145 * t143 - t148, t131, -g(1) * (t139 * t138 + t140 * t146) - g(2) * (-t141 * t138 + t140 * t147) - t140 * t148, -g(1) * (-t138 * t146 + t139 * t140) - g(2) * (-t138 * t147 - t141 * t140) + t138 * t148, -t131, -g(3) * (t142 * pkin(2) - t143 * qJ(3) + qJ(1)) + (g(2) * pkin(5) - g(1) * t144) * t141 + (-g(1) * pkin(5) - g(2) * t144) * t139, 0, 0, 0, 0, 0, -g(1) * (t139 * t134 + t135 * t146) - g(2) * (-t141 * t134 + t135 * t147) - t135 * t148, -g(1) * (-t134 * t146 + t139 * t135) - g(2) * (-t134 * t147 - t141 * t135) + t134 * t148, 0, 0, 0, 0, 0, -g(1) * (t139 * t132 + t133 * t146) - g(2) * (-t141 * t132 + t133 * t147) - t133 * t148, -g(1) * (-t132 * t146 + t139 * t133) - g(2) * (-t132 * t147 - t141 * t133) + t132 * t148;];
U_reg = t1;
