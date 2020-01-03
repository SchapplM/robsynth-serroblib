% Calculate minimal parameter regressor of potential energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:28
% EndTime: 2020-01-03 12:09:28
% DurationCPUTime: 0.03s
% Computational Cost: add. (55->22), mult. (42->27), div. (0->0), fcn. (39->8), ass. (0->15)
t141 = qJ(1) + qJ(2);
t139 = sin(t141);
t140 = cos(t141);
t147 = g(2) * t139 - g(3) * t140;
t146 = cos(qJ(1));
t145 = cos(qJ(3));
t144 = sin(qJ(1));
t143 = sin(qJ(3));
t142 = -qJ(4) - pkin(7);
t138 = qJ(3) + pkin(9) + qJ(5);
t137 = t145 * pkin(3) + pkin(2);
t136 = cos(t138);
t135 = sin(t138);
t134 = g(2) * t140 + g(3) * t139;
t1 = [0, -g(2) * t144 + g(3) * t146, -g(2) * t146 - g(3) * t144, 0, -t147, -t134, 0, 0, 0, 0, 0, -g(1) * t143 - t147 * t145, -g(1) * t145 + t147 * t143, t134, -g(1) * (t143 * pkin(3) + pkin(5) + pkin(6)) - g(2) * (t144 * pkin(1) + t139 * t137 + t140 * t142) - g(3) * (-t146 * pkin(1) - t140 * t137 + t139 * t142), 0, 0, 0, 0, 0, -g(1) * t135 - t147 * t136, -g(1) * t136 + t147 * t135;];
U_reg = t1;
