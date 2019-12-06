% Calculate minimal parameter regressor of potential energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:35
% EndTime: 2019-12-05 19:00:35
% DurationCPUTime: 0.03s
% Computational Cost: add. (50->13), mult. (38->20), div. (0->0), fcn. (38->10), ass. (0->15)
t140 = qJ(3) + qJ(4);
t141 = qJ(1) + qJ(2);
t136 = sin(t141);
t138 = cos(t141);
t146 = g(2) * t136 - g(3) * t138;
t145 = cos(qJ(1));
t144 = cos(qJ(3));
t143 = sin(qJ(1));
t142 = sin(qJ(3));
t139 = qJ(5) + t140;
t137 = cos(t140);
t135 = sin(t140);
t134 = cos(t139);
t133 = sin(t139);
t1 = [0, g(2) * t143 - g(3) * t145, g(2) * t145 + g(3) * t143, 0, t146, g(2) * t138 + g(3) * t136, 0, 0, 0, 0, 0, -g(1) * t142 + t146 * t144, -g(1) * t144 - t146 * t142, 0, 0, 0, 0, 0, -g(1) * t135 + t146 * t137, -g(1) * t137 - t146 * t135, 0, 0, 0, 0, 0, -g(1) * t133 + t146 * t134, -g(1) * t134 - t146 * t133;];
U_reg = t1;
