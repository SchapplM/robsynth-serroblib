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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:08:21
% EndTime: 2022-01-20 12:08:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->14), mult. (38->20), div. (0->0), fcn. (38->10), ass. (0->15)
t141 = qJ(3) + qJ(4);
t142 = qJ(1) + qJ(2);
t137 = sin(t142);
t139 = cos(t142);
t147 = g(1) * t139 + g(2) * t137;
t146 = cos(qJ(1));
t145 = cos(qJ(3));
t144 = sin(qJ(1));
t143 = sin(qJ(3));
t140 = qJ(5) + t141;
t138 = cos(t141);
t136 = sin(t141);
t135 = cos(t140);
t134 = sin(t140);
t1 = [0, -g(1) * t146 - g(2) * t144, g(1) * t144 - g(2) * t146, 0, -t147, g(1) * t137 - g(2) * t139, 0, 0, 0, 0, 0, -g(3) * t143 - t147 * t145, -g(3) * t145 + t147 * t143, 0, 0, 0, 0, 0, -g(3) * t136 - t147 * t138, -g(3) * t138 + t147 * t136, 0, 0, 0, 0, 0, -g(3) * t134 - t147 * t135, -g(3) * t135 + t147 * t134;];
U_reg = t1;
