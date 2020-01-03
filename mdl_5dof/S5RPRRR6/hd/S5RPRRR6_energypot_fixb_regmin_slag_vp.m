% Calculate minimal parameter regressor of potential energy for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:39
% EndTime: 2019-12-31 19:01:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (49->19), mult. (49->32), div. (0->0), fcn. (50->10), ass. (0->18)
t141 = qJ(3) + qJ(4);
t138 = sin(t141);
t152 = g(3) * t138;
t139 = cos(t141);
t142 = sin(qJ(5));
t151 = t139 * t142;
t145 = cos(qJ(5));
t150 = t139 * t145;
t140 = qJ(1) + pkin(9);
t136 = sin(t140);
t137 = cos(t140);
t149 = g(1) * t137 + g(2) * t136;
t144 = sin(qJ(1));
t147 = cos(qJ(1));
t148 = -g(1) * t147 - g(2) * t144;
t146 = cos(qJ(3));
t143 = sin(qJ(3));
t1 = [0, t148, g(1) * t144 - g(2) * t147, -g(3) * (qJ(2) + pkin(5)) + t148 * pkin(1), 0, 0, 0, 0, 0, -g(3) * t143 - t149 * t146, -g(3) * t146 + t149 * t143, 0, 0, 0, 0, 0, -t149 * t139 - t152, -g(3) * t139 + t149 * t138, 0, 0, 0, 0, 0, -g(1) * (t136 * t142 + t137 * t150) - g(2) * (t136 * t150 - t137 * t142) - t145 * t152, -g(1) * (t136 * t145 - t137 * t151) - g(2) * (-t136 * t151 - t137 * t145) + t142 * t152;];
U_reg = t1;
