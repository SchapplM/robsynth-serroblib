% Calculate minimal parameter regressor of potential energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:28
% EndTime: 2019-12-05 18:56:28
% DurationCPUTime: 0.05s
% Computational Cost: add. (54->25), mult. (64->40), div. (0->0), fcn. (72->10), ass. (0->23)
t141 = qJ(2) + qJ(3);
t137 = sin(t141);
t157 = g(3) * t137;
t140 = qJ(4) + qJ(5);
t136 = sin(t140);
t144 = sin(qJ(1));
t156 = t144 * t136;
t138 = cos(t140);
t155 = t144 * t138;
t142 = sin(qJ(4));
t154 = t144 * t142;
t145 = cos(qJ(4));
t153 = t144 * t145;
t147 = cos(qJ(1));
t152 = t147 * t136;
t151 = t147 * t138;
t150 = t147 * t142;
t149 = t147 * t145;
t148 = g(1) * t147 + g(2) * t144;
t146 = cos(qJ(2));
t143 = sin(qJ(2));
t139 = cos(t141);
t1 = [0, -t148, g(1) * t144 - g(2) * t147, 0, 0, 0, 0, 0, -g(3) * t143 - t148 * t146, -g(3) * t146 + t148 * t143, 0, 0, 0, 0, 0, -t148 * t139 - t157, -g(3) * t139 + t148 * t137, 0, 0, 0, 0, 0, -g(1) * (t139 * t149 + t154) - g(2) * (t139 * t153 - t150) - t145 * t157, -g(1) * (-t139 * t150 + t153) - g(2) * (-t139 * t154 - t149) + t142 * t157, 0, 0, 0, 0, 0, -g(1) * (t139 * t151 + t156) - g(2) * (t139 * t155 - t152) - t138 * t157, -g(1) * (-t139 * t152 + t155) - g(2) * (-t139 * t156 - t151) + t136 * t157;];
U_reg = t1;
