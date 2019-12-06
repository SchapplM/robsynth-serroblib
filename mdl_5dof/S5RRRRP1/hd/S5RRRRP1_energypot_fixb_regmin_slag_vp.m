% Calculate minimal parameter regressor of potential energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:59
% EndTime: 2019-12-05 18:45:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (61->23), mult. (52->29), div. (0->0), fcn. (49->8), ass. (0->15)
t147 = qJ(2) + qJ(3);
t149 = sin(qJ(1));
t151 = cos(qJ(1));
t152 = g(1) * t151 + g(2) * t149;
t150 = cos(qJ(2));
t148 = sin(qJ(2));
t146 = qJ(4) + t147;
t145 = -qJ(5) - pkin(8) - pkin(7) - pkin(6);
t144 = cos(t147);
t143 = sin(t147);
t142 = cos(t146);
t141 = sin(t146);
t140 = g(1) * t149 - g(2) * t151;
t139 = t150 * pkin(2) + pkin(3) * t144 + pkin(4) * t142 + pkin(1);
t1 = [0, -t152, t140, 0, 0, 0, 0, 0, -g(3) * t148 - t152 * t150, -g(3) * t150 + t152 * t148, 0, 0, 0, 0, 0, -g(3) * t143 - t152 * t144, -g(3) * t144 + t152 * t143, 0, 0, 0, 0, 0, -g(3) * t141 - t152 * t142, -g(3) * t142 + t152 * t141, -t140, -g(1) * (t151 * t139 - t149 * t145) - g(2) * (t149 * t139 + t151 * t145) - g(3) * (t148 * pkin(2) + pkin(3) * t143 + pkin(4) * t141 + pkin(5));];
U_reg = t1;
