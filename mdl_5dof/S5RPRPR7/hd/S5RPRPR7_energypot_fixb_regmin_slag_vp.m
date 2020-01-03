% Calculate minimal parameter regressor of potential energy for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:50
% EndTime: 2019-12-31 18:19:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (54->27), mult. (53->40), div. (0->0), fcn. (51->10), ass. (0->22)
t150 = qJ(3) + pkin(9);
t166 = g(3) * sin(t150);
t165 = qJ(2) + pkin(5);
t151 = qJ(1) + pkin(8);
t147 = sin(t151);
t153 = sin(qJ(5));
t164 = t147 * t153;
t156 = cos(qJ(5));
t163 = t147 * t156;
t149 = cos(t151);
t162 = t149 * t153;
t161 = t149 * t156;
t160 = g(1) * t149 + g(2) * t147;
t155 = sin(qJ(1));
t158 = cos(qJ(1));
t159 = -g(1) * t158 - g(2) * t155;
t157 = cos(qJ(3));
t154 = sin(qJ(3));
t152 = -qJ(4) - pkin(6);
t148 = cos(t150);
t145 = t157 * pkin(3) + pkin(2);
t1 = [0, t159, g(1) * t155 - g(2) * t158, t159 * pkin(1) - g(3) * t165, 0, 0, 0, 0, 0, -g(3) * t154 - t160 * t157, -g(3) * t157 + t160 * t154, -g(1) * t147 + g(2) * t149, -g(1) * (t158 * pkin(1) + t149 * t145 - t147 * t152) - g(2) * (t155 * pkin(1) + t147 * t145 + t149 * t152) - g(3) * (t154 * pkin(3) + t165), 0, 0, 0, 0, 0, -g(1) * (t148 * t161 + t164) - g(2) * (t148 * t163 - t162) - t156 * t166, -g(1) * (-t148 * t162 + t163) - g(2) * (-t148 * t164 - t161) + t153 * t166;];
U_reg = t1;
