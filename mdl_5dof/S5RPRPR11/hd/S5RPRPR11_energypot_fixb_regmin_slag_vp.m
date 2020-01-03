% Calculate minimal parameter regressor of potential energy for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:54
% EndTime: 2019-12-31 18:27:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (75->30), mult. (91->40), div. (0->0), fcn. (91->8), ass. (0->18)
t159 = sin(qJ(1));
t161 = cos(qJ(1));
t165 = g(1) * t161 + g(2) * t159;
t154 = pkin(8) + qJ(3);
t151 = sin(t154);
t152 = cos(t154);
t158 = sin(qJ(5));
t160 = cos(qJ(5));
t164 = t151 * t160 - t152 * t158;
t163 = t151 * t158 + t152 * t160;
t156 = cos(pkin(8));
t162 = t156 * pkin(2) + pkin(3) * t152 + qJ(4) * t151 + pkin(1);
t157 = -pkin(6) - qJ(2);
t155 = sin(pkin(8));
t149 = g(1) * t159 - g(2) * t161;
t148 = -g(3) * t151 - t165 * t152;
t147 = -g(3) * t152 + t165 * t151;
t1 = [0, -t165, t149, -g(3) * t155 - t165 * t156, -g(3) * t156 + t165 * t155, -t149, -g(1) * (t161 * pkin(1) + t159 * qJ(2)) - g(2) * (t159 * pkin(1) - t161 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, t148, t147, t148, -t149, -t147, -g(3) * (t155 * pkin(2) + t151 * pkin(3) - t152 * qJ(4) + pkin(5)) + (-g(1) * t162 - g(2) * t157) * t161 + (g(1) * t157 - g(2) * t162) * t159, 0, 0, 0, 0, 0, -g(3) * t164 - t165 * t163, g(3) * t163 - t165 * t164;];
U_reg = t1;
