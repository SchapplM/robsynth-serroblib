% Calculate minimal parameter regressor of potential energy for
% S5RPRRR7
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:09
% EndTime: 2019-12-31 19:04:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (55->25), mult. (59->44), div. (0->0), fcn. (64->10), ass. (0->20)
t158 = sin(qJ(3));
t169 = g(3) * t158;
t156 = qJ(4) + qJ(5);
t153 = sin(t156);
t161 = cos(qJ(3));
t168 = t153 * t161;
t154 = cos(t156);
t167 = t154 * t161;
t157 = sin(qJ(4));
t166 = t157 * t161;
t160 = cos(qJ(4));
t165 = t160 * t161;
t155 = qJ(1) + pkin(9);
t151 = sin(t155);
t152 = cos(t155);
t164 = g(1) * t152 + g(2) * t151;
t159 = sin(qJ(1));
t162 = cos(qJ(1));
t163 = -g(1) * t162 - g(2) * t159;
t1 = [0, t163, g(1) * t159 - g(2) * t162, -g(3) * (qJ(2) + pkin(5)) + t163 * pkin(1), 0, 0, 0, 0, 0, -t164 * t161 - t169, -g(3) * t161 + t164 * t158, 0, 0, 0, 0, 0, -g(1) * (t151 * t157 + t152 * t165) - g(2) * (t151 * t165 - t152 * t157) - t160 * t169, -g(1) * (t151 * t160 - t152 * t166) - g(2) * (-t151 * t166 - t152 * t160) + t157 * t169, 0, 0, 0, 0, 0, -g(1) * (t151 * t153 + t152 * t167) - g(2) * (t151 * t167 - t152 * t153) - t154 * t169, -g(1) * (t151 * t154 - t152 * t168) - g(2) * (-t151 * t168 - t152 * t154) + t153 * t169;];
U_reg = t1;
