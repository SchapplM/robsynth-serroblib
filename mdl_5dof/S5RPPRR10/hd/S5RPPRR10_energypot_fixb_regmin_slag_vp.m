% Calculate minimal parameter regressor of potential energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:17
% DurationCPUTime: 0.07s
% Computational Cost: add. (58->28), mult. (104->40), div. (0->0), fcn. (110->8), ass. (0->21)
t162 = sin(qJ(1));
t164 = cos(qJ(1));
t170 = g(1) * t164 + g(2) * t162;
t172 = t164 * pkin(1) + t162 * qJ(2);
t171 = t162 * pkin(1) - t164 * qJ(2);
t159 = sin(pkin(8));
t160 = cos(pkin(8));
t169 = pkin(2) * t160 + qJ(3) * t159;
t158 = qJ(4) + qJ(5);
t152 = sin(t158);
t153 = cos(t158);
t168 = t160 * t152 - t159 * t153;
t167 = t159 * t152 + t160 * t153;
t161 = sin(qJ(4));
t163 = cos(qJ(4));
t166 = t159 * t163 - t160 * t161;
t165 = t159 * t161 + t160 * t163;
t151 = g(1) * t162 - g(2) * t164;
t150 = -g(3) * t159 - t170 * t160;
t149 = -g(3) * t160 + t170 * t159;
t1 = [0, -t170, t151, t150, t149, -t151, -g(3) * pkin(5) - g(1) * t172 - g(2) * t171, t150, -t151, -t149, -g(1) * (t169 * t164 + t172) - g(2) * (t169 * t162 + t171) - g(3) * (t159 * pkin(2) - t160 * qJ(3) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t166 - t170 * t165, g(3) * t165 - t170 * t166, 0, 0, 0, 0, 0, g(3) * t168 - t170 * t167, g(3) * t167 + t170 * t168;];
U_reg = t1;
