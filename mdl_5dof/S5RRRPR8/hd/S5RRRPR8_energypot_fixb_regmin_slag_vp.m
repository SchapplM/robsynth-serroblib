% Calculate minimal parameter regressor of potential energy for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:25
% EndTime: 2019-12-31 21:20:25
% DurationCPUTime: 0.05s
% Computational Cost: add. (62->30), mult. (76->39), div. (0->0), fcn. (77->8), ass. (0->21)
t155 = qJ(2) + qJ(3);
t154 = cos(t155);
t169 = g(3) * t154;
t156 = sin(qJ(5));
t158 = sin(qJ(1));
t168 = t158 * t156;
t159 = cos(qJ(5));
t167 = t158 * t159;
t161 = cos(qJ(1));
t166 = t161 * t156;
t165 = t161 * t159;
t164 = g(1) * t161 + g(2) * t158;
t153 = sin(t155);
t160 = cos(qJ(2));
t163 = t160 * pkin(2) + pkin(3) * t154 + qJ(4) * t153 + pkin(1);
t162 = -pkin(7) - pkin(6);
t157 = sin(qJ(2));
t151 = g(1) * t158 - g(2) * t161;
t150 = g(3) * t153 + t164 * t154;
t149 = t164 * t153 - t169;
t1 = [0, -t164, t151, 0, 0, 0, 0, 0, -g(3) * t157 - t164 * t160, -g(3) * t160 + t164 * t157, 0, 0, 0, 0, 0, -t150, t149, -t151, t150, -t149, -g(3) * (t157 * pkin(2) + t153 * pkin(3) - t154 * qJ(4) + pkin(5)) + (-g(1) * t163 - g(2) * t162) * t161 + (g(1) * t162 - g(2) * t163) * t158, 0, 0, 0, 0, 0, -g(1) * (t153 * t166 + t167) - g(2) * (t153 * t168 - t165) + t156 * t169, -g(1) * (t153 * t165 - t168) - g(2) * (t153 * t167 + t166) + t159 * t169;];
U_reg = t1;
