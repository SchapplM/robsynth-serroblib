% Calculate minimal parameter regressor of potential energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:50
% EndTime: 2019-12-31 17:30:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (48->30), mult. (122->59), div. (0->0), fcn. (156->10), ass. (0->26)
t162 = sin(pkin(4));
t166 = sin(qJ(2));
t179 = t162 * t166;
t169 = cos(qJ(3));
t178 = t162 * t169;
t170 = cos(qJ(2));
t177 = t162 * t170;
t171 = cos(qJ(1));
t176 = t162 * t171;
t167 = sin(qJ(1));
t175 = t167 * t166;
t174 = t167 * t170;
t173 = t171 * t166;
t172 = t171 * t170;
t168 = cos(qJ(4));
t165 = sin(qJ(3));
t164 = sin(qJ(4));
t163 = cos(pkin(4));
t161 = -t163 * t175 + t172;
t160 = t163 * t174 + t173;
t159 = t163 * t173 + t174;
t158 = -t163 * t172 + t175;
t157 = t163 * t165 + t166 * t178;
t156 = t167 * t162 * t165 + t161 * t169;
t155 = t159 * t169 - t165 * t176;
t1 = [0, -g(1) * t171 - g(2) * t167, g(1) * t167 - g(2) * t171, 0, 0, 0, 0, 0, -g(1) * t161 - g(2) * t159 - g(3) * t179, g(1) * t160 + g(2) * t158 - g(3) * t177, 0, 0, 0, 0, 0, -g(1) * t156 - g(2) * t155 - g(3) * t157, -g(1) * (-t161 * t165 + t167 * t178) - g(2) * (-t159 * t165 - t169 * t176) - g(3) * (t163 * t169 - t165 * t179), 0, 0, 0, 0, 0, -g(1) * (t156 * t168 + t160 * t164) - g(2) * (t155 * t168 + t158 * t164) - g(3) * (t157 * t168 - t164 * t177), -g(1) * (-t156 * t164 + t160 * t168) - g(2) * (-t155 * t164 + t158 * t168) - g(3) * (-t157 * t164 - t168 * t177);];
U_reg = t1;
