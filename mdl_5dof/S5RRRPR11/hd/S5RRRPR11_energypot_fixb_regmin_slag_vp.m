% Calculate minimal parameter regressor of potential energy for
% S5RRRPR11
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:13
% EndTime: 2019-12-31 21:35:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (62->40), mult. (142->60), div. (0->0), fcn. (165->8), ass. (0->21)
t169 = sin(qJ(2));
t180 = g(3) * t169;
t170 = sin(qJ(1));
t173 = cos(qJ(2));
t179 = t170 * t173;
t168 = sin(qJ(3));
t174 = cos(qJ(1));
t178 = t174 * t168;
t172 = cos(qJ(3));
t177 = t174 * t172;
t176 = g(1) * t174 + g(2) * t170;
t162 = t168 * t179 + t177;
t164 = -t170 * t172 + t173 * t178;
t175 = g(1) * t164 + g(2) * t162 + t168 * t180;
t171 = cos(qJ(5));
t167 = sin(qJ(5));
t165 = t170 * t168 + t173 * t177;
t163 = t172 * t179 - t178;
t161 = -g(3) * t173 + t176 * t169;
t160 = -g(1) * t165 - g(2) * t163 - t172 * t180;
t1 = [0, -t176, g(1) * t170 - g(2) * t174, 0, 0, 0, 0, 0, -t176 * t173 - t180, t161, 0, 0, 0, 0, 0, t160, t175, t160, -t161, -t175, -g(1) * (t165 * pkin(3) + t170 * pkin(6) + t164 * qJ(4) + (pkin(2) * t173 + pkin(1)) * t174) - g(2) * (t170 * pkin(1) + pkin(2) * t179 + t163 * pkin(3) - t174 * pkin(6) + t162 * qJ(4)) - g(3) * (-t173 * pkin(7) + pkin(5)) + (-g(3) * (pkin(3) * t172 + qJ(4) * t168 + pkin(2)) - t176 * pkin(7)) * t169, 0, 0, 0, 0, 0, -g(1) * (t164 * t167 + t165 * t171) - g(2) * (t162 * t167 + t163 * t171) - (t167 * t168 + t171 * t172) * t180, -g(1) * (t164 * t171 - t165 * t167) - g(2) * (t162 * t171 - t163 * t167) - (-t167 * t172 + t168 * t171) * t180;];
U_reg = t1;
