% Calculate minimal parameter regressor of potential energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:36
% EndTime: 2019-12-31 20:30:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->28), mult. (105->44), div. (0->0), fcn. (118->8), ass. (0->19)
t169 = sin(qJ(1));
t173 = cos(qJ(1));
t177 = g(1) * t173 + g(2) * t169;
t167 = sin(qJ(4));
t168 = sin(qJ(2));
t171 = cos(qJ(4));
t172 = cos(qJ(2));
t176 = t172 * t167 - t168 * t171;
t178 = g(3) * t176;
t175 = t168 * t167 + t172 * t171;
t174 = pkin(2) * t172 + qJ(3) * t168 + pkin(1);
t170 = cos(qJ(5));
t166 = sin(qJ(5));
t165 = g(1) * t169 - g(2) * t173;
t163 = t175 * t173;
t162 = t175 * t169;
t161 = -g(3) * t168 - t177 * t172;
t160 = -g(3) * t172 + t177 * t168;
t1 = [0, -t177, t165, 0, 0, 0, 0, 0, t161, t160, t161, -t165, -t160, -g(3) * (t168 * pkin(2) - t172 * qJ(3) + pkin(5)) + (g(2) * pkin(6) - g(1) * t174) * t173 + (-g(1) * pkin(6) - g(2) * t174) * t169, 0, 0, 0, 0, 0, -g(1) * t163 - g(2) * t162 + t178, g(3) * t175 + t177 * t176, 0, 0, 0, 0, 0, -g(1) * (t163 * t170 - t169 * t166) - g(2) * (t162 * t170 + t173 * t166) + t170 * t178, -g(1) * (-t163 * t166 - t169 * t170) - g(2) * (-t162 * t166 + t173 * t170) - t166 * t178;];
U_reg = t1;
