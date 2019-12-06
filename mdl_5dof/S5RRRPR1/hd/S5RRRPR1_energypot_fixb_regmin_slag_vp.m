% Calculate minimal parameter regressor of potential energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:50
% EndTime: 2019-12-05 18:38:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (56->21), mult. (49->27), div. (0->0), fcn. (46->8), ass. (0->15)
t164 = qJ(2) + qJ(3);
t166 = sin(qJ(1));
t168 = cos(qJ(1));
t169 = g(1) * t168 + g(2) * t166;
t167 = cos(qJ(2));
t165 = sin(qJ(2));
t163 = -qJ(4) - pkin(7) - pkin(6);
t162 = cos(t164);
t161 = sin(t164);
t160 = pkin(9) + qJ(5) + t164;
t159 = cos(t160);
t158 = sin(t160);
t157 = g(1) * t166 - g(2) * t168;
t156 = t167 * pkin(2) + pkin(3) * t162 + pkin(1);
t1 = [0, -t169, t157, 0, 0, 0, 0, 0, -g(3) * t165 - t169 * t167, -g(3) * t167 + t169 * t165, 0, 0, 0, 0, 0, -g(3) * t161 - t169 * t162, -g(3) * t162 + t169 * t161, -t157, -g(1) * (t168 * t156 - t166 * t163) - g(2) * (t166 * t156 + t168 * t163) - g(3) * (t165 * pkin(2) + pkin(3) * t161 + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t158 - t169 * t159, -g(3) * t159 + t169 * t158;];
U_reg = t1;
