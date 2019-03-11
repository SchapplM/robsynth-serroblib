% Calculate minimal parameter regressor of potential energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:01
% EndTime: 2019-03-09 01:32:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (85->32), mult. (75->45), div. (0->0), fcn. (70->10), ass. (0->25)
t161 = pkin(10) + qJ(5);
t157 = cos(t161);
t179 = g(3) * t157;
t165 = qJ(2) + pkin(6);
t178 = g(3) * t165;
t162 = qJ(1) + pkin(9);
t156 = sin(t162);
t166 = sin(qJ(6));
t177 = t156 * t166;
t168 = cos(qJ(6));
t176 = t156 * t168;
t158 = cos(t162);
t175 = t158 * t166;
t174 = t158 * t168;
t169 = cos(qJ(1));
t173 = t169 * pkin(1) + t158 * pkin(2) + t156 * qJ(3);
t172 = -g(1) * t156 + g(2) * t158;
t167 = sin(qJ(1));
t171 = -g(1) * t169 - g(2) * t167;
t170 = t167 * pkin(1) + t156 * pkin(2) - t158 * qJ(3);
t164 = cos(pkin(10));
t163 = sin(pkin(10));
t155 = sin(t161);
t151 = g(1) * t158 + g(2) * t156;
t1 = [0, t171, g(1) * t167 - g(2) * t169, t171 * pkin(1) - t178, t151, t172, -g(1) * t173 - g(2) * t170 - t178, -g(3) * t164 + t172 * t163, g(3) * t163 + t172 * t164, -t151, -g(1) * (t158 * qJ(4) + t173) - g(2) * (t156 * qJ(4) + t170) - g(3) * (pkin(3) + t165) 0, 0, 0, 0, 0, t172 * t155 - t179, g(3) * t155 + t172 * t157, 0, 0, 0, 0, 0, -g(1) * (t155 * t176 + t175) - g(2) * (-t155 * t174 + t177) - t168 * t179, -g(1) * (-t155 * t177 + t174) - g(2) * (t155 * t175 + t176) + t166 * t179;];
U_reg  = t1;
