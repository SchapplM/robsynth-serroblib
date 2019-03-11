% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:28
% EndTime: 2019-03-09 02:31:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (47->34), mult. (78->48), div. (0->0), fcn. (80->8), ass. (0->23)
t159 = cos(qJ(4));
t171 = g(3) * t159;
t154 = qJ(5) + qJ(6);
t148 = sin(t154);
t157 = sin(qJ(1));
t170 = t157 * t148;
t149 = cos(t154);
t169 = t157 * t149;
t155 = sin(qJ(5));
t168 = t157 * t155;
t158 = cos(qJ(5));
t167 = t157 * t158;
t160 = cos(qJ(1));
t166 = t160 * t148;
t165 = t160 * t149;
t164 = t160 * t155;
t163 = t160 * t158;
t162 = t160 * pkin(1) + t157 * qJ(2);
t161 = t157 * pkin(1) - t160 * qJ(2);
t147 = g(1) * t160 + g(2) * t157;
t156 = sin(qJ(4));
t146 = g(1) * t157 - g(2) * t160;
t1 = [0, -t147, t146, t147, -t146, -g(3) * pkin(6) - g(1) * t162 - g(2) * t161, -t146, -t147, -g(1) * (t160 * qJ(3) + t162) - g(2) * (t157 * qJ(3) + t161) - g(3) * (pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -t147 * t156 - t171, g(3) * t156 - t147 * t159, 0, 0, 0, 0, 0, -g(1) * (t156 * t163 - t168) - g(2) * (t156 * t167 + t164) - t158 * t171, -g(1) * (-t156 * t164 - t167) - g(2) * (-t156 * t168 + t163) + t155 * t171, 0, 0, 0, 0, 0, -g(1) * (t156 * t165 - t170) - g(2) * (t156 * t169 + t166) - t149 * t171, -g(1) * (-t156 * t166 - t169) - g(2) * (-t156 * t170 + t165) + t148 * t171;];
U_reg  = t1;
