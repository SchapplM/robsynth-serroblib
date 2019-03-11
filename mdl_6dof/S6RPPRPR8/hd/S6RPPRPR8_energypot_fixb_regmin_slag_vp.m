% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:10
% EndTime: 2019-03-09 01:56:10
% DurationCPUTime: 0.07s
% Computational Cost: add. (79->43), mult. (102->52), div. (0->0), fcn. (97->8), ass. (0->26)
t183 = pkin(2) + pkin(6);
t166 = pkin(9) + qJ(4);
t161 = sin(t166);
t182 = g(3) * t161;
t170 = sin(qJ(6));
t171 = sin(qJ(1));
t181 = t171 * t170;
t172 = cos(qJ(6));
t180 = t171 * t172;
t173 = cos(qJ(1));
t179 = t173 * t170;
t178 = t173 * t172;
t177 = t173 * pkin(1) + t171 * qJ(2);
t176 = g(1) * t177;
t164 = t171 * pkin(1);
t175 = -t173 * qJ(2) + t164;
t159 = g(1) * t171 - g(2) * t173;
t162 = cos(t166);
t167 = sin(pkin(9));
t174 = pkin(3) * t167 + pkin(4) * t161 - qJ(5) * t162;
t169 = -pkin(7) - qJ(3);
t168 = cos(pkin(9));
t160 = g(1) * t173 + g(2) * t171;
t158 = t159 * t162 - t182;
t157 = g(3) * t162 + t159 * t161;
t1 = [0, -t160, t159, t160, -t159, -g(3) * pkin(6) - g(2) * t175 - t176, -g(3) * t168 - t159 * t167, g(3) * t167 - t159 * t168, -t160, -g(1) * (t173 * qJ(3) + t177) - g(2) * (t171 * qJ(3) + t175) - g(3) * t183, 0, 0, 0, 0, 0, -t157, -t158, -t160, t157, t158, -t176 - g(2) * t164 - g(3) * (t168 * pkin(3) + t162 * pkin(4) + t161 * qJ(5) + t183) + (-g(1) * t174 + g(2) * t169) * t171 + (g(1) * t169 - g(2) * (-qJ(2) - t174)) * t173, 0, 0, 0, 0, 0, -g(1) * (-t162 * t181 + t178) - g(2) * (t162 * t179 + t180) - t170 * t182, -g(1) * (-t162 * t180 - t179) - g(2) * (t162 * t178 - t181) - t172 * t182;];
U_reg  = t1;
