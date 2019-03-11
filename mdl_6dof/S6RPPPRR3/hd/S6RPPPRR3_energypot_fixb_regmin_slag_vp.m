% Calculate minimal parameter regressor of potential energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:58
% EndTime: 2019-03-09 01:33:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->36), mult. (135->54), div. (0->0), fcn. (156->10), ass. (0->26)
t192 = cos(qJ(1));
t191 = sin(qJ(1));
t173 = pkin(10) + qJ(5);
t166 = sin(t173);
t190 = g(3) * t166;
t189 = g(3) * (-qJ(3) + pkin(6));
t188 = cos(pkin(9));
t187 = sin(pkin(9));
t167 = cos(t173);
t177 = sin(qJ(6));
t186 = t167 * t177;
t178 = cos(qJ(6));
t185 = t167 * t178;
t184 = t192 * pkin(1) + t191 * qJ(2);
t183 = t192 * pkin(2) + t184;
t159 = -t191 * t187 - t192 * t188;
t160 = t192 * t187 - t191 * t188;
t182 = g(1) * t160 - g(2) * t159;
t181 = g(1) * t159 + g(2) * t160;
t180 = t191 * pkin(1) - t192 * qJ(2);
t179 = t191 * pkin(2) + t180;
t175 = cos(pkin(10));
t174 = sin(pkin(10));
t162 = -g(1) * t192 - g(2) * t191;
t161 = g(1) * t191 - g(2) * t192;
t1 = [0, t162, t161, t162, -t161, -g(3) * pkin(6) - g(1) * t184 - g(2) * t180, t181, t182, -g(1) * t183 - g(2) * t179 - t189, g(3) * t174 + t181 * t175, g(3) * t175 - t181 * t174, -t182, -g(1) * (-t159 * pkin(3) + t160 * qJ(4) + t183) - g(2) * (-t160 * pkin(3) - t159 * qJ(4) + t179) - t189, 0, 0, 0, 0, 0, t181 * t167 + t190, g(3) * t167 - t181 * t166, 0, 0, 0, 0, 0, -g(1) * (-t159 * t185 + t160 * t177) - g(2) * (-t159 * t177 - t160 * t185) + t178 * t190, -g(1) * (t159 * t186 + t160 * t178) - g(2) * (-t159 * t178 + t160 * t186) - t177 * t190;];
U_reg  = t1;
