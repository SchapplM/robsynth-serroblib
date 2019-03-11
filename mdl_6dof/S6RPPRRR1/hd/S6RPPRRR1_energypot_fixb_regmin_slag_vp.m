% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:38
% EndTime: 2019-03-09 02:18:38
% DurationCPUTime: 0.06s
% Computational Cost: add. (89->29), mult. (70->44), div. (0->0), fcn. (68->12), ass. (0->24)
t187 = pkin(11) + qJ(4);
t186 = qJ(5) + t187;
t180 = sin(t186);
t203 = g(3) * t180;
t202 = g(3) * (qJ(2) + pkin(6));
t188 = qJ(1) + pkin(10);
t183 = sin(t188);
t192 = sin(qJ(6));
t201 = t183 * t192;
t194 = cos(qJ(6));
t200 = t183 * t194;
t185 = cos(t188);
t199 = t185 * t192;
t198 = t185 * t194;
t197 = g(1) * t185 + g(2) * t183;
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t196 = -g(1) * t195 - g(2) * t193;
t190 = cos(pkin(11));
t189 = sin(pkin(11));
t184 = cos(t187);
t182 = sin(t187);
t181 = cos(t186);
t1 = [0, t196, g(1) * t193 - g(2) * t195, t196 * pkin(1) - t202, -g(3) * t189 - t197 * t190, -g(3) * t190 + t197 * t189, -g(1) * t183 + g(2) * t185, -g(1) * (t195 * pkin(1) + t185 * pkin(2) + t183 * qJ(3)) - g(2) * (t193 * pkin(1) + t183 * pkin(2) - t185 * qJ(3)) - t202, 0, 0, 0, 0, 0, -g(3) * t182 - t197 * t184, -g(3) * t184 + t197 * t182, 0, 0, 0, 0, 0, -t197 * t181 - t203, -g(3) * t181 + t197 * t180, 0, 0, 0, 0, 0, -g(1) * (t181 * t198 + t201) - g(2) * (t181 * t200 - t199) - t194 * t203, -g(1) * (-t181 * t199 + t200) - g(2) * (-t181 * t201 - t198) + t192 * t203;];
U_reg  = t1;
