% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR5
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:00
% EndTime: 2019-03-09 02:29:00
% DurationCPUTime: 0.05s
% Computational Cost: add. (45->28), mult. (68->38), div. (0->0), fcn. (66->8), ass. (0->19)
t148 = qJ(4) + qJ(5);
t143 = cos(t148);
t161 = g(3) * t143;
t149 = sin(qJ(6));
t151 = sin(qJ(1));
t160 = t151 * t149;
t152 = cos(qJ(6));
t159 = t151 * t152;
t154 = cos(qJ(1));
t158 = t154 * t149;
t157 = t154 * t152;
t156 = t154 * pkin(1) + t151 * qJ(2);
t155 = t151 * pkin(1) - t154 * qJ(2);
t141 = g(1) * t154 + g(2) * t151;
t153 = cos(qJ(4));
t150 = sin(qJ(4));
t142 = sin(t148);
t140 = g(1) * t151 - g(2) * t154;
t1 = [0, -t141, t140, t141, -t140, -g(3) * pkin(6) - g(1) * t156 - g(2) * t155, -t140, -t141, -g(1) * (t154 * qJ(3) + t156) - g(2) * (t151 * qJ(3) + t155) - g(3) * (pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -g(3) * t153 - t141 * t150, g(3) * t150 - t141 * t153, 0, 0, 0, 0, 0, -t141 * t142 - t161, g(3) * t142 - t141 * t143, 0, 0, 0, 0, 0, -g(1) * (t142 * t157 - t160) - g(2) * (t142 * t159 + t158) - t152 * t161, -g(1) * (-t142 * t158 - t159) - g(2) * (-t142 * t160 + t157) + t149 * t161;];
U_reg  = t1;
