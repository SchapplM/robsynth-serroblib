% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:51:24
% EndTime: 2019-03-10 03:51:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->40), mult. (94->66), div. (0->0), fcn. (110->12), ass. (0->21)
t233 = sin(qJ(2));
t242 = g(3) * t233;
t234 = sin(qJ(1));
t236 = cos(qJ(2));
t241 = t234 * t236;
t237 = cos(qJ(1));
t240 = t236 * t237;
t231 = qJ(3) + qJ(4);
t230 = qJ(5) + t231;
t227 = qJ(6) + t230;
t224 = cos(t227);
t239 = t237 * t224;
t238 = g(1) * t237 + g(2) * t234;
t235 = cos(qJ(3));
t232 = sin(qJ(3));
t229 = cos(t231);
t228 = sin(t231);
t226 = cos(t230);
t225 = sin(t230);
t223 = sin(t227);
t1 = [0, -t238, g(1) * t234 - g(2) * t237, 0, 0, 0, 0, 0, -t238 * t236 - t242, -g(3) * t236 + t238 * t233, 0, 0, 0, 0, 0, -g(1) * (t234 * t232 + t235 * t240) - g(2) * (-t232 * t237 + t235 * t241) - t235 * t242, -g(1) * (-t232 * t240 + t234 * t235) - g(2) * (-t232 * t241 - t235 * t237) + t232 * t242, 0, 0, 0, 0, 0, -g(1) * (t234 * t228 + t229 * t240) - g(2) * (-t228 * t237 + t229 * t241) - t229 * t242, -g(1) * (-t228 * t240 + t234 * t229) - g(2) * (-t228 * t241 - t229 * t237) + t228 * t242, 0, 0, 0, 0, 0, -g(1) * (t234 * t225 + t226 * t240) - g(2) * (-t225 * t237 + t226 * t241) - t226 * t242, -g(1) * (-t225 * t240 + t234 * t226) - g(2) * (-t225 * t241 - t226 * t237) + t225 * t242, 0, 0, 0, 0, 0, -g(1) * (t234 * t223 + t236 * t239) - g(2) * (-t223 * t237 + t224 * t241) - t224 * t242, -g(1) * (-t223 * t240 + t234 * t224) - g(2) * (-t223 * t241 - t239) + t223 * t242;];
U_reg  = t1;
