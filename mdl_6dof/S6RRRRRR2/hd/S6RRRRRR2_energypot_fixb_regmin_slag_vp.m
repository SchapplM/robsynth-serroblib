% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR2
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:47
% EndTime: 2019-03-10 03:35:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->28), mult. (74->44), div. (0->0), fcn. (82->12), ass. (0->26)
t221 = qJ(2) + qJ(3);
t219 = qJ(4) + t221;
t213 = sin(t219);
t237 = g(3) * t213;
t220 = qJ(5) + qJ(6);
t215 = sin(t220);
t224 = sin(qJ(1));
t236 = t224 * t215;
t217 = cos(t220);
t235 = t224 * t217;
t222 = sin(qJ(5));
t234 = t224 * t222;
t225 = cos(qJ(5));
t233 = t224 * t225;
t227 = cos(qJ(1));
t232 = t227 * t215;
t231 = t227 * t217;
t230 = t227 * t222;
t229 = t227 * t225;
t228 = g(1) * t227 + g(2) * t224;
t226 = cos(qJ(2));
t223 = sin(qJ(2));
t218 = cos(t221);
t216 = sin(t221);
t214 = cos(t219);
t1 = [0, -t228, g(1) * t224 - g(2) * t227, 0, 0, 0, 0, 0, -g(3) * t223 - t228 * t226, -g(3) * t226 + t228 * t223, 0, 0, 0, 0, 0, -g(3) * t216 - t228 * t218, -g(3) * t218 + t228 * t216, 0, 0, 0, 0, 0, -t228 * t214 - t237, -g(3) * t214 + t228 * t213, 0, 0, 0, 0, 0, -g(1) * (t214 * t229 + t234) - g(2) * (t214 * t233 - t230) - t225 * t237, -g(1) * (-t214 * t230 + t233) - g(2) * (-t214 * t234 - t229) + t222 * t237, 0, 0, 0, 0, 0, -g(1) * (t214 * t231 + t236) - g(2) * (t214 * t235 - t232) - t217 * t237, -g(1) * (-t214 * t232 + t235) - g(2) * (-t214 * t236 - t231) + t215 * t237;];
U_reg  = t1;
