% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:48
% EndTime: 2019-03-09 14:32:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (75->42), mult. (103->61), div. (0->0), fcn. (112->10), ass. (0->31)
t229 = cos(qJ(2));
t245 = g(3) * t229;
t224 = qJ(4) + qJ(5);
t223 = qJ(6) + t224;
t219 = sin(t223);
t227 = sin(qJ(1));
t244 = t227 * t219;
t220 = cos(t223);
t243 = t227 * t220;
t221 = sin(t224);
t242 = t227 * t221;
t222 = cos(t224);
t241 = t227 * t222;
t225 = sin(qJ(4));
t240 = t227 * t225;
t228 = cos(qJ(4));
t239 = t227 * t228;
t230 = cos(qJ(1));
t238 = t230 * t219;
t237 = t230 * t220;
t236 = t230 * t221;
t235 = t230 * t222;
t234 = t230 * t225;
t233 = t230 * t228;
t232 = g(1) * t230 + g(2) * t227;
t226 = sin(qJ(2));
t231 = pkin(2) * t229 + qJ(3) * t226 + pkin(1);
t218 = g(1) * t227 - g(2) * t230;
t217 = g(3) * t226 + t232 * t229;
t216 = t232 * t226 - t245;
t1 = [0, -t232, t218, 0, 0, 0, 0, 0, -t217, t216, -t218, t217, -t216, -g(3) * (t226 * pkin(2) - t229 * qJ(3) + pkin(6)) + (g(2) * pkin(7) - g(1) * t231) * t230 + (-g(1) * pkin(7) - g(2) * t231) * t227, 0, 0, 0, 0, 0, -g(1) * (t226 * t234 + t239) - g(2) * (t226 * t240 - t233) + t225 * t245, -g(1) * (t226 * t233 - t240) - g(2) * (t226 * t239 + t234) + t228 * t245, 0, 0, 0, 0, 0, -g(1) * (t226 * t236 + t241) - g(2) * (t226 * t242 - t235) + t221 * t245, -g(1) * (t226 * t235 - t242) - g(2) * (t226 * t241 + t236) + t222 * t245, 0, 0, 0, 0, 0, -g(1) * (t226 * t238 + t243) - g(2) * (t226 * t244 - t237) + t219 * t245, -g(1) * (t226 * t237 - t244) - g(2) * (t226 * t243 + t238) + t220 * t245;];
U_reg  = t1;
