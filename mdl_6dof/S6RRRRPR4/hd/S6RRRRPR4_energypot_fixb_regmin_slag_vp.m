% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:19
% EndTime: 2019-03-09 22:09:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (98->40), mult. (96->55), div. (0->0), fcn. (101->10), ass. (0->28)
t221 = qJ(2) + qJ(3);
t219 = sin(t221);
t241 = g(3) * t219;
t218 = qJ(4) + pkin(11) + qJ(6);
t214 = sin(t218);
t225 = sin(qJ(1));
t240 = t225 * t214;
t215 = cos(t218);
t239 = t225 * t215;
t223 = sin(qJ(4));
t238 = t225 * t223;
t226 = cos(qJ(4));
t237 = t225 * t226;
t228 = cos(qJ(1));
t236 = t228 * t214;
t235 = t228 * t215;
t234 = t228 * t223;
t233 = t228 * t226;
t232 = pkin(4) * t223 + pkin(7) + pkin(8);
t231 = g(1) * t228 + g(2) * t225;
t216 = t226 * pkin(4) + pkin(3);
t220 = cos(t221);
t222 = -qJ(5) - pkin(9);
t227 = cos(qJ(2));
t230 = t227 * pkin(2) + t216 * t220 - t219 * t222 + pkin(1);
t224 = sin(qJ(2));
t213 = -g(3) * t220 + t231 * t219;
t1 = [0, -t231, g(1) * t225 - g(2) * t228, 0, 0, 0, 0, 0, -g(3) * t224 - t231 * t227, -g(3) * t227 + t231 * t224, 0, 0, 0, 0, 0, -t231 * t220 - t241, t213, 0, 0, 0, 0, 0, -g(1) * (t220 * t233 + t238) - g(2) * (t220 * t237 - t234) - t226 * t241, -g(1) * (-t220 * t234 + t237) - g(2) * (-t220 * t238 - t233) + t223 * t241, -t213, -g(3) * (t224 * pkin(2) + t219 * t216 + t220 * t222 + pkin(6)) + (-g(1) * t230 + g(2) * t232) * t228 + (-g(1) * t232 - g(2) * t230) * t225, 0, 0, 0, 0, 0, -g(1) * (t220 * t235 + t240) - g(2) * (t220 * t239 - t236) - t215 * t241, -g(1) * (-t220 * t236 + t239) - g(2) * (-t220 * t240 - t235) + t214 * t241;];
U_reg  = t1;
