% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR2
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
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:13
% EndTime: 2019-03-09 21:59:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (125->42), mult. (102->59), div. (0->0), fcn. (107->12), ass. (0->29)
t223 = qJ(2) + qJ(3);
t220 = qJ(4) + t223;
t214 = sin(t220);
t240 = g(3) * t214;
t221 = pkin(11) + qJ(6);
t216 = sin(t221);
t227 = sin(qJ(1));
t239 = t227 * t216;
t217 = cos(t221);
t238 = t227 * t217;
t224 = sin(pkin(11));
t237 = t227 * t224;
t225 = cos(pkin(11));
t236 = t227 * t225;
t229 = cos(qJ(1));
t235 = t229 * t216;
t234 = t229 * t217;
t233 = t229 * t224;
t232 = t229 * t225;
t231 = g(1) * t229 + g(2) * t227;
t215 = cos(t220);
t219 = cos(t223);
t228 = cos(qJ(2));
t230 = t228 * pkin(2) + pkin(3) * t219 + pkin(4) * t215 + qJ(5) * t214 + pkin(1);
t226 = sin(qJ(2));
t222 = -pkin(9) - pkin(8) - pkin(7);
t218 = sin(t223);
t212 = -g(3) * t215 + t231 * t214;
t1 = [0, -t231, g(1) * t227 - g(2) * t229, 0, 0, 0, 0, 0, -g(3) * t226 - t231 * t228, -g(3) * t228 + t231 * t226, 0, 0, 0, 0, 0, -g(3) * t218 - t231 * t219, -g(3) * t219 + t231 * t218, 0, 0, 0, 0, 0, -t231 * t215 - t240, t212, -g(1) * (t215 * t232 + t237) - g(2) * (t215 * t236 - t233) - t225 * t240, -g(1) * (-t215 * t233 + t236) - g(2) * (-t215 * t237 - t232) + t224 * t240, -t212, -g(3) * (t226 * pkin(2) + pkin(3) * t218 + t214 * pkin(4) - t215 * qJ(5) + pkin(6)) + (-g(1) * t230 - g(2) * t222) * t229 + (g(1) * t222 - g(2) * t230) * t227, 0, 0, 0, 0, 0, -g(1) * (t215 * t234 + t239) - g(2) * (t215 * t238 - t235) - t217 * t240, -g(1) * (-t215 * t235 + t238) - g(2) * (-t215 * t239 - t234) + t216 * t240;];
U_reg  = t1;
