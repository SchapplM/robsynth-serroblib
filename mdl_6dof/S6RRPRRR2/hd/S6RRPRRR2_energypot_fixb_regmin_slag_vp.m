% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:09
% EndTime: 2019-03-09 13:19:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->33), mult. (76->49), div. (0->0), fcn. (81->10), ass. (0->26)
t220 = qJ(2) + pkin(11) + qJ(4);
t217 = sin(t220);
t240 = g(3) * t217;
t223 = qJ(5) + qJ(6);
t221 = sin(t223);
t227 = sin(qJ(1));
t239 = t227 * t221;
t222 = cos(t223);
t238 = t227 * t222;
t225 = sin(qJ(5));
t237 = t227 * t225;
t228 = cos(qJ(5));
t236 = t227 * t228;
t230 = cos(qJ(1));
t235 = t230 * t221;
t234 = t230 * t222;
t233 = t230 * t225;
t232 = t230 * t228;
t231 = g(1) * t230 + g(2) * t227;
t229 = cos(qJ(2));
t226 = sin(qJ(2));
t224 = -pkin(7) - qJ(3);
t219 = t229 * pkin(2) + pkin(1);
t218 = cos(t220);
t216 = g(1) * t227 - g(2) * t230;
t1 = [0, -t231, t216, 0, 0, 0, 0, 0, -g(3) * t226 - t231 * t229, -g(3) * t229 + t231 * t226, -t216, -g(1) * (t230 * t219 - t227 * t224) - g(2) * (t227 * t219 + t230 * t224) - g(3) * (t226 * pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -t231 * t218 - t240, -g(3) * t218 + t231 * t217, 0, 0, 0, 0, 0, -g(1) * (t218 * t232 + t237) - g(2) * (t218 * t236 - t233) - t228 * t240, -g(1) * (-t218 * t233 + t236) - g(2) * (-t218 * t237 - t232) + t225 * t240, 0, 0, 0, 0, 0, -g(1) * (t218 * t234 + t239) - g(2) * (t218 * t238 - t235) - t222 * t240, -g(1) * (-t218 * t235 + t238) - g(2) * (-t218 * t239 - t234) + t221 * t240;];
U_reg  = t1;
