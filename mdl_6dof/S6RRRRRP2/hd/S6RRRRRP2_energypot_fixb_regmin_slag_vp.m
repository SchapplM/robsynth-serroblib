% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:10
% EndTime: 2019-03-10 01:03:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (137->40), mult. (122->52), div. (0->0), fcn. (131->10), ass. (0->28)
t213 = qJ(2) + qJ(3);
t211 = qJ(4) + t213;
t208 = cos(t211);
t210 = cos(t213);
t218 = cos(qJ(2));
t228 = t218 * pkin(2) + pkin(3) * t210 + pkin(4) * t208 + pkin(1);
t207 = sin(t211);
t226 = g(3) * t207;
t214 = sin(qJ(5));
t216 = sin(qJ(1));
t225 = t216 * t214;
t217 = cos(qJ(5));
t224 = t216 * t217;
t219 = cos(qJ(1));
t223 = t219 * t214;
t222 = t219 * t217;
t221 = g(1) * t219 + g(2) * t216;
t201 = t208 * t225 + t222;
t203 = t208 * t223 - t224;
t220 = g(1) * t203 + g(2) * t201 + t214 * t226;
t215 = sin(qJ(2));
t212 = -pkin(9) - pkin(8) - pkin(7);
t209 = sin(t213);
t204 = t208 * t222 + t225;
t202 = t208 * t224 - t223;
t200 = -g(3) * t208 + t221 * t207;
t199 = -g(1) * t204 - g(2) * t202 - t217 * t226;
t1 = [0, -t221, g(1) * t216 - g(2) * t219, 0, 0, 0, 0, 0, -g(3) * t215 - t221 * t218, -g(3) * t218 + t221 * t215, 0, 0, 0, 0, 0, -g(3) * t209 - t221 * t210, -g(3) * t210 + t221 * t209, 0, 0, 0, 0, 0, -t221 * t208 - t226, t200, 0, 0, 0, 0, 0, t199, t220, t199, -t200, -t220, -g(1) * (t204 * pkin(5) + t203 * qJ(6) - t216 * t212 + t228 * t219) - g(2) * (t202 * pkin(5) + t201 * qJ(6) + t219 * t212 + t228 * t216) - g(3) * (t215 * pkin(2) + pkin(3) * t209 - t208 * pkin(10) + pkin(6)) + (-g(3) * (pkin(5) * t217 + qJ(6) * t214 + pkin(4)) - t221 * pkin(10)) * t207;];
U_reg  = t1;
