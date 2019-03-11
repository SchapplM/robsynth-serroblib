% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR6
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:24
% EndTime: 2019-03-09 13:53:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (79->33), mult. (131->52), div. (0->0), fcn. (150->10), ass. (0->24)
t221 = sin(qJ(1));
t225 = cos(qJ(1));
t231 = g(1) * t225 + g(2) * t221;
t217 = qJ(4) + qJ(5);
t215 = sin(t217);
t216 = cos(t217);
t220 = sin(qJ(2));
t224 = cos(qJ(2));
t230 = t224 * t215 - t220 * t216;
t232 = g(3) * t230;
t229 = t220 * t215 + t224 * t216;
t219 = sin(qJ(4));
t223 = cos(qJ(4));
t228 = t224 * t219 - t220 * t223;
t227 = t220 * t219 + t224 * t223;
t226 = pkin(2) * t224 + qJ(3) * t220 + pkin(1);
t222 = cos(qJ(6));
t218 = sin(qJ(6));
t214 = g(1) * t221 - g(2) * t225;
t212 = -g(3) * t220 - t231 * t224;
t211 = -g(3) * t224 + t231 * t220;
t210 = t229 * t225;
t209 = t229 * t221;
t1 = [0, -t231, t214, 0, 0, 0, 0, 0, t212, t211, t212, -t214, -t211, -g(3) * (t220 * pkin(2) - t224 * qJ(3) + pkin(6)) + (g(2) * pkin(7) - g(1) * t226) * t225 + (-g(1) * pkin(7) - g(2) * t226) * t221, 0, 0, 0, 0, 0, g(3) * t228 - t231 * t227, g(3) * t227 + t231 * t228, 0, 0, 0, 0, 0, -g(1) * t210 - g(2) * t209 + t232, g(3) * t229 + t231 * t230, 0, 0, 0, 0, 0, -g(1) * (t210 * t222 - t221 * t218) - g(2) * (t209 * t222 + t225 * t218) + t222 * t232, -g(1) * (-t210 * t218 - t221 * t222) - g(2) * (-t209 * t218 + t225 * t222) - t218 * t232;];
U_reg  = t1;
