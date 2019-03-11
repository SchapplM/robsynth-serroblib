% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:56
% EndTime: 2019-03-09 19:21:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (92->51), mult. (188->78), div. (0->0), fcn. (225->10), ass. (0->24)
t223 = sin(qJ(2));
t234 = g(3) * t223;
t224 = sin(qJ(1));
t227 = cos(qJ(2));
t233 = t224 * t227;
t222 = sin(qJ(3));
t228 = cos(qJ(1));
t232 = t228 * t222;
t226 = cos(qJ(3));
t231 = t228 * t226;
t230 = g(1) * t228 + g(2) * t224;
t213 = t222 * t233 + t231;
t215 = -t224 * t226 + t227 * t232;
t229 = g(1) * t215 + g(2) * t213 + t222 * t234;
t225 = cos(qJ(5));
t221 = sin(qJ(5));
t220 = qJ(5) + qJ(6);
t219 = cos(t220);
t218 = sin(t220);
t216 = t224 * t222 + t227 * t231;
t214 = t226 * t233 - t232;
t212 = -g(3) * t227 + t230 * t223;
t211 = -g(1) * t216 - g(2) * t214 - t226 * t234;
t1 = [0, -t230, g(1) * t224 - g(2) * t228, 0, 0, 0, 0, 0, -t230 * t227 - t234, t212, 0, 0, 0, 0, 0, t211, t229, t211, -t212, -t229, -g(1) * (t216 * pkin(3) + t224 * pkin(7) + t215 * qJ(4) + (pkin(2) * t227 + pkin(1)) * t228) - g(2) * (t224 * pkin(1) + pkin(2) * t233 + t214 * pkin(3) - t228 * pkin(7) + t213 * qJ(4)) - g(3) * (-t227 * pkin(8) + pkin(6)) + (-g(3) * (pkin(3) * t226 + qJ(4) * t222 + pkin(2)) - t230 * pkin(8)) * t223, 0, 0, 0, 0, 0, -g(1) * (t215 * t221 + t216 * t225) - g(2) * (t213 * t221 + t214 * t225) - (t221 * t222 + t225 * t226) * t234, -g(1) * (t215 * t225 - t216 * t221) - g(2) * (t213 * t225 - t214 * t221) - (-t221 * t226 + t222 * t225) * t234, 0, 0, 0, 0, 0, -g(1) * (t215 * t218 + t216 * t219) - g(2) * (t213 * t218 + t214 * t219) - (t218 * t222 + t219 * t226) * t234, -g(1) * (t215 * t219 - t216 * t218) - g(2) * (t213 * t219 - t214 * t218) - (-t218 * t226 + t219 * t222) * t234;];
U_reg  = t1;
