% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:07
% EndTime: 2019-03-09 21:09:07
% DurationCPUTime: 0.11s
% Computational Cost: add. (161->55), mult. (207->74), div. (0->0), fcn. (227->8), ass. (0->33)
t218 = sin(qJ(1));
t221 = cos(qJ(1));
t227 = g(1) * t221 + g(2) * t218;
t217 = sin(qJ(2));
t238 = g(3) * t217;
t215 = qJ(3) + qJ(4);
t210 = sin(t215);
t237 = t210 * t217;
t211 = cos(t215);
t236 = t211 * t217;
t222 = -pkin(9) - pkin(8);
t235 = t217 * t222;
t216 = sin(qJ(3));
t234 = t218 * t216;
t220 = cos(qJ(2));
t233 = t218 * t220;
t232 = t221 * t210;
t231 = t221 * t211;
t230 = t221 * t216;
t219 = cos(qJ(3));
t229 = t221 * t219;
t208 = t219 * pkin(3) + pkin(2);
t226 = pkin(4) * t236 + qJ(5) * t237 + t217 * t208 + t220 * t222 + pkin(6);
t196 = t210 * t233 + t231;
t197 = t211 * t233 - t232;
t225 = t218 * pkin(1) + t197 * pkin(4) + t196 * qJ(5) + t208 * t233;
t198 = -t218 * t211 + t220 * t232;
t199 = t218 * t210 + t220 * t231;
t224 = pkin(3) * t234 + t199 * pkin(4) + t218 * pkin(7) + t198 * qJ(5) + (t208 * t220 + pkin(1)) * t221;
t223 = g(1) * t198 + g(2) * t196 + g(3) * t237;
t200 = -g(3) * t220 + t227 * t217;
t193 = -g(1) * t199 - g(2) * t197 - g(3) * t236;
t1 = [0, -t227, g(1) * t218 - g(2) * t221, 0, 0, 0, 0, 0, -t227 * t220 - t238, t200, 0, 0, 0, 0, 0, -g(1) * (t220 * t229 + t234) - g(2) * (t219 * t233 - t230) - t219 * t238, -g(1) * (t218 * t219 - t220 * t230) - g(2) * (-t216 * t233 - t229) + t216 * t238, 0, 0, 0, 0, 0, t193, t223, t193, -t200, -t223, -g(1) * (-t221 * t235 + t224) - g(2) * (-t218 * t235 + (-pkin(3) * t216 - pkin(7)) * t221 + t225) - g(3) * t226, t193, -t223, t200, -g(1) * (t199 * pkin(5) + t224) - g(2) * (-pkin(3) * t230 + t197 * pkin(5) - t221 * pkin(7) + t225) - g(3) * (t220 * qJ(6) + t226) + (-g(3) * pkin(5) * t211 + t227 * (qJ(6) + t222)) * t217;];
U_reg  = t1;
