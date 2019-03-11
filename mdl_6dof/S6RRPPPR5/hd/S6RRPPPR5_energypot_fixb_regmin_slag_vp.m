% Calculate minimal parameter regressor of potential energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:44
% EndTime: 2019-03-09 08:23:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (112->57), mult. (241->76), div. (0->0), fcn. (268->8), ass. (0->32)
t217 = cos(pkin(9));
t219 = sin(qJ(2));
t237 = t217 * t219;
t216 = sin(pkin(9));
t238 = t216 * t219;
t241 = pkin(3) * t237 + qJ(4) * t238;
t240 = g(3) * t219;
t239 = t219 * pkin(2) + pkin(6);
t220 = sin(qJ(1));
t236 = t219 * t220;
t223 = cos(qJ(1));
t235 = t219 * t223;
t222 = cos(qJ(2));
t234 = t220 * t222;
t233 = t223 * t216;
t232 = t223 * t217;
t231 = t220 * pkin(7) + qJ(3) * t235 + (pkin(2) * t222 + pkin(1)) * t223;
t230 = -t222 * qJ(3) + t239;
t229 = g(1) * t223 + g(2) * t220;
t228 = t220 * pkin(1) + pkin(2) * t234 - t223 * pkin(7) + qJ(3) * t236;
t201 = -t220 * t217 + t222 * t233;
t202 = t220 * t216 + t222 * t232;
t227 = t202 * pkin(3) + t201 * qJ(4) + t231;
t199 = t216 * t234 + t232;
t226 = g(1) * t201 + g(2) * t199 + g(3) * t238;
t200 = t217 * t234 - t233;
t225 = g(1) * t202 + g(2) * t200 + g(3) * t237;
t224 = t200 * pkin(3) + t199 * qJ(4) + t228;
t221 = cos(qJ(6));
t218 = sin(qJ(6));
t196 = -g(3) * t222 + t229 * t219;
t1 = [0, -t229, g(1) * t220 - g(2) * t223, 0, 0, 0, 0, 0, -t229 * t222 - t240, t196, -t225, t226, -t196, -g(1) * t231 - g(2) * t228 - g(3) * t230, -t196, t225, -t226, -g(1) * t227 - g(2) * t224 - g(3) * (t230 + t241) -t226, t196, -t225, -g(1) * (pkin(4) * t235 + t202 * qJ(5) + t227) - g(2) * (pkin(4) * t236 + t200 * qJ(5) + t224) - g(3) * (qJ(5) * t237 + (-pkin(4) - qJ(3)) * t222 + t239 + t241) 0, 0, 0, 0, 0, -g(1) * (t201 * t221 + t202 * t218) - g(2) * (t199 * t221 + t200 * t218) - (t216 * t221 + t217 * t218) * t240, -g(1) * (-t201 * t218 + t202 * t221) - g(2) * (-t199 * t218 + t200 * t221) - (-t216 * t218 + t217 * t221) * t240;];
U_reg  = t1;
