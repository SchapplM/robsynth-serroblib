% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:38
% EndTime: 2019-03-08 19:37:38
% DurationCPUTime: 0.19s
% Computational Cost: add. (182->62), mult. (432->103), div. (0->0), fcn. (550->12), ass. (0->42)
t238 = pkin(7) + qJ(3);
t215 = sin(pkin(6));
t220 = sin(qJ(4));
t237 = t215 * t220;
t221 = sin(qJ(2));
t236 = t215 * t221;
t223 = cos(qJ(4));
t235 = t215 * t223;
t218 = cos(pkin(6));
t234 = t218 * t221;
t224 = cos(qJ(2));
t233 = t218 * t224;
t202 = pkin(2) * t234 - t238 * t215;
t210 = t224 * pkin(2) + pkin(1);
t214 = sin(pkin(10));
t217 = cos(pkin(10));
t232 = t217 * t202 + t214 * t210;
t231 = -t214 * t202 + t217 * t210;
t230 = pkin(2) * t236 + t238 * t218 + qJ(1);
t213 = sin(pkin(11));
t216 = cos(pkin(11));
t229 = t224 * t213 + t221 * t216;
t228 = t221 * t213 - t224 * t216;
t227 = t228 * t218;
t201 = t229 * t218;
t192 = t217 * t201 - t214 * t228;
t187 = t192 * t220 + t217 * t235;
t194 = -t214 * t201 - t217 * t228;
t189 = t194 * t220 - t214 * t235;
t200 = t229 * t215;
t195 = t200 * t220 - t218 * t223;
t226 = g(1) * t189 + g(2) * t187 + g(3) * t195;
t188 = t192 * t223 - t217 * t237;
t190 = t194 * t223 + t214 * t237;
t196 = t200 * t223 + t218 * t220;
t225 = g(1) * t190 + g(2) * t188 + g(3) * t196;
t222 = cos(qJ(6));
t219 = sin(qJ(6));
t199 = t228 * t215;
t193 = t214 * t227 - t217 * t229;
t191 = -t214 * t229 - t217 * t227;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t214 * t234 + t217 * t224) - g(2) * (t214 * t224 + t217 * t234) - g(3) * t236, -g(1) * (-t214 * t233 - t217 * t221) - g(2) * (-t214 * t221 + t217 * t233) - g(3) * t215 * t224, -g(1) * t231 - g(2) * t232 - g(3) * t230, 0, 0, 0, 0, 0, -t225, t226, g(1) * t193 + g(2) * t191 - g(3) * t199, t225, -t226, -g(1) * (t194 * pkin(3) + t190 * pkin(4) - t193 * pkin(8) + t189 * qJ(5) + t231) - g(2) * (t192 * pkin(3) + t188 * pkin(4) - t191 * pkin(8) + t187 * qJ(5) + t232) - g(3) * (t200 * pkin(3) + t196 * pkin(4) + t199 * pkin(8) + t195 * qJ(5) + t230) 0, 0, 0, 0, 0, -g(1) * (t189 * t219 - t193 * t222) - g(2) * (t187 * t219 - t191 * t222) - g(3) * (t195 * t219 + t199 * t222) -g(1) * (t189 * t222 + t193 * t219) - g(2) * (t187 * t222 + t191 * t219) - g(3) * (t195 * t222 - t199 * t219);];
U_reg  = t1;
