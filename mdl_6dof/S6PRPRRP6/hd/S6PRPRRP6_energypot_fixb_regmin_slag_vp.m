% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:55
% EndTime: 2019-03-08 20:20:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (166->66), mult. (392->98), div. (0->0), fcn. (488->10), ass. (0->42)
t214 = sin(pkin(6));
t236 = pkin(7) * t214;
t218 = sin(qJ(4));
t235 = t214 * t218;
t219 = sin(qJ(2));
t234 = t214 * t219;
t221 = cos(qJ(4));
t233 = t214 * t221;
t222 = cos(qJ(2));
t232 = t214 * t222;
t216 = cos(pkin(6));
t231 = t216 * t219;
t230 = t216 * t222;
t229 = pkin(2) * t234 + t216 * pkin(7) + qJ(1);
t213 = sin(pkin(10));
t215 = cos(pkin(10));
t199 = t213 * t219 - t215 * t230;
t200 = t213 * t222 + t215 * t231;
t228 = t213 * pkin(1) + t200 * pkin(2) + t199 * qJ(3);
t201 = t213 * t230 + t215 * t219;
t202 = -t213 * t231 + t215 * t222;
t227 = t215 * pkin(1) + t202 * pkin(2) + t201 * qJ(3) + t213 * t236;
t189 = t201 * t218 + t213 * t233;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t184 = t189 * t217 - t202 * t220;
t191 = t199 * t218 - t215 * t233;
t186 = t191 * t217 - t200 * t220;
t204 = t216 * t221 - t218 * t232;
t192 = t204 * t217 - t220 * t234;
t226 = g(1) * t184 + g(2) * t186 + g(3) * t192;
t188 = -t201 * t221 + t213 * t235;
t190 = t199 * t221 + t215 * t235;
t203 = t216 * t218 + t221 * t232;
t225 = g(1) * t188 - g(2) * t190 + g(3) * t203;
t224 = -g(1) * t201 - g(2) * t199 + g(3) * t232;
t223 = g(1) * t202 + g(2) * t200 + g(3) * t234;
t193 = t204 * t220 + t217 * t234;
t187 = t191 * t220 + t200 * t217;
t185 = t189 * t220 + t202 * t217;
t183 = -g(1) * t185 - g(2) * t187 - g(3) * t193;
t1 = [-g(3) * qJ(1), 0, -t223, -t224, t223, t224, -g(1) * t227 - g(2) * (-t215 * t236 + t228) - g(3) * (-qJ(3) * t232 + t229) 0, 0, 0, 0, 0, -g(1) * t189 - g(2) * t191 - g(3) * t204, t225, 0, 0, 0, 0, 0, t183, t226, t183, -t225, -t226, -g(1) * (t189 * pkin(4) + t185 * pkin(5) + t202 * pkin(8) + t188 * pkin(9) + t184 * qJ(6) + t227) - g(2) * (t191 * pkin(4) + t187 * pkin(5) + t200 * pkin(8) - t190 * pkin(9) + t186 * qJ(6) + t228) - g(3) * (t216 * pkin(3) + t204 * pkin(4) + t193 * pkin(5) + t203 * pkin(9) + t192 * qJ(6) + t229) + (-g(1) * pkin(3) * t213 - g(3) * (pkin(8) * t219 - qJ(3) * t222) - g(2) * (-pkin(3) - pkin(7)) * t215) * t214;];
U_reg  = t1;
