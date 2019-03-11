% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:12
% EndTime: 2019-03-09 11:45:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (137->44), mult. (124->55), div. (0->0), fcn. (130->10), ass. (0->30)
t218 = qJ(2) + pkin(10);
t214 = qJ(4) + t218;
t211 = cos(t214);
t224 = cos(qJ(2));
t213 = t224 * pkin(2) + pkin(1);
t235 = pkin(3) * cos(t218) + t213 + pkin(4) * t211;
t210 = sin(t214);
t233 = g(3) * t210;
t219 = -pkin(7) - qJ(3);
t221 = sin(qJ(2));
t232 = t221 * pkin(2) + pkin(6);
t220 = sin(qJ(5));
t222 = sin(qJ(1));
t231 = t222 * t220;
t223 = cos(qJ(5));
t230 = t222 * t223;
t225 = cos(qJ(1));
t229 = t225 * t220;
t228 = t225 * t223;
t227 = g(1) * t225 + g(2) * t222;
t204 = t211 * t231 + t228;
t206 = t211 * t229 - t230;
t226 = g(1) * t206 + g(2) * t204 + t220 * t233;
t217 = -pkin(8) + t219;
t209 = g(1) * t222 - g(2) * t225;
t207 = t211 * t228 + t231;
t205 = t211 * t230 - t229;
t203 = -g(3) * t211 + t227 * t210;
t202 = -g(1) * t207 - g(2) * t205 - t223 * t233;
t1 = [0, -t227, t209, 0, 0, 0, 0, 0, -g(3) * t221 - t227 * t224, -g(3) * t224 + t227 * t221, -t209, -g(1) * (t225 * t213 - t222 * t219) - g(2) * (t222 * t213 + t225 * t219) - g(3) * t232, 0, 0, 0, 0, 0, -t227 * t211 - t233, t203, 0, 0, 0, 0, 0, t202, t226, t202, -t203, -t226, -g(1) * (t207 * pkin(5) + t206 * qJ(6) - t222 * t217 + t235 * t225) - g(2) * (t205 * pkin(5) + t204 * qJ(6) + t225 * t217 + t235 * t222) - g(3) * (pkin(3) * sin(t218) - t211 * pkin(9) + t232) + (-g(3) * (pkin(5) * t223 + qJ(6) * t220 + pkin(4)) - t227 * pkin(9)) * t210;];
U_reg  = t1;
