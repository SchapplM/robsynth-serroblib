% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:03
% EndTime: 2019-03-09 08:48:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (96->37), mult. (120->55), div. (0->0), fcn. (130->10), ass. (0->25)
t228 = sin(qJ(1));
t232 = cos(qJ(1));
t235 = g(1) * t232 + g(2) * t228;
t223 = qJ(2) + pkin(10);
t219 = sin(t223);
t220 = cos(t223);
t226 = sin(qJ(5));
t230 = cos(qJ(5));
t213 = t219 * t230 - t220 * t226;
t239 = g(3) * t213;
t227 = sin(qJ(2));
t238 = t227 * pkin(2) + pkin(6);
t231 = cos(qJ(2));
t218 = t231 * pkin(2) + pkin(1);
t224 = -pkin(7) - qJ(3);
t237 = t228 * t218 + t232 * t224;
t236 = t232 * t218 - t228 * t224;
t234 = pkin(3) * t220 + qJ(4) * t219;
t233 = t219 * t226 + t220 * t230;
t229 = cos(qJ(6));
t225 = sin(qJ(6));
t214 = g(1) * t228 - g(2) * t232;
t212 = t233 * t232;
t211 = t233 * t228;
t1 = [0, -t235, t214, 0, 0, 0, 0, 0, -g(3) * t227 - t235 * t231, -g(3) * t231 + t235 * t227, -t214, -g(1) * t236 - g(2) * t237 - g(3) * t238, -g(3) * t219 - t235 * t220, -t214, g(3) * t220 - t235 * t219, -g(1) * (t234 * t232 + t236) - g(2) * (t234 * t228 + t237) - g(3) * (t219 * pkin(3) - t220 * qJ(4) + t238) 0, 0, 0, 0, 0, -g(1) * t212 - g(2) * t211 - t239, g(3) * t233 - t235 * t213, 0, 0, 0, 0, 0, -g(1) * (t212 * t229 - t228 * t225) - g(2) * (t211 * t229 + t232 * t225) - t229 * t239, -g(1) * (-t212 * t225 - t228 * t229) - g(2) * (-t211 * t225 + t232 * t229) + t225 * t239;];
U_reg  = t1;
