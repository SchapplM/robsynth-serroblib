% Calculate minimal parameter regressor of potential energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:29
% EndTime: 2019-03-08 22:58:29
% DurationCPUTime: 0.12s
% Computational Cost: add. (244->67), mult. (584->93), div. (0->0), fcn. (744->10), ass. (0->41)
t247 = pkin(5) + pkin(9);
t225 = sin(pkin(6));
t246 = pkin(7) * t225;
t229 = sin(qJ(3));
t245 = t225 * t229;
t230 = sin(qJ(2));
t244 = t225 * t230;
t232 = cos(qJ(3));
t243 = t225 * t232;
t233 = cos(qJ(2));
t242 = t225 * t233;
t227 = cos(pkin(6));
t241 = t227 * t230;
t240 = t227 * t233;
t224 = sin(pkin(10));
t226 = cos(pkin(10));
t210 = t224 * t233 + t226 * t241;
t199 = t210 * t232 - t226 * t245;
t209 = t224 * t230 - t226 * t240;
t228 = sin(qJ(4));
t231 = cos(qJ(4));
t191 = t199 * t228 - t209 * t231;
t212 = -t224 * t241 + t226 * t233;
t201 = t212 * t232 + t224 * t245;
t211 = t224 * t240 + t226 * t230;
t193 = t201 * t228 - t211 * t231;
t214 = t227 * t229 + t230 * t243;
t202 = t214 * t228 + t231 * t242;
t239 = g(1) * t193 + g(2) * t191 + g(3) * t202;
t192 = t199 * t231 + t209 * t228;
t194 = t201 * t231 + t211 * t228;
t203 = t214 * t231 - t228 * t242;
t238 = g(1) * t194 + g(2) * t192 + g(3) * t203;
t198 = t210 * t229 + t226 * t243;
t200 = t212 * t229 - t224 * t243;
t213 = -t227 * t232 + t229 * t244;
t237 = g(1) * t200 + g(2) * t198 + g(3) * t213;
t236 = t226 * pkin(1) + t212 * pkin(2) + t201 * pkin(3) + t194 * pkin(4) + t211 * pkin(8) + t193 * qJ(5) + t224 * t246;
t235 = pkin(2) * t244 + t214 * pkin(3) + t203 * pkin(4) + t227 * pkin(7) - pkin(8) * t242 + t202 * qJ(5) + qJ(1);
t234 = t224 * pkin(1) + t210 * pkin(2) + t199 * pkin(3) + t192 * pkin(4) + t209 * pkin(8) + t191 * qJ(5) - t226 * t246;
t1 = [-g(3) * qJ(1), 0, -g(1) * t212 - g(2) * t210 - g(3) * t244, g(1) * t211 + g(2) * t209 - g(3) * t242, 0, 0, 0, 0, 0, -g(1) * t201 - g(2) * t199 - g(3) * t214, t237, 0, 0, 0, 0, 0, -t238, t239, -t237, t238, -t239, -g(1) * (t200 * pkin(9) + t236) - g(2) * (t198 * pkin(9) + t234) - g(3) * (t213 * pkin(9) + t235) -t237, -t239, -t238, -g(1) * (t194 * qJ(6) + t247 * t200 + t236) - g(2) * (t192 * qJ(6) + t247 * t198 + t234) - g(3) * (t203 * qJ(6) + t247 * t213 + t235);];
U_reg  = t1;
