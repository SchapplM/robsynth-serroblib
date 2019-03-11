% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:50
% EndTime: 2019-03-09 04:29:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (153->50), mult. (141->71), div. (0->0), fcn. (144->10), ass. (0->33)
t208 = -qJ(5) - pkin(8);
t210 = sin(qJ(3));
t228 = -t208 * t210 + pkin(2);
t227 = g(3) * t210;
t226 = qJ(2) + pkin(6);
t207 = qJ(1) + pkin(9);
t200 = sin(t207);
t209 = sin(qJ(4));
t225 = t200 * t209;
t213 = cos(qJ(3));
t224 = t200 * t213;
t202 = cos(t207);
t223 = t202 * t213;
t221 = t209 * t213;
t212 = cos(qJ(4));
t220 = t212 * t213;
t198 = t212 * pkin(4) + pkin(3);
t219 = t210 * t198 + t213 * t208 + t226;
t218 = g(1) * t202 + g(2) * t200;
t211 = sin(qJ(1));
t214 = cos(qJ(1));
t217 = -g(1) * t214 - g(2) * t211;
t216 = t214 * pkin(1) + pkin(4) * t225 + t200 * pkin(7) + t198 * t223 + t228 * t202;
t215 = t198 * t224 + t211 * pkin(1) + (-pkin(4) * t209 - pkin(7)) * t202 + t228 * t200;
t206 = qJ(4) + pkin(10);
t201 = cos(t206);
t199 = sin(t206);
t189 = -g(3) * t213 + t218 * t210;
t188 = t200 * t199 + t201 * t223;
t187 = t199 * t223 - t200 * t201;
t186 = -t202 * t199 + t201 * t224;
t185 = t199 * t224 + t202 * t201;
t1 = [0, t217, g(1) * t211 - g(2) * t214, t217 * pkin(1) - g(3) * t226, 0, 0, 0, 0, 0, -t218 * t213 - t227, t189, 0, 0, 0, 0, 0, -g(1) * (t202 * t220 + t225) - g(2) * (t200 * t220 - t202 * t209) - t212 * t227, -g(1) * (t200 * t212 - t202 * t221) - g(2) * (-t200 * t221 - t202 * t212) + t209 * t227, -t189, -g(1) * t216 - g(2) * t215 - g(3) * t219, -g(1) * t188 - g(2) * t186 - t201 * t227, -t189, -g(1) * t187 - g(2) * t185 - t199 * t227, -g(1) * (t188 * pkin(5) + t187 * qJ(6) + t216) - g(2) * (t186 * pkin(5) + t185 * qJ(6) + t215) - g(3) * ((pkin(5) * t201 + qJ(6) * t199) * t210 + t219);];
U_reg  = t1;
