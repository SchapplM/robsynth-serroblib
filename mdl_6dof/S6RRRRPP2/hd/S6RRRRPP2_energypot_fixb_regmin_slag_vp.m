% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:39
% EndTime: 2019-03-09 20:51:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (156->46), mult. (189->59), div. (0->0), fcn. (205->8), ass. (0->31)
t210 = qJ(2) + qJ(3);
t208 = cos(t210);
t215 = cos(qJ(2));
t232 = t215 * pkin(2) + pkin(3) * t208 + pkin(1);
t207 = sin(t210);
t211 = sin(qJ(4));
t230 = t207 * t211;
t213 = sin(qJ(1));
t229 = t207 * t213;
t214 = cos(qJ(4));
t228 = t207 * t214;
t216 = cos(qJ(1));
t227 = t207 * t216;
t226 = t213 * t211;
t225 = t213 * t214;
t224 = t216 * t211;
t223 = t216 * t214;
t222 = g(1) * t216 + g(2) * t213;
t212 = sin(qJ(2));
t221 = t212 * pkin(2) + t207 * pkin(3) + pkin(4) * t228 + qJ(5) * t230 + pkin(6);
t191 = t208 * t226 + t223;
t192 = t208 * t225 - t224;
t217 = -pkin(8) - pkin(7);
t220 = t192 * pkin(4) + pkin(9) * t229 + t191 * qJ(5) + t232 * t213 + t216 * t217;
t193 = t208 * t224 - t225;
t219 = g(1) * t193 + g(2) * t191 + g(3) * t230;
t194 = t208 * t223 + t226;
t218 = t194 * pkin(4) + pkin(9) * t227 + t193 * qJ(5) - t213 * t217 + t232 * t216;
t188 = -g(3) * t208 + t222 * t207;
t187 = -g(1) * t194 - g(2) * t192 - g(3) * t228;
t1 = [0, -t222, g(1) * t213 - g(2) * t216, 0, 0, 0, 0, 0, -g(3) * t212 - t222 * t215, -g(3) * t215 + t222 * t212, 0, 0, 0, 0, 0, -g(3) * t207 - t222 * t208, t188, 0, 0, 0, 0, 0, t187, t219, t187, -t188, -t219, -g(1) * t218 - g(2) * t220 - g(3) * (-t208 * pkin(9) + t221) t187, -t219, t188, -g(1) * (t194 * pkin(5) - qJ(6) * t227 + t218) - g(2) * (t192 * pkin(5) - qJ(6) * t229 + t220) - g(3) * (pkin(5) * t228 + (-pkin(9) + qJ(6)) * t208 + t221);];
U_reg  = t1;
