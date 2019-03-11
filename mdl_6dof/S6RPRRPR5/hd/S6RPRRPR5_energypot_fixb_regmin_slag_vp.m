% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:33
% EndTime: 2019-03-09 05:13:33
% DurationCPUTime: 0.08s
% Computational Cost: add. (109->40), mult. (98->52), div. (0->0), fcn. (96->10), ass. (0->24)
t200 = pkin(10) + qJ(3);
t197 = qJ(4) + t200;
t194 = cos(t197);
t213 = g(3) * t194;
t203 = sin(qJ(6));
t204 = sin(qJ(1));
t212 = t204 * t203;
t205 = cos(qJ(6));
t211 = t204 * t205;
t206 = cos(qJ(1));
t210 = t206 * t203;
t209 = t206 * t205;
t208 = g(1) * t206 + g(2) * t204;
t193 = sin(t197);
t196 = cos(t200);
t202 = cos(pkin(10));
t207 = t202 * pkin(2) + pkin(3) * t196 + pkin(4) * t194 + qJ(5) * t193 + pkin(1);
t201 = sin(pkin(10));
t199 = -pkin(8) - pkin(7) - qJ(2);
t195 = sin(t200);
t192 = g(1) * t204 - g(2) * t206;
t190 = g(3) * t193 + t208 * t194;
t189 = t208 * t193 - t213;
t1 = [0, -t208, t192, -g(3) * t201 - t208 * t202, -g(3) * t202 + t208 * t201, -t192, -g(1) * (t206 * pkin(1) + t204 * qJ(2)) - g(2) * (t204 * pkin(1) - t206 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t195 - t208 * t196, -g(3) * t196 + t208 * t195, 0, 0, 0, 0, 0, -t190, t189, -t192, t190, -t189, -g(3) * (t201 * pkin(2) + pkin(3) * t195 + t193 * pkin(4) - t194 * qJ(5) + pkin(6)) + (-g(1) * t207 - g(2) * t199) * t206 + (g(1) * t199 - g(2) * t207) * t204, 0, 0, 0, 0, 0, -g(1) * (t193 * t210 + t211) - g(2) * (t193 * t212 - t209) + t203 * t213, -g(1) * (t193 * t209 - t212) - g(2) * (t193 * t211 + t210) + t205 * t213;];
U_reg  = t1;
