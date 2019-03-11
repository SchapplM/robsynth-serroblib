% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:16
% EndTime: 2019-03-09 02:01:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (148->47), mult. (127->60), div. (0->0), fcn. (130->10), ass. (0->33)
t186 = pkin(10) + qJ(4);
t182 = cos(t186);
t189 = cos(pkin(10));
t206 = t189 * pkin(3) + pkin(4) * t182 + pkin(2);
t180 = sin(t186);
t204 = g(3) * t180;
t190 = qJ(2) + pkin(6);
t203 = g(3) * t190;
t187 = qJ(1) + pkin(9);
t181 = sin(t187);
t192 = sin(qJ(5));
t202 = t181 * t192;
t194 = cos(qJ(5));
t201 = t181 * t194;
t183 = cos(t187);
t200 = t183 * t192;
t199 = t183 * t194;
t198 = g(1) * t183 + g(2) * t181;
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t197 = -g(1) * t195 - g(2) * t193;
t174 = t182 * t202 + t199;
t176 = t182 * t200 - t201;
t196 = g(1) * t176 + g(2) * t174 + t192 * t204;
t191 = -pkin(7) - qJ(3);
t188 = sin(pkin(10));
t185 = t195 * pkin(1);
t184 = t193 * pkin(1);
t177 = t182 * t199 + t202;
t175 = t182 * t201 - t200;
t173 = -g(3) * t182 + t198 * t180;
t172 = -g(1) * t177 - g(2) * t175 - t194 * t204;
t1 = [0, t197, g(1) * t193 - g(2) * t195, t197 * pkin(1) - t203, -g(3) * t188 - t198 * t189, -g(3) * t189 + t198 * t188, -g(1) * t181 + g(2) * t183, -g(1) * (t183 * pkin(2) + t181 * qJ(3) + t185) - g(2) * (t181 * pkin(2) - t183 * qJ(3) + t184) - t203, 0, 0, 0, 0, 0, -t198 * t182 - t204, t173, 0, 0, 0, 0, 0, t172, t196, t172, -t173, -t196, -g(1) * (t177 * pkin(5) + t176 * qJ(6) - t181 * t191 + t206 * t183 + t185) - g(2) * (t175 * pkin(5) + t174 * qJ(6) + t206 * t181 + t183 * t191 + t184) - g(3) * (t188 * pkin(3) - t182 * pkin(8) + t190) + (-g(3) * (pkin(5) * t194 + qJ(6) * t192 + pkin(4)) - t198 * pkin(8)) * t180;];
U_reg  = t1;
