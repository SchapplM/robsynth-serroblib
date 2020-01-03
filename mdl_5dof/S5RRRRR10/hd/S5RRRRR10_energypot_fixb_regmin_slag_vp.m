% Calculate minimal parameter regressor of potential energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:57
% EndTime: 2019-12-31 22:35:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (86->41), mult. (158->77), div. (0->0), fcn. (202->12), ass. (0->30)
t217 = sin(pkin(5));
t221 = sin(qJ(2));
t235 = t217 * t221;
t222 = sin(qJ(1));
t234 = t217 * t222;
t224 = cos(qJ(3));
t233 = t217 * t224;
t225 = cos(qJ(2));
t232 = t217 * t225;
t226 = cos(qJ(1));
t231 = t217 * t226;
t230 = t222 * t221;
t229 = t222 * t225;
t228 = t226 * t221;
t227 = t226 * t225;
t223 = cos(qJ(5));
t220 = sin(qJ(3));
t219 = sin(qJ(5));
t218 = cos(pkin(5));
t216 = qJ(3) + qJ(4);
t215 = cos(t216);
t214 = sin(t216);
t213 = -t218 * t230 + t227;
t212 = t218 * t229 + t228;
t211 = t218 * t228 + t229;
t210 = -t218 * t227 + t230;
t209 = t218 * t214 + t215 * t235;
t208 = t213 * t215 + t214 * t234;
t207 = t211 * t215 - t214 * t231;
t1 = [0, -g(1) * t226 - g(2) * t222, g(1) * t222 - g(2) * t226, 0, 0, 0, 0, 0, -g(1) * t213 - g(2) * t211 - g(3) * t235, g(1) * t212 + g(2) * t210 - g(3) * t232, 0, 0, 0, 0, 0, -g(1) * (t213 * t224 + t220 * t234) - g(2) * (t211 * t224 - t220 * t231) - g(3) * (t218 * t220 + t221 * t233), -g(1) * (-t213 * t220 + t222 * t233) - g(2) * (-t211 * t220 - t224 * t231) - g(3) * (t218 * t224 - t220 * t235), 0, 0, 0, 0, 0, -g(1) * t208 - g(2) * t207 - g(3) * t209, -g(1) * (-t213 * t214 + t215 * t234) - g(2) * (-t211 * t214 - t215 * t231) - g(3) * (-t214 * t235 + t218 * t215), 0, 0, 0, 0, 0, -g(1) * (t208 * t223 + t212 * t219) - g(2) * (t207 * t223 + t210 * t219) - g(3) * (t209 * t223 - t219 * t232), -g(1) * (-t208 * t219 + t212 * t223) - g(2) * (-t207 * t219 + t210 * t223) - g(3) * (-t209 * t219 - t223 * t232);];
U_reg = t1;
