% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:32
% EndTime: 2019-02-26 21:13:32
% DurationCPUTime: 0.30s
% Computational Cost: add. (237->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
t229 = cos(pkin(6));
t224 = sin(pkin(12));
t237 = cos(qJ(1));
t245 = t237 * t224;
t227 = cos(pkin(12));
t233 = sin(qJ(1));
t247 = t233 * t227;
t219 = t229 * t245 + t247;
t232 = sin(qJ(3));
t236 = cos(qJ(3));
t244 = t237 * t227;
t248 = t233 * t224;
t218 = -t229 * t244 + t248;
t228 = cos(pkin(7));
t225 = sin(pkin(7));
t226 = sin(pkin(6));
t251 = t226 * t237;
t243 = t225 * t251;
t241 = t218 * t228 + t243;
t207 = -t219 * t236 + t241 * t232;
t213 = -t218 * t225 + t228 * t251;
t231 = sin(qJ(4));
t235 = cos(qJ(4));
t199 = t207 * t235 + t213 * t231;
t230 = sin(qJ(5));
t261 = t199 * t230;
t234 = cos(qJ(5));
t260 = t199 * t234;
t197 = t207 * t231 - t213 * t235;
t240 = t229 * t247 + t245;
t252 = t226 * t233;
t257 = -t225 * t252 + t240 * t228;
t256 = t219 * t232;
t254 = t225 * t229;
t253 = t226 * t227;
t250 = t228 * t236;
t249 = t230 * t235;
t246 = t234 * t235;
t238 = t240 * t225 + t228 * t252;
t220 = -t229 * t248 + t244;
t217 = -t225 * t253 + t229 * t228;
t212 = t232 * t254 + (t227 * t228 * t232 + t224 * t236) * t226;
t211 = t226 * t224 * t232 - t236 * t254 - t250 * t253;
t209 = t220 * t236 - t257 * t232;
t208 = t220 * t232 + t257 * t236;
t206 = -t241 * t236 - t256;
t204 = t218 * t250 + t236 * t243 + t256;
t203 = t212 * t235 + t217 * t231;
t202 = -t212 * t231 + t217 * t235;
t201 = t209 * t235 + t238 * t231;
t200 = t209 * t231 - t238 * t235;
t196 = t201 * t234 + t208 * t230;
t195 = -t201 * t230 + t208 * t234;
t1 = [t206 * t230 + t260, 0, -t208 * t246 + t209 * t230, -t200 * t234, t195, 0; t196, 0, -t204 * t246 - t207 * t230, t197 * t234, t204 * t234 + t261, 0; 0, 0, -t211 * t246 + t212 * t230, t202 * t234, -t203 * t230 + t211 * t234, 0; t206 * t234 - t261, 0, t208 * t249 + t209 * t234, t200 * t230, -t196, 0; t195, 0, t204 * t249 - t207 * t234, -t197 * t230, -t204 * t230 + t260, 0; 0, 0, t211 * t249 + t212 * t234, -t202 * t230, -t203 * t234 - t211 * t230, 0; t197, 0, -t208 * t231, t201, 0, 0; t200, 0, -t204 * t231, -t199, 0, 0; 0, 0, -t211 * t231, t203, 0, 0;];
JR_rot  = t1;
