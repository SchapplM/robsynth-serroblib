% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (296->62), mult. (885->138), div. (0->0), fcn. (1202->16), ass. (0->72)
t236 = sin(pkin(14));
t244 = cos(pkin(7));
t272 = t236 * t244;
t238 = sin(pkin(8));
t239 = sin(pkin(7));
t271 = t238 * t239;
t240 = sin(pkin(6));
t270 = t239 * t240;
t243 = cos(pkin(8));
t269 = t239 * t243;
t245 = cos(pkin(6));
t268 = t239 * t245;
t267 = t240 * t244;
t248 = sin(qJ(2));
t266 = t240 * t248;
t241 = cos(pkin(14));
t265 = t241 * t244;
t264 = t244 * t248;
t251 = cos(qJ(2));
t263 = t244 * t251;
t262 = t245 * t248;
t261 = t245 * t251;
t260 = t239 * t266;
t237 = sin(pkin(13));
t242 = cos(pkin(13));
t232 = t237 * t251 + t242 * t262;
t231 = -t237 * t248 + t242 * t261;
t254 = t231 * t244 - t242 * t270;
t215 = -t232 * t236 + t254 * t241;
t226 = -t231 * t239 - t242 * t267;
t259 = t215 * t243 + t226 * t238;
t234 = -t237 * t262 + t242 * t251;
t233 = -t237 * t261 - t242 * t248;
t253 = t233 * t244 + t237 * t270;
t217 = -t234 * t236 + t253 * t241;
t227 = -t233 * t239 + t237 * t267;
t258 = t217 * t243 + t227 * t238;
t224 = t241 * t268 + (-t236 * t248 + t241 * t263) * t240;
t230 = t245 * t244 - t251 * t270;
t257 = t224 * t243 + t230 * t238;
t219 = -t231 * t236 - t232 * t265;
t256 = t219 * t243 + t232 * t271;
t221 = -t233 * t236 - t234 * t265;
t255 = t221 * t243 + t234 * t271;
t228 = (-t236 * t251 - t241 * t264) * t240;
t252 = t228 * t243 + t238 * t260;
t250 = cos(qJ(4));
t249 = cos(qJ(5));
t247 = sin(qJ(4));
t246 = sin(qJ(5));
t229 = (-t236 * t264 + t241 * t251) * t240;
t225 = t241 * t266 + (t240 * t263 + t268) * t236;
t223 = -t228 * t238 + t243 * t260;
t222 = t233 * t241 - t234 * t272;
t220 = t231 * t241 - t232 * t272;
t218 = t234 * t241 + t253 * t236;
t216 = t232 * t241 + t254 * t236;
t214 = -t224 * t238 + t230 * t243;
t213 = t229 * t250 + t252 * t247;
t212 = -t221 * t238 + t234 * t269;
t211 = -t219 * t238 + t232 * t269;
t210 = -t217 * t238 + t227 * t243;
t209 = -t215 * t238 + t226 * t243;
t208 = t225 * t250 + t257 * t247;
t207 = -t225 * t247 + t257 * t250;
t206 = t222 * t250 + t255 * t247;
t205 = t220 * t250 + t256 * t247;
t204 = t218 * t250 + t258 * t247;
t203 = -t218 * t247 + t258 * t250;
t202 = t216 * t250 + t259 * t247;
t201 = -t216 * t247 + t259 * t250;
t1 = [0, t206 * t249 + t212 * t246, 0, t203 * t249, -t204 * t246 + t210 * t249, 0; 0, t205 * t249 + t211 * t246, 0, t201 * t249, -t202 * t246 + t209 * t249, 0; 0, t213 * t249 + t223 * t246, 0, t207 * t249, -t208 * t246 + t214 * t249, 0; 0, -t206 * t246 + t212 * t249, 0, -t203 * t246, -t204 * t249 - t210 * t246, 0; 0, -t205 * t246 + t211 * t249, 0, -t201 * t246, -t202 * t249 - t209 * t246, 0; 0, -t213 * t246 + t223 * t249, 0, -t207 * t246, -t208 * t249 - t214 * t246, 0; 0, t222 * t247 - t255 * t250, 0, t204, 0, 0; 0, t220 * t247 - t256 * t250, 0, t202, 0, 0; 0, t229 * t247 - t252 * t250, 0, t208, 0, 0;];
JR_rot  = t1;
