% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPPRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:45
% DurationCPUTime: 0.24s
% Computational Cost: add. (492->53), mult. (1433->111), div. (0->0), fcn. (1946->18), ass. (0->60)
t242 = sin(pkin(8));
t248 = cos(pkin(8));
t239 = sin(pkin(14));
t240 = sin(pkin(13));
t244 = sin(pkin(6));
t245 = cos(pkin(14));
t246 = cos(pkin(13));
t249 = cos(pkin(7));
t273 = t246 * t249;
t243 = sin(pkin(7));
t250 = cos(pkin(6));
t275 = t243 * t250;
t259 = t245 * t275 + (-t239 * t240 + t245 * t273) * t244;
t276 = t243 * t244;
t265 = -t246 * t276 + t250 * t249;
t281 = t265 * t242 + t259 * t248;
t247 = cos(pkin(12));
t241 = sin(pkin(12));
t277 = t241 * t250;
t238 = -t240 * t277 + t247 * t246;
t237 = -t247 * t240 - t246 * t277;
t266 = t237 * t249 + t241 * t276;
t260 = t238 * t239 - t266 * t245;
t274 = t244 * t249;
t267 = -t237 * t243 + t241 * t274;
t280 = -t267 * t242 + t260 * t248;
t272 = t247 * t250;
t236 = t240 * t272 + t241 * t246;
t235 = -t241 * t240 + t246 * t272;
t268 = t235 * t249 - t247 * t276;
t261 = t236 * t239 - t268 * t245;
t269 = -t235 * t243 - t247 * t274;
t279 = -t269 * t242 + t261 * t248;
t278 = cos(qJ(4));
t251 = sin(qJ(6));
t255 = cos(qJ(5));
t271 = t251 * t255;
t254 = cos(qJ(6));
t270 = t254 * t255;
t253 = sin(qJ(4));
t252 = sin(qJ(5));
t233 = t244 * t240 * t245 + (t244 * t273 + t275) * t239;
t229 = -t259 * t242 + t265 * t248;
t228 = t238 * t245 + t266 * t239;
t227 = t236 * t245 + t268 * t239;
t224 = t260 * t242 + t267 * t248;
t223 = t261 * t242 + t269 * t248;
t222 = t233 * t278 + t281 * t253;
t221 = t233 * t253 - t281 * t278;
t220 = t228 * t278 - t280 * t253;
t219 = t228 * t253 + t280 * t278;
t218 = t227 * t278 - t279 * t253;
t217 = t227 * t253 + t279 * t278;
t216 = t222 * t255 + t229 * t252;
t215 = -t222 * t252 + t229 * t255;
t214 = t220 * t255 + t224 * t252;
t213 = -t220 * t252 + t224 * t255;
t212 = t218 * t255 + t223 * t252;
t211 = -t218 * t252 + t223 * t255;
t1 = [0, 0, 0, -t219 * t270 + t220 * t251, t213 * t254, -t214 * t251 + t219 * t254; 0, 0, 0, -t217 * t270 + t218 * t251, t211 * t254, -t212 * t251 + t217 * t254; 0, 0, 0, -t221 * t270 + t222 * t251, t215 * t254, -t216 * t251 + t221 * t254; 0, 0, 0, t219 * t271 + t220 * t254, -t213 * t251, -t214 * t254 - t219 * t251; 0, 0, 0, t217 * t271 + t218 * t254, -t211 * t251, -t212 * t254 - t217 * t251; 0, 0, 0, t221 * t271 + t222 * t254, -t215 * t251, -t216 * t254 - t221 * t251; 0, 0, 0, -t219 * t252, t214, 0; 0, 0, 0, -t217 * t252, t212, 0; 0, 0, 0, -t221 * t252, t216, 0;];
JR_rot  = t1;
