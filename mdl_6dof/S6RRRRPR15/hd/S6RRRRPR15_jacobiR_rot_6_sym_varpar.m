% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR15_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:48
% EndTime: 2019-02-26 22:38:48
% DurationCPUTime: 0.41s
% Computational Cost: add. (290->66), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
t252 = sin(qJ(2));
t253 = sin(qJ(1));
t257 = cos(qJ(2));
t258 = cos(qJ(1));
t281 = cos(pkin(6));
t261 = t258 * t281;
t238 = t252 * t261 + t253 * t257;
t251 = sin(qJ(3));
t256 = cos(qJ(3));
t237 = t253 * t252 - t257 * t261;
t248 = cos(pkin(7));
t246 = sin(pkin(7));
t247 = sin(pkin(6));
t274 = t247 * t258;
t263 = t246 * t274;
t259 = t237 * t248 + t263;
t219 = -t238 * t256 + t251 * t259;
t229 = -t237 * t246 + t248 * t274;
t250 = sin(qJ(4));
t255 = cos(qJ(4));
t209 = t219 * t250 - t229 * t255;
t249 = sin(qJ(6));
t286 = t209 * t249;
t254 = cos(qJ(6));
t285 = t209 * t254;
t284 = t219 * t255 + t229 * t250;
t280 = t238 * t251;
t278 = t246 * t247;
t277 = t246 * t250;
t276 = t246 * t255;
t275 = t247 * t253;
t273 = t248 * t251;
t272 = t248 * t256;
t271 = t249 * t250;
t270 = t250 * t254;
t269 = t251 * t252;
t268 = t251 * t257;
t267 = t252 * t256;
t266 = t256 * t257;
t265 = t252 * t278;
t264 = t246 * t275;
t262 = t253 * t281;
t260 = t281 * t246;
t240 = -t252 * t262 + t258 * t257;
t239 = -t258 * t252 - t257 * t262;
t236 = t248 * t281 - t257 * t278;
t235 = (-t248 * t269 + t266) * t247;
t234 = (t248 * t267 + t268) * t247;
t231 = -t239 * t246 + t248 * t275;
t228 = t251 * t260 + (t248 * t268 + t267) * t247;
t227 = -t256 * t260 + (-t248 * t266 + t269) * t247;
t226 = t235 * t250 - t255 * t265;
t225 = t239 * t256 - t240 * t273;
t224 = t239 * t251 + t240 * t272;
t223 = -t237 * t256 - t238 * t273;
t222 = -t237 * t251 + t238 * t272;
t221 = t240 * t256 + (t239 * t248 + t264) * t251;
t220 = -t239 * t272 + t240 * t251 - t256 * t264;
t218 = -t256 * t259 - t280;
t216 = t237 * t272 + t256 * t263 + t280;
t215 = t228 * t255 + t236 * t250;
t214 = t228 * t250 - t236 * t255;
t213 = t225 * t250 - t240 * t276;
t212 = t223 * t250 - t238 * t276;
t211 = t221 * t255 + t231 * t250;
t210 = t221 * t250 - t231 * t255;
t206 = t210 * t249 + t220 * t254;
t205 = t210 * t254 - t220 * t249;
t1 = [t218 * t254 + t286, t213 * t249 + t224 * t254, -t220 * t271 + t221 * t254, t211 * t249, 0, t205; t206, t212 * t249 + t222 * t254, -t216 * t271 - t219 * t254, -t284 * t249, 0, -t216 * t249 - t285; 0, t226 * t249 + t234 * t254, -t227 * t271 + t228 * t254, t215 * t249, 0, t214 * t254 - t227 * t249; -t218 * t249 + t285, t213 * t254 - t224 * t249, -t220 * t270 - t221 * t249, t211 * t254, 0, -t206; t205, t212 * t254 - t222 * t249, -t216 * t270 + t219 * t249, -t284 * t254, 0, -t216 * t254 + t286; 0, t226 * t254 - t234 * t249, -t227 * t270 - t228 * t249, t215 * t254, 0, -t214 * t249 - t227 * t254; t284, t225 * t255 + t240 * t277, -t220 * t255, -t210, 0, 0; t211, t223 * t255 + t238 * t277, -t216 * t255, t209, 0, 0; 0, t235 * t255 + t250 * t265, -t227 * t255, -t214, 0, 0;];
JR_rot  = t1;
