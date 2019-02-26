% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR14_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.40s
% Computational Cost: add. (343->67), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->70)
t260 = sin(qJ(2));
t261 = sin(qJ(1));
t264 = cos(qJ(2));
t265 = cos(qJ(1));
t289 = cos(pkin(6));
t269 = t265 * t289;
t244 = t260 * t269 + t261 * t264;
t259 = sin(qJ(3));
t263 = cos(qJ(3));
t243 = t261 * t260 - t264 * t269;
t257 = cos(pkin(7));
t255 = sin(pkin(7));
t256 = sin(pkin(6));
t280 = t256 * t265;
t271 = t255 * t280;
t267 = t243 * t257 + t271;
t225 = -t244 * t263 + t259 * t267;
t236 = -t243 * t255 + t257 * t280;
t258 = sin(qJ(4));
t262 = cos(qJ(4));
t215 = t225 * t262 + t236 * t258;
t254 = pkin(13) + qJ(6);
t252 = sin(t254);
t293 = t215 * t252;
t253 = cos(t254);
t292 = t215 * t253;
t213 = t225 * t258 - t236 * t262;
t288 = t244 * t259;
t286 = t252 * t262;
t285 = t253 * t262;
t284 = t255 * t256;
t283 = t255 * t258;
t282 = t255 * t262;
t281 = t256 * t261;
t279 = t257 * t259;
t278 = t257 * t263;
t277 = t259 * t260;
t276 = t259 * t264;
t275 = t260 * t263;
t274 = t263 * t264;
t273 = t260 * t284;
t272 = t255 * t281;
t270 = t261 * t289;
t268 = t289 * t255;
t245 = -t265 * t260 - t264 * t270;
t266 = -t245 * t255 + t257 * t281;
t246 = -t260 * t270 + t265 * t264;
t242 = t257 * t289 - t264 * t284;
t241 = (-t257 * t277 + t274) * t256;
t240 = (t257 * t275 + t276) * t256;
t235 = t259 * t268 + (t257 * t276 + t275) * t256;
t234 = -t263 * t268 + (-t257 * t274 + t277) * t256;
t232 = t241 * t262 + t258 * t273;
t231 = t245 * t263 - t246 * t279;
t230 = t245 * t259 + t246 * t278;
t229 = -t243 * t263 - t244 * t279;
t228 = -t243 * t259 + t244 * t278;
t227 = t246 * t263 + (t245 * t257 + t272) * t259;
t226 = -t245 * t278 + t246 * t259 - t263 * t272;
t224 = -t263 * t267 - t288;
t222 = t243 * t278 + t263 * t271 + t288;
t221 = t235 * t262 + t242 * t258;
t220 = -t235 * t258 + t242 * t262;
t219 = t231 * t262 + t246 * t283;
t218 = t229 * t262 + t244 * t283;
t217 = t227 * t262 + t258 * t266;
t216 = t227 * t258 - t262 * t266;
t212 = t217 * t253 + t226 * t252;
t211 = -t217 * t252 + t226 * t253;
t1 = [t224 * t252 + t292, t219 * t253 + t230 * t252, -t226 * t285 + t227 * t252, -t216 * t253, 0, t211; t212, t218 * t253 + t228 * t252, -t222 * t285 - t225 * t252, t213 * t253, 0, t222 * t253 + t293; 0, t232 * t253 + t240 * t252, -t234 * t285 + t235 * t252, t220 * t253, 0, -t221 * t252 + t234 * t253; t224 * t253 - t293, -t219 * t252 + t230 * t253, t226 * t286 + t227 * t253, t216 * t252, 0, -t212; t211, -t218 * t252 + t228 * t253, t222 * t286 - t225 * t253, -t213 * t252, 0, -t222 * t252 + t292; 0, -t232 * t252 + t240 * t253, t234 * t286 + t235 * t253, -t220 * t252, 0, -t221 * t253 - t234 * t252; t213, t231 * t258 - t246 * t282, -t226 * t258, t217, 0, 0; t216, t229 * t258 - t244 * t282, -t222 * t258, -t215, 0, 0; 0, t241 * t258 - t262 * t273, -t234 * t258, t221, 0, 0;];
JR_rot  = t1;
