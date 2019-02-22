% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:12
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR15_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:11:53
% EndTime: 2019-02-22 12:11:54
% DurationCPUTime: 0.34s
% Computational Cost: add. (296->66), mult. (867->137), div. (0->0), fcn. (1196->14), ass. (0->72)
t258 = sin(qJ(2));
t259 = sin(qJ(1));
t263 = cos(qJ(2));
t264 = cos(qJ(1));
t287 = cos(pkin(6));
t267 = t264 * t287;
t244 = t258 * t267 + t259 * t263;
t257 = sin(qJ(3));
t262 = cos(qJ(3));
t243 = t259 * t258 - t263 * t267;
t254 = cos(pkin(7));
t252 = sin(pkin(7));
t253 = sin(pkin(6));
t280 = t253 * t264;
t269 = t252 * t280;
t265 = t243 * t254 + t269;
t225 = -t244 * t262 + t257 * t265;
t255 = sin(qJ(6));
t291 = t225 * t255;
t260 = cos(qJ(6));
t290 = t225 * t260;
t235 = -t243 * t252 + t254 * t280;
t256 = sin(qJ(5));
t289 = t235 * t256;
t261 = cos(qJ(5));
t288 = t235 * t261;
t286 = t244 * t257;
t284 = t252 * t253;
t283 = t252 * t256;
t282 = t252 * t261;
t281 = t253 * t259;
t279 = t254 * t257;
t278 = t254 * t262;
t277 = t255 * t256;
t276 = t256 * t260;
t275 = t257 * t258;
t274 = t257 * t263;
t273 = t258 * t262;
t272 = t262 * t263;
t271 = t258 * t284;
t270 = t252 * t281;
t268 = t259 * t287;
t266 = t287 * t252;
t246 = -t258 * t268 + t264 * t263;
t245 = -t264 * t258 - t263 * t268;
t242 = t254 * t287 - t263 * t284;
t241 = (-t254 * t275 + t272) * t253;
t240 = (t254 * t273 + t274) * t253;
t237 = -t245 * t252 + t254 * t281;
t234 = t257 * t266 + (t254 * t274 + t273) * t253;
t233 = -t262 * t266 + (-t254 * t272 + t275) * t253;
t232 = t240 * t256 + t261 * t271;
t231 = t245 * t262 - t246 * t279;
t230 = t245 * t257 + t246 * t278;
t229 = -t243 * t262 - t244 * t279;
t228 = -t243 * t257 + t244 * t278;
t227 = t246 * t262 + (t245 * t254 + t270) * t257;
t226 = -t245 * t278 + t246 * t257 - t262 * t270;
t224 = -t262 * t265 - t286;
t222 = t243 * t278 + t262 * t269 + t286;
t221 = t233 * t256 + t242 * t261;
t220 = t233 * t261 - t242 * t256;
t218 = t230 * t256 + t246 * t282;
t217 = t228 * t256 + t244 * t282;
t216 = t226 * t256 + t237 * t261;
t215 = -t226 * t261 + t237 * t256;
t214 = t222 * t256 - t288;
t213 = t222 * t261 + t289;
t212 = t224 * t256 + t288;
t211 = t216 * t260 + t227 * t255;
t210 = -t216 * t255 + t227 * t260;
t1 = [t212 * t260 + t291, t218 * t260 + t231 * t255, -t226 * t255 + t227 * t276, 0, -t215 * t260, t210; t211, t217 * t260 + t229 * t255, -t222 * t255 - t225 * t276, 0, t213 * t260, -t214 * t255 - t290; 0, t232 * t260 + t241 * t255, -t233 * t255 + t234 * t276, 0, t220 * t260, -t221 * t255 + t234 * t260; -t212 * t255 + t290, -t218 * t255 + t231 * t260, -t226 * t260 - t227 * t277, 0, t215 * t255, -t211; t210, -t217 * t255 + t229 * t260, -t222 * t260 + t225 * t277, 0, -t213 * t255, -t214 * t260 + t291; 0, -t232 * t255 + t241 * t260, -t233 * t260 - t234 * t277, 0, -t220 * t255, -t221 * t260 - t234 * t255; -t224 * t261 + t289, -t230 * t261 + t246 * t283, -t227 * t261, 0, t216, 0; t215, -t228 * t261 + t244 * t283, t225 * t261, 0, t214, 0; 0, -t240 * t261 + t256 * t271, -t234 * t261, 0, t221, 0;];
JR_rot  = t1;
