% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:32
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:32:42
% EndTime: 2019-02-22 12:32:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (293->66), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
t256 = sin(qJ(2));
t257 = sin(qJ(1));
t261 = cos(qJ(2));
t262 = cos(qJ(1));
t286 = cos(pkin(6));
t266 = t262 * t286;
t242 = t256 * t266 + t257 * t261;
t255 = sin(qJ(3));
t260 = cos(qJ(3));
t241 = t257 * t256 - t261 * t266;
t252 = cos(pkin(7));
t250 = sin(pkin(7));
t251 = sin(pkin(6));
t279 = t251 * t262;
t268 = t250 * t279;
t264 = t241 * t252 + t268;
t223 = -t242 * t260 + t264 * t255;
t234 = -t241 * t250 + t252 * t279;
t254 = sin(qJ(4));
t259 = cos(qJ(4));
t213 = t223 * t259 + t234 * t254;
t253 = sin(qJ(5));
t290 = t213 * t253;
t258 = cos(qJ(5));
t289 = t213 * t258;
t211 = t223 * t254 - t234 * t259;
t285 = t242 * t255;
t283 = t250 * t251;
t282 = t250 * t254;
t281 = t250 * t259;
t280 = t251 * t257;
t278 = t252 * t255;
t277 = t252 * t260;
t276 = t253 * t259;
t275 = t255 * t256;
t274 = t255 * t261;
t273 = t256 * t260;
t272 = t258 * t259;
t271 = t260 * t261;
t270 = t256 * t283;
t269 = t250 * t280;
t267 = t257 * t286;
t265 = t286 * t250;
t243 = -t262 * t256 - t261 * t267;
t263 = -t243 * t250 + t252 * t280;
t244 = -t256 * t267 + t262 * t261;
t240 = t286 * t252 - t261 * t283;
t239 = (-t252 * t275 + t271) * t251;
t238 = (t252 * t273 + t274) * t251;
t233 = t255 * t265 + (t252 * t274 + t273) * t251;
t232 = -t260 * t265 + (-t252 * t271 + t275) * t251;
t230 = t239 * t259 + t254 * t270;
t229 = t243 * t260 - t244 * t278;
t228 = t243 * t255 + t244 * t277;
t227 = -t241 * t260 - t242 * t278;
t226 = -t241 * t255 + t242 * t277;
t225 = t244 * t260 + (t243 * t252 + t269) * t255;
t224 = -t243 * t277 + t244 * t255 - t260 * t269;
t222 = -t264 * t260 - t285;
t220 = t241 * t277 + t260 * t268 + t285;
t219 = t233 * t259 + t240 * t254;
t218 = -t233 * t254 + t240 * t259;
t217 = t229 * t259 + t244 * t282;
t216 = t227 * t259 + t242 * t282;
t215 = t225 * t259 + t263 * t254;
t214 = t225 * t254 - t263 * t259;
t210 = t215 * t258 + t224 * t253;
t209 = -t215 * t253 + t224 * t258;
t1 = [t222 * t253 + t289, t217 * t258 + t228 * t253, -t224 * t272 + t225 * t253, -t214 * t258, t209, 0; t210, t216 * t258 + t226 * t253, -t220 * t272 - t223 * t253, t211 * t258, t220 * t258 + t290, 0; 0, t230 * t258 + t238 * t253, -t232 * t272 + t233 * t253, t218 * t258, -t219 * t253 + t232 * t258, 0; t222 * t258 - t290, -t217 * t253 + t228 * t258, t224 * t276 + t225 * t258, t214 * t253, -t210, 0; t209, -t216 * t253 + t226 * t258, t220 * t276 - t223 * t258, -t211 * t253, -t220 * t253 + t289, 0; 0, -t230 * t253 + t238 * t258, t232 * t276 + t233 * t258, -t218 * t253, -t219 * t258 - t232 * t253, 0; t211, t229 * t254 - t244 * t281, -t224 * t254, t215, 0, 0; t214, t227 * t254 - t242 * t281, -t220 * t254, -t213, 0, 0; 0, t239 * t254 - t259 * t270, -t232 * t254, t219, 0, 0;];
JR_rot  = t1;
