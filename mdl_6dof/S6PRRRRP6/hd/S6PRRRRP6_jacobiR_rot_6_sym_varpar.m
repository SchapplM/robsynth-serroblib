% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:03
% EndTime: 2019-02-26 20:18:03
% DurationCPUTime: 0.21s
% Computational Cost: add. (228->58), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->60)
t253 = sin(pkin(7));
t254 = sin(pkin(6));
t282 = t253 * t254;
t259 = sin(qJ(4));
t281 = t253 * t259;
t263 = cos(qJ(4));
t280 = t253 * t263;
t264 = cos(qJ(3));
t279 = t253 * t264;
t256 = cos(pkin(7));
t278 = t254 * t256;
t260 = sin(qJ(3));
t277 = t256 * t260;
t276 = t256 * t264;
t257 = cos(pkin(6));
t261 = sin(qJ(2));
t275 = t257 * t261;
t265 = cos(qJ(2));
t274 = t257 * t265;
t258 = sin(qJ(5));
t273 = t258 * t263;
t272 = t260 * t261;
t271 = t260 * t265;
t270 = t261 * t264;
t262 = cos(qJ(5));
t269 = t262 * t263;
t268 = t264 * t265;
t267 = t261 * t282;
t266 = t254 * t279;
t255 = cos(pkin(12));
t252 = sin(pkin(12));
t247 = -t252 * t275 + t255 * t265;
t246 = -t252 * t274 - t255 * t261;
t245 = t252 * t265 + t255 * t275;
t244 = -t252 * t261 + t255 * t274;
t243 = t257 * t256 - t265 * t282;
t242 = (-t256 * t272 + t268) * t254;
t241 = (t256 * t270 + t271) * t254;
t238 = -t246 * t253 + t252 * t278;
t237 = -t244 * t253 - t255 * t278;
t236 = t257 * t253 * t260 + (t256 * t271 + t270) * t254;
t235 = t254 * t272 - t257 * t279 - t268 * t278;
t234 = t242 * t263 + t259 * t267;
t233 = t246 * t264 - t247 * t277;
t232 = t246 * t260 + t247 * t276;
t231 = t244 * t264 - t245 * t277;
t230 = t244 * t260 + t245 * t276;
t229 = t236 * t263 + t243 * t259;
t228 = -t236 * t259 + t243 * t263;
t227 = t247 * t264 + (t246 * t256 + t252 * t282) * t260;
t226 = -t246 * t276 + t247 * t260 - t252 * t266;
t225 = t245 * t264 + (t244 * t256 - t255 * t282) * t260;
t224 = -t244 * t276 + t245 * t260 + t255 * t266;
t223 = t233 * t263 + t247 * t281;
t222 = t231 * t263 + t245 * t281;
t221 = t227 * t263 + t238 * t259;
t220 = -t227 * t259 + t238 * t263;
t219 = t225 * t263 + t237 * t259;
t218 = -t225 * t259 + t237 * t263;
t1 = [0, t223 * t262 + t232 * t258, -t226 * t269 + t227 * t258, t220 * t262, -t221 * t258 + t226 * t262, 0; 0, t222 * t262 + t230 * t258, -t224 * t269 + t225 * t258, t218 * t262, -t219 * t258 + t224 * t262, 0; 0, t234 * t262 + t241 * t258, -t235 * t269 + t236 * t258, t228 * t262, -t229 * t258 + t235 * t262, 0; 0, t233 * t259 - t247 * t280, -t226 * t259, t221, 0, 0; 0, t231 * t259 - t245 * t280, -t224 * t259, t219, 0, 0; 0, t242 * t259 - t263 * t267, -t235 * t259, t229, 0, 0; 0, t223 * t258 - t232 * t262, -t226 * t273 - t227 * t262, t220 * t258, t221 * t262 + t226 * t258, 0; 0, t222 * t258 - t230 * t262, -t224 * t273 - t225 * t262, t218 * t258, t219 * t262 + t224 * t258, 0; 0, t234 * t258 - t241 * t262, -t235 * t273 - t236 * t262, t228 * t258, t229 * t262 + t235 * t258, 0;];
JR_rot  = t1;
