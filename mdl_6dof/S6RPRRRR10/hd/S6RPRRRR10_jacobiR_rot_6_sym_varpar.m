% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:51
% DurationCPUTime: 0.30s
% Computational Cost: add. (363->52), mult. (854->91), div. (0->0), fcn. (1181->14), ass. (0->61)
t251 = cos(pkin(6));
t246 = sin(pkin(13));
t257 = cos(qJ(1));
t265 = t257 * t246;
t249 = cos(pkin(13));
t254 = sin(qJ(1));
t266 = t254 * t249;
t238 = t251 * t265 + t266;
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t264 = t257 * t249;
t267 = t254 * t246;
t237 = -t251 * t264 + t267;
t250 = cos(pkin(7));
t247 = sin(pkin(7));
t248 = sin(pkin(6));
t269 = t248 * t257;
t263 = t247 * t269;
t261 = t237 * t250 + t263;
t226 = -t238 * t256 + t253 * t261;
t232 = -t237 * t247 + t250 * t269;
t245 = qJ(4) + qJ(5);
t243 = sin(t245);
t244 = cos(t245);
t217 = t226 * t244 + t232 * t243;
t252 = sin(qJ(6));
t284 = t217 * t252;
t255 = cos(qJ(6));
t283 = t217 * t255;
t215 = t226 * t243 - t232 * t244;
t260 = t251 * t266 + t265;
t270 = t248 * t254;
t280 = -t247 * t270 + t260 * t250;
t279 = t215 * t252;
t239 = -t251 * t267 + t264;
t228 = t239 * t256 - t280 * t253;
t258 = t247 * t260 + t250 * t270;
t218 = t228 * t243 - t244 * t258;
t278 = t218 * t252;
t272 = t247 * t251;
t231 = t253 * t272 + (t249 * t250 * t253 + t246 * t256) * t248;
t271 = t248 * t249;
t236 = -t247 * t271 + t251 * t250;
t221 = -t231 * t243 + t236 * t244;
t277 = t221 * t252;
t276 = t238 * t253;
t274 = t244 * t252;
t273 = t244 * t255;
t268 = t250 * t256;
t230 = t248 * t246 * t253 - t256 * t272 - t268 * t271;
t227 = t239 * t253 + t280 * t256;
t225 = -t256 * t261 - t276;
t223 = t237 * t268 + t256 * t263 + t276;
t222 = t231 * t244 + t236 * t243;
t220 = t221 * t255;
t219 = t228 * t244 + t243 * t258;
t214 = t218 * t255;
t213 = t215 * t255;
t212 = t219 * t255 + t227 * t252;
t211 = -t219 * t252 + t227 * t255;
t1 = [t225 * t252 + t283, 0, -t227 * t273 + t228 * t252, -t214, -t214, t211; t212, 0, -t223 * t273 - t226 * t252, t213, t213, t223 * t255 + t284; 0, 0, -t230 * t273 + t231 * t252, t220, t220, -t222 * t252 + t230 * t255; t225 * t255 - t284, 0, t227 * t274 + t228 * t255, t278, t278, -t212; t211, 0, t223 * t274 - t226 * t255, -t279, -t279, -t223 * t252 + t283; 0, 0, t230 * t274 + t231 * t255, -t277, -t277, -t222 * t255 - t230 * t252; t215, 0, -t227 * t243, t219, t219, 0; t218, 0, -t223 * t243, -t217, -t217, 0; 0, 0, -t230 * t243, t222, t222, 0;];
JR_rot  = t1;
