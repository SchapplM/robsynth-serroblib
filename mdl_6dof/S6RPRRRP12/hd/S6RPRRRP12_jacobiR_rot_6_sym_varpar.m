% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP12
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
% Datum: 2019-02-22 10:57
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:57:45
% EndTime: 2019-02-22 10:57:45
% DurationCPUTime: 0.28s
% Computational Cost: add. (234->46), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
t248 = cos(pkin(6));
t243 = sin(pkin(12));
t256 = cos(qJ(1));
t264 = t256 * t243;
t246 = cos(pkin(12));
t252 = sin(qJ(1));
t266 = t252 * t246;
t238 = t248 * t264 + t266;
t251 = sin(qJ(3));
t255 = cos(qJ(3));
t263 = t256 * t246;
t267 = t252 * t243;
t237 = -t248 * t263 + t267;
t247 = cos(pkin(7));
t244 = sin(pkin(7));
t245 = sin(pkin(6));
t270 = t245 * t256;
t262 = t244 * t270;
t260 = t237 * t247 + t262;
t226 = -t238 * t255 + t260 * t251;
t232 = -t237 * t244 + t247 * t270;
t250 = sin(qJ(4));
t254 = cos(qJ(4));
t218 = t226 * t254 + t232 * t250;
t249 = sin(qJ(5));
t280 = t218 * t249;
t253 = cos(qJ(5));
t279 = t218 * t253;
t216 = t226 * t250 - t232 * t254;
t259 = t248 * t266 + t264;
t271 = t245 * t252;
t276 = -t244 * t271 + t259 * t247;
t275 = t238 * t251;
t273 = t244 * t248;
t272 = t245 * t246;
t269 = t247 * t255;
t268 = t249 * t254;
t265 = t253 * t254;
t257 = t259 * t244 + t247 * t271;
t239 = -t248 * t267 + t263;
t236 = -t244 * t272 + t248 * t247;
t231 = t251 * t273 + (t246 * t247 * t251 + t243 * t255) * t245;
t230 = t245 * t243 * t251 - t255 * t273 - t269 * t272;
t228 = t239 * t255 - t276 * t251;
t227 = t239 * t251 + t276 * t255;
t225 = -t260 * t255 - t275;
t223 = t237 * t269 + t255 * t262 + t275;
t222 = t231 * t254 + t236 * t250;
t221 = -t231 * t250 + t236 * t254;
t220 = t228 * t254 + t257 * t250;
t219 = t228 * t250 - t257 * t254;
t215 = t220 * t253 + t227 * t249;
t214 = t220 * t249 - t227 * t253;
t1 = [t225 * t249 + t279, 0, -t227 * t265 + t228 * t249, -t219 * t253, -t214, 0; t215, 0, -t223 * t265 - t226 * t249, t216 * t253, t223 * t253 + t280, 0; 0, 0, -t230 * t265 + t231 * t249, t221 * t253, -t222 * t249 + t230 * t253, 0; t216, 0, -t227 * t250, t220, 0, 0; t219, 0, -t223 * t250, -t218, 0, 0; 0, 0, -t230 * t250, t222, 0, 0; -t225 * t253 + t280, 0, -t227 * t268 - t228 * t253, -t219 * t249, t215, 0; t214, 0, -t223 * t268 + t226 * t253, t216 * t249, t223 * t249 - t279, 0; 0, 0, -t230 * t268 - t231 * t253, t221 * t249, t222 * t253 + t230 * t249, 0;];
JR_rot  = t1;
