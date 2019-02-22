% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:03
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:03:49
% EndTime: 2019-02-22 10:03:49
% DurationCPUTime: 0.21s
% Computational Cost: add. (363->65), mult. (853->133), div. (0->0), fcn. (1177->14), ass. (0->67)
t249 = sin(pkin(13));
t252 = cos(pkin(13));
t257 = sin(qJ(2));
t254 = cos(pkin(6));
t260 = cos(qJ(2));
t267 = t254 * t260;
t238 = -t249 * t257 + t252 * t267;
t268 = t254 * t257;
t239 = t249 * t260 + t252 * t268;
t253 = cos(pkin(7));
t256 = sin(qJ(3));
t259 = cos(qJ(3));
t250 = sin(pkin(7));
t251 = sin(pkin(6));
t273 = t250 * t251;
t221 = t239 * t259 + (t238 * t253 - t252 * t273) * t256;
t271 = t251 * t253;
t231 = -t238 * t250 - t252 * t271;
t248 = qJ(4) + qJ(5);
t246 = sin(t248);
t247 = cos(t248);
t211 = -t221 * t246 + t231 * t247;
t255 = sin(qJ(6));
t280 = t211 * t255;
t240 = -t249 * t267 - t252 * t257;
t241 = -t249 * t268 + t252 * t260;
t223 = t241 * t259 + (t240 * t253 + t249 * t273) * t256;
t232 = -t240 * t250 + t249 * t271;
t213 = -t223 * t246 + t232 * t247;
t279 = t213 * t255;
t264 = t257 * t259;
t265 = t256 * t260;
t230 = t254 * t250 * t256 + (t253 * t265 + t264) * t251;
t237 = t254 * t253 - t260 * t273;
t218 = -t230 * t246 + t237 * t247;
t278 = t218 * t255;
t277 = t246 * t250;
t276 = t247 * t250;
t275 = t247 * t255;
t258 = cos(qJ(6));
t274 = t247 * t258;
t272 = t250 * t259;
t270 = t253 * t256;
t269 = t253 * t259;
t266 = t256 * t257;
t263 = t259 * t260;
t262 = t257 * t273;
t261 = t251 * t272;
t236 = (-t253 * t266 + t263) * t251;
t235 = (t253 * t264 + t265) * t251;
t229 = t251 * t266 - t254 * t272 - t263 * t271;
t228 = t236 * t247 + t246 * t262;
t227 = t240 * t259 - t241 * t270;
t226 = t240 * t256 + t241 * t269;
t225 = t238 * t259 - t239 * t270;
t224 = t238 * t256 + t239 * t269;
t222 = -t240 * t269 + t241 * t256 - t249 * t261;
t220 = -t238 * t269 + t239 * t256 + t252 * t261;
t219 = t230 * t247 + t237 * t246;
t217 = t218 * t258;
t216 = t227 * t247 + t241 * t277;
t215 = t225 * t247 + t239 * t277;
t214 = t223 * t247 + t232 * t246;
t212 = t221 * t247 + t231 * t246;
t210 = t213 * t258;
t209 = t211 * t258;
t1 = [0, t216 * t258 + t226 * t255, -t222 * t274 + t223 * t255, t210, t210, -t214 * t255 + t222 * t258; 0, t215 * t258 + t224 * t255, -t220 * t274 + t221 * t255, t209, t209, -t212 * t255 + t220 * t258; 0, t228 * t258 + t235 * t255, -t229 * t274 + t230 * t255, t217, t217, -t219 * t255 + t229 * t258; 0, -t216 * t255 + t226 * t258, t222 * t275 + t223 * t258, -t279, -t279, -t214 * t258 - t222 * t255; 0, -t215 * t255 + t224 * t258, t220 * t275 + t221 * t258, -t280, -t280, -t212 * t258 - t220 * t255; 0, -t228 * t255 + t235 * t258, t229 * t275 + t230 * t258, -t278, -t278, -t219 * t258 - t229 * t255; 0, t227 * t246 - t241 * t276, -t222 * t246, t214, t214, 0; 0, t225 * t246 - t239 * t276, -t220 * t246, t212, t212, 0; 0, t236 * t246 - t247 * t262, -t229 * t246, t219, t219, 0;];
JR_rot  = t1;
