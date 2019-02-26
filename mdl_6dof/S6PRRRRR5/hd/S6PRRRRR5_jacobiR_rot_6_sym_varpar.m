% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:11
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.21s
% Computational Cost: add. (347->62), mult. (869->133), div. (0->0), fcn. (1198->14), ass. (0->67)
t243 = qJ(5) + qJ(6);
t241 = sin(t243);
t253 = cos(qJ(4));
t272 = t241 * t253;
t242 = cos(t243);
t271 = t242 * t253;
t245 = sin(pkin(7));
t246 = sin(pkin(6));
t270 = t245 * t246;
t250 = sin(qJ(4));
t269 = t245 * t250;
t268 = t245 * t253;
t254 = cos(qJ(3));
t267 = t245 * t254;
t248 = cos(pkin(7));
t266 = t246 * t248;
t251 = sin(qJ(3));
t265 = t248 * t251;
t264 = t248 * t254;
t249 = cos(pkin(6));
t252 = sin(qJ(2));
t263 = t249 * t252;
t255 = cos(qJ(2));
t262 = t249 * t255;
t261 = t251 * t252;
t260 = t251 * t255;
t259 = t252 * t254;
t258 = t254 * t255;
t257 = t252 * t270;
t256 = t246 * t267;
t247 = cos(pkin(13));
t244 = sin(pkin(13));
t236 = -t244 * t263 + t247 * t255;
t235 = -t244 * t262 - t247 * t252;
t234 = t244 * t255 + t247 * t263;
t233 = -t244 * t252 + t247 * t262;
t232 = t249 * t248 - t255 * t270;
t231 = (-t248 * t261 + t258) * t246;
t230 = (t248 * t259 + t260) * t246;
t227 = -t235 * t245 + t244 * t266;
t226 = -t233 * t245 - t247 * t266;
t225 = t249 * t245 * t251 + (t248 * t260 + t259) * t246;
t224 = t246 * t261 - t249 * t267 - t258 * t266;
t223 = t231 * t253 + t250 * t257;
t222 = t235 * t254 - t236 * t265;
t221 = t235 * t251 + t236 * t264;
t220 = t233 * t254 - t234 * t265;
t219 = t233 * t251 + t234 * t264;
t218 = t225 * t253 + t232 * t250;
t217 = -t225 * t250 + t232 * t253;
t216 = t236 * t254 + (t235 * t248 + t244 * t270) * t251;
t215 = -t235 * t264 + t236 * t251 - t244 * t256;
t214 = t234 * t254 + (t233 * t248 - t247 * t270) * t251;
t213 = -t233 * t264 + t234 * t251 + t247 * t256;
t212 = t222 * t253 + t236 * t269;
t211 = t220 * t253 + t234 * t269;
t210 = t216 * t253 + t227 * t250;
t209 = -t216 * t250 + t227 * t253;
t208 = t214 * t253 + t226 * t250;
t207 = -t214 * t250 + t226 * t253;
t206 = -t218 * t242 - t224 * t241;
t205 = -t218 * t241 + t224 * t242;
t204 = -t210 * t242 - t215 * t241;
t203 = -t210 * t241 + t215 * t242;
t202 = -t208 * t242 - t213 * t241;
t201 = -t208 * t241 + t213 * t242;
t1 = [0, t212 * t242 + t221 * t241, -t215 * t271 + t216 * t241, t209 * t242, t203, t203; 0, t211 * t242 + t219 * t241, -t213 * t271 + t214 * t241, t207 * t242, t201, t201; 0, t223 * t242 + t230 * t241, -t224 * t271 + t225 * t241, t217 * t242, t205, t205; 0, -t212 * t241 + t221 * t242, t215 * t272 + t216 * t242, -t209 * t241, t204, t204; 0, -t211 * t241 + t219 * t242, t213 * t272 + t214 * t242, -t207 * t241, t202, t202; 0, -t223 * t241 + t230 * t242, t224 * t272 + t225 * t242, -t217 * t241, t206, t206; 0, t222 * t250 - t236 * t268, -t215 * t250, t210, 0, 0; 0, t220 * t250 - t234 * t268, -t213 * t250, t208, 0, 0; 0, t231 * t250 - t253 * t257, -t224 * t250, t218, 0, 0;];
JR_rot  = t1;
