% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:14
% EndTime: 2019-02-26 20:05:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (369->60), mult. (1054->136), div. (0->0), fcn. (1453->16), ass. (0->62)
t270 = cos(qJ(3));
t269 = cos(pkin(13));
t240 = sin(pkin(12));
t242 = sin(pkin(6));
t268 = t240 * t242;
t241 = sin(pkin(7));
t267 = t241 * t242;
t247 = sin(qJ(5));
t266 = t241 * t247;
t251 = cos(qJ(5));
t265 = t241 * t251;
t243 = cos(pkin(12));
t264 = t242 * t243;
t244 = cos(pkin(7));
t263 = t242 * t244;
t245 = cos(pkin(6));
t249 = sin(qJ(2));
t262 = t245 * t249;
t252 = cos(qJ(2));
t261 = t245 * t252;
t246 = sin(qJ(6));
t260 = t246 * t251;
t250 = cos(qJ(6));
t259 = t250 * t251;
t258 = t249 * t267;
t239 = sin(pkin(13));
t248 = sin(qJ(3));
t257 = t270 * t239 + t248 * t269;
t256 = -t248 * t239 + t270 * t269;
t227 = t257 * t241;
t229 = t257 * t244;
t231 = -t240 * t249 + t243 * t261;
t232 = t240 * t252 + t243 * t262;
t255 = -t227 * t264 + t231 * t229 + t232 * t256;
t233 = -t240 * t261 - t243 * t249;
t234 = -t240 * t262 + t243 * t252;
t254 = t227 * t268 + t233 * t229 + t234 * t256;
t253 = t245 * t227 + (t229 * t252 + t249 * t256) * t242;
t230 = t245 * t244 - t252 * t267;
t228 = t256 * t244;
t226 = t256 * t241;
t224 = -t233 * t241 + t240 * t263;
t223 = -t231 * t241 - t243 * t263;
t222 = (-t229 * t249 + t252 * t256) * t242;
t221 = (t228 * t249 + t252 * t257) * t242;
t220 = t222 * t251 + t247 * t258;
t219 = -t234 * t229 + t233 * t256;
t218 = t234 * t228 + t233 * t257;
t217 = -t232 * t229 + t231 * t256;
t216 = t232 * t228 + t231 * t257;
t214 = t245 * t226 + (t228 * t252 - t249 * t257) * t242;
t212 = t230 * t247 + t251 * t253;
t211 = t230 * t251 - t247 * t253;
t209 = t226 * t268 + t233 * t228 - t234 * t257;
t206 = -t226 * t264 + t231 * t228 - t232 * t257;
t204 = t219 * t251 + t234 * t266;
t203 = t217 * t251 + t232 * t266;
t202 = t224 * t247 + t251 * t254;
t201 = t224 * t251 - t247 * t254;
t200 = t223 * t247 + t251 * t255;
t199 = t223 * t251 - t247 * t255;
t1 = [0, t204 * t250 + t218 * t246, t209 * t259 + t246 * t254, 0, t201 * t250, -t202 * t246 - t209 * t250; 0, t203 * t250 + t216 * t246, t206 * t259 + t246 * t255, 0, t199 * t250, -t200 * t246 - t206 * t250; 0, t220 * t250 + t221 * t246, t214 * t259 + t246 * t253, 0, t211 * t250, -t212 * t246 - t214 * t250; 0, -t204 * t246 + t218 * t250, -t209 * t260 + t250 * t254, 0, -t201 * t246, -t202 * t250 + t209 * t246; 0, -t203 * t246 + t216 * t250, -t206 * t260 + t250 * t255, 0, -t199 * t246, -t200 * t250 + t206 * t246; 0, -t220 * t246 + t221 * t250, -t214 * t260 + t250 * t253, 0, -t211 * t246, -t212 * t250 + t214 * t246; 0, t219 * t247 - t234 * t265, t209 * t247, 0, t202, 0; 0, t217 * t247 - t232 * t265, t206 * t247, 0, t200, 0; 0, t222 * t247 - t251 * t258, t214 * t247, 0, t212, 0;];
JR_rot  = t1;
