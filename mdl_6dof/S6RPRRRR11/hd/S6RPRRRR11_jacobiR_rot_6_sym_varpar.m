% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR11
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

function JR_rot = S6RPRRRR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:26
% EndTime: 2019-02-26 21:20:26
% DurationCPUTime: 0.27s
% Computational Cost: add. (349->49), mult. (870->91), div. (0->0), fcn. (1202->14), ass. (0->59)
t241 = cos(pkin(6));
t236 = sin(pkin(13));
t247 = cos(qJ(1));
t255 = t247 * t236;
t239 = cos(pkin(13));
t244 = sin(qJ(1));
t256 = t244 * t239;
t228 = t241 * t255 + t256;
t243 = sin(qJ(3));
t246 = cos(qJ(3));
t254 = t247 * t239;
t257 = t244 * t236;
t227 = -t241 * t254 + t257;
t240 = cos(pkin(7));
t237 = sin(pkin(7));
t238 = sin(pkin(6));
t259 = t238 * t247;
t253 = t237 * t259;
t251 = t227 * t240 + t253;
t216 = -t228 * t246 + t251 * t243;
t222 = -t227 * t237 + t240 * t259;
t242 = sin(qJ(4));
t245 = cos(qJ(4));
t208 = t216 * t245 + t222 * t242;
t235 = qJ(5) + qJ(6);
t233 = sin(t235);
t271 = t208 * t233;
t234 = cos(t235);
t270 = t208 * t234;
t206 = t216 * t242 - t222 * t245;
t250 = t241 * t256 + t255;
t260 = t238 * t244;
t267 = -t237 * t260 + t250 * t240;
t266 = t228 * t243;
t264 = t233 * t245;
t263 = t234 * t245;
t262 = t237 * t241;
t261 = t238 * t239;
t258 = t240 * t246;
t248 = t250 * t237 + t240 * t260;
t229 = -t241 * t257 + t254;
t226 = -t237 * t261 + t241 * t240;
t221 = t243 * t262 + (t239 * t240 * t243 + t236 * t246) * t238;
t220 = t238 * t236 * t243 - t246 * t262 - t258 * t261;
t218 = t229 * t246 - t267 * t243;
t217 = t229 * t243 + t267 * t246;
t215 = -t251 * t246 - t266;
t213 = t227 * t258 + t246 * t253 + t266;
t212 = t221 * t245 + t226 * t242;
t211 = -t221 * t242 + t226 * t245;
t210 = t218 * t245 + t248 * t242;
t209 = t218 * t242 - t248 * t245;
t205 = -t212 * t234 - t220 * t233;
t204 = -t212 * t233 + t220 * t234;
t203 = t210 * t234 + t217 * t233;
t202 = -t210 * t233 + t217 * t234;
t201 = -t213 * t233 + t270;
t200 = t213 * t234 + t271;
t1 = [t215 * t233 + t270, 0, -t217 * t263 + t218 * t233, -t209 * t234, t202, t202; t203, 0, -t213 * t263 - t216 * t233, t206 * t234, t200, t200; 0, 0, -t220 * t263 + t221 * t233, t211 * t234, t204, t204; t215 * t234 - t271, 0, t217 * t264 + t218 * t234, t209 * t233, -t203, -t203; t202, 0, t213 * t264 - t216 * t234, -t206 * t233, t201, t201; 0, 0, t220 * t264 + t221 * t234, -t211 * t233, t205, t205; t206, 0, -t217 * t242, t210, 0, 0; t209, 0, -t213 * t242, -t208, 0, 0; 0, 0, -t220 * t242, t212, 0, 0;];
JR_rot  = t1;
