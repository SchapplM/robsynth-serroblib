% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:32
% EndTime: 2019-02-26 20:54:32
% DurationCPUTime: 0.24s
% Computational Cost: add. (288->48), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->55)
t230 = cos(pkin(6));
t225 = sin(pkin(12));
t236 = cos(qJ(1));
t244 = t236 * t225;
t228 = cos(pkin(12));
t233 = sin(qJ(1));
t245 = t233 * t228;
t217 = t230 * t244 + t245;
t232 = sin(qJ(3));
t235 = cos(qJ(3));
t243 = t236 * t228;
t246 = t233 * t225;
t216 = -t230 * t243 + t246;
t229 = cos(pkin(7));
t226 = sin(pkin(7));
t227 = sin(pkin(6));
t248 = t227 * t236;
t242 = t226 * t248;
t240 = t216 * t229 + t242;
t205 = -t217 * t235 + t240 * t232;
t211 = -t216 * t226 + t229 * t248;
t224 = pkin(13) + qJ(5);
t222 = sin(t224);
t223 = cos(t224);
t197 = t205 * t223 + t211 * t222;
t231 = sin(qJ(6));
t260 = t197 * t231;
t234 = cos(qJ(6));
t259 = t197 * t234;
t195 = t205 * t222 - t211 * t223;
t239 = t230 * t245 + t244;
t249 = t227 * t233;
t256 = -t226 * t249 + t239 * t229;
t255 = t217 * t232;
t253 = t223 * t231;
t252 = t223 * t234;
t251 = t226 * t230;
t250 = t227 * t228;
t247 = t229 * t235;
t237 = t239 * t226 + t229 * t249;
t218 = -t230 * t246 + t243;
t215 = -t226 * t250 + t230 * t229;
t210 = t232 * t251 + (t228 * t229 * t232 + t225 * t235) * t227;
t209 = t227 * t225 * t232 - t235 * t251 - t247 * t250;
t207 = t218 * t235 - t256 * t232;
t206 = t218 * t232 + t256 * t235;
t204 = -t240 * t235 - t255;
t202 = t216 * t247 + t235 * t242 + t255;
t201 = t210 * t223 + t215 * t222;
t200 = -t210 * t222 + t215 * t223;
t199 = t207 * t223 + t237 * t222;
t198 = t207 * t222 - t237 * t223;
t194 = t199 * t234 + t206 * t231;
t193 = -t199 * t231 + t206 * t234;
t1 = [t204 * t231 + t259, 0, -t206 * t252 + t207 * t231, 0, -t198 * t234, t193; t194, 0, -t202 * t252 - t205 * t231, 0, t195 * t234, t202 * t234 + t260; 0, 0, -t209 * t252 + t210 * t231, 0, t200 * t234, -t201 * t231 + t209 * t234; t204 * t234 - t260, 0, t206 * t253 + t207 * t234, 0, t198 * t231, -t194; t193, 0, t202 * t253 - t205 * t234, 0, -t195 * t231, -t202 * t231 + t259; 0, 0, t209 * t253 + t210 * t234, 0, -t200 * t231, -t201 * t234 - t209 * t231; t195, 0, -t206 * t222, 0, t199, 0; t198, 0, -t202 * t222, 0, -t197, 0; 0, 0, -t209 * t222, 0, t201, 0;];
JR_rot  = t1;
