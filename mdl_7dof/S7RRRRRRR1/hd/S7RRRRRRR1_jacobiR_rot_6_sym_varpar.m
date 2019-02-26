% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:21
% EndTime: 2019-02-26 22:54:22
% DurationCPUTime: 0.39s
% Computational Cost: add. (197->62), mult. (593->136), div. (0->0), fcn. (835->12), ass. (0->65)
t233 = sin(qJ(1));
t237 = cos(qJ(3));
t238 = cos(qJ(2));
t245 = t237 * t238;
t231 = sin(qJ(3));
t239 = cos(qJ(1));
t252 = t231 * t239;
t219 = t233 * t245 + t252;
t236 = cos(qJ(4));
t230 = sin(qJ(4));
t232 = sin(qJ(2));
t256 = t230 * t232;
t204 = t219 * t236 + t233 * t256;
t244 = t239 * t237;
t248 = t233 * t231;
t218 = t238 * t248 - t244;
t229 = sin(qJ(5));
t235 = cos(qJ(5));
t193 = t204 * t235 - t218 * t229;
t251 = t232 * t236;
t203 = t219 * t230 - t233 * t251;
t228 = sin(qJ(6));
t234 = cos(qJ(6));
t265 = t193 * t228 - t203 * t234;
t264 = -t193 * t234 - t203 * t228;
t192 = -t204 * t229 - t218 * t235;
t259 = t228 * t230;
t258 = t228 * t235;
t257 = t229 * t236;
t255 = t230 * t234;
t254 = t231 * t232;
t253 = t231 * t238;
t250 = t232 * t237;
t249 = t232 * t239;
t247 = t234 * t235;
t246 = t235 * t236;
t243 = t229 * t254;
t242 = t230 * t254;
t241 = t235 * t254;
t240 = t231 * t249;
t217 = -t230 * t238 + t236 * t250;
t216 = t230 * t250 + t236 * t238;
t223 = t238 * t244 - t248;
t222 = -t233 * t237 - t238 * t252;
t221 = t236 * t245 + t256;
t220 = t230 * t245 - t251;
t214 = t217 * t239;
t213 = t216 * t239;
t212 = t217 * t233;
t211 = t216 * t233;
t210 = (-t229 * t237 - t231 * t246) * t232;
t209 = t223 * t236 + t230 * t249;
t208 = t223 * t230 - t236 * t249;
t207 = t221 * t235 - t229 * t253;
t202 = t217 * t235 - t243;
t201 = -t217 * t229 - t241;
t200 = -t214 * t235 + t229 * t240;
t199 = -t212 * t235 + t233 * t243;
t198 = t222 * t246 - t223 * t229;
t197 = -t218 * t246 - t219 * t229;
t196 = t209 * t235 + t222 * t229;
t195 = t209 * t229 - t222 * t235;
t191 = t196 * t234 + t208 * t228;
t190 = -t196 * t228 + t208 * t234;
t1 = [t264, t200 * t234 - t213 * t228, t198 * t234 + t222 * t259, -t208 * t247 + t209 * t228, -t195 * t234, t190, 0; t191, t199 * t234 - t211 * t228, t197 * t234 - t218 * t259, -t203 * t247 + t204 * t228, t192 * t234, -t265, 0; 0, t207 * t234 + t220 * t228, t210 * t234 - t228 * t242, -t216 * t247 + t217 * t228, t201 * t234, -t202 * t228 + t216 * t234, 0; t265, -t200 * t228 - t213 * t234, -t198 * t228 + t222 * t255, t208 * t258 + t209 * t234, t195 * t228, -t191, 0; t190, -t199 * t228 - t211 * t234, -t197 * t228 - t218 * t255, t203 * t258 + t204 * t234, -t192 * t228, t264, 0; 0, -t207 * t228 + t220 * t234, -t210 * t228 - t234 * t242, t216 * t258 + t217 * t234, -t201 * t228, -t202 * t234 - t216 * t228, 0; t192, -t214 * t229 - t235 * t240, t222 * t257 + t223 * t235, -t208 * t229, t196, 0, 0; t195, -t212 * t229 - t233 * t241, -t218 * t257 + t219 * t235, -t203 * t229, t193, 0, 0; 0, t221 * t229 + t235 * t253 (-t231 * t257 + t235 * t237) * t232, -t216 * t229, t202, 0, 0;];
JR_rot  = t1;
