% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:04
% EndTime: 2019-02-26 20:14:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (273->62), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->61)
t227 = pkin(13) + qJ(6);
t225 = sin(t227);
t237 = cos(qJ(4));
t256 = t225 * t237;
t226 = cos(t227);
t255 = t226 * t237;
t229 = sin(pkin(7));
t230 = sin(pkin(6));
t254 = t229 * t230;
t234 = sin(qJ(4));
t253 = t229 * t234;
t252 = t229 * t237;
t238 = cos(qJ(3));
t251 = t229 * t238;
t232 = cos(pkin(7));
t250 = t230 * t232;
t235 = sin(qJ(3));
t249 = t232 * t235;
t248 = t232 * t238;
t233 = cos(pkin(6));
t236 = sin(qJ(2));
t247 = t233 * t236;
t239 = cos(qJ(2));
t246 = t233 * t239;
t245 = t235 * t236;
t244 = t235 * t239;
t243 = t236 * t238;
t242 = t238 * t239;
t241 = t236 * t254;
t240 = t230 * t251;
t231 = cos(pkin(12));
t228 = sin(pkin(12));
t220 = -t228 * t247 + t231 * t239;
t219 = -t228 * t246 - t231 * t236;
t218 = t228 * t239 + t231 * t247;
t217 = -t228 * t236 + t231 * t246;
t216 = t233 * t232 - t239 * t254;
t215 = (-t232 * t245 + t242) * t230;
t214 = (t232 * t243 + t244) * t230;
t211 = -t219 * t229 + t228 * t250;
t210 = -t217 * t229 - t231 * t250;
t209 = t233 * t229 * t235 + (t232 * t244 + t243) * t230;
t208 = t230 * t245 - t233 * t251 - t242 * t250;
t207 = t215 * t237 + t234 * t241;
t206 = t219 * t238 - t220 * t249;
t205 = t219 * t235 + t220 * t248;
t204 = t217 * t238 - t218 * t249;
t203 = t217 * t235 + t218 * t248;
t202 = t209 * t237 + t216 * t234;
t201 = -t209 * t234 + t216 * t237;
t200 = t220 * t238 + (t219 * t232 + t228 * t254) * t235;
t199 = -t219 * t248 + t220 * t235 - t228 * t240;
t198 = t218 * t238 + (t217 * t232 - t231 * t254) * t235;
t197 = -t217 * t248 + t218 * t235 + t231 * t240;
t196 = t206 * t237 + t220 * t253;
t195 = t204 * t237 + t218 * t253;
t194 = t200 * t237 + t211 * t234;
t193 = -t200 * t234 + t211 * t237;
t192 = t198 * t237 + t210 * t234;
t191 = -t198 * t234 + t210 * t237;
t1 = [0, t196 * t226 + t205 * t225, -t199 * t255 + t200 * t225, t193 * t226, 0, -t194 * t225 + t199 * t226; 0, t195 * t226 + t203 * t225, -t197 * t255 + t198 * t225, t191 * t226, 0, -t192 * t225 + t197 * t226; 0, t207 * t226 + t214 * t225, -t208 * t255 + t209 * t225, t201 * t226, 0, -t202 * t225 + t208 * t226; 0, -t196 * t225 + t205 * t226, t199 * t256 + t200 * t226, -t193 * t225, 0, -t194 * t226 - t199 * t225; 0, -t195 * t225 + t203 * t226, t197 * t256 + t198 * t226, -t191 * t225, 0, -t192 * t226 - t197 * t225; 0, -t207 * t225 + t214 * t226, t208 * t256 + t209 * t226, -t201 * t225, 0, -t202 * t226 - t208 * t225; 0, t206 * t234 - t220 * t252, -t199 * t234, t194, 0, 0; 0, t204 * t234 - t218 * t252, -t197 * t234, t192, 0, 0; 0, t215 * t234 - t237 * t241, -t208 * t234, t202, 0, 0;];
JR_rot  = t1;
