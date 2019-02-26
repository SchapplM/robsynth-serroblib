% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
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
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:29
% EndTime: 2019-02-26 20:17:29
% DurationCPUTime: 0.18s
% Computational Cost: add. (231->61), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->60)
t227 = sin(pkin(7));
t228 = sin(pkin(6));
t256 = t227 * t228;
t233 = sin(qJ(4));
t255 = t227 * t233;
t237 = cos(qJ(4));
t254 = t227 * t237;
t238 = cos(qJ(3));
t253 = t227 * t238;
t230 = cos(pkin(7));
t252 = t228 * t230;
t234 = sin(qJ(3));
t251 = t230 * t234;
t250 = t230 * t238;
t231 = cos(pkin(6));
t235 = sin(qJ(2));
t249 = t231 * t235;
t239 = cos(qJ(2));
t248 = t231 * t239;
t232 = sin(qJ(5));
t247 = t232 * t237;
t246 = t234 * t235;
t245 = t234 * t239;
t244 = t235 * t238;
t236 = cos(qJ(5));
t243 = t236 * t237;
t242 = t238 * t239;
t241 = t235 * t256;
t240 = t228 * t253;
t229 = cos(pkin(12));
t226 = sin(pkin(12));
t221 = -t226 * t249 + t229 * t239;
t220 = -t226 * t248 - t229 * t235;
t219 = t226 * t239 + t229 * t249;
t218 = -t226 * t235 + t229 * t248;
t217 = t231 * t230 - t239 * t256;
t216 = (-t230 * t246 + t242) * t228;
t215 = (t230 * t244 + t245) * t228;
t212 = -t220 * t227 + t226 * t252;
t211 = -t218 * t227 - t229 * t252;
t210 = t231 * t227 * t234 + (t230 * t245 + t244) * t228;
t209 = t228 * t246 - t231 * t253 - t242 * t252;
t208 = t216 * t237 + t233 * t241;
t207 = t220 * t238 - t221 * t251;
t206 = t220 * t234 + t221 * t250;
t205 = t218 * t238 - t219 * t251;
t204 = t218 * t234 + t219 * t250;
t203 = t210 * t237 + t217 * t233;
t202 = -t210 * t233 + t217 * t237;
t201 = t221 * t238 + (t220 * t230 + t226 * t256) * t234;
t200 = -t220 * t250 + t221 * t234 - t226 * t240;
t199 = t219 * t238 + (t218 * t230 - t229 * t256) * t234;
t198 = -t218 * t250 + t219 * t234 + t229 * t240;
t197 = t207 * t237 + t221 * t255;
t196 = t205 * t237 + t219 * t255;
t195 = t201 * t237 + t212 * t233;
t194 = -t201 * t233 + t212 * t237;
t193 = t199 * t237 + t211 * t233;
t192 = -t199 * t233 + t211 * t237;
t1 = [0, t197 * t236 + t206 * t232, -t200 * t243 + t201 * t232, t194 * t236, -t195 * t232 + t200 * t236, 0; 0, t196 * t236 + t204 * t232, -t198 * t243 + t199 * t232, t192 * t236, -t193 * t232 + t198 * t236, 0; 0, t208 * t236 + t215 * t232, -t209 * t243 + t210 * t232, t202 * t236, -t203 * t232 + t209 * t236, 0; 0, -t197 * t232 + t206 * t236, t200 * t247 + t201 * t236, -t194 * t232, -t195 * t236 - t200 * t232, 0; 0, -t196 * t232 + t204 * t236, t198 * t247 + t199 * t236, -t192 * t232, -t193 * t236 - t198 * t232, 0; 0, -t208 * t232 + t215 * t236, t209 * t247 + t210 * t236, -t202 * t232, -t203 * t236 - t209 * t232, 0; 0, t207 * t233 - t221 * t254, -t200 * t233, t195, 0, 0; 0, t205 * t233 - t219 * t254, -t198 * t233, t193, 0, 0; 0, t216 * t233 - t237 * t241, -t209 * t233, t203, 0, 0;];
JR_rot  = t1;
