% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR14_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.28s
% Computational Cost: add. (231->57), mult. (689->127), div. (0->0), fcn. (950->14), ass. (0->61)
t235 = sin(qJ(2));
t236 = sin(qJ(1));
t239 = cos(qJ(2));
t240 = cos(qJ(1));
t263 = cos(pkin(6));
t245 = t240 * t263;
t222 = t235 * t245 + t236 * t239;
t234 = sin(qJ(3));
t238 = cos(qJ(3));
t221 = t236 * t235 - t239 * t245;
t229 = sin(pkin(7));
t232 = cos(pkin(7));
t230 = sin(pkin(6));
t256 = t230 * t240;
t243 = t221 * t232 + t229 * t256;
t204 = -t222 * t238 + t243 * t234;
t215 = -t221 * t229 + t232 * t256;
t233 = sin(qJ(4));
t237 = cos(qJ(4));
t194 = t204 * t233 - t215 * t237;
t195 = t204 * t237 + t215 * t233;
t228 = sin(pkin(13));
t261 = t228 * t237;
t260 = t229 * t230;
t259 = t229 * t233;
t258 = t229 * t237;
t257 = t230 * t236;
t231 = cos(pkin(13));
t255 = t231 * t237;
t254 = t232 * t234;
t253 = t232 * t238;
t252 = t234 * t235;
t251 = t234 * t239;
t250 = t235 * t238;
t249 = t238 * t239;
t248 = t235 * t260;
t247 = t229 * t257;
t246 = t236 * t263;
t244 = t263 * t229;
t223 = -t240 * t235 - t239 * t246;
t242 = -t223 * t229 + t232 * t257;
t241 = -t222 * t234 - t243 * t238;
t224 = -t235 * t246 + t240 * t239;
t220 = t263 * t232 - t239 * t260;
t219 = (-t232 * t252 + t249) * t230;
t218 = (t232 * t250 + t251) * t230;
t214 = t234 * t244 + (t232 * t251 + t250) * t230;
t213 = t238 * t244 + (t232 * t249 - t252) * t230;
t211 = t219 * t237 + t233 * t248;
t210 = t223 * t238 - t224 * t254;
t209 = t223 * t234 + t224 * t253;
t208 = -t221 * t238 - t222 * t254;
t207 = -t221 * t234 + t222 * t253;
t206 = t224 * t238 + (t223 * t232 + t247) * t234;
t205 = -t223 * t253 + t224 * t234 - t238 * t247;
t200 = -t214 * t233 + t220 * t237;
t199 = t210 * t237 + t224 * t259;
t198 = t208 * t237 + t222 * t259;
t197 = t206 * t237 + t242 * t233;
t196 = t206 * t233 - t242 * t237;
t1 = [t195 * t231 + t228 * t241, t199 * t231 + t209 * t228, -t205 * t255 + t206 * t228, -t196 * t231, 0, 0; t197 * t231 + t205 * t228, t198 * t231 + t207 * t228, -t204 * t228 + t241 * t255, t194 * t231, 0, 0; 0, t211 * t231 + t218 * t228, t213 * t255 + t214 * t228, t200 * t231, 0, 0; -t195 * t228 + t231 * t241, -t199 * t228 + t209 * t231, t205 * t261 + t206 * t231, t196 * t228, 0, 0; -t197 * t228 + t205 * t231, -t198 * t228 + t207 * t231, -t204 * t231 - t241 * t261, -t194 * t228, 0, 0; 0, -t211 * t228 + t218 * t231, -t213 * t261 + t214 * t231, -t200 * t228, 0, 0; t194, t210 * t233 - t224 * t258, -t205 * t233, t197, 0, 0; t196, t208 * t233 - t222 * t258, t241 * t233, -t195, 0, 0; 0, t219 * t233 - t237 * t248, t213 * t233, t214 * t237 + t220 * t233, 0, 0;];
JR_rot  = t1;
