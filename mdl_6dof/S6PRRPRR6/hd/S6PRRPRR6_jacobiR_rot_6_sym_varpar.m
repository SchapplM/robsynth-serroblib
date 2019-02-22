% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
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
% Datum: 2019-02-22 09:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:50:54
% EndTime: 2019-02-22 09:50:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (288->62), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->61)
t226 = pkin(13) + qJ(5);
t224 = sin(t226);
t228 = sin(pkin(7));
t255 = t224 * t228;
t225 = cos(t226);
t254 = t225 * t228;
t233 = sin(qJ(6));
t253 = t225 * t233;
t236 = cos(qJ(6));
t252 = t225 * t236;
t229 = sin(pkin(6));
t251 = t228 * t229;
t237 = cos(qJ(3));
t250 = t228 * t237;
t231 = cos(pkin(7));
t249 = t229 * t231;
t234 = sin(qJ(3));
t248 = t231 * t234;
t247 = t231 * t237;
t232 = cos(pkin(6));
t235 = sin(qJ(2));
t246 = t232 * t235;
t238 = cos(qJ(2));
t245 = t232 * t238;
t244 = t234 * t235;
t243 = t234 * t238;
t242 = t235 * t237;
t241 = t237 * t238;
t240 = t235 * t251;
t239 = t229 * t250;
t230 = cos(pkin(12));
t227 = sin(pkin(12));
t219 = -t227 * t246 + t230 * t238;
t218 = -t227 * t245 - t230 * t235;
t217 = t227 * t238 + t230 * t246;
t216 = -t227 * t235 + t230 * t245;
t215 = t232 * t231 - t238 * t251;
t214 = (-t231 * t244 + t241) * t229;
t213 = (t231 * t242 + t243) * t229;
t210 = -t218 * t228 + t227 * t249;
t209 = -t216 * t228 - t230 * t249;
t208 = t232 * t228 * t234 + (t231 * t243 + t242) * t229;
t207 = t229 * t244 - t232 * t250 - t241 * t249;
t206 = t214 * t225 + t224 * t240;
t205 = t218 * t237 - t219 * t248;
t204 = t218 * t234 + t219 * t247;
t203 = t216 * t237 - t217 * t248;
t202 = t216 * t234 + t217 * t247;
t201 = t219 * t237 + (t218 * t231 + t227 * t251) * t234;
t200 = -t218 * t247 + t219 * t234 - t227 * t239;
t199 = t217 * t237 + (t216 * t231 - t230 * t251) * t234;
t198 = -t216 * t247 + t217 * t234 + t230 * t239;
t197 = t208 * t225 + t215 * t224;
t196 = -t208 * t224 + t215 * t225;
t195 = t205 * t225 + t219 * t255;
t194 = t203 * t225 + t217 * t255;
t193 = t201 * t225 + t210 * t224;
t192 = -t201 * t224 + t210 * t225;
t191 = t199 * t225 + t209 * t224;
t190 = -t199 * t224 + t209 * t225;
t1 = [0, t195 * t236 + t204 * t233, -t200 * t252 + t201 * t233, 0, t192 * t236, -t193 * t233 + t200 * t236; 0, t194 * t236 + t202 * t233, -t198 * t252 + t199 * t233, 0, t190 * t236, -t191 * t233 + t198 * t236; 0, t206 * t236 + t213 * t233, -t207 * t252 + t208 * t233, 0, t196 * t236, -t197 * t233 + t207 * t236; 0, -t195 * t233 + t204 * t236, t200 * t253 + t201 * t236, 0, -t192 * t233, -t193 * t236 - t200 * t233; 0, -t194 * t233 + t202 * t236, t198 * t253 + t199 * t236, 0, -t190 * t233, -t191 * t236 - t198 * t233; 0, -t206 * t233 + t213 * t236, t207 * t253 + t208 * t236, 0, -t196 * t233, -t197 * t236 - t207 * t233; 0, t205 * t224 - t219 * t254, -t200 * t224, 0, t193, 0; 0, t203 * t224 - t217 * t254, -t198 * t224, 0, t191, 0; 0, t214 * t224 - t225 * t240, -t207 * t224, 0, t197, 0;];
JR_rot  = t1;
