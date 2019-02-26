% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:30
% DurationCPUTime: 0.23s
% Computational Cost: add. (288->48), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->55)
t229 = cos(pkin(6));
t224 = sin(pkin(12));
t235 = cos(qJ(1));
t243 = t235 * t224;
t227 = cos(pkin(12));
t232 = sin(qJ(1));
t244 = t232 * t227;
t216 = t229 * t243 + t244;
t231 = sin(qJ(3));
t234 = cos(qJ(3));
t242 = t235 * t227;
t245 = t232 * t224;
t215 = -t229 * t242 + t245;
t228 = cos(pkin(7));
t225 = sin(pkin(7));
t226 = sin(pkin(6));
t247 = t226 * t235;
t241 = t225 * t247;
t239 = t215 * t228 + t241;
t204 = -t216 * t234 + t239 * t231;
t210 = -t215 * t225 + t228 * t247;
t223 = qJ(4) + pkin(13);
t221 = sin(t223);
t222 = cos(t223);
t196 = t204 * t222 + t210 * t221;
t230 = sin(qJ(6));
t259 = t196 * t230;
t233 = cos(qJ(6));
t258 = t196 * t233;
t194 = t204 * t221 - t210 * t222;
t238 = t229 * t244 + t243;
t248 = t226 * t232;
t255 = -t225 * t248 + t238 * t228;
t254 = t216 * t231;
t252 = t222 * t230;
t251 = t222 * t233;
t250 = t225 * t229;
t249 = t226 * t227;
t246 = t228 * t234;
t236 = t238 * t225 + t228 * t248;
t217 = -t229 * t245 + t242;
t214 = -t225 * t249 + t229 * t228;
t209 = t231 * t250 + (t227 * t228 * t231 + t224 * t234) * t226;
t208 = t226 * t224 * t231 - t234 * t250 - t246 * t249;
t206 = t217 * t234 - t255 * t231;
t205 = t217 * t231 + t255 * t234;
t203 = -t239 * t234 - t254;
t201 = t215 * t246 + t234 * t241 + t254;
t200 = t209 * t222 + t214 * t221;
t199 = -t209 * t221 + t214 * t222;
t198 = t206 * t222 + t236 * t221;
t197 = t206 * t221 - t236 * t222;
t193 = t198 * t233 + t205 * t230;
t192 = -t198 * t230 + t205 * t233;
t1 = [t203 * t230 + t258, 0, -t205 * t251 + t206 * t230, -t197 * t233, 0, t192; t193, 0, -t201 * t251 - t204 * t230, t194 * t233, 0, t201 * t233 + t259; 0, 0, -t208 * t251 + t209 * t230, t199 * t233, 0, -t200 * t230 + t208 * t233; t203 * t233 - t259, 0, t205 * t252 + t206 * t233, t197 * t230, 0, -t193; t192, 0, t201 * t252 - t204 * t233, -t194 * t230, 0, -t201 * t230 + t258; 0, 0, t208 * t252 + t209 * t233, -t199 * t230, 0, -t200 * t233 - t208 * t230; t194, 0, -t205 * t221, t198, 0, 0; t197, 0, -t201 * t221, -t196, 0, 0; 0, 0, -t208 * t221, t200, 0, 0;];
JR_rot  = t1;
