% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
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
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:44
% EndTime: 2019-02-26 21:06:45
% DurationCPUTime: 0.23s
% Computational Cost: add. (275->48), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->55)
t233 = cos(pkin(6));
t228 = sin(pkin(12));
t239 = cos(qJ(1));
t247 = t239 * t228;
t231 = cos(pkin(12));
t236 = sin(qJ(1));
t248 = t236 * t231;
t220 = t233 * t247 + t248;
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t246 = t239 * t231;
t249 = t236 * t228;
t219 = -t233 * t246 + t249;
t232 = cos(pkin(7));
t229 = sin(pkin(7));
t230 = sin(pkin(6));
t251 = t230 * t239;
t245 = t229 * t251;
t243 = t219 * t232 + t245;
t208 = -t220 * t238 + t243 * t235;
t214 = -t219 * t229 + t232 * t251;
t234 = sin(qJ(4));
t237 = cos(qJ(4));
t200 = t208 * t237 + t214 * t234;
t227 = pkin(13) + qJ(6);
t225 = sin(t227);
t263 = t200 * t225;
t226 = cos(t227);
t262 = t200 * t226;
t198 = t208 * t234 - t214 * t237;
t242 = t233 * t248 + t247;
t252 = t230 * t236;
t259 = -t229 * t252 + t242 * t232;
t258 = t220 * t235;
t256 = t225 * t237;
t255 = t226 * t237;
t254 = t229 * t233;
t253 = t230 * t231;
t250 = t232 * t238;
t240 = t242 * t229 + t232 * t252;
t221 = -t233 * t249 + t246;
t218 = -t229 * t253 + t233 * t232;
t213 = t235 * t254 + (t231 * t232 * t235 + t228 * t238) * t230;
t212 = t230 * t228 * t235 - t238 * t254 - t250 * t253;
t210 = t221 * t238 - t259 * t235;
t209 = t221 * t235 + t259 * t238;
t207 = -t243 * t238 - t258;
t205 = t219 * t250 + t238 * t245 + t258;
t204 = t213 * t237 + t218 * t234;
t203 = -t213 * t234 + t218 * t237;
t202 = t210 * t237 + t240 * t234;
t201 = t210 * t234 - t240 * t237;
t197 = t202 * t226 + t209 * t225;
t196 = -t202 * t225 + t209 * t226;
t1 = [t207 * t225 + t262, 0, -t209 * t255 + t210 * t225, -t201 * t226, 0, t196; t197, 0, -t205 * t255 - t208 * t225, t198 * t226, 0, t205 * t226 + t263; 0, 0, -t212 * t255 + t213 * t225, t203 * t226, 0, -t204 * t225 + t212 * t226; t207 * t226 - t263, 0, t209 * t256 + t210 * t226, t201 * t225, 0, -t197; t196, 0, t205 * t256 - t208 * t226, -t198 * t225, 0, -t205 * t225 + t262; 0, 0, t212 * t256 + t213 * t226, -t203 * t225, 0, -t204 * t226 - t212 * t225; t198, 0, -t209 * t234, t202, 0, 0; t201, 0, -t205 * t234, -t200, 0, 0; 0, 0, -t212 * t234, t204, 0, 0;];
JR_rot  = t1;
