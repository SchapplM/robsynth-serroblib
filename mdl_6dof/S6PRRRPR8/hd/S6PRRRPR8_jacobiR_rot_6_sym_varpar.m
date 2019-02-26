% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:37
% EndTime: 2019-02-26 20:14:37
% DurationCPUTime: 0.21s
% Computational Cost: add. (228->61), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->60)
t220 = sin(pkin(7));
t221 = sin(pkin(6));
t249 = t220 * t221;
t226 = sin(qJ(4));
t248 = t220 * t226;
t230 = cos(qJ(4));
t247 = t220 * t230;
t231 = cos(qJ(3));
t246 = t220 * t231;
t223 = cos(pkin(7));
t245 = t221 * t223;
t227 = sin(qJ(3));
t244 = t223 * t227;
t243 = t223 * t231;
t224 = cos(pkin(6));
t228 = sin(qJ(2));
t242 = t224 * t228;
t232 = cos(qJ(2));
t241 = t224 * t232;
t225 = sin(qJ(6));
t240 = t225 * t226;
t229 = cos(qJ(6));
t239 = t226 * t229;
t238 = t227 * t228;
t237 = t227 * t232;
t236 = t228 * t231;
t235 = t231 * t232;
t234 = t228 * t249;
t233 = t221 * t246;
t222 = cos(pkin(12));
t219 = sin(pkin(12));
t214 = -t219 * t242 + t222 * t232;
t213 = -t219 * t241 - t222 * t228;
t212 = t219 * t232 + t222 * t242;
t211 = -t219 * t228 + t222 * t241;
t210 = t224 * t223 - t232 * t249;
t209 = (-t223 * t238 + t235) * t221;
t208 = (t223 * t236 + t237) * t221;
t205 = -t213 * t220 + t219 * t245;
t204 = -t211 * t220 - t222 * t245;
t203 = t224 * t220 * t227 + (t223 * t237 + t236) * t221;
t202 = t221 * t238 - t224 * t246 - t235 * t245;
t201 = t209 * t226 - t230 * t234;
t200 = t213 * t231 - t214 * t244;
t199 = t213 * t227 + t214 * t243;
t198 = t211 * t231 - t212 * t244;
t197 = t211 * t227 + t212 * t243;
t196 = t203 * t230 + t210 * t226;
t195 = t203 * t226 - t210 * t230;
t194 = t214 * t231 + (t213 * t223 + t219 * t249) * t227;
t193 = -t213 * t243 + t214 * t227 - t219 * t233;
t192 = t212 * t231 + (t211 * t223 - t222 * t249) * t227;
t191 = -t211 * t243 + t212 * t227 + t222 * t233;
t190 = t200 * t226 - t214 * t247;
t189 = t198 * t226 - t212 * t247;
t188 = t194 * t230 + t205 * t226;
t187 = t194 * t226 - t205 * t230;
t186 = t192 * t230 + t204 * t226;
t185 = t192 * t226 - t204 * t230;
t1 = [0, t190 * t225 + t199 * t229, -t193 * t240 + t194 * t229, t188 * t225, 0, t187 * t229 - t193 * t225; 0, t189 * t225 + t197 * t229, -t191 * t240 + t192 * t229, t186 * t225, 0, t185 * t229 - t191 * t225; 0, t201 * t225 + t208 * t229, -t202 * t240 + t203 * t229, t196 * t225, 0, t195 * t229 - t202 * t225; 0, t190 * t229 - t199 * t225, -t193 * t239 - t194 * t225, t188 * t229, 0, -t187 * t225 - t193 * t229; 0, t189 * t229 - t197 * t225, -t191 * t239 - t192 * t225, t186 * t229, 0, -t185 * t225 - t191 * t229; 0, t201 * t229 - t208 * t225, -t202 * t239 - t203 * t225, t196 * t229, 0, -t195 * t225 - t202 * t229; 0, t200 * t230 + t214 * t248, -t193 * t230, -t187, 0, 0; 0, t198 * t230 + t212 * t248, -t191 * t230, -t185, 0, 0; 0, t209 * t230 + t226 * t234, -t202 * t230, -t195, 0, 0;];
JR_rot  = t1;
