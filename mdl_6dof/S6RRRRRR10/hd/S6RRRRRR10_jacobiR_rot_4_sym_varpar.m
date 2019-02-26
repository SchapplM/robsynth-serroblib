% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:57
% EndTime: 2019-02-26 22:52:57
% DurationCPUTime: 0.36s
% Computational Cost: add. (239->56), mult. (721->123), div. (0->0), fcn. (987->14), ass. (0->62)
t225 = sin(qJ(2));
t226 = sin(qJ(1));
t229 = cos(qJ(2));
t230 = cos(qJ(1));
t259 = cos(pkin(6));
t241 = t230 * t259;
t211 = t226 * t225 - t229 * t241;
t212 = t225 * t241 + t226 * t229;
t224 = sin(qJ(3));
t228 = cos(qJ(3));
t219 = sin(pkin(7));
t220 = sin(pkin(6));
t253 = t220 * t230;
t243 = t219 * t253;
t222 = cos(pkin(7));
t250 = t222 * t224;
t196 = t211 * t250 - t212 * t228 + t224 * t243;
t223 = sin(qJ(4));
t227 = cos(qJ(4));
t195 = (t211 * t222 + t243) * t228 + t212 * t224;
t205 = -t211 * t219 + t222 * t253;
t218 = sin(pkin(8));
t221 = cos(pkin(8));
t238 = t195 * t221 + t205 * t218;
t264 = t196 * t227 + t238 * t223;
t263 = -t196 * t223 + t238 * t227;
t257 = t218 * t219;
t256 = t219 * t220;
t255 = t219 * t221;
t254 = t220 * t226;
t252 = t221 * t223;
t251 = t221 * t227;
t249 = t222 * t228;
t248 = t224 * t225;
t247 = t224 * t229;
t246 = t225 * t228;
t245 = t228 * t229;
t244 = t225 * t256;
t242 = t226 * t259;
t240 = t259 * t219;
t214 = -t225 * t242 + t230 * t229;
t213 = -t230 * t225 - t229 * t242;
t232 = t213 * t222 + t219 * t254;
t197 = -t214 * t224 + t232 * t228;
t207 = -t213 * t219 + t222 * t254;
t237 = t197 * t221 + t207 * t218;
t203 = t228 * t240 + (t222 * t245 - t248) * t220;
t236 = t203 * t221 + (t259 * t222 - t229 * t256) * t218;
t199 = t211 * t224 - t212 * t249;
t235 = t199 * t221 + t212 * t257;
t201 = -t213 * t224 - t214 * t249;
t234 = t201 * t221 + t214 * t257;
t208 = (-t222 * t246 - t247) * t220;
t231 = t208 * t221 + t218 * t244;
t209 = (-t222 * t248 + t245) * t220;
t204 = t224 * t240 + (t222 * t247 + t246) * t220;
t202 = t213 * t228 - t214 * t250;
t200 = -t211 * t228 - t212 * t250;
t198 = t214 * t228 + t232 * t224;
t192 = t198 * t227 + t237 * t223;
t191 = -t198 * t223 + t237 * t227;
t1 = [t264, t202 * t227 + t234 * t223, t197 * t227 - t198 * t252, t191, 0, 0; t192, t200 * t227 + t235 * t223, -t195 * t227 + t196 * t252, -t263, 0, 0; 0, t209 * t227 + t231 * t223, t203 * t227 - t204 * t252, -t204 * t223 + t236 * t227, 0, 0; t263, -t202 * t223 + t234 * t227, -t197 * t223 - t198 * t251, -t192, 0, 0; t191, -t200 * t223 + t235 * t227, t195 * t223 + t196 * t251, t264, 0, 0; 0, -t209 * t223 + t231 * t227, -t203 * t223 - t204 * t251, -t204 * t227 - t236 * t223, 0, 0; -t195 * t218 + t205 * t221, -t201 * t218 + t214 * t255, t198 * t218, 0, 0, 0; -t197 * t218 + t207 * t221, -t199 * t218 + t212 * t255, -t196 * t218, 0, 0, 0; 0, -t208 * t218 + t221 * t244, t204 * t218, 0, 0, 0;];
JR_rot  = t1;
