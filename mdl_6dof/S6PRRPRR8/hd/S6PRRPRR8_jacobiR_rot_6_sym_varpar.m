% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (234->61), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->60)
t218 = sin(pkin(7));
t219 = sin(pkin(6));
t247 = t218 * t219;
t224 = sin(qJ(5));
t246 = t218 * t224;
t228 = cos(qJ(5));
t245 = t218 * t228;
t229 = cos(qJ(3));
t244 = t218 * t229;
t221 = cos(pkin(7));
t243 = t219 * t221;
t225 = sin(qJ(3));
t242 = t221 * t225;
t241 = t221 * t229;
t222 = cos(pkin(6));
t226 = sin(qJ(2));
t240 = t222 * t226;
t230 = cos(qJ(2));
t239 = t222 * t230;
t223 = sin(qJ(6));
t238 = t223 * t224;
t227 = cos(qJ(6));
t237 = t224 * t227;
t236 = t225 * t226;
t235 = t225 * t230;
t234 = t226 * t229;
t233 = t229 * t230;
t232 = t226 * t247;
t231 = t219 * t244;
t220 = cos(pkin(12));
t217 = sin(pkin(12));
t212 = -t217 * t240 + t220 * t230;
t211 = -t217 * t239 - t220 * t226;
t210 = t217 * t230 + t220 * t240;
t209 = -t217 * t226 + t220 * t239;
t208 = t222 * t221 - t230 * t247;
t207 = (-t221 * t236 + t233) * t219;
t206 = (t221 * t234 + t235) * t219;
t203 = -t211 * t218 + t217 * t243;
t202 = -t209 * t218 - t220 * t243;
t201 = t222 * t218 * t225 + (t221 * t235 + t234) * t219;
t200 = t219 * t236 - t222 * t244 - t233 * t243;
t199 = t206 * t224 + t228 * t232;
t198 = t211 * t229 - t212 * t242;
t197 = t211 * t225 + t212 * t241;
t196 = t209 * t229 - t210 * t242;
t195 = t209 * t225 + t210 * t241;
t194 = t200 * t224 + t208 * t228;
t193 = t200 * t228 - t208 * t224;
t192 = t212 * t229 + (t211 * t221 + t217 * t247) * t225;
t191 = -t211 * t241 + t212 * t225 - t217 * t231;
t190 = t210 * t229 + (t209 * t221 - t220 * t247) * t225;
t189 = -t209 * t241 + t210 * t225 + t220 * t231;
t188 = t197 * t224 + t212 * t245;
t187 = t195 * t224 + t210 * t245;
t186 = t191 * t224 + t203 * t228;
t185 = t191 * t228 - t203 * t224;
t184 = t189 * t224 + t202 * t228;
t183 = t189 * t228 - t202 * t224;
t1 = [0, t188 * t227 + t198 * t223, -t191 * t223 + t192 * t237, 0, t185 * t227, -t186 * t223 + t192 * t227; 0, t187 * t227 + t196 * t223, -t189 * t223 + t190 * t237, 0, t183 * t227, -t184 * t223 + t190 * t227; 0, t199 * t227 + t207 * t223, -t200 * t223 + t201 * t237, 0, t193 * t227, -t194 * t223 + t201 * t227; 0, -t188 * t223 + t198 * t227, -t191 * t227 - t192 * t238, 0, -t185 * t223, -t186 * t227 - t192 * t223; 0, -t187 * t223 + t196 * t227, -t189 * t227 - t190 * t238, 0, -t183 * t223, -t184 * t227 - t190 * t223; 0, -t199 * t223 + t207 * t227, -t200 * t227 - t201 * t238, 0, -t193 * t223, -t194 * t227 - t201 * t223; 0, -t197 * t228 + t212 * t246, -t192 * t228, 0, t186, 0; 0, -t195 * t228 + t210 * t246, -t190 * t228, 0, t184, 0; 0, -t206 * t228 + t224 * t232, -t201 * t228, 0, t194, 0;];
JR_rot  = t1;
