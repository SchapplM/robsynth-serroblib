% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR13_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:37:56
% EndTime: 2019-02-22 10:37:56
% DurationCPUTime: 0.21s
% Computational Cost: add. (240->48), mult. (692->94), div. (0->0), fcn. (956->14), ass. (0->57)
t218 = cos(pkin(6));
t213 = sin(pkin(12));
t226 = cos(qJ(1));
t231 = t226 * t213;
t216 = cos(pkin(12));
t222 = sin(qJ(1));
t232 = t222 * t216;
t206 = t218 * t231 + t232;
t221 = sin(qJ(3));
t225 = cos(qJ(3));
t230 = t226 * t216;
t233 = t222 * t213;
t205 = -t218 * t230 + t233;
t217 = cos(pkin(7));
t214 = sin(pkin(7));
t215 = sin(pkin(6));
t237 = t215 * t226;
t228 = t214 * t237;
t227 = t205 * t217 + t228;
t194 = -t206 * t225 + t227 * t221;
t219 = sin(qJ(6));
t246 = t194 * t219;
t223 = cos(qJ(6));
t245 = t194 * t223;
t199 = -t205 * t214 + t217 * t237;
t220 = sin(qJ(5));
t244 = t199 * t220;
t224 = cos(qJ(5));
t243 = t199 * t224;
t242 = t206 * t221;
t240 = t214 * t218;
t239 = t215 * t216;
t238 = t215 * t222;
t236 = t217 * t225;
t235 = t219 * t220;
t234 = t220 * t223;
t229 = t214 * t238;
t208 = -t218 * t233 + t230;
t207 = -t218 * t232 - t231;
t204 = -t214 * t239 + t218 * t217;
t201 = -t207 * t214 + t217 * t238;
t198 = t221 * t240 + (t216 * t217 * t221 + t213 * t225) * t215;
t197 = t215 * t213 * t221 - t225 * t240 - t236 * t239;
t196 = t208 * t225 + (t207 * t217 + t229) * t221;
t195 = -t207 * t236 + t208 * t221 - t225 * t229;
t193 = -t227 * t225 - t242;
t191 = t205 * t236 + t225 * t228 + t242;
t190 = t197 * t220 + t204 * t224;
t189 = t197 * t224 - t204 * t220;
t187 = t195 * t220 + t201 * t224;
t186 = -t195 * t224 + t201 * t220;
t185 = t191 * t220 - t243;
t184 = t191 * t224 + t244;
t183 = t193 * t220 + t243;
t182 = t187 * t223 + t196 * t219;
t181 = -t187 * t219 + t196 * t223;
t1 = [t183 * t223 + t246, 0, -t195 * t219 + t196 * t234, 0, -t186 * t223, t181; t182, 0, -t191 * t219 - t194 * t234, 0, t184 * t223, -t185 * t219 - t245; 0, 0, -t197 * t219 + t198 * t234, 0, t189 * t223, -t190 * t219 + t198 * t223; -t183 * t219 + t245, 0, -t195 * t223 - t196 * t235, 0, t186 * t219, -t182; t181, 0, -t191 * t223 + t194 * t235, 0, -t184 * t219, -t185 * t223 + t246; 0, 0, -t197 * t223 - t198 * t235, 0, -t189 * t219, -t190 * t223 - t198 * t219; -t193 * t224 + t244, 0, -t196 * t224, 0, t187, 0; t186, 0, t194 * t224, 0, t185, 0; 0, 0, -t198 * t224, 0, t190, 0;];
JR_rot  = t1;
