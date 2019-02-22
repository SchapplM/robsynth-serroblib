% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:37
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:37:43
% EndTime: 2019-02-22 12:37:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (222->44), mult. (512->90), div. (0->0), fcn. (717->12), ass. (0->53)
t208 = sin(qJ(2));
t209 = sin(qJ(1));
t211 = cos(qJ(2));
t212 = cos(qJ(1));
t231 = cos(pkin(6));
t214 = t212 * t231;
t195 = t208 * t214 + t209 * t211;
t207 = sin(qJ(3));
t210 = cos(qJ(3));
t194 = t208 * t209 - t211 * t214;
t204 = sin(pkin(7));
t206 = cos(pkin(7));
t205 = sin(pkin(6));
t225 = t205 * t212;
t213 = t194 * t206 + t204 * t225;
t181 = -t195 * t210 + t207 * t213;
t188 = -t194 * t204 + t206 * t225;
t203 = qJ(4) + qJ(5);
t201 = sin(t203);
t202 = cos(t203);
t173 = t181 * t201 - t188 * t202;
t174 = t181 * t202 + t188 * t201;
t229 = t201 * t204;
t228 = t202 * t204;
t227 = t204 * t205;
t226 = t205 * t209;
t224 = t206 * t207;
t223 = t206 * t210;
t222 = t207 * t208;
t221 = t207 * t211;
t220 = t208 * t210;
t219 = t210 * t211;
t218 = t208 * t227;
t217 = t204 * t226;
t216 = t204 * t231;
t215 = t209 * t231;
t179 = -t195 * t207 - t210 * t213;
t197 = -t208 * t215 + t211 * t212;
t196 = -t212 * t208 - t211 * t215;
t193 = t231 * t206 - t211 * t227;
t192 = (-t206 * t222 + t219) * t205;
t190 = -t196 * t204 + t206 * t226;
t187 = t207 * t216 + (t206 * t221 + t220) * t205;
t186 = t210 * t216 + (t206 * t219 - t222) * t205;
t185 = t196 * t210 - t197 * t224;
t184 = -t194 * t210 - t195 * t224;
t183 = t197 * t210 + (t196 * t206 + t217) * t207;
t182 = -t196 * t223 + t197 * t207 - t210 * t217;
t178 = -t187 * t202 - t193 * t201;
t177 = -t187 * t201 + t193 * t202;
t176 = t183 * t202 + t190 * t201;
t175 = -t183 * t201 + t190 * t202;
t1 = [t174, t185 * t202 + t197 * t229, -t182 * t202, t175, t175, 0; t176, t184 * t202 + t195 * t229, t179 * t202, t173, t173, 0; 0, t192 * t202 + t201 * t218, t186 * t202, t177, t177, 0; -t173, -t185 * t201 + t197 * t228, t182 * t201, -t176, -t176, 0; t175, -t184 * t201 + t195 * t228, -t179 * t201, t174, t174, 0; 0, -t192 * t201 + t202 * t218, -t186 * t201, t178, t178, 0; t179, t196 * t207 + t197 * t223, t183, 0, 0, 0; t182, -t194 * t207 + t195 * t223, -t181, 0, 0, 0; 0 (t206 * t220 + t221) * t205, t187, 0, 0, 0;];
JR_rot  = t1;
