% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR15_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:25:42
% EndTime: 2019-02-22 12:25:42
% DurationCPUTime: 0.23s
% Computational Cost: add. (136->42), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->50)
t209 = sin(qJ(2));
t210 = sin(qJ(1));
t213 = cos(qJ(2));
t214 = cos(qJ(1));
t233 = cos(pkin(6));
t216 = t214 * t233;
t198 = t209 * t216 + t210 * t213;
t208 = sin(qJ(3));
t212 = cos(qJ(3));
t197 = t209 * t210 - t213 * t216;
t204 = sin(pkin(7));
t206 = cos(pkin(7));
t205 = sin(pkin(6));
t227 = t205 * t214;
t215 = t197 * t206 + t204 * t227;
t184 = -t198 * t212 + t215 * t208;
t191 = -t197 * t204 + t206 * t227;
t207 = sin(qJ(4));
t211 = cos(qJ(4));
t237 = t184 * t207 - t191 * t211;
t236 = -t184 * t211 - t191 * t207;
t231 = t204 * t205;
t230 = t204 * t207;
t229 = t204 * t211;
t228 = t205 * t210;
t226 = t206 * t208;
t225 = t206 * t212;
t224 = t208 * t209;
t223 = t208 * t213;
t222 = t209 * t212;
t221 = t212 * t213;
t220 = t209 * t231;
t219 = t204 * t228;
t218 = t204 * t233;
t217 = t210 * t233;
t182 = -t198 * t208 - t215 * t212;
t200 = -t209 * t217 + t213 * t214;
t199 = -t214 * t209 - t213 * t217;
t196 = t233 * t206 - t213 * t231;
t195 = (-t206 * t224 + t221) * t205;
t193 = -t199 * t204 + t206 * t228;
t190 = t208 * t218 + (t206 * t223 + t222) * t205;
t189 = t212 * t218 + (t206 * t221 - t224) * t205;
t188 = t199 * t212 - t200 * t226;
t187 = -t197 * t212 - t198 * t226;
t186 = t200 * t212 + (t199 * t206 + t219) * t208;
t185 = -t199 * t225 + t200 * t208 - t212 * t219;
t181 = t186 * t211 + t193 * t207;
t180 = t186 * t207 - t193 * t211;
t1 = [t182, t199 * t208 + t200 * t225, t186, 0, 0, 0; t185, -t197 * t208 + t198 * t225, -t184, 0, 0, 0; 0 (t206 * t222 + t223) * t205, t190, 0, 0, 0; t236, -t188 * t211 - t200 * t230, t185 * t211, t180, 0, 0; -t181, -t187 * t211 - t198 * t230, -t182 * t211, -t237, 0, 0; 0, -t195 * t211 - t207 * t220, -t189 * t211, t190 * t207 - t196 * t211, 0, 0; t237, t188 * t207 - t200 * t229, -t185 * t207, t181, 0, 0; t180, t187 * t207 - t198 * t229, t182 * t207, t236, 0, 0; 0, t195 * t207 - t211 * t220, t189 * t207, t190 * t211 + t196 * t207, 0, 0;];
JR_rot  = t1;
