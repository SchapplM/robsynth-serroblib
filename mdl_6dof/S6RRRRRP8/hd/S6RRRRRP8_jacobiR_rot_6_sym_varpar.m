% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:30
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:30:20
% EndTime: 2019-02-22 12:30:20
% DurationCPUTime: 0.15s
% Computational Cost: add. (155->33), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
t211 = cos(pkin(6));
t213 = sin(qJ(2));
t217 = cos(qJ(1));
t219 = t217 * t213;
t214 = sin(qJ(1));
t216 = cos(qJ(2));
t221 = t214 * t216;
t202 = t211 * t219 + t221;
t209 = qJ(3) + qJ(4);
t207 = sin(t209);
t208 = cos(t209);
t210 = sin(pkin(6));
t224 = t210 * t217;
t194 = -t202 * t208 + t207 * t224;
t218 = t217 * t216;
t222 = t214 * t213;
t201 = -t211 * t218 + t222;
t212 = sin(qJ(5));
t215 = cos(qJ(5));
t232 = t194 * t212 + t201 * t215;
t231 = t194 * t215 - t201 * t212;
t228 = t208 * t212;
t227 = t208 * t215;
t226 = t210 * t213;
t225 = t210 * t214;
t223 = t212 * t216;
t220 = t215 * t216;
t192 = -t202 * t207 - t208 * t224;
t204 = -t211 * t222 + t218;
t203 = t211 * t221 + t219;
t200 = t211 * t207 + t208 * t226;
t199 = -t207 * t226 + t211 * t208;
t198 = t199 * t215;
t197 = t199 * t212;
t196 = t204 * t208 + t207 * t225;
t195 = t204 * t207 - t208 * t225;
t191 = t195 * t215;
t190 = t195 * t212;
t189 = t192 * t215;
t188 = t192 * t212;
t187 = t196 * t215 + t203 * t212;
t186 = t196 * t212 - t203 * t215;
t1 = [t231, -t203 * t227 + t204 * t212, -t191, -t191, -t186, 0; t187, -t201 * t227 + t202 * t212, t189, t189, t232, 0; 0 (t208 * t220 + t212 * t213) * t210, t198, t198, -t200 * t212 - t210 * t220, 0; t192, -t203 * t207, t196, t196, 0, 0; t195, -t201 * t207, -t194, -t194, 0, 0; 0, t210 * t216 * t207, t200, t200, 0, 0; t232, -t203 * t228 - t204 * t215, -t190, -t190, t187, 0; t186, -t201 * t228 - t202 * t215, t188, t188, -t231, 0; 0 (t208 * t223 - t213 * t215) * t210, t197, t197, t200 * t215 - t210 * t223, 0;];
JR_rot  = t1;
