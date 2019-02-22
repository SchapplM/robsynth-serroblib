% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:41:31
% EndTime: 2019-02-22 11:41:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (256->36), mult. (515->65), div. (0->0), fcn. (735->12), ass. (0->44)
t199 = qJ(4) + qJ(5);
t197 = sin(t199);
t198 = cos(t199);
t203 = cos(pkin(6));
t200 = sin(pkin(12));
t202 = cos(pkin(12));
t205 = sin(qJ(2));
t208 = cos(qJ(2));
t211 = t208 * t200 + t205 * t202;
t191 = t211 * t203;
t192 = t205 * t200 - t208 * t202;
t206 = sin(qJ(1));
t209 = cos(qJ(1));
t213 = t209 * t191 - t206 * t192;
t201 = sin(pkin(6));
t216 = t201 * t209;
t176 = t197 * t216 - t198 * t213;
t210 = t192 * t203;
t180 = -t206 * t211 - t209 * t210;
t204 = sin(qJ(6));
t207 = cos(qJ(6));
t226 = t176 * t204 - t180 * t207;
t225 = t176 * t207 + t180 * t204;
t174 = -t197 * t213 - t198 * t216;
t224 = t174 * t204;
t212 = -t206 * t191 - t209 * t192;
t217 = t201 * t206;
t177 = t197 * t212 - t198 * t217;
t223 = t177 * t204;
t190 = t211 * t201;
t186 = -t190 * t197 + t203 * t198;
t220 = t186 * t204;
t219 = t198 * t204;
t218 = t198 * t207;
t189 = t192 * t201;
t187 = t190 * t198 + t203 * t197;
t185 = t186 * t207;
t183 = t206 * t210 - t209 * t211;
t178 = t197 * t217 + t198 * t212;
t173 = t177 * t207;
t172 = t174 * t207;
t171 = t178 * t207 - t183 * t204;
t170 = -t178 * t204 - t183 * t207;
t1 = [t225, t183 * t218 + t204 * t212, 0, -t173, -t173, t170; t171, t180 * t218 + t204 * t213, 0, t172, t172, t226; 0, -t189 * t218 + t190 * t204, 0, t185, t185, -t187 * t204 + t189 * t207; -t226, -t183 * t219 + t207 * t212, 0, t223, t223, -t171; t170, -t180 * t219 + t207 * t213, 0, -t224, -t224, t225; 0, t189 * t219 + t190 * t207, 0, -t220, -t220, -t187 * t207 - t189 * t204; t174, t183 * t197, 0, t178, t178, 0; t177, t180 * t197, 0, -t176, -t176, 0; 0, -t189 * t197, 0, t187, t187, 0;];
JR_rot  = t1;
