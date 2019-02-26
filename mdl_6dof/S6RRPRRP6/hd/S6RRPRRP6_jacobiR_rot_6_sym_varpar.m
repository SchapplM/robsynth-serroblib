% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:50
% EndTime: 2019-02-26 21:48:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (151->30), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
t202 = sin(qJ(4));
t206 = cos(qJ(4));
t200 = cos(pkin(6));
t197 = sin(pkin(11));
t199 = cos(pkin(11));
t203 = sin(qJ(2));
t207 = cos(qJ(2));
t210 = t207 * t197 + t203 * t199;
t191 = t210 * t200;
t192 = t203 * t197 - t207 * t199;
t204 = sin(qJ(1));
t208 = cos(qJ(1));
t212 = t208 * t191 - t204 * t192;
t198 = sin(pkin(6));
t217 = t198 * t208;
t177 = t202 * t217 - t206 * t212;
t209 = t192 * t200;
t181 = -t204 * t210 - t208 * t209;
t201 = sin(qJ(5));
t205 = cos(qJ(5));
t222 = t177 * t201 - t181 * t205;
t221 = t177 * t205 + t181 * t201;
t218 = t198 * t204;
t216 = t201 * t206;
t214 = t205 * t206;
t211 = -t204 * t191 - t208 * t192;
t175 = -t202 * t212 - t206 * t217;
t190 = t210 * t198;
t189 = t192 * t198;
t187 = t190 * t206 + t200 * t202;
t186 = -t190 * t202 + t200 * t206;
t184 = t204 * t209 - t208 * t210;
t179 = t202 * t218 + t206 * t211;
t178 = t202 * t211 - t206 * t218;
t174 = t179 * t205 - t184 * t201;
t173 = t179 * t201 + t184 * t205;
t1 = [t221, t184 * t214 + t201 * t211, 0, -t178 * t205, -t173, 0; t174, t181 * t214 + t201 * t212, 0, t175 * t205, t222, 0; 0, -t189 * t214 + t190 * t201, 0, t186 * t205, -t187 * t201 + t189 * t205, 0; t175, t184 * t202, 0, t179, 0, 0; t178, t181 * t202, 0, -t177, 0, 0; 0, -t189 * t202, 0, t187, 0, 0; t222, t184 * t216 - t205 * t211, 0, -t178 * t201, t174, 0; t173, t181 * t216 - t205 * t212, 0, t175 * t201, -t221, 0; 0, -t189 * t216 - t190 * t205, 0, t186 * t201, t187 * t205 + t189 * t201, 0;];
JR_rot  = t1;
