% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
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
% Datum: 2019-02-22 12:31
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:31:29
% EndTime: 2019-02-22 12:31:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (141->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
t196 = cos(pkin(6));
t198 = sin(qJ(2));
t202 = cos(qJ(1));
t204 = t202 * t198;
t199 = sin(qJ(1));
t201 = cos(qJ(2));
t206 = t199 * t201;
t187 = t196 * t204 + t206;
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t195 = sin(pkin(6));
t208 = t195 * t202;
t180 = -t187 * t200 + t197 * t208;
t203 = t202 * t201;
t207 = t199 * t198;
t186 = -t196 * t203 + t207;
t194 = qJ(4) + qJ(5);
t192 = sin(t194);
t193 = cos(t194);
t172 = t180 * t192 + t186 * t193;
t216 = t180 * t193 - t186 * t192;
t213 = t192 * t200;
t212 = t193 * t200;
t211 = t195 * t197;
t210 = t195 * t200;
t209 = t195 * t201;
t205 = t200 * t201;
t178 = -t187 * t197 - t200 * t208;
t189 = -t196 * t207 + t203;
t188 = t196 * t206 + t204;
t185 = t196 * t197 + t198 * t210;
t184 = t196 * t200 - t198 * t211;
t182 = t189 * t200 + t199 * t211;
t181 = t189 * t197 - t199 * t210;
t177 = t185 * t193 - t192 * t209;
t176 = -t185 * t192 - t193 * t209;
t175 = t182 * t193 + t188 * t192;
t174 = t182 * t192 - t188 * t193;
t1 = [t216, -t188 * t212 + t189 * t192, -t181 * t193, -t174, -t174, 0; t175, -t186 * t212 + t187 * t192, t178 * t193, t172, t172, 0; 0 (t192 * t198 + t193 * t205) * t195, t184 * t193, t176, t176, 0; t178, -t188 * t197, t182, 0, 0, 0; t181, -t186 * t197, -t180, 0, 0, 0; 0, t197 * t209, t185, 0, 0, 0; t172, -t188 * t213 - t189 * t193, -t181 * t192, t175, t175, 0; t174, -t186 * t213 - t187 * t193, t178 * t192, -t216, -t216, 0; 0 (t192 * t205 - t193 * t198) * t195, t184 * t192, t177, t177, 0;];
JR_rot  = t1;
