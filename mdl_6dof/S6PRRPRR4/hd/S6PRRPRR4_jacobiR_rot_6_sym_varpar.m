% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:49
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:49:46
% EndTime: 2019-02-22 09:49:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (151->39), mult. (430->78), div. (0->0), fcn. (608->12), ass. (0->46)
t192 = sin(pkin(11));
t194 = cos(pkin(11));
t203 = cos(qJ(2));
t195 = cos(pkin(6));
t199 = sin(qJ(2));
t207 = t195 * t199;
t184 = t192 * t203 + t194 * t207;
t198 = sin(qJ(3));
t193 = sin(pkin(6));
t202 = cos(qJ(3));
t209 = t193 * t202;
t178 = t184 * t198 + t194 * t209;
t211 = t193 * t198;
t179 = t184 * t202 - t194 * t211;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t168 = t178 * t201 - t179 * t197;
t196 = sin(qJ(6));
t217 = t168 * t196;
t200 = cos(qJ(6));
t216 = t168 * t200;
t186 = -t192 * t207 + t194 * t203;
t180 = t186 * t198 - t192 * t209;
t181 = t186 * t202 + t192 * t211;
t171 = t180 * t201 - t181 * t197;
t215 = t171 * t196;
t214 = t171 * t200;
t210 = t193 * t199;
t187 = -t195 * t202 + t198 * t210;
t188 = t195 * t198 + t199 * t209;
t176 = t187 * t201 - t188 * t197;
t213 = t176 * t196;
t212 = t176 * t200;
t208 = t193 * t203;
t206 = t195 * t203;
t169 = t178 * t197 + t179 * t201;
t172 = t180 * t197 + t181 * t201;
t177 = t187 * t197 + t188 * t201;
t205 = t197 * t202 - t198 * t201;
t204 = t197 * t198 + t201 * t202;
t185 = -t192 * t206 - t194 * t199;
t183 = -t192 * t199 + t194 * t206;
t182 = t204 * t208;
t174 = t204 * t185;
t173 = t204 * t183;
t1 = [0, t174 * t200 - t186 * t196, -t214, 0, t214, -t172 * t196 + t185 * t200; 0, t173 * t200 - t184 * t196, -t216, 0, t216, -t169 * t196 + t183 * t200; 0, t182 * t200 - t196 * t210, -t212, 0, t212, -t177 * t196 + t200 * t208; 0, -t174 * t196 - t186 * t200, t215, 0, -t215, -t172 * t200 - t185 * t196; 0, -t173 * t196 - t184 * t200, t217, 0, -t217, -t169 * t200 - t183 * t196; 0, -t182 * t196 - t200 * t210, t213, 0, -t213, -t177 * t200 - t196 * t208; 0, t205 * t185, -t172, 0, t172, 0; 0, t205 * t183, -t169, 0, t169, 0; 0, t205 * t208, -t177, 0, t177, 0;];
JR_rot  = t1;
