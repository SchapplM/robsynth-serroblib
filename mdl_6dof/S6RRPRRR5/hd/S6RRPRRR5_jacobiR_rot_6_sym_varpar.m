% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
% Datum: 2019-02-22 11:42
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:42:09
% EndTime: 2019-02-22 11:42:09
% DurationCPUTime: 0.17s
% Computational Cost: add. (244->33), mult. (533->65), div. (0->0), fcn. (760->12), ass. (0->40)
t194 = sin(qJ(4));
t197 = cos(qJ(4));
t193 = cos(pkin(6));
t190 = sin(pkin(12));
t192 = cos(pkin(12));
t195 = sin(qJ(2));
t198 = cos(qJ(2));
t201 = t198 * t190 + t195 * t192;
t181 = t201 * t193;
t182 = t195 * t190 - t198 * t192;
t196 = sin(qJ(1));
t199 = cos(qJ(1));
t203 = t199 * t181 - t196 * t182;
t191 = sin(pkin(6));
t206 = t191 * t199;
t167 = t194 * t206 - t197 * t203;
t200 = t182 * t193;
t171 = -t196 * t201 - t199 * t200;
t189 = qJ(5) + qJ(6);
t187 = sin(t189);
t188 = cos(t189);
t159 = t167 * t187 - t171 * t188;
t160 = t167 * t188 + t171 * t187;
t209 = t187 * t197;
t208 = t188 * t197;
t207 = t191 * t196;
t202 = -t196 * t181 - t199 * t182;
t165 = -t194 * t203 - t197 * t206;
t180 = t201 * t191;
t179 = t182 * t191;
t177 = t180 * t197 + t193 * t194;
t176 = -t180 * t194 + t193 * t197;
t174 = t196 * t200 - t199 * t201;
t169 = t194 * t207 + t197 * t202;
t168 = t194 * t202 - t197 * t207;
t164 = -t177 * t188 - t179 * t187;
t163 = -t177 * t187 + t179 * t188;
t162 = t169 * t188 - t174 * t187;
t161 = -t169 * t187 - t174 * t188;
t1 = [t160, t174 * t208 + t187 * t202, 0, -t168 * t188, t161, t161; t162, t171 * t208 + t187 * t203, 0, t165 * t188, t159, t159; 0, -t179 * t208 + t180 * t187, 0, t176 * t188, t163, t163; -t159, -t174 * t209 + t188 * t202, 0, t168 * t187, -t162, -t162; t161, -t171 * t209 + t188 * t203, 0, -t165 * t187, t160, t160; 0, t179 * t209 + t180 * t188, 0, -t176 * t187, t164, t164; t165, t174 * t194, 0, t169, 0, 0; t168, t171 * t194, 0, -t167, 0, 0; 0, -t179 * t194, 0, t177, 0, 0;];
JR_rot  = t1;
