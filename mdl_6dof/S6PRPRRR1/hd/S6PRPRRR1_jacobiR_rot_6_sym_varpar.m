% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:37
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:37:49
% EndTime: 2019-02-22 09:37:49
% DurationCPUTime: 0.12s
% Computational Cost: add. (204->32), mult. (409->65), div. (0->0), fcn. (583->12), ass. (0->40)
t180 = sin(pkin(12));
t183 = cos(pkin(12));
t187 = sin(qJ(2));
t189 = cos(qJ(2));
t173 = t187 * t180 - t189 * t183;
t179 = qJ(4) + qJ(5);
t177 = sin(t179);
t178 = cos(t179);
t185 = cos(pkin(6));
t191 = t189 * t180 + t187 * t183;
t172 = t191 * t185;
t181 = sin(pkin(11));
t184 = cos(pkin(11));
t193 = t184 * t172 - t181 * t173;
t182 = sin(pkin(6));
t196 = t182 * t184;
t156 = -t177 * t193 - t178 * t196;
t186 = sin(qJ(6));
t202 = t156 * t186;
t192 = -t181 * t172 - t184 * t173;
t197 = t181 * t182;
t158 = -t177 * t192 + t178 * t197;
t201 = t158 * t186;
t171 = t191 * t182;
t167 = -t171 * t177 + t185 * t178;
t200 = t167 * t186;
t199 = t178 * t186;
t188 = cos(qJ(6));
t198 = t178 * t188;
t190 = t173 * t185;
t170 = t173 * t182;
t168 = t171 * t178 + t185 * t177;
t166 = t167 * t188;
t164 = t181 * t190 - t184 * t191;
t161 = -t181 * t191 - t184 * t190;
t159 = t177 * t197 + t178 * t192;
t157 = -t177 * t196 + t178 * t193;
t155 = t158 * t188;
t154 = t156 * t188;
t1 = [0, t164 * t198 + t186 * t192, 0, t155, t155, -t159 * t186 - t164 * t188; 0, t161 * t198 + t186 * t193, 0, t154, t154, -t157 * t186 - t161 * t188; 0, -t170 * t198 + t171 * t186, 0, t166, t166, -t168 * t186 + t170 * t188; 0, -t164 * t199 + t188 * t192, 0, -t201, -t201, -t159 * t188 + t164 * t186; 0, -t161 * t199 + t188 * t193, 0, -t202, -t202, -t157 * t188 + t161 * t186; 0, t170 * t199 + t171 * t188, 0, -t200, -t200, -t168 * t188 - t170 * t186; 0, t164 * t177, 0, t159, t159, 0; 0, t161 * t177, 0, t157, t157, 0; 0, -t170 * t177, 0, t168, t168, 0;];
JR_rot  = t1;
