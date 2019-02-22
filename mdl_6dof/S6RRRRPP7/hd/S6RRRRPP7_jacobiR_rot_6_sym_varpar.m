% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:15
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:15:46
% EndTime: 2019-02-22 12:15:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (109->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
t174 = cos(pkin(6));
t176 = sin(qJ(2));
t180 = cos(qJ(1));
t182 = t180 * t176;
t177 = sin(qJ(1));
t179 = cos(qJ(2));
t184 = t177 * t179;
t165 = t174 * t182 + t184;
t175 = sin(qJ(3));
t178 = cos(qJ(3));
t173 = sin(pkin(6));
t186 = t173 * t180;
t159 = -t165 * t178 + t175 * t186;
t181 = t180 * t179;
t185 = t177 * t176;
t164 = -t174 * t181 + t185;
t172 = qJ(4) + pkin(11);
t170 = sin(t172);
t171 = cos(t172);
t195 = t159 * t170 + t164 * t171;
t194 = t159 * t171 - t164 * t170;
t191 = t170 * t178;
t190 = t171 * t178;
t189 = t173 * t175;
t188 = t173 * t178;
t187 = t173 * t179;
t183 = t178 * t179;
t157 = -t165 * t175 - t178 * t186;
t167 = -t174 * t185 + t181;
t166 = t174 * t184 + t182;
t163 = t174 * t175 + t176 * t188;
t162 = t174 * t178 - t176 * t189;
t161 = t167 * t178 + t177 * t189;
t160 = t167 * t175 - t177 * t188;
t156 = t161 * t171 + t166 * t170;
t155 = t161 * t170 - t166 * t171;
t1 = [t194, -t166 * t190 + t167 * t170, -t160 * t171, -t155, 0, 0; t156, -t164 * t190 + t165 * t170, t157 * t171, t195, 0, 0; 0 (t170 * t176 + t171 * t183) * t173, t162 * t171, -t163 * t170 - t171 * t187, 0, 0; t157, -t166 * t175, t161, 0, 0, 0; t160, -t164 * t175, -t159, 0, 0, 0; 0, t175 * t187, t163, 0, 0, 0; t195, -t166 * t191 - t167 * t171, -t160 * t170, t156, 0, 0; t155, -t164 * t191 - t165 * t171, t157 * t170, -t194, 0, 0; 0 (t170 * t183 - t171 * t176) * t173, t162 * t170, t163 * t171 - t170 * t187, 0, 0;];
JR_rot  = t1;
