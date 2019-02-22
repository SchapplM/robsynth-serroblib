% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:01
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:01:28
% EndTime: 2019-02-22 12:01:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (109->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
t173 = cos(pkin(6));
t175 = sin(qJ(2));
t179 = cos(qJ(1));
t181 = t179 * t175;
t176 = sin(qJ(1));
t178 = cos(qJ(2));
t183 = t176 * t178;
t164 = t173 * t181 + t183;
t174 = sin(qJ(3));
t177 = cos(qJ(3));
t172 = sin(pkin(6));
t185 = t172 * t179;
t158 = -t164 * t177 + t174 * t185;
t180 = t179 * t178;
t184 = t176 * t175;
t163 = -t173 * t180 + t184;
t171 = pkin(11) + qJ(5);
t169 = sin(t171);
t170 = cos(t171);
t194 = t158 * t169 + t163 * t170;
t193 = t158 * t170 - t163 * t169;
t190 = t169 * t177;
t189 = t170 * t177;
t188 = t172 * t174;
t187 = t172 * t177;
t186 = t172 * t178;
t182 = t177 * t178;
t156 = -t164 * t174 - t177 * t185;
t166 = -t173 * t184 + t180;
t165 = t173 * t183 + t181;
t162 = t173 * t174 + t175 * t187;
t161 = t173 * t177 - t175 * t188;
t160 = t166 * t177 + t176 * t188;
t159 = t166 * t174 - t176 * t187;
t155 = t160 * t170 + t165 * t169;
t154 = t160 * t169 - t165 * t170;
t1 = [t193, -t165 * t189 + t166 * t169, -t159 * t170, 0, -t154, 0; t155, -t163 * t189 + t164 * t169, t156 * t170, 0, t194, 0; 0 (t169 * t175 + t170 * t182) * t172, t161 * t170, 0, -t162 * t169 - t170 * t186, 0; t156, -t165 * t174, t160, 0, 0, 0; t159, -t163 * t174, -t158, 0, 0, 0; 0, t174 * t186, t162, 0, 0, 0; t194, -t165 * t190 - t166 * t170, -t159 * t169, 0, t155, 0; t154, -t163 * t190 - t164 * t170, t156 * t169, 0, -t193, 0; 0 (t169 * t182 - t170 * t175) * t172, t161 * t169, 0, t162 * t170 - t169 * t186, 0;];
JR_rot  = t1;
