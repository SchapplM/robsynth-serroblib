% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:01
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:00:58
% EndTime: 2019-02-22 12:00:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (50->25), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->23)
t165 = sin(qJ(2));
t163 = sin(qJ(5));
t164 = sin(qJ(3));
t167 = cos(qJ(5));
t168 = cos(qJ(3));
t173 = t163 * t168 - t164 * t167;
t180 = t173 * t165;
t170 = cos(qJ(1));
t166 = sin(qJ(1));
t169 = cos(qJ(2));
t177 = t166 * t169;
t157 = t164 * t177 + t170 * t168;
t158 = -t170 * t164 + t168 * t177;
t179 = t157 * t163 + t158 * t167;
t176 = t170 * t169;
t174 = t157 * t167 - t158 * t163;
t159 = t164 * t176 - t166 * t168;
t160 = t166 * t164 + t168 * t176;
t152 = -t159 * t167 + t160 * t163;
t153 = t159 * t163 + t160 * t167;
t172 = t163 * t164 + t167 * t168;
t155 = t172 * t165;
t1 = [-t179, -t170 * t155, t152, 0, -t152, 0; t153, -t166 * t155, -t174, 0, t174, 0; 0, t172 * t169, t180, 0, -t180, 0; t166 * t165, -t176, 0, 0, 0, 0; -t170 * t165, -t177, 0, 0, 0, 0; 0, -t165, 0, 0, 0, 0; t174, -t170 * t180, -t153, 0, t153, 0; t152, -t166 * t180, -t179, 0, t179, 0; 0, t173 * t169, -t155, 0, t155, 0;];
JR_rot  = t1;
