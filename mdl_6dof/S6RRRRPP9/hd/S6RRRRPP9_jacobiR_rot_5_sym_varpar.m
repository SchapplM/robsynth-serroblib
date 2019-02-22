% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:17:00
% EndTime: 2019-02-22 12:17:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (74->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t161 = cos(pkin(6));
t164 = sin(qJ(2));
t169 = cos(qJ(1));
t171 = t169 * t164;
t165 = sin(qJ(1));
t168 = cos(qJ(2));
t174 = t165 * t168;
t155 = t161 * t171 + t174;
t163 = sin(qJ(3));
t167 = cos(qJ(3));
t160 = sin(pkin(6));
t177 = t160 * t169;
t149 = -t155 * t167 + t163 * t177;
t170 = t169 * t168;
t175 = t165 * t164;
t154 = -t161 * t170 + t175;
t162 = sin(qJ(4));
t166 = cos(qJ(4));
t184 = t149 * t162 + t154 * t166;
t183 = -t149 * t166 + t154 * t162;
t180 = t160 * t163;
t179 = t160 * t167;
t178 = t160 * t168;
t176 = t162 * t167;
t173 = t166 * t167;
t172 = t167 * t168;
t147 = -t155 * t163 - t167 * t177;
t157 = -t161 * t175 + t170;
t156 = t161 * t174 + t171;
t153 = t161 * t163 + t164 * t179;
t152 = t161 * t167 - t164 * t180;
t151 = t157 * t167 + t165 * t180;
t150 = t157 * t163 - t165 * t179;
t146 = t151 * t166 + t156 * t162;
t145 = t151 * t162 - t156 * t166;
t1 = [t147, -t156 * t163, t151, 0, 0, 0; t150, -t154 * t163, -t149, 0, 0, 0; 0, t163 * t178, t153, 0, 0, 0; t183, t156 * t173 - t157 * t162, t150 * t166, t145, 0, 0; -t146, t154 * t173 - t155 * t162, -t147 * t166, -t184, 0, 0; 0 (-t162 * t164 - t166 * t172) * t160, -t152 * t166, t153 * t162 + t166 * t178, 0, 0; t184, -t156 * t176 - t157 * t166, -t150 * t162, t146, 0, 0; t145, -t154 * t176 - t155 * t166, t147 * t162, t183, 0, 0; 0 (t162 * t172 - t164 * t166) * t160, t152 * t162, t153 * t166 - t162 * t178, 0, 0;];
JR_rot  = t1;
