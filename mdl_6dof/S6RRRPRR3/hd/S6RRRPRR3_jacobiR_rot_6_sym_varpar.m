% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:08
% EndTime: 2019-02-26 22:17:08
% DurationCPUTime: 0.06s
% Computational Cost: add. (136->27), mult. (170->24), div. (0->0), fcn. (256->8), ass. (0->26)
t151 = qJ(2) + qJ(3);
t150 = cos(t151);
t153 = sin(qJ(5));
t166 = sin(t151);
t167 = cos(qJ(5));
t144 = -t150 * t153 + t166 * t167;
t143 = t150 * t167 + t166 * t153;
t154 = sin(qJ(1));
t138 = t144 * t154;
t152 = sin(qJ(6));
t165 = t138 * t152;
t155 = cos(qJ(6));
t136 = t138 * t155;
t156 = cos(qJ(1));
t141 = t144 * t156;
t164 = t141 * t152;
t137 = t141 * t155;
t163 = t143 * t152;
t142 = t143 * t155;
t139 = t143 * t154;
t158 = -t139 * t155 - t156 * t152;
t157 = t139 * t152 - t156 * t155;
t140 = t143 * t156;
t135 = t140 * t155 - t154 * t152;
t134 = -t140 * t152 - t154 * t155;
t1 = [t158, -t137, -t137, 0, t137, t134; t135, -t136, -t136, 0, t136, -t157; 0, t142, t142, 0, -t142, -t144 * t152; t157, t164, t164, 0, -t164, -t135; t134, t165, t165, 0, -t165, t158; 0, -t163, -t163, 0, t163, -t144 * t155; t138, -t140, -t140, 0, t140, 0; -t141, -t139, -t139, 0, t139, 0; 0, -t144, -t144, 0, t144, 0;];
JR_rot  = t1;
