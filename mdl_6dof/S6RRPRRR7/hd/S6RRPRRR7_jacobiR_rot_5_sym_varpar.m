% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:35
% EndTime: 2019-02-26 21:57:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (46->20), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->25)
t115 = sin(qJ(4));
t118 = cos(qJ(2));
t132 = sin(qJ(2));
t133 = cos(qJ(4));
t109 = -t118 * t115 + t132 * t133;
t108 = t132 * t115 + t118 * t133;
t116 = sin(qJ(1));
t104 = t109 * t116;
t114 = sin(qJ(5));
t131 = t104 * t114;
t117 = cos(qJ(5));
t130 = t104 * t117;
t119 = cos(qJ(1));
t107 = t109 * t119;
t129 = t107 * t114;
t128 = t107 * t117;
t127 = t108 * t114;
t126 = t108 * t117;
t105 = t108 * t116;
t121 = -t105 * t117 - t119 * t114;
t120 = t105 * t114 - t119 * t117;
t106 = t108 * t119;
t103 = t106 * t117 - t116 * t114;
t102 = -t106 * t114 - t116 * t117;
t1 = [t121, -t128, 0, t128, t102, 0; t103, -t130, 0, t130, -t120, 0; 0, t126, 0, -t126, -t109 * t114, 0; t120, t129, 0, -t129, -t103, 0; t102, t131, 0, -t131, t121, 0; 0, -t127, 0, t127, -t109 * t117, 0; t104, -t106, 0, t106, 0, 0; -t107, -t105, 0, t105, 0, 0; 0, -t109, 0, t109, 0, 0;];
JR_rot  = t1;
